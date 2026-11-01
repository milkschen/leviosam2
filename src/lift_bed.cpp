/*
 * lift_bed.cpp
 *
 * Lifting over genomic intervals in the BED format
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "lift_bed.hpp"

#include <getopt.h>

#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iostream>

#include "leviosam_utils.hpp"
#include "version.hpp"

void print_lift_bed_help_msg() {
    std::cerr << "\nLift over a BED file\n";
    std::cerr << "Version: " << VERSION << "\n";
    std::cerr << "Usage:   leviosam2 bed [options] -b <bed> -C <clft> -p "
                 "<prefix>\n\n";
    std::cerr << "Inputs:  -b string   Path to the input BED.\n";
    std::cerr << "         -C string   Path to an indexed ChainMap. See "
                 "`leviosam2 index` for details.\n";
    std::cerr << "         -p string   Prefix to the output files.\n";
    std::cerr << "Options: -h          Print detailed usage.\n";
    std::cerr << "         -G INT      Maximum endpoint clipping, each source or target gap,\n"
                 "                     and net length change (bases), inclusive [500].\n";
    std::cerr << "         -V INT      Verbose level [0].\n";
    std::cerr << "         BED3-6 and BED12+ are supported; BED7-11 are unmapped.\n";
    std::cerr << "\n";
}

namespace {
struct BedProjection {
    std::string target;
    int64_t start = 0, end = 0;
    bool forward = true;
};

bool number(const std::string& text, int64_t& result) {
    if (text.empty() || !std::all_of(text.begin(), text.end(),
        [](char c) { return c >= '0' && c <= '9'; })) return false;
    try {
        size_t used = 0;
        result = std::stoll(text, &used);
        return used == text.size() && result >= 0 &&
               result <= std::numeric_limits<int32_t>::max();
    } catch (...) { return false; }
}

// Boundaries, unlike bases, need no -1 correction on the reverse strand.
int64_t boundary(const chain::Interval& block, int64_t pos) {
    return block.strand ? pos + block.offset :
        int64_t(block.offset) + block.source_start + block.source_end - pos;
}

bool project(const chain::ChainMap& cmap, const std::string& contig,
             int64_t start, int64_t end, int allowed_gaps, BedProjection& out) {
    const auto* blocks = cmap.intervals_for(contig);
    if (!blocks || start > end || allowed_gaps < 0) return false;
    auto first_block = std::lower_bound(blocks->begin(), blocks->end(), start,
        [](const chain::Interval& block, int64_t pos) { return block.source_end <= pos; });
    // A zero-length BED record denotes an insertion point; require it inside a
    // block rather than guessing at an ambiguous gap or chain boundary.
    if (start == end) {
        if (first_block != blocks->end() && first_block->source_start < start &&
            start < first_block->source_end) {
            out = {first_block->target, boundary(*first_block, start),
                   boundary(*first_block, start), first_block->strand};
            return out.start >= 0;
        }
        return false;
    }
    const chain::Interval* prev = nullptr;
    int64_t first = 0, last = 0, lo = 0, hi = 0;
    for (auto it = first_block; it != blocks->end(); ++it) {
        const auto& block = *it;
        if (block.source_start >= end) break;
        int64_t a = std::max(start, int64_t(block.source_start));
        int64_t b = std::min(end, int64_t(block.source_end));
        int64_t x = boundary(block, a), y = boundary(block, b);
        int64_t l = std::min(x,y), h = std::max(x,y);
        if (!prev) {
            out.target = block.target; out.forward = block.strand;
            first = a; lo = l; hi = h;
        } else {
            if (block.target != out.target || block.strand != out.forward) return false;
            int64_t source_gap = int64_t(block.source_start) - prev->source_end;
            int64_t target_gap = out.forward ? l - hi : lo - h;
            if (source_gap < 0 || source_gap > allowed_gaps || target_gap < 0 ||
                target_gap > allowed_gaps) return false;
            lo = std::min(lo,l); hi = std::max(hi,h);
        }
        last = b; prev = &block;
    }
    if (!prev || first - start > allowed_gaps || end - last > allowed_gaps || lo < 0)
        return false;
    // Bound the net length change across the aligned portion as well.
    if (std::llabs((hi - lo) - (last - first)) > allowed_gaps) return false;
    out.start = lo; out.end = hi;
    return true;
}

bool parse_list(const std::string& value, std::vector<int64_t>& result) {
    std::stringstream stream(value);
    std::string token;
    while (std::getline(stream, token, ',')) {
        int64_t n;
        if (!number(token, n)) return false;
        result.push_back(n);
    }
    return !result.empty();
}

std::string list_string(const std::vector<int64_t>& values) {
    std::string result;
    for (auto value : values) result += std::to_string(value) + ",";
    return result;
}
}

int lift_bed_record(chain::ChainMap& cmap, std::string& line,
                    const int& allowed_gaps) {
    auto fields = LevioSamUtils::str_to_vector(line, "\t");
    auto reject = [&](const char* reason) { line += "\t# "; line += reason; line += "\n"; return 0; };
    int64_t start, end;
    if (fields.size() < 3 || !number(fields[1], start) || !number(fields[2], end) || start > end)
        return reject("invalid BED coordinates or shape");
    if (fields.size() > 6 && fields.size() < 12) return reject("BED7-11 unsupported");
    BedProjection outer;
    if (!project(cmap, fields[0], start, end, allowed_gaps, outer)) return reject("unmappable interval or gap exceeds -G");
    if (fields.size() >= 6 && fields[5] != "+" && fields[5] != "-" && fields[5] != ".")
        return reject("invalid or unmappable BED fields");
    if (fields.size() >= 12) {
        int64_t thick_start, thick_end, count;
        std::vector<int64_t> sizes, starts;
        if (!number(fields[6], thick_start) || !number(fields[7], thick_end) ||
            !number(fields[9], count) || count == 0 || thick_start < start ||
            thick_end < thick_start || thick_end > end ||
            !parse_list(fields[10], sizes) || !parse_list(fields[11], starts) ||
            int64_t(sizes.size()) != count || int64_t(starts.size()) != count)
            return reject("invalid or unmappable BED fields");
        std::vector<std::pair<int64_t,int64_t>> mapped;
        int64_t prev_end = start;
        for (size_t i = 0; i < sizes.size(); ++i) {
            int64_t a = start + starts[i], b = a + sizes[i];
            if (sizes[i] == 0 || a < prev_end || b > end) return reject("invalid or unmappable BED fields");
            BedProjection part;
            if (!project(cmap, fields[0], a, b, allowed_gaps, part) ||
                part.target != outer.target || part.forward != outer.forward ||
                part.start < outer.start || part.end > outer.end) return reject("invalid or unmappable BED fields");
            mapped.emplace_back(part.start, part.end); prev_end = b;
        }
        if (starts.front() != 0 || prev_end != end)
            return reject("BED12 blocks must cover the outer boundaries");
        std::sort(mapped.begin(), mapped.end());
        if (mapped.front().first != outer.start || mapped.back().second != outer.end)
            return reject("mapped BED12 blocks do not cover the outer boundaries");
        std::vector<int64_t> new_sizes, new_starts;
        for (const auto& block : mapped) {
            if (!new_starts.empty() && block.first < outer.start + new_starts.back() + new_sizes.back()) return reject("invalid or unmappable BED fields");
            new_starts.push_back(block.first - outer.start);
            new_sizes.push_back(block.second - block.first);
        }
        // Empty coding range denotes noncoding; anchor it at the new start.
        if (thick_start == thick_end) {
            fields[6] = fields[7] = std::to_string(outer.start);
        } else {
            BedProjection thick;
            if (!project(cmap, fields[0], thick_start, thick_end, allowed_gaps, thick) ||
                thick.target != outer.target || thick.forward != outer.forward) return reject("invalid or unmappable BED fields");
            fields[6] = std::to_string(thick.start); fields[7] = std::to_string(thick.end);
        }
        fields[10] = list_string(new_sizes); fields[11] = list_string(new_starts);
    }
    fields[0] = outer.target; fields[1] = std::to_string(outer.start);
    fields[2] = std::to_string(outer.end);
    if (fields.size() >= 6 && !outer.forward) {
        if (fields[5] == "+") fields[5] = "-";
        else if (fields[5] == "-") fields[5] = "+";
    }
    line.clear();
    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) line += "\t";
        line += fields[i];
    }
    line += "\n";
    return 1;
}

void lift_bed(lift_bed_opts& args) {
    chain::ChainMap cmap = [&] {
        if (args.chainmap_fname != "") {
            std::cerr << "Loading levioSAM index...";
            std::ifstream in(args.chainmap_fname, std::ios::binary);
            if (!in) {
                std::cerr << "[E::lift_bed] Cannot open chain index " << args.chainmap_fname << "\n";
                exit(1);
            }
            return chain::ChainMap(in, args.verbose, args.allowed_gaps);
        } else {
            std::cerr << "[E::lift_bed] Failed to load chain index. Exit.\n";
            print_lift_bed_help_msg();
            exit(1);
        }
    }();
    std::ofstream out_f(args.outpre + ".bed");
    std::ofstream unmapped_f(args.outpre + "-unmapped.bed");

    std::ifstream bed_f(args.bed_fname);
    if (!bed_f || !out_f || !unmapped_f) {
        std::cerr << "[E::lift_bed] Cannot open BED input or output files\n";
        exit(1);
    }
    int cnt = 0;
    if (bed_f.is_open()) {
        std::string line;
        while (std::getline(bed_f, line)) {
            if (line.empty() || line[0] == '#' || line.compare(0, 6, "track ") == 0 ||
                line.compare(0, 8, "browser ") == 0) continue;
            cnt += 1;
            // std::cerr << line << "\n";
            int lift = lift_bed_record(cmap, line, args.allowed_gaps);
            if (lift == 1) {
                out_f << line;
            } else if (lift == 0) {
                unmapped_f << line;
            } else {
                std::cerr
                    << "[E::lift_bed] Unexpected liftover outcome at record ";
                std::cerr << cnt << "\n";
                exit(1);
            }
        }
        if (bed_f.bad() || !out_f || !unmapped_f) {
            std::cerr << "[E::lift_bed] Error reading or writing BED records\n";
            exit(1);
        }
        bed_f.close();
        out_f.close();
        unmapped_f.close();
        if (!out_f || !unmapped_f) {
            std::cerr << "[E::lift_bed] Error closing output files\n";
            exit(1);
        }
    }
}

int lift_bed_run(int argc, char** argv) {
    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();
    int c;
    lift_bed_opts args;
    args.cmd = LevioSamUtils::make_cmd(argc, argv);
    static struct option long_options[] {
        {"help", no_argument, 0, 'h'},
        {"bed_fname", required_argument, 0, 'b'},
            {"chainmap_fname", required_argument, 0, 'C'},
            {"allowed_gaps", required_argument, 0, 'G'},
            {"output", required_argument, 0, 'p'},
            {"verbose", required_argument, 0, 'V'},
        {0, 0, 0, 0},
    };
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hb:C:G:p:V:", long_options,
                            &long_index)) != -1) {
        switch (c) {
            case 'h':
                print_lift_bed_help_msg();
                exit(0);
            case 'b':
                args.bed_fname = optarg;
                break;
            case 'C':
                args.chainmap_fname = optarg;
                break;
            case 'G':
                {
                    int64_t value;
                    if (!number(optarg, value)) {
                        std::cerr << "[E::lift_bed] -G must be a nonnegative integer\n";
                        return 1;
                    }
                    args.allowed_gaps = static_cast<int>(value);
                }
                break;
            case 'V': {
                int64_t value;
                if (!number(optarg, value)) {
                    std::cerr << "[E::lift_bed] -V must be a nonnegative integer\n";
                    return 1;
                }
                args.verbose = static_cast<int>(value);
                break;
            }
            case 'p':
                args.outpre = optarg;
                break;
            default:
                fprintf(stderr, "ignoring option %c\n", c);
                exit(1);
        }
    }
    if (args.bed_fname == "") {
        std::cerr << "[E::lift_bed_run] Argument -b/--bed_fname is required\n";
        print_lift_bed_help_msg();
        exit(1);
    }
    if (args.chainmap_fname == "") {
        std::cerr
            << "[E::lift_bed_run] Argument -C/--chainmap_fname is required\n";
        print_lift_bed_help_msg();
        exit(1);
    }

    std::cerr << "Inputs:\n";
    std::cerr << " - BED: " << args.bed_fname << "\n";
    std::cerr << " - clft: " << args.chainmap_fname << "\n";

    lift_bed(args);

    double cpu_duration =
        (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration =
        (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "Finished in " << cpu_duration << " CPU seconds, or "
              << wall_duration.count() << " wall clock seconds\n";

    return 0;
}
