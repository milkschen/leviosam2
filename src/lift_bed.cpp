#include <ctime>
#include <getopt.h>
#include <iostream>
#include "leviosam_utils.hpp"
#include "lift_bed.hpp"
#include "version.hpp"


void print_lift_bed_help_msg() {
    std::cerr << "\nLift over a BED file\n";
    std::cerr << "Version: %s\n", VERSION;
    std::cerr << "Usage:   leviosam bed [options] -b <bed> -C <clft> -p <prefix>\n\n";
    std::cerr << "Inputs:  -b string   Path to the input BED.\n";
    std::cerr << "         -C string   Path to the chain index.\n";
    std::cerr << "         -p string   Prefix to the output files.\n";
    std::cerr << "Options: -h          Print detailed usage.\n";
    std::cerr << "         -G INT      Number of allowed gaps for an interval. [500]\n";
    std::cerr << "         -V INT      Verbose level [0].\n";
    std::cerr << "\n";
}


/* Lift a BED record
 *
 * Returns 1 if it's liftable; 0 if not.
 * The lifted record is updated in the `line` string variable by reference
 */
int lift_bed_record(
    chain::ChainMap& cmap,
    std::string& line, const int& allowed_gaps
) {
    std::vector<std::string> fields = LevioSamUtils::str_to_vector(line, "\t");
    std::string s = fields[0];

    hts_pos_t p1 = std::stoi(fields[1]);
    hts_pos_t op1 = cmap.lift_pos(s, p1, allowed_gaps, true);
    // BED end position is open, so we subtract p2 by one and then
    // add it back after lift-over
    hts_pos_t p2 = std::stoi(fields[2]) - 1;
    hts_pos_t op2 = cmap.lift_pos(s, p2, allowed_gaps, false) + 1;
    
    // Lift contig
    std::string t1 = cmap.lift_contig(s, p1);
    std::string t2 = cmap.lift_contig(s, p2);
    if (op1 < 0 || op2 < 0 || (op2 - op1 <= 0) ||
        t1 == "*" || t2 == "*" || (t1 != t2) ||
        (op2 - op1 - p2 + p1 > allowed_gaps) // gap between (intvl(p1), intvl(p2)) is too large
    ) {
        line += ("\t# " + t1 + ":" + std::to_string(op1) + "-" + t2 + ":" + std::to_string(op2) + "\n");
        return 0;
    } else {
        line = t1 + "\t" + std::to_string(op1) + "\t" + std::to_string(op2);
        if (fields.size() > 3) {
            for (int i = 3; i < fields.size(); i++) {
                line += ("\t" + fields[i]);
            }
        }
        line += "\n";
        return 1;
    }
    // error
    return -1;
}


void lift_bed(lift_bed_opts &args) {
    chain::ChainMap cmap = [&] {
        if (args.chainmap_fname != "") {
            std::cerr << "Loading levioSAM index...";
            std::ifstream in(args.chainmap_fname, std::ios::binary);
            return chain::ChainMap(
                in, args.verbose, args.allowed_gaps
            );
        } else {
            std::cerr << "[E::lift_bed] Failed to load chain index. Exit.\n";
            print_lift_bed_help_msg();
            exit(1);
        }
    }();
    std::ofstream out_f(args.outpre + ".bed");
    std::ofstream unmapped_f(args.outpre + "-unmapped.bed");
    
    std::ifstream bed_f(args.bed_fname);
    int cnt = 0;
    if (bed_f.is_open()) {
        std::string line;
        while (std::getline(bed_f, line)) {
            cnt += 1;
            // std::cerr << line << "\n";
            int lift = lift_bed_record(cmap, line, args.allowed_gaps);
            if (lift == 1) {
                out_f << line;
            } else if (lift == 0) {
                unmapped_f << line;
            } else {
                std::cerr << "[E::lift_bed] Unexpected liftover outcome at record ";
                std::cerr << cnt << "\n";
                exit(1);
            }
        }
        bed_f.close();
    }
}


int lift_bed_run(int argc, char** argv) {
    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();
    int c;
    lift_bed_opts args;
    args.cmd = LevioSamUtils::make_cmd(argc,argv);
    static struct option long_options[] {
        {"bed_fname", required_argument, 0, 'b'},
        {"chainmap_fname", required_argument, 0, 'C'},
        {"allowed_gaps", required_argument, 0, 'G'},
        {"output", required_argument, 0, 'p'},
        {"verbose", required_argument, 0, 'V'},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "h:b:C:G:p:V:", long_options, &long_index)) != -1) {
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
                args.allowed_gaps = atoi(optarg);
                break;
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
        std::cerr << "[E::lift_bed_run] Argument -C/--chainmap_fname is required\n";
        print_lift_bed_help_msg();
        exit(1);
    }

    std::cerr << "Inputs:\n";
    std::cerr << " - BED: " << args.bed_fname << "\n";
    std::cerr << " - clft: " << args.chainmap_fname << "\n";

    lift_bed(args);
    
    double cpu_duration = (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration = (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "Finished in " << cpu_duration << " CPU seconds, or " << 
                                   wall_duration.count() << " wall clock seconds\n";

    return 0;
}
