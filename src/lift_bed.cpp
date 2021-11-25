#include <ctime>
#include <getopt.h>
#include <iostream>
#include "chain.hpp"
#include "leviosam_utils.hpp"
#include "lift_bed.hpp"
#include "version.hpp"


void print_lift_bed_help_msg() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Lift over a BED file\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Usage:   leviosam bed [options] -b <bed> -C <clft> -p <prefix>\n\n");
    fprintf(stderr, "Inputs:  -b string   Path to the input BED.\n");
    fprintf(stderr, "         -C string   Path to the chain index.\n");
    fprintf(stderr, "         -p string   Prefix to the output files.\n");
    fprintf(stderr, "Options: -h          Print detailed usage.\n");
    fprintf(stderr, "         -V INT      Verbose level [0].\n");
    fprintf(stderr, "\n");
}


int lift_bed(lift_bed_opts &args) {
    chain::ChainMap cmap = [&] {
        if (args.chainmap_fname != "") {
            std::cerr << "Loading levioSAM index...";
            std::ifstream in(args.chainmap_fname, std::ios::binary);
            return chain::ChainMap(
                in, args.verbose, args.allowed_cigar_changes
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
            std::vector<std::string> fields = LevioSamUtils::str_to_vector(line, "\t");
            std::cerr << line << "\n";
            size_t p1 = std::stoi(fields[1]);
            size_t p2 = std::stoi(fields[2]);
            size_t op1 = cmap.lift_pos(fields[0], p1);
            size_t op2 = cmap.lift_pos(fields[0], p2);
            out_f << fields[0] << "\t" << op1 << "\t" << op2;
            if (fields.size() > 3) {
                for (int i = 3; i < fields.size(); i++) {
                    out_f << "\t" << fields[i];
                }
            }
            out_f << "\n";
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
        {"allowed_cigar_changes", required_argument, 0, 'G'},
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
                args.allowed_cigar_changes = atoi(optarg);
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
