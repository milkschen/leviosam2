#include <liftover.hpp>
#include <getopt.h>
#include <htslib/vcf.h>

struct lift_opts {
    std::string vcf_fname = "";
    std::string sample = "";
    std::string output_fname = "liftover.lft";
};

int main(int argc, char** argv) {
    // read vcf file name
    // read other options (which ones do we want to include?)
    int c;
    lift_opts args;
    static struct option long_options[] {
        {"vcf", required_argument, 0, 'v'},
        {"sample", required_argument, 0, 's'},
        {"output", required_argument, 0, 'o'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "v:s:o:", long_options, &long_index)) != -1) { 
        switch (c) {
            case 0:
                break;
            case 'v':
                args.vcf_fname = optarg;
                break;
            case 's':
                args.sample = optarg;
                break;
            case 'o':
                args.output_fname = optarg;
                break;
            default:
                fprintf(stderr, "ignoring option %s\n", c);
                exit(1);
                break;
        }
    }

    if (args.vcf_fname == "" && args.sample == "") {
        fprintf(stderr, "vcf file name and sample name are required!! \n");
        exit(1);
    }
    vcfFile* fp = bcf_open(args.vcf_fname.data(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    lift::LiftMap l(fp, hdr, args.sample);
    std::ofstream o(args.output_fname);
    l.serialize(o);
    fprintf(stderr, "liftover file saved to %s\n", args.output_fname.data());
    return 0;
}
