'''
This script includes two functionalities: chain2vcf and chain2bed

chain2vcf: Converts a chain file to a VCF file.

    This conversion is not capble of handling translocations and inversions. Users should specify
    a list of "chain ids" that contain only indels.

Example: 
    # This output VCF is used to convert h38 coordinates to T2T coordinates.
    python chain_utils.py chain2vcf -i hg38.t2t-chm13-v1.0.over.chain -c 1-23 -o main_chain-h38.t2t_chm13.vcf


chain2bed: Converts a chain file to a BED file.

    The conversion only supports one chain per contig (we currently don't support translocations).
    Users need to provide a list of chain ids for conversion.
    The output BED file will contain the "chain block" regions.
'''

import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    #subparsers = parser.add_subparsers(dest='subparser')
    parser_c2v = subparsers.add_parser('chain2vcf')
    parser_c2v.add_argument(
        '-i', '--in_chain', required=True,
        help='Path to the chain file to be converted.'
    )
    parser_c2v.add_argument(
        '-c', '--chain_ids', required=True,
        help='Chain ids to be converted, separated by commas and hyphens. E.g., "1-9", "1,5-8".'
    )
    parser_c2v.add_argument(
        '-o', '--out_vcf',
        help='Path to the output VCF file.'
    )
    parser_c2v.set_defaults(func=chain2vcf)

    parser_c2b = subparsers.add_parser('chain2bed')
    parser_c2b.add_argument(
        '-i', '--in_chain', required=True,
        help='Path to the chain file to be converted.'
    )
    parser_c2b.add_argument(
        '-c', '--chain_ids', required=True,
        help='Chain ids to be converted, separated by commas and hyphens. E.g., "1-9", "1,5-8".'
    )
    parser_c2b.add_argument(
        '-o', '--out_bed',
        help='Path to the output BED file.'
    )
    parser_c2b.set_defaults(func=chain2bed)

    args = parser.parse_args()
    args.func(args)

    #kwargs = vars(parser.parse_args())
    #globals()[kwargs.pop('subparser')](**kwargs)


class Chain_Fields():
    HEADER = 0
    SCORE = 1
    TNAME = 2
    TSIZE = 3
    TSTRAND = 4
    TSTART = 5
    TEND = 6
    QNAME = 7
    QSIZE = 8
    QSTRAND = 9
    QSTART = 10
    QEND = 11
    CID = 12
    ALN_SIZE = 0
    ALN_DT = 1
    ALN_DQ = 2


class Chain2Bed():
    def __init__(self, out_bed=sys.stdout, chain_hdr='', CF=None):
        assert(chain_hdr[CF.HEADER] == 'chain')
        assert(chain_hdr[CF.TNAME] == chain_hdr[CF.QNAME])
        self.out_bed = out_bed
        self.CF = CF
        self.contig = chain_hdr[self.CF.QNAME]
        self.qstart = int(chain_hdr[self.CF.QSTART])
        self.tstart = int(chain_hdr[self.CF.TSTART])
        self.bed_records = []

    def add(self, chain_aln):
        if len(self.bed_records) == 0 or self.tstart != self.bed_records[-1][1]:
            self.bed_records.append([self.tstart, self.tstart + int(chain_aln[self.CF.ALN_SIZE])])
        else:
            self.bed_records[-1][1] = self.tstart + int(chain_aln[self.CF.ALN_SIZE])

        self.tstart += int(chain_aln[self.CF.ALN_SIZE]) + int(chain_aln[self.CF.ALN_DT])
        self.qstart += int(chain_aln[self.CF.ALN_SIZE]) + int(chain_aln[self.CF.ALN_DQ])

    def write(self):
        for record in self.bed_records:
            print(f'{self.contig}\t{record[0]}\t{record[1]}', file=self.out_bed)



class Chain2Vcf():
    def __init__(self, out_vcf=sys.stdout, chain_hdr='', CF=None):
        assert(chain_hdr[CF.HEADER] == 'chain')
        # assert(chain_hdr[CF.TNAME] == chain_hdr[CF.QNAME])

        self.out_vcf = out_vcf
        self.CF = CF

        self.contig = chain_hdr[self.CF.QNAME]
        # Add one b/c VCF format is one-based.
        self.qstart = int(chain_hdr[self.CF.QSTART]) + 1
        self.tstart = int(chain_hdr[self.CF.TSTART]) + 1

        if chain_hdr[self.CF.TNAME] == chain_hdr[self.CF.QNAME]:
            # INS w.r.t the query sequence.
            if self.qstart < self.tstart:
                ref = 'A'
                qry = 'A' * (self.tstart - self.qstart + 1)
                print(f'{self.contig}\t{self.qstart}\t.\t{ref}\t{qry}\t.\tAUTO\tTPOS={self.tstart}',
                      file=self.out_vcf)
            # DEL
            elif self.qstart > self.tstart:
                ref = 'A' * (self.qstart - self.tstart + 1)
                qry = 'A'
                print(f'{self.contig}\t{self.qstart}\t.\t{ref}\t{qry}\t.\tAUTO\tTPOS={self.tstart}',
                      file=self.out_vcf)
        else:
            len_t = int(chain_hdr[self.CF.TEND]) - int(chain_hdr[self.CF.TSTART])
            len_q = int(chain_hdr[self.CF.QEND]) - int(chain_hdr[self.CF.QSTART])
            if len_t > len_q:
                ref = 'A'
                qry = 'A' * (len_t - len_q + 1)
                print(f'{self.contig}\t{self.qstart}\t.\t{ref}\t{qry}\t.\tAUTO\tTPOS={self.tstart}',
                      file=self.out_vcf)
            elif len_q > len_t:
                ref = 'A' * (len_q - len_t + 1)
                qry = 'A'
                print(f'{self.contig}\t{self.qstart - len_q}\t.\t{ref}\t{qry}\t.\tAUTO\tTPOS={self.tstart}',
                      file=self.out_vcf)




    def add(self, chain_aln):
        self.tstart += int(chain_aln[self.CF.ALN_SIZE])
        self.qstart += int(chain_aln[self.CF.ALN_SIZE])

        len_t = int(chain_aln[self.CF.ALN_DT])
        len_q = int(chain_aln[self.CF.ALN_DQ])

        if (len_t > len_q):
            # Insertion wrt the query seq
            len_t -= len_q
            len_q = 0
            seq_t = 'A' * (len_t + 1)
            seq_q = 'A' * (len_q + 1)
            print(f'{self.contig}\t{self.qstart}\t.\t{seq_q}\t{seq_t}\t.\tAUTO\t' +
                  f'TPOS={self.tstart};ALN_DT={int(chain_aln[self.CF.ALN_DT])};ALN_DQ={int(chain_aln[self.CF.ALN_DQ])}',
                  file=self.out_vcf)
        elif (len_t < len_q):
            # Deletion wrt the query seq
            len_q -= len_t
            len_t = 0
            seq_t = 'A' * (len_t + 1)
            seq_q = 'A' * (len_q + 1)
            print(f'{self.contig}\t{self.qstart}\t.\t{seq_q}\t{seq_t}\t.\tAUTO\t' +
                  f'TPOS={self.tstart};ALN_DT={int(chain_aln[self.CF.ALN_DT])};ALN_DQ={int(chain_aln[self.CF.ALN_DQ])}',
                  file=self.out_vcf)
        self.tstart += int(chain_aln[self.CF.ALN_DT])
        self.qstart += int(chain_aln[self.CF.ALN_DQ])



# Parse raw `chain_ids` and return a list of chain ids.
#
# Input:
#   chain_ids: A string of chain_ids that are separated by commas and hyphens.
# Output:
#   A list containing individual chain ids.
def parse_chain_id(chain_ids):
    list_chain = []
    chain_ids = chain_ids.split(',')
    for c in chain_ids:
        c = c.split('-')
        if len(c) == 1:
            list_chain.append(c[0])
        elif len(c) == 2:
            for i in range(int(c[0]), int(c[1]) + 1):
                list_chain.append(str(i))
    return list_chain


def get_contig_name(fn, chain_ids, CF):
    dict_cid_contig = {}
    with open(fn, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) == 13:
                if line[CF.CID] in chain_ids:
                    dict_cid_contig[line[CF.CID]] = line[CF.TNAME]
    return dict_cid_contig


def write_vcf_hdr(f_out, chain_ids, dict_cid_contig):
    print('##fileformat=VCFv4.3', file=f_out)
    print('##FILTER=<ID=AUTO,Description="Generated automatically.">', file=f_out)
    print(dict_cid_contig)
    for cid, contig in dict_cid_contig.items():
        print(f'##contig=<ID={contig}>', file=f_out)
    print('##INFO=<ID=TPOS,Number=A,Type=Integer,Description="Variant position on SEQ_T.">',
          file=f_out)
    print('##INFO=<ID=ALN_DT,Number=A,Type=Integer,Description="Length of gap on SEQ_T.">',
          file=f_out)
    print('##INFO=<ID=ALN_DQ,Number=A,Type=Integer,Description="Length of gap on SEQ_Q.">',
          file=f_out)
    print('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO', file=f_out)


def chain2bed(args):
    print('chain2bed', file=sys.stderr)
    # Constants.
    CF = Chain_Fields()

    chain_ids = parse_chain_id(args.chain_ids)
    print('Chain IDs to be processed', chain_ids, file=sys.stderr)

    if not args.out_bed:
        f_out = sys.stdout
    else:
        f_out = open(args.out_bed, 'w')

    f = open(args.in_chain, 'r')
    current_id = ''
    chain = None
    for line in f:
        line = line.split()
        # Chain header
        if len(line) == 13:
            if line[CF.CID] in chain_ids:
                chain = Chain2Bed(out_bed=f_out, chain_hdr=line, CF=CF)
                if line[CF.TNAME] != line[CF.QNAME]:
                    del chain
                    chain = None
        elif len(line) == 0:
            current_id = ''
            if chain:
                del chain
                chain = None
        elif len(line) == 3:
            if chain:
                chain.add(line)
        elif len(line) == 1:
            if chain:
                chain.write()
                del chain
                chain = None



def chain2vcf(args):
    # Constants.
    CF = Chain_Fields()

    chain_ids = parse_chain_id(args.chain_ids)
    if not args.out_vcf:
        f_out = sys.stdout
    else:
        f_out = open(args.out_vcf, 'w')

    write_vcf_hdr(f_out, chain_ids,
                  get_contig_name(args.in_chain, chain_ids, CF))

    f = open(args.in_chain, 'r')
    current_id = ''
    chain = None
    for line in f:
        line = line.split()
        # Chain header
        if len(line) == 13:
            if line[CF.CID] in chain_ids:
            # if line[CF.QNAME] in ['chr2']:
                chain = Chain2Vcf(out_vcf=f_out, chain_hdr=line, CF=CF)
                if line[CF.TNAME] != line[CF.QNAME]:
                    print('[Error] Contig names mismatch',
                          line[CF.TNAME], line[CF.QNAME], file=sys.stderr)
                    exit()
                    del chain
                    chain = None
        elif len(line) == 0:
            current_id = ''
            if chain:
                del chain
                chain = None
        elif len(line) == 3:
            if chain:
                chain.add(line)
        elif len(line) == 1:
            if chain:
                del chain
                chain = None


if __name__ == '__main__':
    parse_args()

