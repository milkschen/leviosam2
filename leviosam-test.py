import unittest
import subprocess

BWA_SE_GOLD = 'testdata/bwa-single_end-cigar_op-grch38.sam'
BWA_SE_SAM = 'testdata/bwa-single_end-cigar_op-major.sam'
BWA_SE_OUT_PREFIX = 'testdata/bwa-single_end-cigar_op-major-lifted'

BT2_SE_GOLD = 'testdata/bt2-single_end-cigar_op-grch38.sam'
BT2_SE_SAM = 'testdata/bt2-single_end-cigar_op-major.sam'
BT2_SE_OUT_PREFIX = 'testdata/bt2-single_end-cigar_op-major-lifted'

BWA_PE_GOLD = 'testdata/bwa-paired_end-grch38.sam'
BWA_PE_SAM = 'testdata/bwa-paired_end-major.sam'
BWA_PE_OUT_PREFIX = 'testdata/bwa-paired_end-major-lifted'

BT2_PE_GOLD = 'testdata/bt2-paired_end-grch38.sam'
BT2_PE_SAM = 'testdata/bt2-paired_end-major.sam'
BT2_PE_OUT_PREFIX = 'testdata/bt2-paired_end-major-lifted'

params = [
    {'name': 'BWA-SE', 'gold': BWA_SE_GOLD, 'sam': BWA_SE_SAM, 'out_prefix': BWA_SE_OUT_PREFIX, 'mode': 'SE'},
    {'name': 'BT2-PE', 'gold': BT2_SE_GOLD, 'sam': BT2_SE_SAM, 'out_prefix': BT2_SE_OUT_PREFIX, 'mode': 'SE'},
    {'name': 'BWA-PE', 'gold': BWA_PE_GOLD, 'sam': BWA_PE_SAM, 'out_prefix': BWA_PE_OUT_PREFIX, 'mode': 'PE'},
    {'name': 'BT2-PE','gold': BT2_PE_GOLD, 'sam': BT2_PE_SAM, 'out_prefix': BT2_PE_OUT_PREFIX, 'mode': 'PE'}
]


SAM_NAME = 0
SAM_FLAG = 1
SAM_CHROM = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_MCHROM = 6
SAM_MPOS = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10
dict_test_field = {
    'FLAG': SAM_FLAG, 'CHROM': SAM_CHROM, 'POS': SAM_POS,
    'CIGAR': SAM_CIGAR, 'MCHROM': SAM_MCHROM, 'MPOS': SAM_MPOS, 'TLEN': SAM_TLEN,
    'SEQ': SAM_SEQ, 'QUAL': SAM_QUAL}
    # 'MAPQ': SAM_MAPQ,


class SamProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for param in params:
            process = subprocess.Popen(
                ['./leviosam', 
                'lift', '-l', 'testdata/wg-maj.lft', '-a', param['sam'],
                '-p', param['out_prefix'] + '.sam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            print(f'Lifted {param["out_prefix"]}.sam')

    @classmethod
    def tearDownClass(cls):
        for param in params:
            cmd = f'rm {param["out_prefix"]}.sam; rm {param["out_prefix"]}-RG.sam'
            #process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            #stdout, stderr = process_rm.communicate()
            #print(f'Cleaned up {param["out_prefix"]}-lifted.sam and {param["out_prefix"]}-RG.sam')

    """Read a SAM file as a dict."""
    def read_sam_file_as_dict(self, fn):
        d = {}
        with open(fn, 'r') as f:
            for line in f:
                if line[0] != '@':
                    line = line.split()
                    d.setdefault(line[0], []).append(line)
        return d

    """Test if the lifted fields are identical to the alignments against standard reference.

    This is not globally true, since some reads could be mapped to other regions in an augmented
    reference, but we only provide cases where mapping fields are not altered.
    """
    def test_identical_fields(self):
        for i, param in enumerate(params):
            assert param['mode'] in ['SE', 'PE']
            
            # Each subtest here is one parameter set (e.g. BWA-SE, BT2-PE)
            with self.subTest(aln_param_idx=param):
                dict_lifted = self.read_sam_file_as_dict(f'{param["out_prefix"]}.sam')
                dict_gold = self.read_sam_file_as_dict(f'{param["gold"]}')

                for test_idx, test_field in enumerate(dict_test_field.keys()):
                    with self.subTest(test_idx=test_field):
                        for k in dict_lifted.keys():
                            # Single-end files.
                            if param['mode'] == 'SE':
                                aln_lifted = dict_lifted[k][0]
                                aln_gold = dict_gold[k][0]
                                self.assertEqual(
                                    aln_lifted[dict_test_field[test_field]],
                                    aln_gold[dict_test_field[test_field]])
                            # Paired-end files. Get first/second segments and then compare.
                            elif param['mode'] == 'PE':
                                if int(dict_lifted[k][0][SAM_FLAG]) & 64:
                                    aln_lifted_first = dict_lifted[k][0]
                                    aln_lifted_second = dict_lifted[k][1]
                                if int(dict_gold[k][0][SAM_FLAG]) & 64:
                                    aln_gold_first = dict_gold[k][0]
                                    aln_gold_second = dict_gold[k][1]
                                self.assertEqual(
                                    aln_lifted_first[dict_test_field[test_field]],
                                    aln_gold_first[dict_test_field[test_field]])
                                self.assertEqual(
                                    aln_lifted_second[dict_test_field[test_field]],
                                    aln_gold_second[dict_test_field[test_field]])
                        print(f'[PASSED] {len(dict_lifted.keys())} {test_field} for {param["name"]}')

    ''' Test if the lifted SAM file passes `picard ValidateSamFile`. '''
    def test_picard_validate(self):
        for i, param in enumerate(params):
            with self.subTest(aln_param_idx=i):
                process = subprocess.Popen(
                    ['picard',
                        f'AddOrReplaceReadGroups I={param["out_prefix"]}.sam \
                        O={param["out_prefix"]}-RG.sam \
                        RGID=test RGLB=test RGPL=illumina RGSM=test RGPU=TEST'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                process2 = subprocess.Popen(
                    ['picard',
                        f'ValidateSamFile I={param["out_prefix"]}-RG.sam MODE= VERBOSE'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process2.communicate()
                self.assertNotIn(
                    'ERROR', stderr.decode(), msg='Check "No errors found"')
                print('[PASSED] Picard AddOrReplaceReadGroups & ValidateSamFile')


if __name__ == '__main__':
    unittest.main()
