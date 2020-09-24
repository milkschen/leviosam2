import unittest
import subprocess
from parameterized import parameterized, parameterized_class

BWA_OP_GOLD = 'testdata/bwa-single_end-cigar_op-grch38.sam'
BWA_OP_SAM = 'testdata/bwa-single_end-cigar_op-major.sam'
BWA_OP_OUT_PREFIX = 'testdata/bwa-single_end-cigar_op-major-lifted'

BT2_OP_GOLD = 'testdata/bt2-single_end-cigar_op-grch38.sam'
BT2_OP_SAM = 'testdata/bt2-single_end-cigar_op-major.sam'
BT2_OP_OUT_PREFIX = 'testdata/bt2-single_end-cigar_op-major-lifted'

params = [
    {'gold': BWA_OP_GOLD, 'sam': BWA_OP_SAM, 'out_prefix': BWA_OP_OUT_PREFIX},
    {'gold': BT2_OP_GOLD, 'sam': BT2_OP_SAM, 'out_prefix': BT2_OP_OUT_PREFIX}]

class TestBwaLocal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for param in params:
            process = subprocess.Popen(
                ['./leviosam', 
                'lift', '-l', 'testdata/wg-maj.lft', '-a', param['sam'], '-p', param['out_prefix']],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            print(f'Lifted {param["out_prefix"]}.sam')

    '''
    Test if the lifted POSs are identical to the alignments against standard reference.

    This is not globally true, since some reads could be mapped to other regions in an augmented
    reference, but we only provide cases where mapping positions are not altered in the test case.
    '''
    def test_identical_pos(self):
        for i, param in enumerate(params):
            with self.subTest(i=i):
                dict_name_pos_lifted = {}
                with open(f'{param["out_prefix"]}.sam', 'r') as f:
                    for line in f:
                        if line[0] != '@':
                            line = line.split()
                            dict_name_pos_lifted[line[0]] = f'{line[2]}_{line[3]}'
                dict_name_pos_gold = {}
                with open(f'{param["gold"]}', 'r') as f:
                    for line in f:
                        if line[0] != '@':
                            line = line.split()
                            dict_name_pos_gold[line[0]] = f'{line[2]}_{line[3]}'
                for k in dict_name_pos_lifted.keys():
                    self.assertEqual(dict_name_pos_lifted[k], dict_name_pos_gold[k])
                print(f'[PASSED] {len(dict_name_pos_lifted.keys())} contig-position pairs match')

    '''
    Test if the lifted CIGARs are identical to the alignments against standard reference.

    This is not globally true, since some reads could be mapped to other regions in an augmented
    reference, but we only provide cases where CIGARs are not altered in the test case.
    '''
    def test_identical_cigar(self):
        for i, param in enumerate(params):
            with self.subTest(i=i):
                dict_cigar_lifted = {}
                with open(f'{param["out_prefix"]}.sam', 'r') as f:
                    for line in f:
                        if line[0] != '@':
                            line = line.split()
                            dict_cigar_lifted[line[0]] = f'{line[5]}'
                dict_cigar_gold = {}
                with open(f'{param["gold"]}', 'r') as f:
                    for line in f:
                        if line[0] != '@':
                            line = line.split()
                            dict_cigar_gold[line[0]] = f'{line[5]}'
                for k in dict_cigar_lifted.keys():
                    self.assertEqual(dict_cigar_lifted[k], dict_cigar_gold[k])
                print(f'[PASSED] {len(dict_cigar_lifted.keys())} CIGAR pairs match')

    ''' Test if the lifted SAM file passes `picard ValidateSamFile`. '''
    def test_picard_validate(self):
        for i, param in enumerate(params):
            with self.subTest(i=i):
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
