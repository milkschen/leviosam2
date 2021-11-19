import pysam
import unittest
import subprocess
import sys


LEVIOSAM = '../leviosam'
WG_MAJOR_LFT = 'major.lft'

BWA_SE_H38 = 'bwa-se-grch38.bam'
BWA_SE_MAJOR = 'bwa-se-major.bam'
BWA_SE_OUT_PREFIX = 'bwa-se-major-lifted'

BT2_SE_H38 = 'bt2-se-grch38.bam'
BT2_SE_MAJOR = 'bt2-se-major.bam'
BT2_SE_OUT_PREFIX = 'bt2-se-major-lifted'

BWA_PE_H38 = 'bwa-pe-grch38.bam'
BWA_PE_MAJOR = 'bwa-pe-major.bam'
BWA_PE_OUT_PREFIX = 'bwa-pe-major-lifted'

BT2_PE_H38 = 'bt2-pe-grch38.bam'
BT2_PE_MAJOR = 'bt2-pe-major.bam'
BT2_PE_OUT_PREFIX = 'bt2-pe-major-lifted'

OVRLP_LFT = 'overlapping_example.lft'
OVRLP_GOLD = 'overlapping_example-lifted-gold.bam'
OVRLP_SAM = 'overlapping_example.bam'
OVRLP_OUT_PREFIX = 'overlapping_example-lifted'

params = [
    {'name': 'BWA-SE', 'gold': BWA_SE_H38, 'sam': BWA_SE_MAJOR, 'out_prefix': BWA_SE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'SE'},
    {'name': 'BT2-PE', 'gold': BT2_SE_H38, 'sam': BT2_SE_MAJOR, 'out_prefix': BT2_SE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'SE'},
    {'name': 'BWA-PE', 'gold': BWA_PE_H38, 'sam': BWA_PE_MAJOR, 'out_prefix': BWA_PE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'PE'},
    {'name': 'BT2-PE', 'gold': BT2_PE_H38, 'sam': BT2_PE_MAJOR, 'out_prefix': BT2_PE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'PE'},
    {'name': 'OVRLP', 'gold': OVRLP_GOLD, 'sam': OVRLP_SAM, 'out_prefix': OVRLP_OUT_PREFIX, 'lft': OVRLP_LFT, 'mode': 'SE'}
]


class SamProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for param in params:
            process = subprocess.Popen(
                [LEVIOSAM,
                'lift', '-l', param['lft'], '-a', param['sam'],
                '-p', param['out_prefix'], '-O', 'bam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

    @classmethod
    def tearDownClass(cls):
        for param in params:
            cmd = f'rm {param["out_prefix"]}.bam; rm {param["out_prefix"]}-RG.bam'
            process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            stdout, stderr = process_rm.communicate()
            print(f'Cleaned up {param["out_prefix"]}.bam and {param["out_prefix"]}-RG.bam')

    """Read a SAM file as a dict."""
    def read_sam_file_as_dict(self, fn):
        d = {}
        print(fn, file=sys.stderr)
        f = pysam.AlignmentFile(fn)
        for rec in f:
            qname = rec.query_name
            if not d.get(qname):
                d[qname] = [None, None]
            if (not rec.is_paired) or rec.is_read1:
                d[qname][0] = rec
            else:
                d[qname][1] = rec
        return d

    """Test if the lifted fields are identical to the alignments against standard reference.

    This is not globally true, since some reads could be mapped to other regions in an augmented
    reference, but we only provide cases where mapping fields are not altered.
    """
    def test_identical_fields(self):
        def compare_fields(lift, gold):
            self.assertEqual(lift.reference_start, gold.reference_start)
            self.assertEqual(lift.reference_name, gold.reference_name)
            self.assertEqual(lift.cigarstring, gold.cigarstring)
            # self.assertEqual(lift.next_reference_start, gold.next_reference_start)
            self.assertEqual(lift.next_reference_name, gold.next_reference_name)
            self.assertEqual(lift.template_length, gold.template_length)
            self.assertEqual(lift.query_sequence, gold.query_sequence)
            self.assertEqual(lift.query_qualities, gold.query_qualities)

        for i, param in enumerate(params):
            assert param['mode'] in ['SE', 'PE']
            
            # Each subtest here is one parameter set (e.g. BWA-SE, BT2-PE)
            with self.subTest(aln_param_idx=param):
                dict_lifted = self.read_sam_file_as_dict(f'{param["out_prefix"]}.bam')
                dict_gold = self.read_sam_file_as_dict(f'{param["gold"]}')

                for k in dict_lifted.keys():
                    # Single-end
                    if param['mode'] == 'SE':
                        compare_fields(dict_lifted[k][0], dict_gold[k][0])
                    # Paired-end
                    elif param['mode'] == 'PE':
                        compare_fields(dict_lifted[k][0], dict_gold[k][0])
                        compare_fields(dict_lifted[k][1], dict_gold[k][1])


    ''' Test if the lifted SAM file passes `picard ValidateSamFile`. '''
    def test_picard_validate(self):
        for i, param in enumerate(params):
            with self.subTest(aln_param_idx=i):
                process = subprocess.Popen(
                    ['picard',
                        f'AddOrReplaceReadGroups I={param["out_prefix"]}.bam \
                        O={param["out_prefix"]}-RG.bam \
                        RGID=test RGLB=test RGPL=illumina RGSM=test RGPU=TEST'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                process2 = subprocess.Popen(
                    ['picard',
                        f'ValidateSamFile I={param["out_prefix"]}-RG.bam MODE= VERBOSE'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process2.communicate()
                self.assertNotIn(
                    'ERROR', stderr.decode(), msg='Check "No errors found"')
                print(f'[PASSED] Picard AddOrReplaceReadGroups & ValidateSamFile ({param["out_prefix"]}.bam)')


BT2_PE_CHM13 = 'bt2-pe-chm13_v1.1_hg2y.bam'
BT2_SE_CHM13 = 'bt2-se-chm13_v1.1_hg2y.bam'
BWA_PE_CHM13 = 'bwa-pe-chm13_v1.1_hg2y.bam'
BWA_SE_CHM13 = 'bwa-se-chm13_v1.1_hg2y.bam'
MM2_HIFI_1 = 'mm2-test1.sam'
MM2_HIFI_2 = 'mm2-test2.sam'
CLFT_CHM13_TO_H38 = 'chm13_v1.1_hg2Y-grch38.clft'
PE_CHAIN_TESTS = {BT2_PE_CHM13: BT2_PE_H38,
                  BWA_PE_CHM13: BWA_PE_H38}
SE_CHAIN_TESTS = {BT2_SE_CHM13: BT2_SE_H38,
                  BWA_SE_CHM13: BWA_SE_H38}

class Chain(unittest.TestCase):
    def compare_aln(self, gold, result):
        if gold.is_unmapped or result.is_unmapped or gold.is_supplementary or result.is_supplementary:
            return
        self.assertLess(abs(gold.reference_start-result.reference_start), 5)
        self.assertEqual(gold.reference_name, result.reference_name)
        self.assertEqual(gold.next_reference_name, result.next_reference_name)
        self.assertEqual(gold.query_sequence, result.query_sequence)
        self.assertEqual(gold.query_qualities, result.query_qualities)

    def read_single_end(self, fn):
        f = pysam.AlignmentFile(fn)
        reads1 = {}
        for r in f:
            reads1[r.query_name] = r
        return reads1

    def test_single_end_chain(self):
        for i, (orig, gold) in enumerate(SE_CHAIN_TESTS.items()):
            lifted_pre = orig + '-lifted'
            lifted_fn = lifted_pre + '.bam'
            process = subprocess.Popen(
                [LEVIOSAM,
                'lift', '-C', CLFT_CHM13_TO_H38, '-a', orig,
                '-p', lifted_pre, '-O', 'bam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            reads_gold1 = self.read_single_end(gold)
            reads_lift1 = self.read_single_end(lifted_fn)

            for _, (name, v) in enumerate(reads_gold1.items()):
                result = reads_lift1.get(name)
                self.assertFalse(result is None)
                self.compare_aln(v, result)

            cmd = f'rm {lifted_fn}'
            process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            stdout, stderr = process_rm.communicate()

    def read_paired_end(self, fn):
        f = pysam.AlignmentFile(fn)
        reads1 = {}
        reads2 = {}
        for r in f:
            if r.is_read1:
                reads1[r.query_name] = r
            elif r.is_read2:
                reads2[r.query_name] = r
            else:
                print(f'Invalid record in {fn}. Please check')
                print(r)
                exit(1)
        return reads1, reads2

    def test_paired_end_chain(self):
        for i, (orig, gold) in enumerate(PE_CHAIN_TESTS.items()):
            lifted_pre = orig + '-lifted'
            lifted_fn = lifted_pre + '.bam'
            process = subprocess.Popen(
                [LEVIOSAM,
                'lift', '-C', CLFT_CHM13_TO_H38, '-a', orig,
                '-p', lifted_pre, '-O', 'bam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            reads_gold1, reads_gold2 = self.read_paired_end(gold)
            reads_lift1, reads_lift2 = self.read_paired_end(lifted_fn)

            for _, (name, v) in enumerate(reads_gold1.items()):
                result = reads_lift1.get(name)
                self.assertFalse(result is None)
                self.compare_aln(v, result)
            for _, (name, v) in enumerate(reads_gold2.items()):
                result = reads_lift2.get(name)
                self.assertFalse(result is None)
                self.compare_aln(v, result)

            cmd = f'rm {lifted_fn}'
            process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            stdout, stderr = process_rm.communicate()

    def test_validate_long_reads(self):
        for fn in [MM2_HIFI_1, MM2_HIFI_2]:
            process = subprocess.Popen(
                [LEVIOSAM,
                'lift', '-C', CLFT_CHM13_TO_H38, '-a', fn,
                '-p', fn+'-lifted', '-O', 'bam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            process2 = subprocess.Popen(
                ['picard',
                    f'ValidateSamFile I={fn}-lifted.bam MODE= VERBOSE'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process2.communicate()
            self.assertNotIn(
                'ERROR', stderr.decode(), msg='Check "No errors found"')
            cmd = f'rm {fn}-lifted.bam'
            process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            stdout, stderr = process_rm.communicate()
            print(f'[PASSED] Picard ValidateSamFile ({fn})')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        LEVIOSAM = sys.argv.pop()
        print(f'Use the leviosam software at {LEVIOSAM}', file=sys.stderr)

    unittest.main()
