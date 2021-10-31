import pysam
import unittest
import subprocess

WG_MAJOR_LFT = 'testdata/wg-maj.lft'

BWA_SE_GOLD = 'testdata/bwa-single_end-cigar_op-grch38.bam'
BWA_SE_SAM = 'testdata/bwa-single_end-cigar_op-major.bam'
BWA_SE_OUT_PREFIX = 'testdata/bwa-single_end-cigar_op-major-lifted'

BT2_SE_GOLD = 'testdata/bt2-single_end-cigar_op-grch38.bam'
BT2_SE_SAM = 'testdata/bt2-single_end-cigar_op-major.bam'
BT2_SE_OUT_PREFIX = 'testdata/bt2-single_end-cigar_op-major-lifted'

BWA_PE_GOLD = 'testdata/bwa-paired_end-grch38.bam'
BWA_PE_SAM = 'testdata/bwa-paired_end-major.bam'
BWA_PE_OUT_PREFIX = 'testdata/bwa-paired_end-major-lifted'

BT2_PE_GOLD = 'testdata/bt2-paired_end-grch38.bam'
BT2_PE_SAM = 'testdata/bt2-paired_end-major.bam'
BT2_PE_OUT_PREFIX = 'testdata/bt2-paired_end-major-lifted'

OVRLP_LFT = 'testdata/overlapping_example.lft'
OVRLP_GOLD = 'testdata/overlapping_example-lifted-gold.bam'
OVRLP_SAM = 'testdata/overlapping_example.bam'
OVRLP_OUT_PREFIX = 'testdata/overlapping_example-lifted'

params = [
    {'name': 'BWA-SE', 'gold': BWA_SE_GOLD, 'sam': BWA_SE_SAM, 'out_prefix': BWA_SE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'SE'},
    {'name': 'BT2-PE', 'gold': BT2_SE_GOLD, 'sam': BT2_SE_SAM, 'out_prefix': BT2_SE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'SE'},
    {'name': 'BWA-PE', 'gold': BWA_PE_GOLD, 'sam': BWA_PE_SAM, 'out_prefix': BWA_PE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'PE'},
    {'name': 'BT2-PE', 'gold': BT2_PE_GOLD, 'sam': BT2_PE_SAM, 'out_prefix': BT2_PE_OUT_PREFIX, 'lft': WG_MAJOR_LFT, 'mode': 'PE'},
    {'name': 'OVRLP', 'gold': OVRLP_GOLD, 'sam': OVRLP_SAM, 'out_prefix': OVRLP_OUT_PREFIX, 'lft': OVRLP_LFT, 'mode': 'SE'}
]


class SamProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for param in params:
            process = subprocess.Popen(
                ['./leviosam', 
                'lift', '-l', param['lft'], '-a', param['sam'],
                '-p', param['out_prefix'], '-O', 'bam'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            print(f'Lifted {param["out_prefix"]}.bam')

    @classmethod
    def tearDownClass(cls):
        for param in params:
            cmd = f'rm {param["out_prefix"]}.bam; rm {param["out_prefix"]}-RG.bam'
            process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            stdout, stderr = process_rm.communicate()
            print(f'Cleaned up {param["out_prefix"]}-lifted.bam and {param["out_prefix"]}-RG.bam')

    """Read a SAM file as a dict."""
    def read_sam_file_as_dict(self, fn):
        d = {}
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
                print('[PASSED] Picard AddOrReplaceReadGroups & ValidateSamFile')


class Chain(unittest.TestCase):
    def compare_aln(self, gold, result):
        if gold.is_unmapped or result.is_unmapped:
            return
        self.assertEqual(gold.reference_start, result.reference_start)

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

    def test_paired_endh38_to_h37(self):
        INPUT_FN = 'testdata/bt2-paired_end-grch38.bam'
        GOLD_FN = 'testdata/bt2-paired_end-grch37.bam'
        OUTPUT_PREFIX = 'testdata/bt2-paired_end-grch38_to_grch37'
        OUTPUT_FN = OUTPUT_PREFIX + '.bam'
        process = subprocess.Popen(
            ['./leviosam', 
            'lift', '-C', 'testdata/hg38ToHg19.over.clft', '-a', INPUT_FN,
            '-p', OUTPUT_PREFIX, '-O', 'bam'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(f'Lifted {INPUT_FN}')

        reads_gold1, reads_gold2 = self.read_paired_end(GOLD_FN)
        reads_lift1, reads_lift2 = self.read_paired_end(OUTPUT_FN)

        for _, (name, v) in enumerate(reads_gold1.items()):
            result = reads_lift1.get(name)
            self.assertFalse(result is None)
            self.compare_aln(v, result)
        for _, (name, v) in enumerate(reads_gold2.items()):
            result = reads_lift2.get(name)
            self.assertFalse(result is None)
            self.compare_aln(v, result)

        cmd = f'rm {OUTPUT_FN}'
        process_rm = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = process_rm.communicate()
        print(f'Cleaned up {OUTPUT_FN}')


if __name__ == '__main__':
    unittest.main()
