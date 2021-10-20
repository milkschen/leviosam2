'''
Tests for levioSAM utilities

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import unittest
import subprocess
import mask_fasta_with_bed
import leviosam_utils

class TestMakeFastaWithBed(unittest.TestCase):
    def test_chr1_1(self):
        ref = mask_fasta_with_bed.mask_fasta_with_bed(bed='testdata/chr1_test1.bed', fasta='testdata/chr1_1_100000.fa')
        gold_ref = leviosam_utils.read_fasta('testdata/mask_fasta_with_bed-chr1_test1.fa.results')
        self.assertDictEqual(ref, gold_ref)

if __name__ == '__main__':
    unittest.main()
