import compare_sam
import pysam
import unittest

def create_sam_header(name: list, length: list) -> pysam.AlignmentHeader:
    assert len(name) == len(length)
    sq = ''
    for i, n in enumerate(name):
        sq += f'@SQ	SN:{n}	LN:{length[i]}\n'
    pg = f'@PG	ID:bowtie2	PN:bowtie2	VN:2.4.2	CL:"bowtie2-align-s --wrapper basic-0 -x chm13.draft_v1.0 -S chm13_1.0.sam -1 HG002.R1.fq -2 HG002.R2.fq"'
    return pysam.AlignmentHeader.from_text(
        sq + pg
    )

class TestCompareSam(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chm13_v1_0_header = create_sam_header(
            [f'chr{i}' for i in range(1,23)],
            [248387497, 242696747, 201106605, 193575430,
             182045437, 172126870, 160567423, 146259322,
             150617274, 134758122, 135127772, 133324781,
             114240146, 101219177, 100338308, 96330493,
             84277185, 80542536, 61707359, 66210247,
             45827691, 51353906])
        
    def test_calc_identity(self):
        baseline = pysam.AlignedSegment.fromstring(
            'read	161	chr8	8132831	42	151M	=	8133221	541	CCCTGGGATGAATGATATAGGTTCATTCATTGAACCTCATTCATTCAACCTCACTGTTGGCAGTGGGAAGTTCCTACTGGGGTAGCAAAACTGACAATAGGAGTCTAGAGCTGTCATTAATTTCTTTCTCCTCTGTATGGAGAAAGCTTGC	FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFF',
            self.chm13_v1_0_header)
        query = pysam.AlignedSegment.fromstring(
            'read	161	chr8	8132845	24	4M14I133M	=	8133221	3473114	GCAAGCTTTCTCCATACAGAGGAGAAAGAAATTAATGACAGCTCTAGACTCCTATTGTCAGTTTTGCTACCCCAGTAGGAACTTCCCACTGCCAACAGTGAGGTTGAATGAATGAGGTTCAATGAATGAACCTATATCATTCATCCCAGGG	FFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF',
            self.chm13_v1_0_header)

        summ = compare_sam.CompareSamSummary()
        idy = summ._calc_identity(query=query, baseline=baseline)
        self.assertEqual(idy, 133/151)

if __name__ == '__main__':
    unittest.main()
