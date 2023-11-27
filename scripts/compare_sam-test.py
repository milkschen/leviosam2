import unittest

import compare_sam
import pysam


def create_sam_header(name: list, length: list) -> pysam.AlignmentHeader:
    assert len(name) == len(length)
    sq = ""
    for i, n in enumerate(name):
        sq += f"@SQ	SN:{n}	LN:{length[i]}\n"
    pg = f'@PG	ID:bowtie2	PN:bowtie2	VN:2.4.2	CL:"bowtie2-align-s --wrapper basic-0 -x chm13.draft_v1.0 -S chm13_1.0.sam -1 HG002.R1.fq -2 HG002.R2.fq"'
    return pysam.AlignmentHeader.from_text(sq + pg)


class TestCompareSam(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chm13_v1_0_header = create_sam_header(
            [f"chr{i}" for i in range(1, 23)],
            [
                248387497,
                242696747,
                201106605,
                193575430,
                182045437,
                172126870,
                160567423,
                146259322,
                150617274,
                134758122,
                135127772,
                133324781,
                114240146,
                101219177,
                100338308,
                96330493,
                84277185,
                80542536,
                61707359,
                66210247,
                45827691,
                51353906,
            ],
        )

    def test_calc_identity(self):
        def _test_calc_identity(b, q, gold_idy):
            baseline = pysam.AlignedSegment.fromstring(
                b, self.chm13_v1_0_header
            )
            query = pysam.AlignedSegment.fromstring(q, self.chm13_v1_0_header)
            summ = compare_sam.CompareSamSummary()
            idy = summ._calc_identity(query=query, baseline=baseline)
            self.assertEqual(idy, gold_idy)

        test_cases = [
            [
                "A00744:46:HV3C3DSXX:2:1101:13946:8171	145	chr1	164180484	42	9M2I140M	=	164179905	-728	ATATGTATATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATTCAGAGCTCAGATATCTGTGACTGTTTTTTCTCTTTTCCATACTTTTCTTTGGGGCAGCTGTTATCTCTGTGGCCTTGAGCTTGATC	F:FFFFFFFFFFFFFFFFF,F:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFF,FFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF",
                "A00744:46:HV3C3DSXX:2:1101:13946:8171	145	chr1	164180472	42	62M10D89M	=	164179905	652819	ATATGTATATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATTCAGAGCTCAGATATCTGTGACTGTTTTTTCTCTTTTCCATACTTTTCTTTGGGGCAGCTGTTATCTCTGTGGCCTTGAGCTTGATC	F:FFFFFFFFFFFFFFFFF,F:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFF,FFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFF",
                89 / 151,
            ],
            [
                "A00744:46:HV3C3DSXX:2:1101:4589:4930	161	chr8	8132831	42	151M	=	8133221	541	CCCTGGGATGAATGATATAGGTTCATTCATTGAACCTCATTCATTCAACCTCACTGTTGGCAGTGGGAAGTTCCTACTGGGGTAGCAAAACTGACAATAGGAGTCTAGAGCTGTCATTAATTTCTTTCTCCTCTGTATGGAGAAAGCTTGC	FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFF",
                "A00744:46:HV3C3DSXX:2:1101:4589:4930	161	chr8	8132845	24	4M14I133M	=	8133221	3473114	GCAAGCTTTCTCCATACAGAGGAGAAAGAAATTAATGACAGCTCTAGACTCCTATTGTCAGTTTTGCTACCCCAGTAGGAACTTCCCACTGCCAACAGTGAGGTTGAATGAATGAGGTTCAATGAATGAACCTATATCATTCATCCCAGGG	FFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF",
                133 / 151,
            ],
        ]
        for t in test_cases:
            _test_calc_identity(t[0], t[1], t[2])


if __name__ == "__main__":
    unittest.main()
