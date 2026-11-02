"""Exercise the single-RG SE/PE workflow with actual BAM operations.

Run in the Pixi integration environment. The aligner fixture honors BWA's -R;
BAM operations use pysam's samtools and paired selection uses core reconcile.
"""

import os
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam

import leviosam2

ROOT = pathlib.Path(__file__).resolve().parents[1]
GROUP = dict(ID="original", SM="sample name", LB="library", PL="ILLUMINA")


def header(groups):
    return pysam.AlignmentHeader.from_dict(
        dict(HD=dict(VN="1.6", SO="unsorted"), SQ=[dict(SN="chr1", LN=1000)], RG=groups)
    )


def records(name, paired, score=100, group="original"):
    tags = f"\tAS:i:{score}\tNM:i:{max(0, (300 - score) // 100)}" + (
        f"\tRG:Z:{group}" if group else ""
    )
    flags = (99, 147) if paired else (0,)
    return [
        f"{name}\t{flag}\tchr1\t{11 + i * 30}\t60\t10M\t"
        + (
            f"=\t{41 if i == 0 else 11}\t{40 if i == 0 else -40}"
            if paired
            else "*\t0\t0"
        )
        + "\tAAAAAAAAAA\tIIIIIIIIII"
        + tags
        for i, flag in enumerate(flags)
    ]


def write_bam(path, groups, lines):
    hdr = header(groups)
    with pysam.AlignmentFile(path, "wb", header=hdr) as stream:
        for line in lines:
            stream.write(pysam.AlignedSegment.fromstring(line, hdr))


class ReadGroupWorkflowIntegration(unittest.TestCase):
    def test_single_group_survives_single_and_paired_end_workflows(self):
        for paired in (False, True):
            with (
                self.subTest(paired=paired),
                tempfile.TemporaryDirectory() as directory,
            ):
                root = pathlib.Path(directory)
                samtools = root / "samtools.py"
                samtools.write_text(
                    f"#!{sys.executable}\nimport pysam, sys\n"
                    "getattr(pysam, sys.argv[1])(*sys.argv[2:], catch_stdout=False)\n"
                )
                target = root / "aligned.sam"
                # The fixture reproduces the aligner's RG assignment from -R.
                aligner = root / "aligner.py"
                aligner.write_text(
                    f"#!{sys.executable}\nimport pathlib, sys\n"
                    "rg = sys.argv[sys.argv.index('-R') + 1].replace(chr(92) + 't', chr(9))\n"
                    "fields = dict(f.split(':', 1) for f in rg.split('\t')[1:])\n"
                    "print('@SQ\tSN:chr1\tLN:1000')\nprint(rg)\n"
                    f"for line in pathlib.Path({str(target)!r}).read_text().splitlines():\n"
                    "    print(line + '\tRG:Z:' + fields['ID'])\n"
                )
                for tool in (samtools, aligner):
                    tool.chmod(0o755)
                argv = [
                    "leviosam2.py",
                    "-s",
                    "ilmn_pe" if paired else "ilmn_se",
                    "-a",
                    "bwamem",
                    "-i",
                    str(root / "input.bam"),
                    "-o",
                    str(root / "result"),
                    "-C",
                    "unused.clft",
                    "-f",
                    str(root / "target.fa"),
                    "--keep_tmp_files",
                    "--aligner_exe",
                    str(aligner),
                    "--samtools_exe",
                    str(samtools),
                    "--leviosam2_exe",
                    os.environ.get("LEVIOSAM2", str(ROOT / "build" / "leviosam2")),
                ]
                with mock.patch.object(sys, "argv", argv):
                    workflow = leviosam2.Leviosam2Workflow(leviosam2.parse_args())
                workflow._set_filenames()
                original = records("source", paired, score=200) + records(
                    "target", paired, score=100
                )
                write_bam(workflow.input_bam, [GROUP], original)
                workflow._configure_read_group()
                target.write_text(
                    "\n".join(
                        records("source", paired, score=100, group=None)
                        + records("target", paired, score=200, group=None)
                    )
                    + "\n"
                )
                if paired:
                    write_bam(workflow._fn_deferred_pe, [GROUP], original)
                    write_bam(
                        workflow._fn_committed_pe, [GROUP], records("committed", True)
                    )
                    workflow._fn_deferred_pe_fq1.touch()
                    workflow._fn_deferred_pe_fq2.touch()
                else:
                    write_bam(
                        workflow._fn_committed_sorted,
                        [GROUP],
                        records("committed", False),
                    )
                    workflow._fn_deferred_fq_se.touch()
                workflow.run_realign_deferred()
                if paired:
                    workflow.run_refflow_merge_pe()
                workflow.run_merge_and_index()
                with pysam.AlignmentFile(workflow._fn_final) as result:
                    self.assertEqual(result.header.to_dict()["RG"], [GROUP])
                    reads = list(result)
                self.assertEqual(len(reads), 6 if paired else 3)
                self.assertTrue(
                    all(read.get_tag("RG") == GROUP["ID"] for read in reads)
                )
                if paired:
                    selected = {
                        r.query_name: r.get_tag("RF") for r in reads if r.has_tag("RF")
                    }
                    self.assertEqual(selected, dict(source="source", target="target"))


if __name__ == "__main__":
    unittest.main()
