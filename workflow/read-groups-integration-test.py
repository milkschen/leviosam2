"""Exercise single/multiple-RG SE/PE workflows with actual BAM operations.

Run in the Pixi integration environment. The aligner fixture honors BWA's -R;
BAM operations use pysam's samtools and paired selection uses core reconcile.
"""

import itertools
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
    def test_read_groups_survive_single_and_paired_end_workflows(self):
        for paired, multiple, fallback in itertools.product((False, True), repeat=3):
            if fallback and not multiple:
                continue
            groups = (
                [GROUP, dict(ID="second", SM="sample name", LB="second-library")]
                if multiple
                else [GROUP]
            )
            with (
                self.subTest(paired=paired, multiple=multiple, fallback=fallback),
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
                    "rg = sys.argv[sys.argv.index('-R') + 1].replace(chr(92) + 't', chr(9)) if '-R' in sys.argv else ''\n"
                    "fields = dict(f.split(':', 1) for f in rg.split('\t')[1:])\n"
                    "print('@SQ\tSN:chr1\tLN:1000')\n"
                    "if rg: print(rg)\n"
                    f"for line in pathlib.Path({str(target)!r}).read_text().splitlines():\n"
                    "    print(line + ('\tRG:Z:' + fields['ID'] if fields else ''))\n"
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
                if fallback:
                    argv.extend(
                        ["--read_group", r"@RG\tID:fallback\tLB:fallback-library"]
                    )
                with mock.patch.object(sys, "argv", argv):
                    workflow = leviosam2.Leviosam2Workflow(leviosam2.parse_args())
                workflow._set_filenames()
                original = records("source", paired, score=200) + records(
                    "target",
                    paired,
                    score=100,
                    group="second" if multiple else "original",
                )
                if paired and multiple:
                    # Each end can have a different original group. The selected
                    # target pair will deliberately be emitted in reversed order.
                    original[-1] = original[-1].replace("RG:Z:second", "RG:Z:original")
                if multiple:
                    original += records("ungrouped", paired, score=100, group=None)
                expected = {
                    (line.split("\t")[0], int(line.split("\t")[1]) & 192): (
                        line.split("RG:Z:")[1]
                        if "RG:Z:" in line
                        else "fallback" if fallback else None
                    )
                    for line in original
                }
                write_bam(workflow.input_bam, groups, original)
                workflow._configure_read_group()
                target_records = records("source", paired, score=100, group=None)
                target_pair = records("target", paired, score=200, group=None)
                if paired and multiple:
                    target_pair.reverse()
                target_records += target_pair
                if multiple:
                    target_records += records(
                        "ungrouped", paired, score=200, group=None
                    )
                    if not paired:
                        target_records += [
                            target_records[0].replace("\t0\tchr1", "\t256\tchr1"),
                            target_records[1].replace("\t0\tchr1", "\t2048\tchr1"),
                        ]
                target.write_text("\n".join(target_records) + "\n")
                if paired:
                    write_bam(workflow._fn_deferred_pe, groups, original)
                    write_bam(
                        workflow._fn_committed_pe, groups, records("committed", True)
                    )
                    workflow._fn_deferred_pe_fq1.touch()
                    workflow._fn_deferred_pe_fq2.touch()
                else:
                    write_bam(workflow._fn_deferred, groups, original)
                    write_bam(
                        workflow._fn_committed_sorted,
                        groups,
                        records("committed", False),
                    )
                    workflow._fn_deferred_fq_se.touch()
                workflow.run_realign_deferred()
                if paired:
                    workflow.run_refflow_merge_pe()
                workflow.run_merge_and_index()
                with pysam.AlignmentFile(workflow._fn_final) as result:
                    result_groups = {
                        group["ID"]: group for group in result.header.to_dict()["RG"]
                    }
                    for group in groups:
                        self.assertEqual(result_groups[group["ID"]], group)
                    if fallback:
                        self.assertEqual(
                            result_groups["fallback"]["LB"], "fallback-library"
                        )
                    reads = list(result)
                for read in reads:
                    group = read.get_tag("RG") if read.has_tag("RG") else None
                    if read.query_name == "committed":
                        self.assertEqual(group, "original")
                    else:
                        self.assertEqual(
                            group, expected[(read.query_name, read.flag & 192)]
                        )
                self.assertEqual(
                    len(reads),
                    (8 if paired else 6) if multiple else (6 if paired else 3),
                )
                if paired:
                    selected = {
                        r.query_name: r.get_tag("RF") for r in reads if r.has_tag("RF")
                    }
                    expected_selection = dict(source="source", target="target")
                    if multiple:
                        expected_selection["ungrouped"] = "target"
                    self.assertEqual(selected, expected_selection)


class ReadGroupRestorationErrors(unittest.TestCase):
    def test_invalid_source_assignments_and_unmatched_output_fail(self):
        for case, message in (
            ("ambiguous", "Ambiguous source RG"),
            ("undefined", "Invalid or undefined source RG"),
            ("unmatched", "No source record"),
            ("conflicting", "Conflicting read-group ID"),
        ):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                root = pathlib.Path(directory)
                groups = [GROUP, dict(ID="second", LB="second-library")]
                original = records("a", False)
                if case == "ambiguous":
                    original += records("a", False, group="second")
                if case == "undefined":
                    original = records("a", False, group="missing")
                write_bam(root / "original.bam", groups, original)
                target_groups = (
                    [dict(GROUP, LB="different")] if case == "conflicting" else []
                )
                write_bam(
                    root / "target.bam",
                    target_groups,
                    records("other" if case == "unmatched" else "a", False, group=None),
                )
                binary = os.environ.get("LEVIOSAM2", str(ROOT / "build" / "leviosam2"))
                result = subprocess.run(
                    [
                        binary,
                        "reconcile",
                        "--restore-rg",
                        str(root / "original.bam"),
                        "-s",
                        f"target:{root / 'target.bam'}",
                        "-o",
                        str(root / "out.bam"),
                    ],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(message, result.stderr)

    def test_paired_restoration_rejects_ambiguous_names_ends_and_undefined_rg(self):
        for case, message in (
            ("names", "Ambiguous source read name"),
            ("ends", "Ambiguous read end"),
            ("undefined", "Invalid or undefined source RG"),
        ):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                root = pathlib.Path(directory)
                original = records(
                    "a",
                    True,
                    score=100,
                    group="missing" if case == "undefined" else "original",
                )
                target = records("a", True, score=200, group=None)
                if case == "names":
                    original *= 2
                    target *= 2
                if case == "ends":
                    original[1] = original[1].replace("\t147\t", "\t99\t")
                write_bam(root / "original.bam", [GROUP], original)
                write_bam(root / "target.bam", [], target)
                binary = os.environ.get("LEVIOSAM2", str(ROOT / "build" / "leviosam2"))
                result = subprocess.run(
                    [
                        binary,
                        "reconcile",
                        "--preserve-rg",
                        "-m",
                        "-s",
                        f"source:{root / 'original.bam'}",
                        "-s",
                        f"target:{root / 'target.bam'}",
                        "-o",
                        str(root / "out.bam"),
                    ],
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn(message, result.stderr)

    def test_secondary_source_records_do_not_override_primary_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            groups = [GROUP, dict(ID="second", LB="second-library")]
            original = records("a", False) + [
                records("a", False, group="second")[0].replace(
                    "\t0\tchr1", "\t256\tchr1"
                )
            ]
            write_bam(root / "original.bam", groups, original)
            write_bam(root / "target.bam", [], records("a", False, group=None))
            binary = os.environ.get("LEVIOSAM2", str(ROOT / "build" / "leviosam2"))
            subprocess.run(
                [
                    binary,
                    "reconcile",
                    "--restore-rg",
                    str(root / "original.bam"),
                    "-s",
                    f"target:{root / 'target.bam'}",
                    "-o",
                    str(root / "out.bam"),
                ],
                check=True,
            )
            with pysam.AlignmentFile(root / "out.bam") as stream:
                self.assertEqual(next(stream).get_tag("RG"), "original")


if __name__ == "__main__":
    unittest.main()
