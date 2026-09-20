#!/usr/bin/env python3

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam


REPO_ROOT = Path(__file__).resolve().parent
TESTDATA = REPO_ROOT / "testdata"
LEVIOSAM = REPO_ROOT / "build" / "leviosam2"

IDENTITY_FLAG_MASK = 0x1 | 0x40 | 0x80 | 0x100 | 0x200 | 0x400 | 0x800
INVALIDATED_TAGS = {"LO", "MC", "MD", "NM"}


@dataclass(frozen=True)
class LiftCase:
    name: str
    source: str
    gold: str
    index: str


LEGACY_CASES = (
    LiftCase("bwa-se", "bwa-se-major.bam", "bwa-se-grch38.bam", "major.lft"),
    LiftCase("bt2-se", "bt2-se-major.bam", "bt2-se-grch38.bam", "major.lft"),
    LiftCase("bwa-pe", "bwa-pe-major.bam", "bwa-pe-grch38.bam", "major.lft"),
    LiftCase("bt2-pe", "bt2-pe-major.bam", "bt2-pe-grch38.bam", "major.lft"),
    LiftCase(
        "overlap",
        "overlapping_example.bam",
        "overlapping_example-lifted-gold.bam",
        "overlapping_example.lft",
    ),
)

CHAIN_INDEX = "chm13_v1.1_hg2Y-grch38.clft"
CHAIN_CASES = (
    LiftCase(
        "bt2-se-chain",
        "bt2-se-chm13_v1.1_hg2y.bam",
        "bt2-se-grch38.bam",
        CHAIN_INDEX,
    ),
    LiftCase(
        "bwa-se-chain",
        "bwa-se-chm13_v1.1_hg2y.bam",
        "bwa-se-grch38.bam",
        CHAIN_INDEX,
    ),
    LiftCase(
        "bt2-pe-chain",
        "bt2-pe-chm13_v1.1_hg2y.bam",
        "bt2-pe-grch38.bam",
        CHAIN_INDEX,
    ),
    LiftCase(
        "bwa-pe-chain",
        "bwa-pe-chm13_v1.1_hg2y.bam",
        "bwa-pe-grch38.bam",
        CHAIN_INDEX,
    ),
)

LONG_READ_INPUTS = ("mm2-test1.sam", "mm2-test2.sam")


def run_checked(command):
    command = [str(arg) for arg in command]
    try:
        return subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        stdout = getattr(exc, "stdout", "") or ""
        stderr = getattr(exc, "stderr", "") or ""
        rendered = " ".join(command)
        raise AssertionError(
            f"command failed: {rendered}\nstdout:\n{stdout}\nstderr:\n{stderr}"
        ) from exc


def read_records(path):
    with pysam.AlignmentFile(str(path)) as alignment_file:
        return list(alignment_file)


def record_key(record):
    if record.is_read1:
        read_end = 1
    elif record.is_read2:
        read_end = 2
    else:
        read_end = 0
    return (
        record.query_name,
        read_end,
        record.is_secondary,
        record.is_supplementary,
    )


def record_counts(records):
    return Counter(record_key(record) for record in records)


def unique_records(records):
    counts = record_counts(records)
    duplicates = {key: count for key, count in counts.items() if count != 1}
    if duplicates:
        raise AssertionError(f"test fixture has ambiguous record identities: {duplicates}")
    return {record_key(record): record for record in records}


def stable_tags(record):
    return {
        tag: value
        for tag, value in record.get_tags()
        if tag not in INVALIDATED_TAGS
    }


def run_lift(case, output_directory, index_option):
    output_prefix = output_directory / case.name
    run_checked(
        [
            LEVIOSAM,
            "lift",
            index_option,
            TESTDATA / case.index,
            "-a",
            TESTDATA / case.source,
            "-p",
            output_prefix,
            "-O",
            "bam",
        ]
    )
    output = output_prefix.with_suffix(".bam")
    if not output.is_file():
        raise AssertionError(f"lift did not create expected output: {output}")
    return output


class IntegrationTestCase(unittest.TestCase):
    def assert_unsorted_header(self, bam_path):
        with pysam.AlignmentFile(str(bam_path)) as alignment_file:
            header = alignment_file.header.to_dict()
        self.assertEqual(header.get("HD", {}).get("SO"), "unsorted")

    def assert_same_membership(self, expected, actual, context):
        self.assertEqual(
            record_counts(actual),
            record_counts(expected),
            msg=f"record membership differs for {context}",
        )

    def assert_pair_consistency(self, records, context):
        primary = {
            (record.query_name, record.is_read1): record
            for record in records
            if record.is_paired
            and not record.is_secondary
            and not record.is_supplementary
        }
        names = {query_name for query_name, _ in primary}
        for query_name in names:
            with self.subTest(context=context, query_name=query_name):
                read1 = primary.get((query_name, True))
                read2 = primary.get((query_name, False))
                self.assertIsNotNone(read1, "paired record is missing read1")
                self.assertIsNotNone(read2, "paired record is missing read2")
                self.assertEqual(read1.mate_is_unmapped, read2.is_unmapped)
                self.assertEqual(read2.mate_is_unmapped, read1.is_unmapped)
                self.assertEqual(read1.mate_is_reverse, read2.is_reverse)
                self.assertEqual(read2.mate_is_reverse, read1.is_reverse)
                self.assertEqual(read1.next_reference_name, read2.reference_name)
                self.assertEqual(read2.next_reference_name, read1.reference_name)
                self.assertEqual(read1.next_reference_start, read2.reference_start)
                self.assertEqual(read2.next_reference_start, read1.reference_start)
                self.assertEqual(read1.template_length, -read2.template_length)

    def validate_with_picard(self, bam_path, output_directory):
        read_group_bam = output_directory / f"{bam_path.stem}-rg.bam"
        run_checked(
            [
                "picard",
                "AddOrReplaceReadGroups",
                f"I={bam_path}",
                f"O={read_group_bam}",
                "RGID=test",
                "RGLB=test",
                "RGPL=illumina",
                "RGSM=test",
                "RGPU=test",
                "VALIDATION_STRINGENCY=LENIENT",
            ]
        )
        result = run_checked(
            [
                "picard",
                "ValidateSamFile",
                f"I={read_group_bam}",
                "MODE=SUMMARY",
                "IGNORE=MISSING_TAG_NM",
                # LevioSAM2 intentionally retains secondary/supplementary
                # provenance when an alignment becomes unmapped.
                "IGNORE=INVALID_FLAG_NOT_PRIM_ALIGNMENT",
                "IGNORE=INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
                # Secondary/supplementary long-read records may legally omit
                # base qualities in the source fixture.
                "IGNORE=QUALITY_NOT_STORED",
            ]
        )
        self.assertNotIn("ERROR::", result.stdout + result.stderr)


class LegacyLiftIntegrationTest(IntegrationTestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory(
            prefix="leviosam2-legacy-integration-"
        )
        cls.output_directory = Path(cls.temporary_directory.name)
        cls.outputs = {
            case.name: run_lift(case, cls.output_directory, "-l")
            for case in LEGACY_CASES
        }

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_alignment_fields(self):
        for case in LEGACY_CASES:
            with self.subTest(case=case.name):
                source_records = read_records(TESTDATA / case.source)
                lifted_records = read_records(self.outputs[case.name])
                gold_records = read_records(TESTDATA / case.gold)
                self.assert_unsorted_header(self.outputs[case.name])

                self.assert_same_membership(source_records, lifted_records, case.name)
                self.assert_same_membership(gold_records, lifted_records, case.name)

                source = unique_records(source_records)
                lifted = unique_records(lifted_records)
                gold = unique_records(gold_records)
                for key, result in lifted.items():
                    original = source[key]
                    expected = gold[key]
                    self.assertEqual(result.flag, expected.flag)
                    self.assertEqual(result.mapping_quality, original.mapping_quality)
                    self.assertEqual(result.reference_name, expected.reference_name)
                    self.assertEqual(result.reference_start, expected.reference_start)
                    self.assertEqual(result.cigarstring, expected.cigarstring)
                    self.assertEqual(
                        result.next_reference_name, expected.next_reference_name
                    )
                    self.assertEqual(
                        result.next_reference_start, expected.next_reference_start
                    )
                    self.assertEqual(result.template_length, expected.template_length)
                    self.assertEqual(result.query_sequence, expected.query_sequence)
                    self.assertEqual(result.query_qualities, expected.query_qualities)
                    self.assertEqual(stable_tags(result), stable_tags(original))
                    self.assertFalse(result.has_tag("NM"))
                    self.assertFalse(result.has_tag("MD"))

                self.assert_pair_consistency(lifted_records, case.name)

    def test_picard_validation(self):
        for case in LEGACY_CASES:
            with self.subTest(case=case.name):
                self.validate_with_picard(
                    self.outputs[case.name], self.output_directory
                )


class ChainLiftIntegrationTest(IntegrationTestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory(
            prefix="leviosam2-chain-integration-"
        )
        cls.output_directory = Path(cls.temporary_directory.name)
        cls.outputs = {
            case.name: run_lift(case, cls.output_directory, "-C")
            for case in CHAIN_CASES
        }

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_record_preservation_and_lifted_fields(self):
        for case in CHAIN_CASES:
            with self.subTest(case=case.name):
                source_records = read_records(TESTDATA / case.source)
                lifted_records = read_records(self.outputs[case.name])
                gold_records = read_records(TESTDATA / case.gold)
                self.assert_unsorted_header(self.outputs[case.name])

                self.assert_same_membership(source_records, lifted_records, case.name)
                source = unique_records(source_records)
                lifted = unique_records(lifted_records)
                gold = unique_records(gold_records)

                for key, result in lifted.items():
                    original = source[key]
                    self.assertEqual(
                        result.flag & IDENTITY_FLAG_MASK,
                        original.flag & IDENTITY_FLAG_MASK,
                    )
                    self.assertEqual(stable_tags(result), stable_tags(original))
                    self.assertFalse(result.has_tag("NM"))
                    self.assertFalse(result.has_tag("MD"))
                    self.assertTrue(result.has_tag("LO"))

                    lift_status = result.get_tag("LO").split("_")
                    self.assertEqual(lift_status[0] == "L", not result.is_unmapped)
                    if result.is_paired:
                        self.assertEqual(
                            lift_status[1] == "L", not result.mate_is_unmapped
                        )

                    expected = gold.get(key)
                    if result.is_unmapped:
                        self.assertEqual(
                            result.query_sequence, original.get_forward_sequence()
                        )
                        self.assertEqual(
                            result.query_qualities, original.get_forward_qualities()
                        )
                        continue
                    self.assertEqual(result.mapping_quality, original.mapping_quality)
                    if expected is None or expected.is_unmapped or expected.is_supplementary:
                        continue

                    self.assertEqual(result.reference_name, expected.reference_name)
                    position_delta = abs(
                        result.reference_start - expected.reference_start
                    )
                    self.assertLess(position_delta, 5)
                    if position_delta == 0:
                        self.assertEqual(result.cigarstring, expected.cigarstring)
                    self.assertEqual(result.is_reverse, expected.is_reverse)
                    self.assertEqual(result.query_sequence, expected.query_sequence)
                    self.assertEqual(result.query_qualities, expected.query_qualities)

                self.assert_pair_consistency(lifted_records, case.name)

    def test_picard_validation(self):
        for case in CHAIN_CASES:
            with self.subTest(case=case.name):
                self.validate_with_picard(
                    self.outputs[case.name], self.output_directory
                )

    def test_long_read_picard_validation(self):
        for input_name in LONG_READ_INPUTS:
            with self.subTest(input=input_name):
                case = LiftCase(
                    f"{Path(input_name).stem}-chain",
                    input_name,
                    input_name,
                    CHAIN_INDEX,
                )
                output = run_lift(case, self.output_directory, "-C")
                self.assert_unsorted_header(output)
                self.validate_with_picard(output, self.output_directory)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run LevioSAM2's extended integration tests."
    )
    parser.add_argument(
        "leviosam",
        nargs="?",
        default=str(LEVIOSAM),
        help="path to the leviosam2 executable (default: build/leviosam2)",
    )
    return parser.parse_known_args()


if __name__ == "__main__":
    arguments, unittest_arguments = parse_arguments()
    LEVIOSAM = Path(arguments.leviosam).expanduser().resolve()
    if not LEVIOSAM.is_file():
        raise SystemExit(f"leviosam2 executable not found: {LEVIOSAM}")
    unittest.main(argv=[sys.argv[0], *unittest_arguments])
