"""
Nae-Chyun Chen
Johns Hopkins University
2021-2023

Tests of the workflow/leviosam2.py workflow.
"""
import argparse
import copy
import pathlib
import sys
import unittest

import leviosam2

TIME_CMDS = ["", "time -v -ao test.time_log"]


class WorkflowStatic(unittest.TestCase):
    def test_validate_exe(self):
        # `ls` should not fail
        leviosam2.validate_executable(cmd="ls")

    def test_check_input_exists(self):
        leviosam2.check_file_exists(pathlib.Path.cwd() / sys.argv[0])
        with self.assertRaises(FileNotFoundError):
            leviosam2.check_file_exists(
                pathlib.Path.cwd() / f"{sys.argv[0]}-not-a-file"
            )


class Workflow(unittest.TestCase):
    """Tests of the leviosam2 workflow script."""

    def setUp(self) -> None:
        """Set up for every test."""
        args = argparse.Namespace
        args.aligner = "bowtie2"
        args.aligner_exe = "auto"
        args.sequence_type = "ilmn_pe"
        args.out_prefix = "out/prefix"
        args.out_format = "bam"
        args.leviosam2_index = "test.clft"
        args.num_threads = 4
        args.source_label = "source"
        args.target_label = "target"
        args.dryrun = True
        args.forcerun = False
        args.samtools_exe = "samtools"
        args.bgzip_exe = "bgzip"
        args.gnu_time_exe = "gtime"
        args.leviosam2_exe = "leviosam2"
        args.measure_time = False
        args.target_fasta = "target.fasta"
        args.target_aligner_index = "target.idx"
        args.input_bam = "input/aln.bam"
        args.read_group = "rg"
        args.lift_commit_min_mapq = None
        args.lift_commit_min_score = None
        args.lift_commit_max_frac_clipped = None
        args.lift_commit_max_isize = None
        args.lift_commit_max_hdist = None
        args.lift_max_gap = None
        args.lift_bed_commit_source = None
        args.lift_bed_defer_target = None
        args.lift_realign_config = None
        args.keep_tmp_files = True

        self.args = args
        self.workflow = leviosam2.Leviosam2Workflow(args)
        self.workflow._set_filenames()

    # TODO
    # def test_set_filenames_single_end(self):
    #     pass

    # TODO
    # def test_set_filenames_paired_end(self):
    #     pass

    def test_infer_aligner_exe(self):
        for aln in ["bowtie2", "minimap2", "winnowmap2", "strobealign"]:
            self.workflow.aligner = aln
            self.workflow._infer_aligner_exe()
            self.assertEqual(self.workflow.aligner_exe, aln)

        self.workflow.aligner = "bwamem"
        self.workflow._infer_aligner_exe()
        self.assertEqual(self.workflow.aligner_exe, "bwa")

        self.workflow.aligner = "bwamem2"
        self.workflow._infer_aligner_exe()
        self.assertEqual(self.workflow.aligner_exe, "bwa-mem2")

    def test_run_leviosam2_basic(self):
        result = self.workflow.run_leviosam2()
        expected = (
            f"{self.args.leviosam2_exe} lift -C {self.args.leviosam2_index} "
            f"-O bam -a {self.args.input_bam} "
            f"-p {self.args.out_prefix} -t 4 -m "
            f"-f {self.args.target_fasta} "
        )
        self.assertEqual(result, expected)

    def test_run_leviosam2_basic_cram(self):
        self.workflow.out_format = "cram"
        result = self.workflow.run_leviosam2()
        # TODO
        # When using `-O cram` in the workflow, we still write to the BAM format
        # in leviosam2 at this moment. We need to look into the way to make a
        # right cram file in leviosam2.
        expected = (
            f"{self.args.leviosam2_exe} lift -C {self.args.leviosam2_index} "
            f"-O bam -a {self.args.input_bam} "
            f"-p {self.args.out_prefix} -t 4 -m "
            f"-f {self.args.target_fasta} "
        )
        self.assertEqual(result, expected)

    def test_run_leviosam2_defer(self):
        new_workflow = copy.deepcopy(self.workflow)
        new_workflow.lift_commit_min_mapq = 10
        new_workflow.lift_commit_min_score = -10
        result = new_workflow.run_leviosam2()
        expected = (
            f"{self.args.leviosam2_exe} lift -C {self.args.leviosam2_index} "
            f"-O bam "
            f"-a {self.args.input_bam} "
            f"-p {self.args.out_prefix} -t 4 -m -f {self.args.target_fasta} "
            f"-S mapq:10 -S aln_score:-10 "
        )
        self.assertEqual(result, expected)

    def test_run_leviosam2_defer_comprehensive(self):
        new_workflow = copy.deepcopy(self.workflow)
        new_workflow.lift_commit_min_mapq = 30
        new_workflow.lift_commit_min_score = 100
        new_workflow.lift_commit_max_frac_clipped = 0.95
        new_workflow.lift_commit_max_isize = 1000
        new_workflow.lift_commit_max_hdist = 5
        new_workflow.lift_max_gap = 20
        new_workflow.lift_bed_commit_source = "commit/source.bed"
        new_workflow.lift_bed_defer_target = "defer/target.bed"
        new_workflow.lift_realign_config = "configs/ilmn.yaml"

        result = new_workflow.run_leviosam2()
        expected = (
            f"{self.args.leviosam2_exe} lift -C {self.args.leviosam2_index} "
            f"-O bam -a {self.args.input_bam} "
            f"-p {self.args.out_prefix} -t 4 -m -f {self.args.target_fasta} "
            f"-S mapq:30 -S aln_score:100 -S clipped_frac:0.95 "
            f"-S isize:1000 -S hdist:5 -G 20 "
            f"-r commit/source.bed -D defer/target.bed "
            f"-x configs/ilmn.yaml "
        )
        self.assertEqual(result, expected)

    def test_run_sort_committed(self):
        self.workflow.time_cmd = "gtime -v -ao test.time_log "

        result = self.workflow.run_sort_committed()
        expected = (
            f"gtime -v -ao test.time_log samtools sort -@ 4 "
            f"-o {self.args.out_prefix}-committed-sorted.bam "
            f"{self.args.out_prefix}-committed.bam"
        )
        self.assertEqual(result, expected)

    def test_run_collate_pe(self):
        result = self.workflow.run_collate_pe()
        expected = (
            f"{self.args.leviosam2_exe} collate "
            f"-a {self.args.out_prefix}-committed-sorted.bam "
            f"-b {self.args.out_prefix}-deferred.bam "
            f"-p {self.args.out_prefix}-paired"
        )
        self.assertEqual(result, expected)

    def test_run_realign_deferred_bt2(self):
        """Tests realign deferred (bt2, paired-end)"""
        result = self.workflow.run_realign_deferred()
        expected = (
            f"bowtie2 rg -p 3 -x target.idx "
            f"-1 {self.args.out_prefix}-paired-deferred-R1.fq.gz "
            f"-2 {self.args.out_prefix}-paired-deferred-R2.fq.gz |  "
            f"samtools view -hbo {self.args.out_prefix}-paired-realigned.bam"
        )
        self.assertEqual(result, expected)

    def test_run_realign_deferred_bt2_single_end(self):
        """Tests realign deferred (bt2, single-end)"""
        self.workflow.is_single_end = True
        result = self.workflow.run_realign_deferred()
        expected = (
            f"bowtie2 rg -p 3 -x target.idx "
            f"-U {self.args.out_prefix}-deferred.fq.gz | "
            f"samtools sort -@ 1 -o {self.args.out_prefix}-realigned.bam"
        )
        self.assertEqual(result, expected)

    def test_run_realign_deferred_bwa(self):
        """Tests realign deferred (bwa, paired-end)"""
        self.workflow.aligner = "bwamem"
        self.workflow.aligner_exe = "bwa"
        result = self.workflow.run_realign_deferred()
        expected = (
            f"bwa mem -R rg -t 3 target.idx "
            f"{self.args.out_prefix}-paired-deferred-R1.fq.gz "
            f"{self.args.out_prefix}-paired-deferred-R2.fq.gz |  "
            f"samtools view -hbo {self.args.out_prefix}-paired-realigned.bam"
        )
        self.assertEqual(result, expected)

    def test_run_realign_deferred_bwa_single_end(self):
        """Tests realign deferred (bwa, single-end)"""
        self.workflow.is_single_end = True
        self.workflow.aligner = "bwamem"
        self.workflow.aligner_exe = "bwa"
        result = self.workflow.run_realign_deferred()
        expected = (
            f"bwa mem -R rg -t 3 target.idx "
            f"{self.args.out_prefix}-deferred.fq.gz | "
            f"samtools sort -@ 1 -o {self.args.out_prefix}-realigned.bam"
        )
        self.assertEqual(result, expected)

    def test_run_bam_to_fastq_se(self):
        result = self.workflow.run_bam_to_fastq_se()
        expected = (
            f"samtools fastq {self.args.out_prefix}-deferred.bam | "
            f"bgzip > {self.args.out_prefix}-deferred.fq.gz"
        )
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
