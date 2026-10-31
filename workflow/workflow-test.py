"""
Nae-Chyun Chen
Johns Hopkins University
2021-2023

Tests of the workflow/leviosam2.py workflow.
"""
import argparse
import copy
import json
import os
import shlex
import pathlib
import sys
import subprocess
import tempfile
from unittest import mock
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

    def test_workflow_output_format_choices(self):
        argv = ["leviosam2.py", "-s", "ilmn_se", "-a", "bowtie2",
                "-i", "input.bam", "-o", "result", "-C", "index.clft",
                "-f", "target.fa"]
        for output_format in ("bam", "cram"):
            with self.subTest(output_format=output_format), mock.patch.object(
                    sys, "argv", argv + ["-O", output_format]):
                self.assertEqual(leviosam2.parse_args().out_format, output_format)
        with mock.patch.object(sys, "argv", argv + ["-O", "sam"]), \
                mock.patch("sys.stderr"), self.assertRaises(SystemExit) as error:
            leviosam2.parse_args()
        self.assertEqual(error.exception.code, 2)

    def failure_workflow(self, directory):
        workflow = copy.deepcopy(self.workflow)
        workflow.dryrun = False
        workflow.forcerun = True
        workflow.keep_tmp_files = False
        workflow.out_prefix = pathlib.Path(directory) / "result"
        workflow.input_bam = pathlib.Path(directory) / "input.bam"
        workflow.input_bam.touch()
        workflow._set_filenames()
        return workflow

    def test_failed_lift_stops_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = self.failure_workflow(directory)
            workflow.leviosam2 = "false"
            with mock.patch.object(workflow, "run_sort_committed") as sort, \
                    mock.patch.object(workflow, "run_clean") as clean:
                with self.assertLogs(leviosam2.logger, level="ERROR") as logs:
                    with self.assertRaises(subprocess.CalledProcessError):
                        workflow.run_workflow()
                self.assertIn("run_leviosam2 failed", logs.output[0])
                sort.assert_not_called()
                clean.assert_not_called()

    def test_failed_pipeline_producer_stops_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = self.failure_workflow(directory)
            workflow.is_single_end = True
            workflow._fn_deferred.touch()
            workflow.samtools = "false"
            workflow.bgzip = "cat"
            with mock.patch.object(workflow, "run_leviosam2"), \
                    mock.patch.object(workflow, "run_sort_committed"), \
                    mock.patch.object(workflow, "run_realign_deferred") as realign, \
                    mock.patch.object(workflow, "run_clean") as clean:
                with self.assertRaises(subprocess.CalledProcessError):
                    workflow.run_workflow()
                realign.assert_not_called()
                clean.assert_not_called()

    def test_failed_aligner_pipeline_stops_before_merge(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = self.failure_workflow(directory)
            workflow._fn_deferred_pe_fq1.touch()
            workflow._fn_deferred_pe_fq2.touch()
            workflow.aligner_exe = "false"
            # The consumer succeeds even though the aligner fails.
            workflow.samtools = "true"
            with mock.patch.object(workflow, "run_leviosam2"), \
                    mock.patch.object(workflow, "run_sort_committed"), \
                    mock.patch.object(workflow, "run_collate_pe"), \
                    mock.patch.object(workflow, "run_refflow_merge_pe") as merge, \
                    mock.patch.object(workflow, "run_clean") as clean:
                with self.assertRaises(subprocess.CalledProcessError):
                    workflow.run_workflow()
                merge.assert_not_called()
                clean.assert_not_called()

    def test_failed_cram_conversion_or_index_retains_bam(self):
        for failed_stage in ("view", "index"):
            with self.subTest(failed_stage=failed_stage), \
                    tempfile.TemporaryDirectory() as directory:
                workflow = self.failure_workflow(directory)
                workflow.out_format = "cram"
                workflow._fn_final.write_text("original BAM")
                bam_index = workflow._fn_final.with_suffix(".bam.bai")
                bam_index.touch()
                marker = pathlib.Path(directory) / "index-ran"
                tool = pathlib.Path(directory) / "samtools.sh"
                tool.write_text(
                    '#!/bin/sh\n'
                    'if [ "$1" = "index" ]; then touch ' + str(marker) + '; fi\n'
                    'if [ "$1" = "' + failed_stage + '" ]; then exit 23; fi\n'
                )
                tool.chmod(0o755)
                workflow.samtools = str(tool)
                # Isolate the finalization path; all preceding stages succeeded.
                stages = ["run_leviosam2", "run_sort_committed", "run_collate_pe",
                          "run_realign_deferred", "run_refflow_merge_pe",
                          "run_merge_and_index"]
                patches = [mock.patch.object(workflow, stage) for stage in stages]
                for patch in patches:
                    patch.start()
                try:
                    with mock.patch.object(workflow, "run_clean") as clean:
                        with self.assertRaises(subprocess.CalledProcessError) as error:
                            workflow.run_workflow()
                        self.assertEqual(error.exception.returncode, 23)
                        clean.assert_not_called()
                    self.assertEqual(workflow._fn_final.read_text(), "original BAM")
                    self.assertTrue(bam_index.is_file())
                    self.assertEqual(marker.exists(), failed_stage == "index")
                finally:
                    for patch in patches:
                        patch.stop()

    # TODO
    # def test_set_filenames_single_end(self):
    #     pass

    # TODO
    # def test_set_filenames_paired_end(self):
    #     pass

    def test_shell_arguments_with_metacharacters(self):
        with tempfile.TemporaryDirectory(prefix="workflow space ; ") as directory:
            root = pathlib.Path(directory)
            log = root / "argv log"
            tool = root / "fake tool ; name"
            tool.write_text(
                "#!/usr/bin/env python3\n"
                "import json, os, sys\n"
                "with open(os.environ['WORKFLOW_ARG_LOG'], 'a') as f: f.write(json.dumps(sys.argv[1:]) + '\\n')\n"
                "print('first line\\nsecond line')\n"
            )
            tool.chmod(0o755)
            self.args.measure_time = True
            self.args.gnu_time_exe = str(tool)
            self.args.out_prefix = str(root / "timed ; result")
            timed = leviosam2.Leviosam2Workflow(self.args)
            self.assertEqual(shlex.split(timed.time_cmd), [str(tool), "-v", "-ao", self.args.out_prefix + ".time_log"])
            self.args.measure_time = False
            workflow = self.failure_workflow(directory)
            workflow.out_prefix = root / "result ; $(touch INJECTED)"
            workflow.target_fasta = root / "target ; ref"
            workflow.leviosam2_index = root / "index ; file"
            workflow.target_aligner_index = root / "aligner index ; file"
            workflow.lift_bed_commit_source = str(root / "bed ; file")
            workflow.source_label = "source ; label"
            workflow.target_label = "target ; label"
            workflow.read_group = "ID:foo SM:bar ; $(touch INJECTED)"
            workflow.leviosam2 = workflow.samtools = workflow.bgzip = workflow.aligner_exe = str(tool)
            workflow._set_filenames()
            for name in ("_fn_committed", "_fn_committed_sorted", "_fn_deferred", "_fn_deferred_pe_fq1", "_fn_deferred_pe_fq2", "_fn_deferred_realigned_pe", "_fn_deferred_pe", "_fn_committed_pe", "_fn_deferred_reconciled", "_fn_final"):
                getattr(workflow, name).touch()
            with mock.patch.dict(os.environ, {"WORKFLOW_ARG_LOG": str(log)}):
                workflow.validate_executables()
                workflow.run_leviosam2()
                workflow.run_sort_committed()
                workflow.run_collate_pe()
                for aligner, flag in (("bowtie2", "--rg"), ("bwamem", "-R"), ("bwamem2", "-R"), ("strobealign", "--rg"), ("minimap2", "-R"), ("winnowmap2", "-R")):
                    workflow.aligner = aligner
                    workflow.read_group = ("--rg-id foo --rg 'SM:sample name' ; $(touch INJECTED)"
                                           if aligner == "bowtie2" else
                                           "ID:foo SM:bar ; $(touch INJECTED)")
                    read_group = workflow.read_group
                    workflow.is_single_end = aligner in ("minimap2", "winnowmap2")
                    workflow._fn_deferred_fq_se.touch()
                    workflow.dryrun = True
                    first = workflow.run_realign_deferred()
                    self.assertEqual(first, workflow.run_realign_deferred())
                    self.assertEqual(workflow.read_group, read_group)
                    workflow.dryrun = False
                    workflow.run_realign_deferred()
                    calls = [json.loads(line) for line in log.read_text().splitlines()]
                    align_call = next(call for call in reversed(calls) if flag in call)
                    if aligner == "bowtie2":
                        self.assertEqual(align_call[align_call.index("--rg-id") + 1], "foo")
                        self.assertEqual(align_call[align_call.index("--rg") + 1], "SM:sample name")
                        self.assertIn(";", align_call)
                    else:
                        self.assertEqual(align_call[align_call.index(flag) + 1], read_group)
                workflow.is_single_end = False
                workflow.run_refflow_merge_pe()
                workflow.run_merge_and_index()
                workflow.run_bam_to_fastq_se()
                workflow.run_bam_to_cram()
            calls = [json.loads(line) for line in log.read_text().splitlines()]
            lift = next(call for call in calls if call and call[0] == "lift")
            self.assertEqual(lift[lift.index("-p") + 1], str(workflow.out_prefix))
            self.assertEqual(lift[lift.index("-r") + 1], workflow.lift_bed_commit_source)
            reconcile = next(call for call in calls if call and call[0] == "reconcile")
            self.assertIn(workflow.source_label + ":" + str(workflow._fn_deferred_pe_sortn), reconcile)
            self.assertFalse((pathlib.Path.cwd() / "INJECTED").exists())

    def test_preset_combination_matrix(self):
        supported = {
            **{(sequence, aligner): leviosam2.Leviosam2WorkflowIlmnPreset
               for sequence in ("ilmn_se", "ilmn_pe")
               for aligner in ("bowtie2", "bwamem", "bwamem2")},
            **{("pb_hifi", aligner): leviosam2.Leviosam2WorkflowPacbioHifiPreset
               for aligner in ("minimap2", "winnowmap2")},
            **{("ont", aligner): leviosam2.Leviosam2WorkflowOntPreset
               for aligner in ("minimap2", "winnowmap2")},
        }
        self.args.use_preset = True
        for sequence in ("ilmn_pe", "ilmn_se", "pb_hifi", "ont"):
            for aligner in ("bowtie2", "bwamem", "bwamem2", "minimap2", "winnowmap2", "strobealign"):
                with self.subTest(sequence=sequence, aligner=aligner):
                    self.args.sequence_type, self.args.aligner = sequence, aligner
                    with mock.patch.object(leviosam2.Leviosam2Workflow, "check_inputs") as check, \
                         mock.patch.object(leviosam2.Leviosam2Workflow, "validate_executables"), \
                         mock.patch.object(leviosam2.Leviosam2Workflow, "run_workflow"):
                        if (sequence, aligner) in supported:
                            with mock.patch.object(supported[(sequence, aligner)], "__init__", return_value=None) as init:
                                leviosam2.run_workflow(self.args)
                                init.assert_called_once_with(self.args)
                            check.assert_called_once()
                        else:
                            with self.assertRaisesRegex(ValueError, f"{sequence} and {aligner}"):
                                leviosam2.run_workflow(self.args)
                            check.assert_not_called()

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
        new_workflow.lift_commit_max_frac_clipped = 0.05
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
            f"-S mapq:30 -S aln_score:100 -S clipped_frac:0.05 "
            f"-S isize:1000 -S hdist:5 -G 20 "
            f"-r commit/source.bed -D defer/target.bed "
            f"-x configs/ilmn.yaml "
        )
        self.assertEqual(result, expected)

    def test_zero_lift_options_from_cli(self):
        argv = ["leviosam2.py", "-s", "ilmn_se", "-a", "bowtie2",
                "-i", "input.bam", "-o", "result", "-C", "index.clft",
                "-f", "target.fa", "--dryrun", "-g", "0",
                "--lift_commit_min_mapq", "0", "--lift_commit_min_score", "0",
                "--lift_commit_max_frac_clipped", "0",
                "--lift_commit_max_isize", "0", "--lift_commit_max_hdist", "0"]
        with mock.patch.object(sys, "argv", argv):
            workflow = leviosam2.Leviosam2Workflow(leviosam2.parse_args())
        command = workflow.run_leviosam2()
        for option in ("-S mapq:0", "-S aln_score:0", "-S clipped_frac:0.0",
                       "-S isize:0", "-S hdist:0", "-G 0"):
            with self.subTest(option=option):
                self.assertIn(option + " ", command)

    def test_illumina_preset_forwards_zero_gap(self):
        args = copy.copy(self.args)
        args.use_preset = True
        workflow = leviosam2.Leviosam2WorkflowIlmnPreset(args)
        self.assertIn("-G 0 ", workflow.run_leviosam2())

    def test_measure_time_initializes_with_output_prefix(self):
        argv = ["leviosam2.py", "-s", "ilmn_se", "-a", "bowtie2",
                "-i", "input.bam", "-o", "out/result", "-C", "index.clft",
                "-f", "target.fa", "--dryrun", "--measure_time",
                "--gnu_time_exe", "my-time"]
        with mock.patch.object(sys, "argv", argv):
            workflow = leviosam2.Leviosam2Workflow(leviosam2.parse_args())
        workflow._set_filenames()
        self.assertEqual(workflow.time_cmd, "my-time -v -ao out/result.time_log ")
        self.assertTrue(workflow.run_leviosam2().startswith(workflow.time_cmd))
        self.assertTrue(workflow.run_sort_committed().startswith(workflow.time_cmd))

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

    def test_long_read_realign_uses_sequence_preset(self):
        for aligner in ("minimap2", "winnowmap2"):
            for sequence_type, preset in (("pb_hifi", "map-hifi"), ("ont", "map-ont")):
                with self.subTest(aligner=aligner, sequence_type=sequence_type):
                    workflow = copy.deepcopy(self.workflow)
                    workflow.aligner = aligner
                    workflow.aligner_exe = aligner
                    workflow.sequence_type = sequence_type
                    workflow.is_single_end = True
                    self.assertIn(f"-a -x {preset} --MD", workflow.run_realign_deferred())

    def test_run_bam_to_fastq_se(self):
        result = self.workflow.run_bam_to_fastq_se()
        expected = (
            f"samtools fastq {self.args.out_prefix}-deferred.bam | "
            f"bgzip > {self.args.out_prefix}-deferred.fq.gz"
        )
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
