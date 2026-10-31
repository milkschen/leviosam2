"""
Nae-Chyun Chen
Johns Hopkins University
2021-2023

Tests of the workflow/leviosam2.py workflow.
"""
import argparse
import copy
import json
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

    def test_sam_conversion_publishes_only_on_success(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = self.failure_workflow(directory)
            workflow.out_format = "sam"
            workflow._fn_final.write_text("recovery BAM")
            sam = workflow._fn_final.with_suffix(".sam")
            sam.write_text("previous SAM")
            temp = sam.with_name(sam.name + ".tmp")

            def fail(cmd, stage):
                temp.write_text("partial")
                raise subprocess.CalledProcessError(1, cmd)

            with mock.patch.object(workflow, "_run_command", side_effect=fail):
                with self.assertRaises(subprocess.CalledProcessError):
                    workflow.run_bam_to_sam()
            self.assertEqual(sam.read_text(), "previous SAM")
            self.assertEqual(workflow._fn_final.read_text(), "recovery BAM")
            self.assertFalse(temp.exists())

            def succeed(cmd, stage):
                self.assertIn("view -h -O sam", cmd)
                temp.write_text("@HD\tVN:1.6\n")

            with mock.patch.object(workflow, "_run_command", side_effect=succeed):
                workflow.run_bam_to_sam()
            self.assertEqual(sam.read_text(), "@HD\tVN:1.6\n")
            self.assertTrue(workflow._sam_completed)

    def test_sam_conversion_quotes_executable_and_paths(self):
        with tempfile.TemporaryDirectory(prefix="sam space ; ") as directory:
            root = pathlib.Path(directory)
            workflow = self.failure_workflow(directory)
            workflow.out_prefix = root / "result ; $(touch INJECTED)"
            workflow._set_filenames()
            workflow._fn_final.write_text("BAM")
            log = root / "argument log"
            tool = root / "samtools ; fake"
            tool.write_text(
                "#!/usr/bin/env python3\n"
                "import json, pathlib, sys\n"
                f"pathlib.Path({str(log)!r}).write_text(json.dumps(sys.argv[1:]))\n"
                "pathlib.Path(sys.argv[sys.argv.index('-o') + 1]).write_text('@HD\\tVN:1.6\\n')\n"
            )
            tool.chmod(0o755)
            workflow.samtools = str(tool)
            workflow.run_bam_to_sam()
            sam = workflow._fn_final.with_suffix(".sam")
            self.assertEqual(json.loads(log.read_text()),
                             ["view", "-h", "-O", "sam", "-o", str(sam) + ".tmp",
                              str(workflow._fn_final)])
            self.assertTrue(sam.is_file())
            self.assertFalse((pathlib.Path.cwd() / "INJECTED").exists())

    def test_workflow_dispatches_sam_conversion(self):
        workflow = self.workflow
        workflow.out_format = "sam"
        with mock.patch.object(workflow, "run_leviosam2"), \
                mock.patch.object(workflow, "run_sort_committed"), \
                mock.patch.object(workflow, "run_collate_pe"), \
                mock.patch.object(workflow, "run_realign_deferred"), \
                mock.patch.object(workflow, "run_refflow_merge_pe"), \
                mock.patch.object(workflow, "run_merge_and_index"), \
                mock.patch.object(workflow, "run_bam_to_sam") as convert:
            workflow.run_workflow()
        convert.assert_called_once_with()

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
                    'if [ "$1" = "index" ]; then touch ' + str(marker) + '; fi\n'
                    'if [ "$1" = "' + failed_stage + '" ]; then exit 23; fi\n'
                )
                workflow.samtools = "bash " + str(tool)
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
