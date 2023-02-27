'''
Tests of the workflow/leviosam2.py workflow.
'''
import argparse
import pathlib
import sys
import unittest

import leviosam2

TIME_CMDS = ['', 'time -v -ao test.time_log']


class WorkflowStatic(unittest.TestCase):

    def test_validate_exe(self):
        # `ls` should not fail
        leviosam2.Leviosam2Workflow.validate_exe(cmd='ls')

    def test_check_input_exists(self):
        leviosam2.Leviosam2Workflow._check_input_exists(pathlib.Path.cwd() /
                                                        sys.argv[0])
        with self.assertRaises(FileNotFoundError):
            leviosam2.Leviosam2Workflow._check_input_exists(
                pathlib.Path.cwd() / f'{sys.argv[0]}-not-a-file')


class Workflow(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        args = argparse.Namespace
        args.aligner = 'TBD'
        args.aligner_exe = 'TBD'
        args.sequence_type = 'TBD'
        args.out_prefix = 'out/prefix'
        args.num_threads = 4
        args.source_label = 'source'
        args.target_label = 'target'
        args.dryrun = True
        args.forcerun = False
        args.samtools_binary = 'samtools'
        args.bgzip_binary = 'bgzip'
        args.gnu_time_binary = 'gtime'
        args.leviosam2_binary = 'leviosam2'
        args.measure_time = False
        args.target_fasta = 'target.fasta'
        args.input_alignment = 'input/aln.bam'

        cls.workflow = leviosam2.Leviosam2Workflow(args)
        cls.workflow._set_filenames()

    def test_infer_aligner_exe(self):
        for aln in ['bowtie2', 'minimap2', 'winnowmap2', 'strobealign']:
            self.workflow.aligner = aln
            self.workflow._infer_aligner_exe()
            self.assertEqual(self.workflow.aligner_exe, aln)

        self.workflow.aligner = 'bwamem'
        self.workflow._infer_aligner_exe()
        self.assertEqual(self.workflow.aligner_exe, 'bwa')

        self.workflow.aligner = 'bwamem2'
        self.workflow._infer_aligner_exe()
        self.assertEqual(self.workflow.aligner_exe, 'bwa-mem2')

    def test_run_leviosam2_basic(self):
        result = self.workflow.run_leviosam2(clft='test.clft')
        expected = (f' leviosam2 lift -C test.clft -a input/aln.bam '
                    f'-p out/prefix -t 4 -m -f target.fasta ')
        self.assertEqual(result, expected)

    def test_run_leviosam2_defer(self):
        result = self.workflow.run_leviosam2(clft='test.clft',
                                             lift_commit_min_mapq=10,
                                             lift_commit_min_score=-10)
        expected = (f' leviosam2 lift -C test.clft -a input/aln.bam '
                    f'-p out/prefix -t 4 -m -f target.fasta '
                    f'-S mapq:10 -S aln_score:-10 ')
        self.assertEqual(result, expected)

    def test_run_leviosam2_defer_comprehensive(self):
        result = self.workflow.run_leviosam2(
            clft='test.clft',
            lift_commit_min_mapq=30,
            lift_commit_min_score=100,
            lift_commit_max_frac_clipped=0.95,
            lift_commit_max_isize=1000,
            lift_commit_max_hdist=5,
            lift_max_gap=20,
            lift_bed_commit_source='commit/source.bed',
            lift_bed_defer_target='defer/target.bed',
            lift_realign_config='configs/ilmn.yaml')
        expected = (f' leviosam2 lift -C test.clft -a input/aln.bam '
                    f'-p out/prefix -t 4 -m -f target.fasta '
                    f'-S mapq:30 -S aln_score:100 -S clipped_frac:0.95 '
                    f'-S isize:1000 -S hdist:5 -G 20 '
                    f'-r commit/source.bed -D defer/target.bed '
                    f'-x configs/ilmn.yaml ')
        self.assertEqual(result, expected)

    def test_run_sort_committed(self):
        self.workflow.time_cmd = 'gtime -v -ao test.time_log'

        result = self.workflow.run_sort_committed()
        expected = (
            f'gtime -v -ao test.time_log samtools sort -@ 4 '
            '-o out/prefix-committed-sorted.bam out/prefix-committed.bam')
        self.assertEqual(result, expected)

    def test_run_collate_pe(self):
        result = self.workflow.run_collate_pe()
        expected = (f' leviosam2 collate '
                    '-a out/prefix-committed-sorted.bam '
                    '-b out/prefix-deferred.bam -p out/prefix-paired')
        self.assertEqual(result, expected)

    # def test_run_realign_deferred_pe(self):
    #     raise NotImplementedError
    #     # leviosam2.run_realign_deferred_pe()


if __name__ == '__main__':
    unittest.main()