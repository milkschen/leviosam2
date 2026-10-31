"""Exercise RG preservation using actual BAM operations and core reconcile.

The aligner is a deterministic fixture emitter; samtools operations run through
pysam's bundled implementation. This tests the workflow's metadata path, not
the mapping accuracy of external aligners. Run in the Pixi integration env.
"""

import gzip
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam

import leviosam2


ROOT = pathlib.Path(__file__).resolve().parents[1]
GROUPS = [dict(ID='one', SM='sample', LB='library-one', PL='ILLUMINA'),
          dict(ID='two', SM='sample', LB='library-two', PL='ILLUMINA')]
FALLBACK = dict(ID='fallback', SM='sample', LB='new-library', PL='ILLUMINA')


def header(groups):
    return pysam.AlignmentHeader.from_dict(dict(HD=dict(VN='1.6', SO='unsorted'),
        SQ=[dict(SN='chr1', LN=1000)], RG=groups))


def alignment(name, group, paired=False, score=100, start=101):
    tags = f'\tAS:i:{score}\tNM:i:{max(0, (300 - score) // 100)}' + (f'\tRG:Z:{group}' if group else '')
    first = f'{name}\t{99 if paired else 0}\tchr1\t{start}\t60\t10M\t'
    if not paired:
        return [first + '*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII' + tags]
    return [first + f'=\t{start + 30}\t40\tAAAAAAAAAA\tIIIIIIIIII' + tags,
        f'{name}\t147\tchr1\t{start + 30}\t60\t10M\t=\t{start}\t-40\tTTTTTTTTTT\tIIIIIIIIII' + tags]


def write_bam(path, groups, lines):
    bam_header = header(groups)
    with pysam.AlignmentFile(path, 'wb', header=bam_header) as stream:
        for line in lines:
            stream.write(pysam.AlignedSegment.fromstring(line, bam_header))


class ReadGroupWorkflowIntegration(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = pathlib.Path(self.temporary.name)
        self.bridge = self.root / 'samtools.py'
        self.bridge.write_text(f'#!{sys.executable}\nimport pysam, sys\n'
                               'getattr(pysam, sys.argv[1])(*sys.argv[2:], catch_stdout=False)\n')
        self.target = self.root / 'aligner-output.sam'
        self.aligner = self.root / 'aligner.py'
        self.aligner.write_text(f'#!{sys.executable}\nimport pathlib, sys\n'
                                f'sys.stdout.write(pathlib.Path({str(self.target)!r}).read_text())\n')
        self.bgzip = self.root / 'bgzip.py'
        self.bgzip.write_text(f'#!{sys.executable}\nimport gzip, sys\n'
                             'sys.stdout.buffer.write(gzip.compress(sys.stdin.buffer.read()))\n')
        for tool in (self.bridge, self.aligner, self.bgzip):
            tool.chmod(0o755)

    def workflow(self, paired=False):
        argv = ['leviosam2.py', '-s', 'ilmn_pe' if paired else 'ilmn_se',
                '-a', 'bwamem', '-i', str(self.root / 'input.bam'),
                '-o', str(self.root / 'result'), '-C', 'unused.clft',
                '-f', str(self.root / 'target.fa'), '--forcerun', '--keep_tmp_files',
                '--aligner_exe', str(self.aligner),
                '--samtools_exe', str(self.bridge),
                '--leviosam2_exe', str(ROOT / 'build' / 'leviosam2'),
                '--read_group', r'@RG\tID:fallback\tSM:sample\tLB:new-library\tPL:ILLUMINA']
        with mock.patch.object(sys, 'argv', argv):
            workflow = leviosam2.Leviosam2Workflow(leviosam2.parse_args())
        workflow._set_filenames()
        workflow.bgzip = str(self.bgzip)
        return workflow

    def assert_final_groups(self, workflow, expected):
        with pysam.AlignmentFile(workflow._fn_final) as stream:
            actual_header = stream.header.to_dict()
            reads = list(stream)
        actual = {}
        for read in reads:
            actual.setdefault(read.query_name, set()).add(read.get_tag('RG') if read.has_tag('RG') else None)
        self.assertEqual(actual, {name: {group} for name, group in expected.items()})
        groups = {group['ID']: group for group in actual_header['RG']}
        self.assertEqual(groups['one']['LB'], 'library-one')
        self.assertEqual(groups['two']['LB'], 'library-two')
        self.assertEqual(groups['fallback']['LB'], 'new-library')
        self.assertEqual(len(groups), 3)
        self.assertTrue(workflow._fn_final.with_suffix('.bam.bai').is_file())
        return reads

    def test_single_end_pipeline_preserves_libraries_and_uses_fallback(self):
        workflow = self.workflow()
        originals = (alignment('a', 'one') + alignment('b', 'two')
                     + alignment('c', None))
        write_bam(workflow._fn_deferred, GROUPS, originals)
        write_bam(workflow._fn_committed_sorted, GROUPS, alignment('committed', 'one', start=11))
        target = (alignment('a', 'fallback') + alignment('b', 'fallback')
                  + alignment('c', 'fallback'))
        # Secondary and supplementary output for an existing read inherits RG too.
        target += [target[0].replace('\t0\tchr1', '\t256\tchr1', 1),
                   target[1].replace('\t0\tchr1', '\t2048\tchr1', 1)]
        self.target.write_text(str(header([FALLBACK])) + '\n'.join(target) + '\n')
        workflow.run_bam_to_fastq_se()
        with gzip.open(workflow._fn_deferred_fq_se, 'rt') as stream:
            self.assertNotIn('RG:', stream.read())
        workflow.run_realign_deferred()
        workflow.run_merge_and_index()
        self.assert_final_groups(workflow, dict(a='one', b='two', c='fallback', committed='one'))

    def test_paired_reconcile_source_and_target_keep_original_groups(self):
        workflow = self.workflow(paired=True)
        original_lines, target_lines = [], []
        for name, group, source_score, target_score in [
            ('a', 'one', 200, 100), ('b', 'two', 100, 200),
            ('c', None, 50, 200), ('d', None, 300, 100),
        ]:
            original_lines.extend(alignment(name, group, paired=True, score=source_score))
            target_lines.extend(alignment(name, 'fallback', paired=True, score=target_score))
        write_bam(workflow._fn_deferred, GROUPS, original_lines)
        write_bam(workflow._fn_committed_sorted, GROUPS, alignment('committed', 'one', paired=True, start=11))
        self.target.write_text(str(header([FALLBACK])) + '\n'.join(target_lines) + '\n')
        workflow.run_collate_pe()
        workflow.run_realign_deferred()
        workflow.run_refflow_merge_pe()
        workflow.run_merge_and_index()
        reads = self.assert_final_groups(workflow,
            dict(a='one', b='two', c='fallback', d=None, committed='one'))
        selected = {read.query_name: read.get_tag('RF') for read in reads if read.has_tag('RF')}
        self.assertEqual(selected, dict(a='source', b='target', c='target', d='source'))

    def test_header_collision_stops_before_final_merge(self):
        workflow = self.workflow()
        write_bam(workflow._fn_committed_sorted, GROUPS, alignment('a', 'one'))
        conflicting = [dict(GROUPS[0], LB='another-library')]
        write_bam(workflow._fn_deferred_realigned_se, conflicting, alignment('b', 'one'))
        with self.assertRaises(subprocess.CalledProcessError):
            workflow.run_merge_and_index()
        self.assertFalse(workflow._fn_final.exists())


if __name__ == '__main__':
    unittest.main()
