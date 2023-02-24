'''
Tests of the workflow/leviosam2.py workflow.
'''
import unittest

import leviosam2


class Workflow(unittest.TestCase):

    def test_run_sort_committed(self):
        for time_cmd in ['', 'time -v -o test.time_log']:
            result = leviosam2.run_sort_committed(time_cmd=time_cmd,
                                                  samtools='samtools',
                                                  num_threads=4,
                                                  out_prefix='test',
                                                  dryrun=True,
                                                  forcerun=False)
            expected = (f'{time_cmd} samtools sort -@ 4 '
                        '-o test-committed-sorted.bam test-committed.bam')
            self.assertEqual(result, expected)

    def test_run_collate_pe(self):
        for time_cmd in ['', 'time -v -o test.time_log ']:
            result = leviosam2.run_collate_pe(time_cmd=time_cmd,
                                              leviosam2='leviosam2',
                                              out_prefix='test',
                                              dryrun=True,
                                              forcerun=False)
            expected = (f'{time_cmd} leviosam2 collate '
                        '-a test-committed-sorted.bam '
                        '-b test-deferred.bam -p test-paired')
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()