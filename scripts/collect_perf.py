'''
Collect computational performance from a collection of GNU time reports.

Usage: python collect_perf.py -a bt2_all.time_log -l lift.time_log -l aln.time_log

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--aln-log',
        help='Paths to the GNU time log for full alignment.')
    parser.add_argument(
        '-l', '--leviosam-logs',
        action='append', required=True,
        help='Paths to GNU time logs for the levioSAM pipeline.')
    parser.add_argument(
        '-o', '--output',
        help='Path to the output TSV file.'
    )
    args = parser.parse_args()
    return args


def collect_perf_core(logs_list, cols):
    print(logs_list)
    ls_perf = []
    for log in logs_list:
        f = open(log, 'r')
        task = os.path.basename(log).split('.')[0]
        for line in f:
            line = line.rstrip()
            if line.count('User time (seconds):') > 0:
                usr_time = float(line.split('):')[1])
            elif line.count('System time (seconds):') > 0:
                sys_time = float(line.split('):')[1])
            elif line.count('Elapsed (wall clock) time (h:mm:ss or m:ss):') > 0:
                wt = line.split('):')[1].split(':')
                if len(wt) == 3:
                    wall_time = 60 * 60 * float(wt[0]) + 60 * float(wt[1]) + float(wt[2])
                elif len(wt) == 2:
                    wall_time = 60 * float(wt[0]) + float(wt[1])
                else:
                    print('error - invalid wall time format')
                    exit(1)
            elif line.count('Maximum resident set size (kbytes):') > 0:
                max_rss = int(line.split('):')[1])
        f.close()
        cpu_time = usr_time + sys_time
        ls_perf.append([
            task, usr_time, sys_time,
            cpu_time, wall_time, max_rss])
    df = pd.DataFrame(ls_perf, columns=cols)
    return df


def collect_perf(args):
    COLS = ['Task', 'User time (s)', 'System time (s)',
            'CPU time (s)', 'Wall time (s)', 'Max RSS (KB)']

    df_l = collect_perf_core(args.leviosam_logs, cols=COLS)
    l_sum = ['summary']
    for i in range(1, 5):
        l_sum.append(sum(df_l.iloc[:, i]))
    l_sum.append(max(df_l.iloc[:, 5]))
    df_l_sum = pd.DataFrame([l_sum], columns=COLS)

    df_l = df_l.append(df_l_sum, ignore_index=True)
    print(df_l)

    df_a = collect_perf_core([args.aln_log], cols=COLS)
    print(df_a)


if __name__ == '__main__':
    args = parse_args()
    collect_perf(args)
