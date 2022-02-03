'''
Collect computational performance from a collection of GNU time reports.

Usage: 
```
python collect_perf.py -a bt2_all.time_log -l lift.time_log -l collate.time_log \
-l to_fastq.time_log -l aln_paired.time_log -l aln_unpaired.time_log \
-l merge.time_log -l sort_all.time_log
```

Nae-Chyun Chen
Johns Hopkins University
2021-2022
'''
import argparse
import os
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--aln-log',
        help='Paths to the GNU time log for full alignment.')
    parser.add_argument(
        '-an', '--aln-name', default='aln',
        help='Label of the full alignment experiment [aln].')
    parser.add_argument(
        '-l', '--leviosam-logs',
        action='append', required=True,
        help='Paths to GNU time logs for the levioSAM pipeline.')
    parser.add_argument(
        '-ln', '--leviosam-name', default='leviosam',
        help='Label of the levioSAM experiment [leviosam].')
    parser.add_argument(
        '-o', '--output',
        help='Path to the output TSV file.'
    )
    args = parser.parse_args()
    return args


def collect_perf_core(f, ls_perf) -> None:
    for line in f:
        line = line.rstrip()
        if line.count('Command being timed:') > 0:
            cmd = line.split('"')[1].split()
            program = cmd[0].split('/')[-1]
            task = program + '_' + cmd[1]
        elif line.count('User time (seconds):') > 0:
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
                print('error - invalid wall time format', file=sys.stderr)
                exit(1)
        elif line.count('Maximum resident set size (kbytes):') > 0:
            max_rss = int(line.split('):')[1])
            cpu_time = usr_time + sys_time
            ls_perf.append([
                task, usr_time, sys_time,
                cpu_time, wall_time, max_rss])
    return


def collect_perf_list(logs_list, cols):
    print(logs_list, file=sys.stderr)
    ls_perf = []
    for log in logs_list:
        f = open(log, 'r')
        # task = os.path.basename(log).split('.')[0]
        collect_perf_core(f, ls_perf)
        f.close()
    df = pd.DataFrame(ls_perf, columns=cols)
    return df


def summarize_df(df, cols):
    l_sum = ['summary']
    for i in range(1, 5):
        l_sum.append(sum(df.iloc[:, i]))
    l_sum.append(max(df.iloc[:, 5]))
    df_sum = pd.DataFrame([l_sum], columns=cols)
    df = df.append(df_sum, ignore_index=True)
    return df


def collect_perf(args):
    COLS = ['Task', 'User time (s)', 'System time (s)',
            'CPU time (s)', 'Wall time (s)', 'Max RSS (KB)']

    df_l = collect_perf_list(args.leviosam_logs, cols=COLS)
    df_l['Method'] = args.aln_name
    # df_l = summarize_df(df=df_l, cols=COLS)
    # print(df_l)

    df_a = collect_perf_list([args.aln_log], cols=COLS)
    df_a['Method'] = args.leviosam_name
    # df_a = summarize_df(df=df_a, cols=COLS)
    # print(df_a)

    df = pd.concat([df_l, df_a])
    if args.output:
        df.to_csv(args.output, sep='\t', index=False)
    else:
        print(df)

if __name__ == '__main__':
    args = parse_args()
    collect_perf(args)
