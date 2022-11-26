'''
Get the count and fraction of each leviosam2 alignment categories
'''
import sys
from collections import defaultdict

import pandas as pd
import pysam

with pysam.AlignmentFile(sys.argv[1]) as f:
    comprehensive_cnt = defaultdict(int)
    old_contig = ''
    print(f.mapped, file=sys.stderr)
    for r in f:
        if r.reference_name != old_contig:
            print(r.reference_name, dict(comprehensive_cnt), file=sys.stderr)
            old_contig = r.reference_name
        if r.is_secondary:
            comprehensive_cnt['Secondary'] += 1
            continue
        if r.is_supplementary:
            comprehensive_cnt['Supplementary'] += 1
            continue
        if r.has_tag('LO') and r.get_tag('LO') == 'L_L':
            comprehensive_cnt['Committed'] += 1
            continue
        if r.is_unmapped and (not r.has_tag('RF')):
            comprehensive_cnt['Suppressed'] += 1
            continue
        if r.has_tag('RF'):
            comprehensive_cnt['Deferred'] += 1
            if r.get_tag('RF') == 'source':
                comprehensive_cnt['Deferred-source'] += 1
            elif r.get_tag('RF') == 'target':
                comprehensive_cnt['Deferred-target'] += 1

total_cnt = 0
cats = ['Deferred', 'Committed', 'Suppress']
cat_cnt = {}
for cat in cats:
    total_cnt += comprehensive_cnt[cat]
    cat_cnt[cat] = comprehensive_cnt[cat]

frac_cnt = {}
for k, v in comprehensive_cnt.items():
    if k in cats:
        frac_cnt[k] = v / total_cnt

df = pd.DataFrame({'Frac': frac_cnt, 'Count': cat_cnt})
df = df.reset_index()
df.columns = ['Category', 'Frac', 'Count']
df.to_csv(f'{sys.argv[1]}.csv', sep='\t', index=None)
