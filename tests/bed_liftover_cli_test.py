#!/usr/bin/env python3
"""CLI regressions for half-open BED intervals on tiny known chains."""
import pathlib, subprocess, sys
exe=pathlib.Path(sys.argv[1]).resolve(); root=pathlib.Path(sys.argv[2]); root.mkdir(parents=True,exist_ok=True)
chain=root/'bed.chain'; fai=root/'target.fai'
chain.write_text('''chain 100 f 2000 + 100 200 d 2000 + 300 400 1
100

chain 100 r 2000 + 100 200 rd 2000 - 600 700 2
100

chain 100 continuous 2000 + 100 200 d 2000 + 1400 1500 13
50 0 0
50

chain 100 sg 2000 + 100 700 d 2000 + 500 600 3
50 500 0
50

chain 100 tg 2000 + 100 200 d 2000 + 700 1300 4
50 0 500
50

chain 100 sg501 2000 + 100 701 d 2000 + 0 100 5
50 501 0
50

chain 100 tg501 2000 + 100 200 d 2000 + 0 601 6
50 0 501
50

chain 100 mix 2000 + 100 150 d 2000 + 100 150 7
50

chain 100 mix 2000 + 150 200 rd 2000 - 600 650 8
50

chain 100 mix 2000 + 200 250 d 2000 + 200 250 9
50

chain 100 ms 2000 + 100 150 d 2000 + 100 150 10
50

chain 100 ms 2000 + 150 200 d 2000 - 1800 1850 11
50

chain 100 ms 2000 + 200 250 d 2000 + 200 250 12
50

''')
fai.write_text('d\t2000\nrd\t2000\n')
def run(args,ok=True):
 p=subprocess.run([str(exe)]+list(map(str,args)),text=True,capture_output=True)
 if ok and p.returncode: raise AssertionError(f'{args}: {p.returncode}: {p.stderr}')
 return p
run(['index','-c',chain,'-F',fai,'-p',root/'idx'])
def bedrun(name,lines,g=500):
 p=root/(name+'.in'); p.write_text(''.join(x+'\n' for x in lines)); prefix=root/name
 result=run(['bed','-C',root/'idx.clft','-b',p,'-p',prefix,'-G',str(g)])
 def rows(path): return [l.split('\t') for l in path.read_text().splitlines() if l and not l.startswith('#')]
 return rows(root/(name+'.bed')),rows(root/(name+'-unmapped.bed')),result
# Half-open block limits, single bases, reverse orientation and BED6 strand.
mapped,unmapped,_=bedrun('basic',[
 'f\t100\t101\tfirst\t0\t+', 'f\t199\t200\tlast\t0\t-',
 'r\t100\t101\trev_first\t0\t+', 'r\t199\t200\trev_last\t0\t-',
 'r\t110\t120\trev_multi\t0\t+', 'f\t110\t120', 'continuous\t140\t160\tadjacent',
 'missing\t1\t2\tunknown\t0\t+', 'f\t200\t201\toutside\t0\t+'],0)
expected={('d','300','301','first','0','+'),('d','399','400','last','0','-'),('rd','1399','1400','rev_first','0','-'),('rd','1300','1301','rev_last','0','+'),('rd','1380','1390','rev_multi','0','-'),('d','310','320'),('d','1440','1460','adjacent')}
assert {tuple(r) for r in mapped}==expected,(mapped,expected)
assert {r[3] for r in unmapped}=={'unknown','outside'},unmapped
# Partially clipped endpoints are projected to the actual aligned boundary.
mapped,unmapped,_=bedrun('clip',['f\t95\t105\tleft','f\t195\t205\tright','r\t95\t105\trleft','r\t195\t205\trright'],5)
assert {tuple(r) for r in mapped}=={('d','300','305','left'),('d','395','400','right'),('rd','1395','1400','rleft'),('rd','1300','1305','rright')},mapped
assert not unmapped
_,unmapped,_=bedrun('clipzero',['f\t95\t105\tleft','f\t195\t205\tright'],0); assert len(unmapped)==2
# Inclusive source and target gap threshold, with one-bp larger controls.
lines=['sg\t140\t660\tsg500','tg\t140\t160\ttg500','sg501\t140\t661\tsg501','tg501\t140\t160\ttg501']
m,u,_=bedrun('g500',lines,500); assert {r[3] for r in m}=={'sg500','tg500'},(m,u); assert {r[3] for r in u}=={'sg501','tg501'},(m,u)
m,u,_=bedrun('g0',lines,0); assert not m and len(u)==4,(m,u)
# Inconsistent intermediate mapping must not be hidden by matching ends.
m,u,_=bedrun('mixed',['mix\t140\t210\tmixed','ms\t140\t210\tstrand'],500); assert not m and len(u)==2,(m,u)
# BED12 blocks and thick region must move with the interval and flip order.
lines=['f\t110\t140\tfwd12\t0\t+\t115\t135\t0\t2\t5,10,\t0,20,','r\t110\t140\trev12\t0\t+\t115\t135\t0\t2\t5,10,\t0,20,']
m,u,_=bedrun('bed12',lines,0); assert len(m)==2 and not u,(m,u)
for row in m:
 a=row[:10]; sizes=[int(x) for x in row[10].rstrip(',').split(',')]; starts=[int(x) for x in row[11].rstrip(',').split(',')]
 if row[3]=='fwd12': assert a==['d','310','340','fwd12','0','+','315','335','0','2'] and sizes==[5,10] and starts==[0,20],row
 elif row[3]=='rev12': assert a==['rd','1360','1390','rev12','0','-','1365','1385','0','2'] and sizes==[10,5] and starts==[0,25],row
 else: raise AssertionError(row)
# Invalid values must be diagnosed and cannot appear in mapped output.
for i,line in enumerate(['f\t-1\t10\tnegative','f\t10\t10\tempty','f\t30\t20\tbackwards','f\tx\t20\tnonnumeric','f\t10','f\t110\t140\tbad12\t0\t+\t115\t135\t0\t2\t5,\t0,20,','f\t999999999999999999999999\t9999999999999999999999999\toverflow']):
 p=root/f'bad{i}.in'; p.write_text(line+'\n'); prefix=root/f'bad{i}'
 result=run(['bed','-C',root/'idx.clft','-b',p,'-p',prefix],ok=False)
 if result.returncode==0:
  assert not (root/f'bad{i}.bed').read_text().strip(),(line,result.stderr)
  assert (root/f'bad{i}-unmapped.bed').read_text().strip() or result.stderr.strip(),(line,'no diagnostic')
 else: assert result.stderr.strip(),(line,'no diagnostic')
# Missing input and index are errors, not successful empty outputs.
for opts in [['-C',root/'idx.clft','-b',root/'absent.bed'],['-C',root/'absent.clft','-b',root/'basic.in']]:
 result=run(['bed']+opts+['-p',root/'missing'],ok=False); assert result.returncode!=0 and result.stderr.strip(),(opts,result.returncode,result.stderr)
# CLI options and BED metadata lines are accepted or rejected explicitly.
help_result=run(['bed','-h'],ok=False)
assert help_result.returncode==0,(help_result.returncode,help_result.stderr)
assert run(['bed','--help'],ok=False).returncode==0
for bad in ('-1','garbage','2x','99999999999999999999999'):
 result=run(['bed','-C',root/'idx.clft','-b',root/'basic.in','-p',root/'bad-g','-G',bad],ok=False)
 assert result.returncode!=0 and result.stderr.strip(),(bad,result.returncode,result.stderr)
meta=root/'metadata.in'; meta.write_text('track name=example\nbrowser position f:1-100\n# comment\n\nf\t100\t101\tkept\n')
prefix=root/'metadata'
run(['bed','-C',root/'idx.clft','-b',meta,'-p',prefix,'-V','0'])
assert (root/'metadata.bed').read_text().strip()=='d\t300\t301\tkept'
assert not (root/'metadata-unmapped.bed').read_text().strip()
