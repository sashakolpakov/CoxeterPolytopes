#!/usr/bin/env python3
"""Native bounded modular-rank certificate, not a compactness search cutoff."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import argparse,hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json';CPP=ROOT/'scripts/fresh_h9_incidence_surgery_20261013_modular.cpp';OUT=ROOT/'results/fresh_h9_incidence_surgery_20261013_modular.json'
SELF=Path(__file__).resolve()
def run(case=None):
 for p in (241,1009):
  assert all(p%d for d in range(2,isqrt(p)+1))
  assert all(any(x*x%p==d for x in range(1,p))for d in (2,3,5))
  assert 2*p*p<2**31
 d=json.loads(SOURCE.read_text());ans=[]
 with tempfile.TemporaryDirectory(prefix='h9_surgery_modular_')as tmp:
  exe=Path(tmp)/'rank';subprocess.run(['c++','-O3','-std=c++17',str(CPP),'-o',str(exe)],check=True)
  for r in sorted(d['cases'],key=lambda r:r['facet_count']):
   if case and r['id']!=case:continue
   G=[[F(x)for x in row]for row in r['source_positive_gram']];n=len(G);P={tuple(p)for p in r['parallel_pairs']};S=r['finite_edge_cover']['minimum_supports'];assert r['finite_edge_cover']['enumeration_complete']
   assert 11<=n<=16 and all(G[i][i]>0 for i in range(n))
   labels=[[1 if i==j else {F(0):2,F(1,4):3,F(1,2):4,F(1):0}[G[i][j]*G[i][j]/(G[i][i]*G[j][j])]for j in range(n)]for i in range(n)]
   B=[b for b in combinations(range(n),11)if not any(i in b and j in b for i,j in P)]
   text=f'{n} {len(S)} {len(B)}\n'+'\n'.join(' '.join(map(str,row))for row in labels)+'\n'
   text+='\n'.join(str(len(s))+' '+' '.join(str(x)for p in s for x in p)for s in S)+'\n'
   text+='\n'.join(' '.join(map(str,b))for b in B)+'\n'
   proc=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True)
   U=[];counts=None
   for line in proc.stdout.splitlines():
    z=line.split()
    if z[0]=='U':U.append(list(map(int,z[1:])))
    else:assert z[0]=='COUNTS';counts=list(map(int,z[1:]))
   assert counts and counts[0]==sum(counts[1:]) and counts[-1]==len(U)
   assert counts[0]==sum(4**len(s)for s in S)
   row={'id':r['id'],'fixed_principal_11_count':len(B),'input_sha256':hashlib.sha256(text.encode()).hexdigest(),'counts':dict(zip(['assignments','certified_rank11_mod241','certified_rank11_mod1009','unresolved'],counts)),'unresolved_support_code':U}
   ans.append(row);print(row['id'],row['counts'],flush=True)
 return {'status':'EXACT_MODULAR_RANK_GATE_WITH_EXPLICIT_UNRESOLVED_LIST','source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[SOURCE,CPP,SELF]},'cases':ans,'proof':'Both primes split2,3,5; norm-two finite Gram entries map to a finite field. Every tested11-node principal block avoids free distances. Nonzero modular determinant implies exact rank>=11. Integer products/updates are bounded by2*p*p<2^31. Zero determinants give no exclusion.','scope':'All minimum supports consisting only of originally nonright finite pairs; every selected label independently changes to a different member of {2,3,4,5,6}. All original parallel pairs have arbitrary real entries. No claim about larger supports, changed original right pairs, other labels, or other sources.'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--case');p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.case)
 out=OUT if not a.case else OUT.with_name(OUT.stem+'_'+a.case+'.json')
 if a.verify:assert d==json.loads(out.read_text());print('PASS exact modular replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n')
if __name__=='__main__':main()
