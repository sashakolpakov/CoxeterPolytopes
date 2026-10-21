#!/usr/bin/env python3
"""Serialize all2552 rank57 low8 profiles for independent Lean checking.

The saved profiles are reordered into the kernel enumerator's canonical
low-rank-first order; no polynomial, cost or bound is recomputed here.
Numerical proofs are split into20 modules to retain successful kernel work.
"""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json'
EXPECTED='7a01260324b2dff2c81bece2f953f9602e1c38537a389f2da2caba439d5e7a43'
DEST=ROOT/'formal/CoxeterProofs/Generated'
def ns(xs):return '['+','.join(map(str,xs))+']'
def frac(x):
 q=Fraction(x);return f'⟨{q.numerator},{q.denominator}⟩'
def high_parts(left,least=8):
 if left==0:yield ()
 for r in range(least,left+1):
  for tail in high_parts(left-r,r):yield (r,)+tail
def render():
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
 record=json.loads(SOURCE.read_text())
 for p,h in record['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 x=record['ideal'];by_rank={tuple(p['ranks']):p for p in x['rows']}
 domain=[(low,)+rs if low else rs for low in range(8)for rs in high_parts(57-low)]
 assert len(domain)==len(by_rank)==len(x['rows'])==2552 and set(domain)==set(by_rank)
 rows=[by_rank[rs]for rs in domain]
 env=['⟨0,1⟩']+[frac(x['current_envelopes'][str(r)])for r in range(1,58)]
 common=['set_option maxRecDepth 100000','set_option maxHeartbeats 0',
   'namespace CoxeterProofs.Generated.Ideal58','open CoxeterProofs.IdealProfile']
 outputs={}
 outputs['Ideal58Data.lean']='\n'.join(['import CoxeterProofs.IdealProfile58Arithmetic',*common,
   'def envelopes : List FractionData := ['+','.join(env)+']',
   'end CoxeterProofs.Generated.Ideal58',''])
 names=[];modules=[]
 for start in range(0,len(rows),128):
  module=f'Ideal58Chunk{start//128:02}';modules.append(module)
  lines=['import CoxeterProofs.Generated.Ideal58Data',*common]
  for k in range(start,min(start+128,len(rows)),32):
   name=f'chunk{k//32:02}';names.append(name);out=[]
   for p in rows[k:k+32]:
    out.append('⟨'+ns(p['ranks'])+','+ns(p['h_polynomial'])+f', {p["link_vertices"]}, {p["G"]}, '+frac(p['normalized_cost'])+'⟩')
   lines += [f'def {name} : List Profile := [\n'+',\n'.join(out)+']',
     f'theorem {name}_checked : {name}.all (CoxeterProofs.Ideal58.check envelopes) = true := by decide']
  lines += ['end CoxeterProofs.Generated.Ideal58',''];outputs[module+'.lean']='\n'.join(lines)
 lines=[*[f'import CoxeterProofs.Generated.{m}'for m in modules],*common,
  'def profiles : List Profile := '+' ++ '.join(names),
  'theorem all_checked : profiles.all (CoxeterProofs.Ideal58.check envelopes) = true := by',
  '  simp only [profiles, List.all_append, '+', '.join(n+'_checked'for n in names)+', Bool.and_self]',
  'theorem every_profile_checked (p : Profile) (member : p ∈ profiles) :',
  '    CoxeterProofs.Ideal58.check envelopes p = true := (List.all_eq_true.mp all_checked) p member',
  'theorem every_profile_gap (p : Profile) (member : p ∈ profiles) :',
  '    72*p.cost.den+57*p.vertices*p.cost.den+2*p.vertices*p.cost.num ≤ 5*p.indexTotal*p.cost.den :=',
  '  CoxeterProofs.Ideal58.absolute_gap_of_check envelopes p (every_profile_checked p member)',
  'theorem domain_complete : profiles.map Profile.ranks = CoxeterProofs.Ideal58.terminalDomain := by decide',
  'theorem profile_count : profiles.length = 2552 := by decide',
  'end CoxeterProofs.Generated.Ideal58','']
 outputs['Ideal58.lean']='\n'.join(lines)
 return outputs
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();files=render()
 for name,s in files.items():
  if a.write:(DEST/name).write_text(s)
  else:assert (DEST/name).read_text()==s,name
 print('PASS exact export: all2552 rank57 low8 profiles, canonical order,20 numerical modules; geometry separate')
