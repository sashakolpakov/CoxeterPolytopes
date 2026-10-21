#!/usr/bin/env python3
"""Small numerical selection diagnostic for the saved bounded two-anchor bank.

No new search, exact angle assertion, or global obstruction is claimed.
"""
import argparse,json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_joint_route_20261015_selection'
SOURCE='results/fresh_h20_a1a1_h18_coupling_20261015.json'
BANK='results/fresh_h20_joint_route_20261015_two_anchor.json'
EXACT='results/fresh_h20_joint_route_20261015_exact.json'
def run():
 d=json.loads((ROOT/SOURCE).read_text());bank=json.loads((ROOT/BANK).read_text());exact=json.loads((ROOT/EXACT).read_text())
 dec=lambda x:float(Fraction(x['a']))+np.sqrt(2)*float(Fraction(x['b']))
 dec4=lambda x:sum(float(Fraction(a))*np.sqrt(b)for a,b in zip(x,[1,2,3,6]))
 Q=np.array([[dec(x)for x in row]for row in d['source_form']])
 rows=[];seen=set();original=[]
 for k,r in enumerate(bank['accepted']):
  key=tuple(np.round(np.r_[r['source_projection'],r['support_contacts']],7))
  if key not in seen:seen.add(key);rows.append(r);original.append(k)
 Y=np.array([r['source_projection']for r in rows]);ab=np.array([r['support_contacts']for r in rows])
 cutters=[i for i,r in enumerate(rows)if r['cuts_escape']]
 assert len(rows)==1844 and cutters==[270,271]
 alphabet=2*np.cos(np.pi/np.arange(2,61));diagnostics=[]
 for i in cutters:
  pairs=Y@Q@Y[i]+ab@ab[i]/2
  compatible=[]
  for j,g in enumerate(pairs):
   if i==j:continue
   if g<-2-1e-7:label='dotted'
   elif abs(g+2)<1e-7:label='parallel'
   else:
    k=int(np.argmin(abs(alphabet+g)))
    if abs(alphabet[k]+g)>=1e-7:continue
    label=k+2
   compatible.append(dict(index=j,relation=label,pairing=float(g)))
  diagnostics.append(dict(index=i,original_bank_index=original[i],support_contacts=rows[i]['support_contacts'],compatible=compatible))
 E=np.array([dec4(x)for x in exact['strict_escape']])
 qE=float(E[:19]@Q@E[:19]+2*(E[19:]@E[19:]))
 products=Y@Q@E[:19]-ab@E[19:]
 assert qE>0 and max(products)<0
 selected_pair=float(Y[270]@Q@Y[271]+ab[270]@ab[271]/2)
 assert selected_pair>0
 paths=[f'scripts/{TAG}.py',SOURCE,BANK,EXACT]
 return dict(status='NUMERICAL_BOUNDED_BANK_SELECTION_AND_COMMON_ESCAPE',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  raw_accepted=len(bank['accepted']),deduplicated_count=len(rows),rounding_digits=7,cutters=cutters,cutter_pairing=selected_pair,
  cutter_compatibility=diagnostics,renewed_escape=exact['strict_escape'],escape_norm=qE,
  bank_escape_pairing_min=float(min(products)),bank_escape_pairing_max=float(max(products)),
  maximum_attainer=int(np.argmax(products)),
  scope='Only the saved bounded bank: numerical deduplication and angle fits with tolerance1e-7 and finite labels2..60. The two cutters are numerically incompatible; the same exact-primary escape is numerically strict for every bank proposal even ignoring mutual incompatibility. This is not an exact bank certificate, all-profile theorem, or finite-volume exclusion for further caps.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['deduplicated_count'],d['cutter_pairing'],d['bank_escape_pairing_max'])
