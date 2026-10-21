#!/usr/bin/env python3
"""Independent literal star-equality replay; no proposal imports."""
import argparse,json
from itertools import product
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_structured_incidence_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return [i for i in range(m.bit_length())if m>>i&1]
def run():
 primary=f'results/{TAG}.json';source=f'results/{TAG}_twoflip.json'
 d=json.loads((ROOT/primary).read_text());z=json.loads((ROOT/source).read_text())
 for p,h in z['sources'].items():assert sha(p)==h,p
 recorded={(r['parent'],r['carrier']):r for r in z['first_flips']};assert len(recorded)==len(z['first_flips'])
 seen=set();attempts=0;witnesses=0
 for case in d['cases']:
  co=case['coefficients'];n=sum(map(len,co));off=[];k=0
  for c in co:off.append(k);k+=len(c)
  allmask=(1<<n)-1
  carriers={allmask^sum(1<<(o+i)for o,i in zip(off,v))for v in product(*[range(len(c))for c in co])}
  V={v['incident']for v in case['vertex_records']};M=set(case['minimal_nonfaces'])
  for A in carriers:
   omitted={i for i in bits(A)if A^(1<<i)in V};P=sum(1<<i for i in omitted);N=A^P
   if len(omitted)not in(2,3,4)or P not in M:continue
   star={v for v in V if v&N==N}
   if star!={A^(1<<i)for i in omitted}:continue
   key=case['case'],A;seen.add(key);r=recorded[key]
   assert r['old_missing']==P and r['temporary_missing']==N and N.bit_count()>5
   V1=(V-star)|{A^(1<<i)for i in bits(N)}
   # Directly certify that the temporary N is a minimal nonface.
   assert not any(v&N==N for v in V1)
   assert all(any(v&(N^(1<<i))==N^(1<<i)for v in V1)for i in bits(N))
   failures={q['carrier']:q for q in r['failed_second_carriers']}
   targets={B for B in carriers if B!=A and B&N==N};assert set(failures)==targets
   for B in targets:
    attempts+=1;Q=B^N;old={v for v in V1 if v&Q==Q};expected={B^(1<<i)for i in bits(N)}
    assert old!=expected
    f=failures[B];w=f['witness'];assert w in(old^expected) and f['unexpected']==(w in old);witnesses+=1
 assert set(recorded)==seen and z['records']==[]
 assert len(seen)==z['counts']['first_exact_large_face_flips']and attempts==z['counts']['candidate_second_carriers']
 paths=[primary,source,f'scripts/{TAG}_twoflip_independent.py']
 return dict(status='PASS_INDEPENDENT_CONTROLLED_TWO_FLIP_STAR_REPLAY',sources={p:sha(p)for p in paths},
  valid_first_flips=len(seen),second_carriers=attempts,explicit_failure_witnesses=witnesses,
  returning_two_flip_recipes=0,scope='Literal completeness for the saved16 parents and product-source10-facet carriers only. Every temporary large minimal nonface and every failed second star is rechecked. No claim about three or more flips, arbitrary carriers, or realization of new abstract spheres.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();z=run();p=ROOT/f'results/{TAG}_twoflip_independent.json'
 if a.verify:assert z==json.loads(p.read_text())
 else:p.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status'],z['valid_first_flips'],z['second_carriers'])
