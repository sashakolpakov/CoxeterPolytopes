#!/usr/bin/env python3
"""Independent exact checks supporting the retained-H21 extension strategy.

No producer imports.  The cusp list is frozen from the completed corner walk;
the complete eighteen-type target also comes from the cited Niemeier
classification, not from this finite computation.  Object integers avoid
native overflow.  This script does not duplicate the edge-walk proof.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from sympy import Matrix

P=Path(__file__).resolve().parents[1]
TAG='fresh_h22_extension_feasibility_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
INPUT=f'data/enumerations/{TAG}_cusp_input.json'
OUT=P/'results'/f'{TAG}.json'
TYPES=[['A5']*4,['D4']*5,['D5','A7','A7'],['A1','A1','A9','A9'],
       ['A1','A1','D6','D6','D6'],['D7','A11'],['A3','E6','A11'],['E6']*3,
       ['D4','D8','D8'],['D5','A15'],['A1','A1','A1','A17'],
       ['A1','A1','A1','E7','D10'],['D6','E7','E7'],['D8','D12'],
       ['D4','D16'],['D12','E8'],['D4','E8','E8'],['D20']]
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def enc(v):return list(map(int,v))
def run():
 raw=json.loads((P/RAW).read_text())['coxeter_input']
 Q=-np.array([[int(x)for x in r]for r in raw['basis_gram']],dtype=object)
 W=np.array([[int(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']],dtype=object)
 assert Matrix(Q.tolist()).det()==-4
 h=np.sum(W,axis=0);h//=math.gcd(*map(int,h));assert h@Q@h==-14
 assert all(w@Q@h<0 for w in W)
 inp=json.loads((P/INPUT).read_text());rows=[]
 for row in inp['cusps']:
  delta=np.array(row['primitive_ray'],dtype=object);assert delta@Q@delta==0
  H=-h@Q@delta;assert H>0
  pair=W@Q@delta;assert max(pair)==0
  ids=[i for i,x in enumerate(pair)if not x];assert ids==row['active_walls']
  ws=W[ids];nq=[int(w@Q@w)for w in ws];assert set(nq)<={2,4}
  # K(r)=H*r+(r,h)*delta projects to h-perp without denominators.
  ks=[H*w+(w@Q@h)*delta for w in ws];cov=[Q@w for w in ws]
  roots=[tuple(k)for k,q in zip(ks,nq)if q==2];seen=set(roots);cursor=0
  assert roots and len(roots)==len(seen)
  while cursor<len(roots):
   y=np.array(roots[cursor],dtype=object)
   assert y@Q@y==2*H*H and y@Q@h==0 and y@Q@delta==0
   for k,c,q in zip(ks,cov,nq):
    num=2*(y@c);den=H*q;assert num%den==0
    t=tuple(y-(num//den)*k)
    if t not in seen:seen.add(t);roots.append(t)
   cursor+=1;assert len(roots)<=840
  assert all(tuple(-x for x in r)in seen for r in roots)
  R=np.array(roots,dtype=object);GG=R@Q@R.T;unseen=set(range(len(roots)))
  components=[]
  while unseen:
   todo=[unseen.pop()];cc=[]
   while todo:
    i=todo.pop();cc.append(i);near={j for j in unseen if GG[i,j]}
    unseen-=near;todo.extend(near)
   rank=int(Matrix(R[cc,:].tolist()).rank());N=len(cc)
   if N==rank*(rank+1):label=f'A{rank}'
   elif rank>=4 and N==2*rank*(rank-1):label=f'D{rank}'
   else:
    assert(rank,N)in[(6,72),(7,126),(8,240)];label=f'E{rank}'
   components.append({'type':label,'rank':rank,'root_count':N})
  types=sorted(c['type']for c in components)
  matches=[j+1 for j,t in enumerate(TYPES)if sorted(t)==types];assert len(matches)==1
  rows.append({'snapshot_orbit_id':row['id'],'orbit_size':row['orbit_size'],
    'norm_two_quotient_root_types':types,'norm_two_root_count':len(roots),
    'components':sorted(components,key=lambda c:c['type']),
    'elkies_schuett_table_6_1_row':matches[0]})
 matched=[r['elkies_schuett_table_6_1_row']for r in rows]
 assert len(set(matched))==len(matched)
 # Three nonzero discriminant classes of D20: e1, half-all-ones,
 # and half-all-ones minus e1.  Their norms are odd modulo 2.
 classes=[[F(1)]+[F(0)]*19,[F(1,2)]*20,[F(-1,2)]+[F(1,2)]*19]
 q=[sum(x*x for x in v)for v in classes];assert q==[1,5,5]
 # Elementary proof in the report reduces candidates to divisors of 120.
 # Real-cyclotomic Galois exponent <=2 means every unit square is +/-1.
 divisors=[m for m in range(1,121)if 120%m==0]
 labels=[m for m in divisors if m>=2 and all((a*a)%m in(1%m,(-1)%m)
          for a in range(1,m)if math.gcd(a,m)==1)]
 assert labels==[2,3,4,5,6,8,10,12,24]
 # The pi/7 two-cycle x=4cos²(pi/7) has this irreducible cubic.
 coeff=[1,-6,5,-1] # increasing powers for -x³+5x²-6x+1
 assert all(sum(c*F(a)**i for i,c in enumerate(coeff))!=0 for a in[-1,1])
 return {'status':'PASS_INDEPENDENT_CUSP_TYPES_AND_CONDITIONAL_FIELD_GATES',
   'sources':{p:sha(p)for p in[RAW,INPUT,'scripts/'+Path(__file__).name]},
   'snapshot_source':inp['source'],'snapshot_source_sha256':inp['source_sha256'],
   'snapshot_status':inp['snapshot_status'],'cusp_type_checks':rows,
   'published_cusp_types':TYPES,
   'published_types_missing_from_snapshot':[{'row':j+1,'types':t}for j,t in enumerate(TYPES)if j+1 not in matched],
   'd20_nonzero_discriminant_norms':list(map(str,q)),
   'pseudo_arithmetic_nonuniform_finite_labels':labels,
   'pi7_cycle_irreducible_polynomial_in_increasing_order':coeff,
   'scope':'The field/label restriction is conditional on pseudo-arithmeticity. The eighteen cusp types use the cited published Niemeier classification; this script independently classifies only the frozen snapshot and does not close its corner queue.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run()
 if args.verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print('classified',len(out['cusp_type_checks']),'missing',out['published_types_missing_from_snapshot'])
