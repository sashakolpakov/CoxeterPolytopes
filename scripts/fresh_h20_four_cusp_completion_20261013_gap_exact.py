#!/usr/bin/env python3
"""Rational changed-profile cap and all-profile fixed17-block obstruction."""
import argparse,json
from fractions import Fraction as F
from collections import Counter
from fresh_h20_four_cusp_completion_20261013 import (
 state,OUT as SOURCE,ROOT,read,sha,B,E,Z,O,vec,dot,mv,encoded,
 pairing,norm,e,label,strict_point,bind_dependencies)
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013_gap_exact.json'
SHORT=[0,1,2,3,4,5,8,9,13,14,15,20,21,22,24,25,26]
R6=[72,-12,-18,9,9,18,44,46,54,-5,-13,-6,-9,-9,-18,16,20,36,2,-2,-8]
POLAR=[48,-24,-6,6,6,12,26,34,36,-14,-22,-12,-6,-3,-3,16,20,24,-4,-8,-5]
POS={18:F(9,31),49:F(85,62),55:F(73,62),181:F(3,31)}
FREE={0:F(18,31),2:F(-339,31),3:F(-219,62),4:F(-267,62),8:F(24,31),9:F(12,31),14:F(-753,62),15:F(-603,62),20:F(6,31),22:F(-3,31),24:F(-228,31),25:F(-264,31),26:F(-6,31)}
def construct():
 Q,h,hc,roots=state();roots.extend((B,vec(c['normal']))for c in read(SOURCE)['new_caps']);assert len(roots)==224
 r=tuple(E(F(x,6))for x in R6)
 return Q,h,hc,roots,r

def run():
 dependencies=bind_dependencies(SOURCE);Q,h,hc,roots,r=construct();rc=mv(Q,r);assert dot(r,rc)==O and dot(r,hc).sign()<0
 pairs=[]
 for j,(f,v)in enumerate(roots):
  H,p=pairing(B,rc,f,v);pairs.append({'old':j,'pairing':H.encode(p),'kind':label(H,p,e(H,O),norm(f,v,Q))})
 assert all(dot(r,mv(Q,roots[i][1]))==Z for i in SHORT)
 target=vec(read(SOURCE)['strict_escape']['vector']);cut=dot(r,mv(Q,target));assert cut.sign()>0
 roots.append((B,r));cusp_signs=[]
 for row in read(SOURCE)['four_completed_links']:
  z=dot(r,mv(Q,vec(row['null_ray'])));assert z.sign()<0;cusp_signs.append(z.enc())
 p=tuple(E(x)for x in POLAR);pc=mv(Q,p);assert dot(p,pc)==E(63)and dot(p,hc).sign()<0
 signs=[]
 for f,v in roots:
  H,z=pairing(B,pc,f,v);assert H.sign(z)<=0;signs.append(H.sign(z))
 assert all(x>0 for x in POS.values())and set(FREE)<=set(SHORT)
 coeffs=dict(FREE)
 for i,w in POS.items():coeffs[i]=coeffs.get(i,F(0))+w
 assert all(roots[i][0]is B for i in coeffs)
 representation=tuple(sum((E(w)*roots[i][1][j]for i,w in coeffs.items()),Z)for j in range(21));assert representation==p
 escape=None
 for k in range(12,65):
  x=tuple(a+E(F(1,2**k))*b for a,b in zip(p,h,strict=True))
  try:escape=strict_point(x,Q,hc,roots)
  except AssertionError:continue
  escape['center_perturbation_power']=k;break
 assert escape is not None
 return {'status':'PASS_EXACT_225_WALL_RATIONAL_CAP_AND_FIXED17_BLOCK_OBSTRUCTION','source':SOURCE,'source_sha256':sha(SOURCE),'bound_prior_dependencies':dependencies,'wall_count':225,'new_cap':{'wall':224,'normal':encoded(r),'norm':['1','0'],'strict_center_pairing':dot(r,hc).enc(),'all224_pair_checks':pairs,'pair_counts':dict(sorted(Counter(x['kind']for x in pairs).items()))},'previous_escape_cut':cut.enc(),'four_completed_cusp_strict_pairings':cusp_signs,'strict_escape':escape,
 'fixed17_block_obstruction':{'short_indices':SHORT,'feasible_spacelike_point':encoded(p),'point_norm':['63','0'],'point_center_pairing':dot(p,hc).enc(),'all225_pairing_signs':signs,'positive_old_coefficients':[[i,str(w)]for i,w in POS.items()],'free_short_coefficients':[[i,str(w)]for i,w in FREE.items()],
 'argument':'The displayed vector p is a positive combination of retained walls18,49,55,181 modulo the span of17 fixed short walls. Any new normal orthogonal to these17 walls and acute to the four retained walls satisfies(r,p)<=0, regardless of norm or angle labels. Thus p survives every such addition. For any finite full-dimensional completion, perturb p slightly toward its own strict interior point to obtain a strict spacelike escape. Therefore this whole fixed17-block mechanism cannot produce finite volume.'},
 'scope':'One rational cap cuts224escape and preserves four bounded cusp links, but225stillhas infinitevolume. The all-profile obstruction retains the specified oldwalls and requires their acute compatibility; it does not rule out changing the orthogonal block, old-facet replacement, or arbitrary redundant inequalities.'}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');args=a.parse_args();d=run()
 if args.verify:assert d==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['new_cap']['pair_counts']);print('escape',d['strict_escape']['norm'])
