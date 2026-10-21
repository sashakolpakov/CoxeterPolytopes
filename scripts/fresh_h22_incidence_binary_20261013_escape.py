#!/usr/bin/env python3
"""Exact strict escape for the incidence-profile H22 chamber."""
import argparse,json
from fractions import Fraction as F
import numpy as np
from fresh_h22_incidence_core_20261013 import ROOT,read,sha
from fresh_h22_new_cusp_completion_20261013_lower_escape import run as lower_check
OUT='results/fresh_h22_incidence_binary_20261013_escape.json'
SOURCE='results/fresh_h22_incidence_binary_20261013.json'
def check(nums,den):
    source=read(SOURCE)
    assert all(sha(p)==v for p,v in source['sources'].items())
    inherited=lower_check(nums,den)
    raw=read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input']
    Q=[[-int(x) for x in row] for row in raw['basis_gram']]
    cov=[sum(x*y for x,y in zip(row,nums[:22])) for row in Q]
    core=read('results/fresh_h22_incidence_core_20261013.json')
    def values(rows): return [sum(int(x)*y for x,y in zip(r,cov))+2*abs(nums[22]) for r in rows]
    a=values(core['old_projections']); b=values(source['accepted_old_projections'])
    assert max(a)<0 and max(b)<0
    q=F(sum(x*y for x,y in zip(nums[:22],cov))+nums[22]**2,den**2)
    assert q>0 and q==F(inherited['square'])
    return {'status':'PASS_EXACT_INCIDENCE_PROFILE_CHAMBER_STRICT_ESCAPE','sources':{SOURCE:sha(SOURCE)},
        'facets':source['facets'],'numerator':nums,'denominator':den,'norm':str(q),
        'four_core_max_pair_numerator':max(a),'binary560_max_pair_numerator':max(b),
        'inherited_bank_escape_replay':inherited,
        'scope':'Strict future spacelike witness on all5,398,110 facets. The added560 walls cut the earlier escape, but the new chamber has infinite volume.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.verify:
        old=read(OUT); out=check(old['numerator'],old['denominator']); assert out==old
    else:
        scout=read('data/enumerations/fresh_h22_incidence_binary_20261013_volume_scout.json')
        x=np.array(scout['runs'][0]['point']);h=np.array(read('results/fresh_h22_incidence_core_20261013.json')['center'])
        x*=.999;x[:22]+=.001*h/np.sqrt(14);den=10**6;nums=np.rint(x*den).astype(np.int64).tolist()
        out=check(nums,den);(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['norm'])
