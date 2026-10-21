#!/usr/bin/env python3
"""Reproduce the numerical relaxed union test; compatibility is not imposed."""
import argparse,json
from fractions import Fraction
import numpy as np
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P
from fresh_h9_makarov_joint_20261013_selection import setup as oldsetup
import fresh_h9_makarov_joint_20261013 as volume
def run(objectives):
    s,Q,W,c,E=setup();_,_,_,R,_,_,_=oldsetup()
    old=json.loads((ROOT/f'results/{P}_inventory.json').read_text())['compatible_indices']
    caps=[{'normal':r.tolist()}for r in W[11:]]+[{'normal':R[i].tolist()}for i in old]
    for suffix in['profiles_all','remote_edges','faces','faces_latest','nonright_faces','coupled_frontier']:
        caps+=json.loads((ROOT/f'results/{P}_{suffix}.json').read_text())['caps']
    print('relaxed union entries',len(caps),flush=True)
    volume.OUT=ROOT/f'results/{P}_nonright_union.json'
    out=volume.volume_scout({'caps':caps},objectives)
    exact=json.loads((ROOT/f'results/{P}_coupled_exact.json').read_text())
    f=lambda x:float(Fraction(x[0]))+np.sqrt(5)*float(Fraction(x[1]))
    new=np.array([[f(x)*np.sqrt(f(row['radicand']))for x in row['projective_normal']]for row in exact['new_caps']])
    pairings=new@Q@np.array(out['best_point']);assert max(pairings)<0
    out['two_exact_added_cap_pairings']=pairings.tolist()
    out['scope']+=' The same numerical point is also strict for both exact added PD7 caps; dictionary entries are not deduplicated.'
    (ROOT/f'results/{P}_nonright_union_volume.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--objectives',type=int,default=100);a=p.parse_args();run(a.objectives)
