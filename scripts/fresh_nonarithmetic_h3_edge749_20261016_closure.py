#!/usr/bin/env python3
"""Try a third side closing the facet of an H4+H4 edge-end cap."""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import numpy as np
from scipy.linalg import null_space
import fresh_nonarithmetic_h3_a2_profiles_20261016 as A
C=A.C;B=A.B;L=A.L;ROOT=A.ROOT;ETA=A.ETA
BASE='fresh_nonarithmetic_h3_edge749_20261016';TAG=BASE+'_closure'

def run():
    d=json.loads((ROOT/f'results/{BASE}.json').read_text());geo=d['recentered_geometry'];N=np.array(geo['old_projected_normals'])
    inventory=json.loads((ROOT/'results/fresh_h9_persistent_cusp_screen_20261016.json').read_text())['pair_matrix']
    support=[0,1,2,3,5,6,7,30]
    assert [(i,j,inventory[i][j])for i,j in itertools.combinations(support,2)if inventory[i][j]!=2]==[(0,1,5),(1,2,3),(2,3,3),(5,6,3),(6,7,3),(7,30,5)]
    rays=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in geo['face_rays']]
    zero=next(p for p in d['profiles']if p['node']==0)
    candidates=[(i,c)for i,c in enumerate(zero['candidates'])if 0 in c['open_end_coverage']]
    counts=Counter();accepted=[];closest=[]
    for index,cap in candidates:
        r=np.array(cap['vector']);assert abs(r@ETA@N[7])<2e-5 and abs(r@ETA@N[30])<2e-5
        frame=np.array([r,N[7],N[30]]);G=frame@ETA@frame.T;assert min(np.linalg.eigvalsh(G))>0
        Ji=np.linalg.inv(G);k=null_space(frame@ETA)[:,0];qk=k@ETA@k;assert qk<0
        intermediate=B.clip(rays,r,55)
        for profile in d['profiles']:
            q=profile['residual_square'];off=np.array(profile['native_offsets'])
            for labels in itertools.product(range(2,11),repeat=3):
                counts['recipes']+=1
                c=np.array([L.product(labels[0]),L.product(labels[1])-off[7],L.product(labels[2])-off[30]])
                if c[0]**2/G[0,0]+c[1]**2/G[1,1]>=q-1e-7 or c[0]**2/G[0,0]+c[2]**2/G[2,2]>=q-1e-7:
                    counts['nonfinite_triangle_corner']+=1;continue
                p=c@Ji@frame;rho=(q-p@ETA@p)/qk
                if rho<=1e-9:counts['no_real_branch']+=1;continue
                for sign in(-1,1):
                    v=p+sign*math.sqrt(rho)*k;counts['real_branches']+=1
                    if v[0]<=profile['center_lower_bound']+1e-7:continue
                    g=N@ETA@v+off;bad=[j for j,x in enumerate(g)if not C.compatible(x,5e-5)]
                    rec=dict(first_candidate=index,first_recipe=cap['recipe'],A2_profile=[profile['node'],profile['label']],
                        labels_to_r_7_30=labels,sign=sign,vector=v,failed_old_walls=bad)
                    if bad:
                        closest.append((len(bad),rec));continue
                    assert abs(v@ETA@v-q)<5e-5 and abs(v@ETA@r-c[0])<5e-5
                    rr=B.clip(intermediate,v,56);face=[x for x in rr if x['mask']>>55&1]
                    rec.update(first_vector=r,full_H3_counts=dict(Counter(x['kind']for x in rr)),first_cap_face_counts=dict(Counter(x['kind']for x in face)),
                        first_cap_face_rays=[dict(active=[j for j in range(x['mask'].bit_length())if x['mask']>>j&1],point=x['v'],kind=x['kind'])for x in face])
                    accepted.append(rec)
        print('TRIANGLE_CLOSURE_SEED',index,'accepted',len(accepted),flush=True)
    closest.sort(key=lambda z:z[0]);paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json','results/fresh_h9_persistent_cusp_screen_20261016.json']
    return C.clean(dict(status='BOUNDED_NUMERICAL_H4_SQUARED_ENDPOINT_FACET_TRIANGLE_SEARCH',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        finite_edge_support_components=[[0,1,2,3],[30,7,6,5]],first_cap_candidates=len(candidates),counts=dict(counts),accepted=accepted,closest=[r for _,r in closest[:20]],
        scope='Seven sampled first caps on the H4+H4 edge, all seven exact A2 source profiles, and finite labels2..10 to firstcap/old7/old30. Both additional triangle corners must have positive projected Gram; all55 native gates retained. Numerical discovery, not an all-parameter or exact obstruction.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved endpoint triangle search')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['accepted']),flush=True)
