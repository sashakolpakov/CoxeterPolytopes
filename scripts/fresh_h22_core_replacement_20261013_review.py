#!/usr/bin/env python3
"""Separate radical replay of the813-wall replacement's new geometry.

Does not import the construction producer. The original812 complete pair
census is covered by its separately replayed standalone verifier.
"""
import argparse,json
from fractions import Fraction as F
from collections import Counter
import fresh_h22_multiquadratic_20261013_review as R
PRE='fresh_h22_core_replacement_20261013'
SOURCE='results/'+PRE+'.json';OUT='results/'+PRE+'_review.json'

def run(verify=False):
    d=R.read(SOURCE)
    for p,s in d['sources'].items():assert R.sha(p)==s,p
    raw=R.read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input']
    Q=[[-int(x)for x in row]for row in raw['basis_gram']]
    M=[row+[0]for row in Q]+[[0]*22+[2]]
    W=[list(map(int,row))for row in raw['wall_coordinates_as_columns_transposed']]
    core=R.read('results/fresh_h22_incidence_core_20261013.json')
    part=core['incidence_parts'][0];v0=core['old_projections'][1];C=[]
    for row in R.read('data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json')['candidates']:
        if row['twice_height_square']!=4:continue
        z=[row['A']*v0[j]-sum(W[part[i]][j]for i in range(21)if row['mask']>>i&1)for j in range(22)]
        assert all(x%2==0 for x in z);C.append([x//2 for x in z])
    assert len(C)==301
    old=[w+[0]for w in W]+[c+[1]for c in C]+[c+[-1]for c in C]
    assert d['deleted_old_wall']==1
    retained=[i for i in range(812)if i!=1];assert retained==d['retained_base_indices']
    oldr=list(map(R.ratvec,old));oc=[R.mv(M,w)for w in oldr]
    nn=[R.pair(w,c)for w,c in zip(oldr,oc,strict=True)]
    h=R.ratvec(core['center']+[0]);hc=R.mv(M,h)
    u=list(map(F,d['old_projection_numerator_u']));U=F(d['u_square'])
    p=list(map(F,d['projection_p']));delta=d['cusp_ray'];t=F(d['parameter_t'])
    E=d['twenty_retained_incident'];coef=list(map(F,d['projection_coefficients']))
    assert p==[F(core['center'][j])-sum(a*W[i][j]for a,i in zip(coef,E,strict=True))for j in range(22)]
    assert u==[-a+t*b for a,b in zip(p,delta,strict=True)]
    rplus=R.vm(R.ratvec(u+[0]),R.sqrtq(1/U));rplus[-1]=R.sqrtq(F(1,2))
    rminus=R.flip(rplus);new=[rplus,rminus];nc=[R.mv(M,r)for r in new]
    assert all(R.pair(r,c)==R.C(2)for r,c in zip(new,nc,strict=True))
    assert R.pair(new[0],nc[1])=={}
    counts=Counter(right=1)
    for r in new:
        assert R.sign(R.pair(r,hc))<0
        for i in retained:
            a=R.pair(r,oc[i])
            if i in E:assert a=={};counts['right']+=1
            else:
                assert R.sign(a)<0 and R.sign(R.add(R.mul(a,a),R.scale(nn[i],-2)))>0
                counts['dotted']+=1
    assert counts=={'right':41,'dotted':1582}
    assert all(R.sign(R.pair(h,oc[i]))<0 for i in retained)
    v=R.ratvec(list(map(F,d['finite_vertex']['old_coordinates']))+[0]);vc=R.mv(M,v)
    assert R.pair(v,vc)==R.C(-U) and R.sign(R.pair(v,hc))<0
    for i in retained:
        a=R.pair(v,oc[i]);assert (a=={})==(i in E)
        assert R.sign(a)<=0
    assert all(R.pair(v,c)=={}for c in nc)
    GI=[]
    for i in E:
        rr=[]
        for j in E:
            a=R.pair(oldr[i],oc[j]);assert set(a)<={1};rr.append(a.get(1,F(0)))
        GI.append(rr)
    R.pd(GI)
    # Together with the orthogonal new pair, this is a positive rank22 Gram.
    x=R.ratvec(d['strict_escape']['coordinates']);xc=R.mv(M,x)
    assert R.pair(x,xc)==R.C(1586) and R.pair(x,hc)==R.C(-814)
    assert all(R.sign(R.pair(x,oc[i]))<0 for i in retained)
    assert all(R.sign(R.pair(x,c))<0 for c in nc)
    dd=R.ratvec(delta+[0]);cuts=[R.pair(dd,c)for c in nc]
    assert all(R.sign(a)>0 for a in cuts)
    cq=d['conditional_nonquasi_arithmetic'];a=R.pair(rplus,oc[cq['cycle_parent_index']])
    cycle=R.mul(a,a);c0,c1=map(F,cq['cycle_coefficients_in_1_sqrt_radicand'])
    assert cycle==R.add(R.C(c0),R.scale(R.sqrtq(F(cq['radicand'])),c1))
    assert c1 and len(cycle)==2 and 1 in cycle
    ii=cq['rational_indefinite_short_block']
    assert all(nn[i]==R.C(2)for i in ii)
    s=R.ratvec([sum(old[i][j]for i in ii)for j in range(23)])
    assert R.pair(s,R.mv(M,s))==R.C(F(cq['ones_square']))
    assert R.sign(R.pair(s,R.mv(M,s)))<0
    out={'status':'PASS_SEPARATE_RADICAL_REPLAY_AND_MATHEMATICAL_SOURCE_REVIEW',
         'sources':{p:R.sha(p)for p in(SOURCE,'scripts/'+PRE+'.py','scripts/'+PRE+'_review.py','scripts/fresh_h22_multiquadratic_20261013_review.py')},
         'new_pair_counts':dict(counts),'retained_normals':811,'new_normals':2,
         'actual_finite_vertex_rank':22,'vertex_square':str(-U),
         'strict_escape_square':'1586','strict_escape_center_pair':'-814',
         'new_cusp_pairings':list(map(R.enc,cuts)),'irrational_cycle':R.enc(cycle),
         'rational_short_block_negative_value':cq['ones_square'],
         'scope':'Source/theorem review plus fresh standalone producer replay; separate sparse-squarefree-radical replay of all1623 new pairings, bothnorms, full813 center/finite-vertex/strict-escape conditions, PD20+2I2 incident Gram, cusp cuts and nonquasi witness. No second full812-base pair census, no finite-volume claim.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],dict(counts),flush=True);return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
