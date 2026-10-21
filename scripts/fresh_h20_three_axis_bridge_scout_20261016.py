#!/usr/bin/env python3
"""Bounded eight-pattern bridge diagnostic; no new exact wall is promoted."""
import argparse,hashlib,itertools,json
from pathlib import Path
import numpy as np
from fresh_h20_release144_cap_20261016_residual_scout import geometry156,gates,constraint,census,directions,radial,score,recipe,close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_three_axis_bridge_scout_20261016'
def run():
    G=geometry156();r156=recipe(G,3);ar=constraint(G,r156)
    G.update(N=np.r_[G['N'],[r156]],norms=np.r_[G['norms'],2],A=np.r_[G['A'],[ar]],b=np.r_[G['b'],1],labels=G['labels']+[156])
    A,b,ids=G['A'],G['b'],G['labels'];Z=directions();rbase=radial(A,b,Z)
    boosts=G['boosts'];cs=np.sqrt(1+boosts*boosts);beta=boosts/cs
    E6=np.array([1,-.01,-.98,-.28,-.05]);records=[];accepted=[]
    for ls in itertools.product([2,3],repeat=3):
        gamma=np.array([0 if m==2 else .5 for m in ls])/cs
        aa=beta@beta-1;bb=-2*beta@gamma;cc=gamma@gamma-1;assert aa>0
        s=(-bb+np.sqrt(bb*bb-4*aa*cc))/(2*aa);coeff=s*beta-gamma
        r=np.sqrt(2)*(s*G['e0']-coeff@G['ee'][:3]);bad=gates(G,r)
        assert abs(r@G['g']@r-2)<1e-8
        for j,m in zip([154,153,155],ls):assert abs(r@G['g']@G['N'][j]+(0 if m==2 else 1))<1e-8
        row=dict(labels_to154_153_155=list(ls),s=float(s),axis_coefficients=coeff.tolist(),chart=r.tolist(),failed_old_gates=bad,E6_pair=float(r@G['g']@E6),center_pair=float(r@G['g']@G['H']))
        if not bad:
            assert row['center_pair']<0 and row['E6_pair']<0
            a=constraint(G,r);res=census(np.r_[A,[a]],np.r_[b,1],ids+[157]);rr=np.minimum(rbase,1/np.maximum(Z@a,1e-30))
            row.update(residual=res,score=score(rr),new_facet_vertices=[v for v in res['vertices']if 157 in v['active']])
            accepted.append((ls,r,a))
        records.append(row)
    assert len(accepted)==6
    pairs=[]
    for (li,ri,ai),(lj,rj,aj)in itertools.combinations(accepted,2):
        pair=ri@G['g']@rj;co=-pair/2
        kind='obtuse'if pair>1e-7 else 'dotted'if co>1+1e-7 else 'parallel'if abs(co-1)<1e-7 else 'right'if abs(co)<1e-7 else 'finite_or_gap'
        data=dict(left=list(li),right=list(lj),pair=float(pair),normalized_pair=float(pair/2),relation=kind)
        if kind=='finite_or_gap':data['label']=float(np.pi/np.arccos(co))
        pairs.append(data)
    target=[next(t for t in accepted if t[0]==ls)for ls in [(2,2,3),(3,2,2)]]
    (li,ri,ai),(lj,rj,aj)=target
    identity=1-1/(4*cs[0]**2);assert identity>0 and abs(ri@G['g']@rj/2-identity)<1e-9
    unionres=census(np.r_[A,[ai,aj]],np.r_[b,1,1],ids+[157,158])
    unionrad=np.minimum(rbase,np.minimum(1/np.maximum(Z@ai,1e-30),1/np.maximum(Z@aj,1e-30)))
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_release144_cap_20261016_residual_scout.py','results/fresh_h20_release144_cap_20261016_residual_scout.json','results/fresh_h20_release144_cap_20261016.json']
    return dict(status='PASS_BOUNDED_EIGHT_THREE_AXIS_BRIDGE_PATTERNS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},base157_score=score(rbase),patterns=records,accepted_patterns=6,pair_relations=pairs,strongest_pair=dict(patterns=[[2,2,3],[3,2,2]],exact_normalized_pair_formula='1-1/(4*c_star^2)>0',normalized_pair=float(identity),compatible=False,geometric_union_only=dict(residual=unionres,score=score(unionrad))),scope='Eight prescribed right/pi3 patterns only; six numerical all157-pair passes. Numerical residual Qhull/1e-6 sphere-band and common2^20-direction volume scores, no exact new-cap or proper-cusp certificate. The strongest two alternatives are mutually obtuse; their geometric union is not a coherent Coxeter extension. All branches remain right144 and miss E6.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['accepted_patterns'],d['strongest_pair']['normalized_pair'],flush=True)
