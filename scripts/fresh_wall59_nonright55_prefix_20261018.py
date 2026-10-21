#!/usr/bin/env python3
"""Bounded new wall67 duty: a nonright finite contact with retained wall55.

The old59/60/61 deletion debit is exact inherited support arithmetic.
New normals, native tests and target-face cones are numerical discovery only.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import numpy as np
import fresh_wall59_bond_replacement_20261018 as R

T=R.T; ROOT=R.ROOT; M=T.M
TAG='fresh_wall59_nonright55_prefix_20261018'
FRAME=[55,30,58]
RELEASE=[59,60,61]


def run():
    old=json.loads((ROOT/f'results/{R.TAG}.json').read_text())
    for p,h in old['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(p,h)
    s=T.Scout()
    s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    s.counts=Counter();s.blockers=Counter();s.near=[];s.improvements=[]
    debit=T.H.deletion(s.primary,s.cur,RELEASE)
    assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    nd={i:v for i,v in s.old.items() if i not in RELEASE}
    ids=[i for i in s.cur['actual_sides'] if i not in RELEASE]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    deletion=s.score(nd,{},ids,'delete59_60_61')
    assert control is not None and deletion is not None
    assert all(s.ps(i,{})==0 for i in FRAME)
    old55_30=T.G.relation(s.dot(nd[55],0,nd[30],0))
    old55_58=T.G.relation(s.dot(nd[55],0,nd[58],0))
    assert old55_30==old55_58==2
    frame=np.array([nd[i] for i in FRAME])
    recipes=[]
    for pi,m55,m30,m58 in product(range(7),(3,4,5),range(2,6),range(2,6)):
        c1=T.P.corner3(0,0,pi,old55_30,m55,m30)
        c2=T.P.corner3(0,0,pi,old55_58,m55,m58)
        if c1 is not None and c2 is not None:
            recipes.append(dict(profile=pi,labels=[m55,m30,m58],corners=[c1,c2]))
    equations=[];accepted=[];hist=Counter();blockers=Counter();unresolved=[]
    for ri,recipe in enumerate(recipes):
        pi=recipe['profile'];ms=recipe['labels']
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]] for i,m in zip(FRAME,ms)]
        roots,info=T.N.norm_roots(frame,rhs,s.profiles[pi]['square'])
        row=dict(index=ri,**recipe,frame=FRAME,norm=info,branches=[])
        equations.append(row)
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(recipe=ri,stage='norm',info=info));continue
        for sign,v in roots:
            err=max(abs(float(v@M@v)-s.profiles[pi]['square']),
                    max(abs(float(v@M@nd[i])-r) for i,r in zip(FRAME,rhs)))
            if err>2e-5:
                unresolved.append(dict(recipe=ri,stage='residual',branch=sign,error=err));continue
            vals,bad,cp=s.native(v,pi,nd,{})
            branch=dict(branch=sign,normal=v.tolist(),norm_residual=err,
                        center_product=cp,native_products=sorted(vals.items()),blockers=bad)
            row['branches'].append(branch)
            blockers[tuple(bad)]+=1
            if cp>=-1e-8:
                hist['wrong_center']+=1;continue
            if bad:
                hist['native_rejection']+=1;continue
            nn={**nd,67:v};types={67:pi}
            assert s.corner([55,30,67],nn,types) is not None
            assert s.corner([55,58,67],nn,types) is not None
            sid=s.score(nn,types,ids+[67],'nonright55_prefix')
            if sid is None:
                hist['cone_incidence_unresolved']+=1;continue
            bd=s.ledgers[sid]
            rays55=[r for r in bd['rays'] if 55 in r['active']]
            whole55=(all(r['kind']=='finite' for r in rays55) and
                {tuple(sorted(r['active'])) for r in rays55}==
                {(7,30,55),(7,55,58),(30,55,67),(55,58,67)})
            accepted.append(dict(index=len(accepted),recipe=ri,profile=pi,labels=ms,
                **branch,score=sid,whole55_quadrilateral=whole55,
                actual_open_edges=[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')],
                exterior_supports=[r['active'] for r in bd['rays'] if r['kind']=='exterior']))
            hist['native_cone_accepted']+=1
        print('NONRIGHT55',ri,pi,ms,'accepted',len(accepted),'unresolved',len(unresolved)+len(s.unresolved),flush=True)
    paths=[Path(__file__),Path(R.__file__),Path(T.__file__),Path(T.H.__file__),
           Path(T.N.__file__),Path(T.P.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY]
    paths += [ROOT/f'{folder}/{R.TAG}.{ext}' for folder,ext in [('results','json'),('results','md')]]
    return dict(status='BOUNDED_NUMERICAL_NONRIGHT55_WALL59_REPLACEMENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        deleted_walls=RELEASE,exact_deletion_debit=debit,control_score=control,deletion_score=deletion,
        frame=FRAME,required_corner_duties=[[55,30,67],[55,58,67]],
        raw_low_domain_size=7*3*4*4,label55_domain=[3,4,5],other_label_domain=[2,3,4,5],
        recipes=recipes,equations=equations,accepted=accepted,
        counts=dict(raw_low_domain_size=336,finite_recipes=len(recipes),
                    real_norm_branches=sum(len(r['branches']) for r in equations),accepted=len(accepted),**hist),
        blocker_histogram=[dict(walls=list(k),count=n) for k,n in blockers.most_common()],
        boundary_ledgers=s.ledgers,unresolved=unresolved+s.unresolved,
        scope='Complete declared336 low profile/contact tuples and both separated numerical norm branches. New nonright55 duty; all60 retained old native products and strict center precede fullS6 cone/independentPD incidence. Inherited finite classifier and exact7/6/1 deletion debit. No high-label, continuous or exact-normal/ambient-volume claim. Any smaller deletion open-end count is not control improvement.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'unresolved',len(out['unresolved']),flush=True)
