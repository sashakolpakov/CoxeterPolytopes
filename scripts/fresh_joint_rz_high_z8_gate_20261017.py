#!/usr/bin/env python3
"""Continuous high-Z8 necessary-contact gate on the complete R circle domain.

Every unclosed interval is exported explicitly.  No numerical sampling is
used as interval coverage; the frozen numeric prefix selects exact signs.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_old8_triangle_double_flip_20261017 as SC
import fresh_old30_connector_native_gate_20261017 as E
Q,P,U,I,G,ROOT=E.Q,E.P,E.U,E.I,E.G,E.ROOT
TAG='fresh_joint_rz_high_z8_gate_20261017'

def run():
    cp=ROOT/'results/fresh_r_circle_admissibility_20261017.json'
    pp=ROOT/'results'/f'{SC.TAG}.json'
    circle=json.loads(cp.read_text());parent=json.loads(pp.read_text())
    assert circle['exact_circle']['necessary_labels']==list(range(2,16))
    a,_,c,x,_,_=map(np.array,parent['prefixes'][0]['normals'])
    _,_,ns,_,_,_=SC.H.setup();base=Q.Engine()
    counts=Counter();records=[];unresolved=[];limit=512;maxdepth=18
    domain=I(Q.theta(6).l,F(2))
    rdata={r['label']:r for r in circle['numerical_R_candidates']}
    for mr in circle['exact_circle']['necessary_labels']:
        g=E.Gram(base,ns)
        g.add('C',[7,30,10],0,[2,2,9],c)
        g.add('A',[7,8,'C'],0,[2,2,3],a)
        g.add('R',[8,44,63],0,[2,mr,2],np.array(rdata[mr]['normal']))
        # Bind the exact chosen branch to the unique permitted circle
        # sheet, independently of the numerical sign-selection match.
        assert g.pair('R',10).h<0
        g.add('X',[8,'A','C'],0,[2,3,2],x)
        assert all(F(r['remaining_square'][0])>0 for r in g.recipes)
        cases=[]
        for mR in [3,4,5]:
            local=Counter()
            def attempt(box):
                rows,info=g.norm_frame([8,'R','X'],0,[-box,-Q.theta(mR),I()])
                if not rows:return dict(reason='negative_remaining_square',frame_certificate=info)
                branches=[]
                for sign,row in rows:
                    bad=[]
                    for wall in E.TARGETS:
                        value=U.rnd(row[wall]) # zero A2 profile
                        witness=Q.exclusion(value)
                        if witness:bad.append(dict(wall=wall,**witness))
                    if not bad:return None
                    branches.append(dict(sign=sign,rejection=bad[0]))
                return dict(reason='both_native_norm_branches_excluded',frame_certificate=info,branches=branches)
            def cover(box,depth=0):
                local['boxes']+=1;counts['boxes']+=1
                local['maximum_depth']=max(local['maximum_depth'],depth)
                proof=attempt(box)
                if proof is not None:
                    local['certified_leaves']+=1;counts['certified_leaves']+=1
                    return dict(box=box.enc(),certificate=proof)
                if depth>=maxdepth or local['boxes']>=limit:
                    local['unresolved_leaves']+=1;counts['unresolved_leaves']+=1
                    record=dict(R44=mr,RZ=mR,box=box.enc(),depth=depth,
                        reason='depth_or_family_box_budget_reached_without_certificate')
                    unresolved.append(record)
                    return dict(box=box.enc(),unresolved=True)
                mid=(box.l+box.h)/2
                return dict(split=str(mid),children=[cover(I(box.l,mid),depth+1),cover(I(mid,box.h),depth+1)])
            tree=cover(domain)
            cases.append(dict(RZ=mR,Z_profile=0,ZX=2,Z8_domain=domain.enc(),counts=dict(local),tree=tree))
            print('HIGH_Z8',mr,mR,dict(local),flush=True)
        records.append(dict(R44=mr,exact_R10=g.pair('R',10).enc(),exact_prefix_recipes=g.recipes,cases=cases))
    assert len(records)==14 and sum(len(r['cases'])for r in records)==42
    counts['families']=42;counts['R_labels']=14
    counts['maximum_depth']=max(c['counts']['maximum_depth']for r in records for c in r['cases'])
    paths=[Path(__file__),Path(SC.__file__),Path(E.__file__),Path(Q.__file__),Path(G.__file__),Path(G.E.X.__file__),
        pp,cp,ROOT/'scripts/fresh_r_circle_admissibility_20261017.py',G.CONTROL,G.PARENT,SC.H.PARENT,SC.H.BOUNDARY]
    return dict(status='PASS_ALL_42_HIGH_Z8_FAMILIES_EXCLUDED'if not unresolved else'PARTIAL_HIGH_Z8_INTERVAL_GATE_WITH_EXPLICIT_UNRESOLVED_BOXES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},counts=dict(counts),records=records,unresolved=unresolved,
        domain='Fixed C10=9,AC=3,X8/XA/XC=2/3/2, all profile0; R rightS6+8+63, R44 in complete necessary2..15. High Z8 forces profile0 and ZX2 in the finite8/X/Z corner. Y right8/R/Z requires that span Lorentzian, leaving RZ3/4/5 in the low-RZ domain. For RZ3 the physical domain has theta(Z8)>sqrt3; for4/5 it includes theta(Z8)=sqrt3. All42 cases are relaxed to the closed interval[sqrt3,2], including the null and infinite-label endpoints.',
        completeness='Every node either has outward-rational necessary native Coxeter exclusions on both norm branches, a negative norm remainder after a positive complementary-line certificate, or two children covering its entire closed interval. Budget/depth failures remain explicit unresolved leaves. The finite-angle witness search through64 is only a rejection witness search; no candidate-label cutoff is imposed on Z8.',
        scope='Exact necessary-condition gate for42 specified high-Z8 families only when every leaf is certified. It does not certify the parent numerical chamber or its full native gates, vary A/C/X, address high RZ or different incidences, or claim finite volume. Shared physical intervals and finite-corner classification are inherited. No unresolved box is accepted or declared impossible.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
