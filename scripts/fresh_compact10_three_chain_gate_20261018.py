#!/usr/bin/env python3
"""Necessary exact gates for a longer compact10 chain: only old44 survives.

All original physical roots come from the frozen exact certificate intervals.
The first cap J has finite8/10/J,8/65/J,63/65/J duties. No later cap is built.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_compact10_cross_bridge_high_interval_20261018 as I
import fresh_moving33_compact10_bridge_20261018 as B

ROOT=I.ROOT;TAG='fresh_compact10_three_chain_gate_20261018'
WHITE=I.WHITE;SIDES=[4,7,8,10,44,63,65];IDS=WHITE+SIDES


def physical(primary):
    metric=[list(map(I.iv,r)) for r in primary['computational_gram']]
    assert len(metric)==10 and all(len(r)==10 for r in metric)
    center=list(map(I.iv,primary['common_center']));normals={}
    for w in primary['facet_witnesses']:
        if w['wall'] not in IDS:continue
        hp=I.iv(primary['center_products'][str(w['wall'])]);assert hp[1]<0
        normals[w['wall']]=[I.div(I.mul(I.iv(2),I.sub(x,y)),hp) for x,y in zip(center,map(I.iv,w['point']))]
    assert set(normals)==set(IDS)
    g=[[I.dot(normals[a],metric,normals[b]) for b in IDS] for a in IDS]
    assert all(g[i][i][0]<=2<=g[i][i][1] for i in range(len(IDS)))
    wi,wp=I.inverse([r[:6] for r in g[:6]]);assert all(p[0]>0 for p in wp)
    offsets={i:g[IDS.index(i)][:6] for i in SIDES}
    expected={8:[0]*6,10:[0,0,0,0,-1,0],44:[0]*6,63:[0]*6,65:[0]*6}
    for i,values in expected.items():
        assert all(q[0]<=x<=q[1] for q,x in zip(offsets[i],values))
    res={a:{b:I.sub(g[IDS.index(a)][IDS.index(b)],I.dot(offsets[a],wi,offsets[b])) for b in SIDES} for a in SIDES}
    # Exact native equalities are inherited from the verified control. These
    # interval containments check their physical interface, not prove equality.
    for a,b in ((8,10),(8,65),(63,65),(8,44),(10,44),(63,44)):
        assert g[IDS.index(a)][IDS.index(b)][0]<=0<=g[IDS.index(a)][IDS.index(b)][1]
    return dict(walls=IDS,physical_gram=g,white_inverse=wi,white_pivots=wp,old_white_covectors=offsets,residual=res)


def frame(data,walls,gate):
    g=[[data['residual'][a][b] for b in walls] for a in walls]
    inverse,pivots=I.inverse(g);assert pivots[0][0]>0 and pivots[1][0]>0 and pivots[2][1]<0
    c=[data['residual'][a][gate] for a in walls]
    v=[I.total(I.mul(inverse[i][j],c[j]) for j in range(3)) for i in range(3)]
    s=I.sub(data['residual'][gate][gate],I.dot(c,inverse,c));assert s[0]>0
    return dict(walls=walls,gate=gate,gram=g,inverse=inverse,pivots=pivots,gate_coefficients=v,gate_residual_square=s)


def products(data,fr,wc,magnitudes,continuous=False):
    wi=data['white_inverse'];square=I.sub(I.iv(2),I.dot(wc,wi,wc))
    offs=[I.dot(wc,wi,data['old_white_covectors'][a]) for a in fr['walls']+[fr['gate']]]
    rhs=[I.sub(I.neg(t),o) for t,o in zip(magnitudes,offs)]
    lin=I.add(I.total(I.mul(q,v) for q,v in zip(rhs,fr['gate_coefficients'])),offs[-1])
    rem=I.sub(square,I.dot(rhs,fr['inverse'],rhs))
    assert rem[1]>=0
    if continuous:
        # A real norm root necessarily has nonnegative remainder. Intersection
        # also covers the zero-remainder endpoint, where the signs coincide.
        rem=(max(F(0),rem[0]),rem[1])
    else:assert rem[0]>0
    rad=I.sqrt(I.mul(rem,fr['gate_residual_square']))
    return dict(white_covector=wc,projected_square=square,offsets=offs,rhs=rhs,
                linear=lin,norm_remainder=rem,radical=rad,branches=[I.sub(lin,rad),I.add(lin,rad)])


def run():
    old=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    primary=json.loads((ROOT/f'results/{I.PRIMARY}.json').read_text())
    assert primary['status']=='PASS_EXACT_JOINT49_62_REPLACEMENT_63_WALLS_199_FINITE'
    data=physical(primary)
    frames={str(g):frame(data,[8,10,65],g) for g in (4,63)}
    c4=[data['residual'][a][4] for a in (8,10,65)]
    c63=[data['residual'][a][63] for a in (8,10,65)]
    gate_cross=I.sub(data['residual'][4][63],I.dot(c4,frames['4']['inverse'],c63))
    assert gate_cross[1]<0
    # The original six-white complement has dimension4; the Lorentzian
    # three-frame leaves a positive line. Negative residual cross product
    # therefore reverses the4/63 radical signs for the same physical normal.
    rt2,rt3,rt5,rt6=[I.sqrt(I.iv(n)) for n in (2,3,5,6)]
    phi=I.div(I.add(I.iv(1),rt5),I.iv(2))
    mag={2:I.iv(0),3:I.iv(1),4:rt2,5:phi}
    profiles=[[I.iv(0)]*6]
    for pos in (4,5):
        for m in (3,4,5):
            w=[I.iv(0)]*6;w[pos]=I.neg(mag[m]);profiles.append(w)
    symbolic=B.T.E.recipes((0,1,0),tuple([2]*9),((0,1),(0,2)))
    low=[(pi,list(ms)) for pi,ms,hi in symbolic if not hi]
    high=[(pi,list(ms),list(hi)) for pi,ms,hi in symbolic if hi]
    expected=[(0,[2,m,6],[2]) for m in range(2,6)]+[(0,[6,2,2],[0])]
    assert high==expected and len(low)==46
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
        if all(B.T.P.corner3((0,1,0)[a],(0,1,0)[b],pi,2,ms[a],ms[b]) is not None for a,b in ((0,1),(0,2)))]
    assert sorted(low)==sorted(direct)
    rows=[];counts={'positive':0,'native_gap':0}
    for pi,labels in low:
        evaluations={g:products(data,fr,profiles[pi],[mag[m] for m in labels]) for g,fr in frames.items()}
        proofs=[]
        for sign in range(2):
            candidates=[]
            for g,ev in evaluations.items():
                si=sign if g=='63' else 1-sign;value=ev['branches'][si]
                if value[0]>0:proof=dict(kind='positive',margin=value[0])
                elif value[0]>-1 and value[1]<0:proof=dict(kind='native_gap',labels=[2,3],margin=min(value[0]+1,-value[1]))
                else:continue
                candidates.append(dict(gate=int(g),gate_sign=si,canonical63_sign=sign,**proof))
            assert candidates,(pi,labels,sign,I.enc(evaluations))
            proofs.append(candidates[0]);counts[candidates[0]['kind']]+=1
        rows.append(dict(index=len(rows),profile=pi,labels=labels,evaluations=evaluations,exclusions=proofs))
    # The isolated high8 factor has profile0 and J10=J65=0. The gate63
    # enclosure retains the full continuous parameter, including its limit2.
    high8=products(data,frames['63'],profiles[0],[(rt3[0],F(2)),I.iv(0),I.iv(0)],True)
    assert all(v[0]>0 for v in high8['branches'])
    # For high65, finite63/65/J forces J63=0: the bond>=6 is an isolated
    # I2 factor. Solve first against8/10/63, leaving J65 as a measured gate.
    fr65=frame(data,[8,10,63],65);high65=[]
    theta12=I.div(I.add(rt2,rt6),I.iv(2))
    # Elementary rational Taylor bounds and 3<pi<22/7 imply these bounds.
    theta11_upper=F(2)-F(3,11)**2+F(3,11)**4/F(12)
    theta13_lower=F(2)-F(22,91)**2
    assert theta11_upper<theta12[0]<theta12[1]<theta13_lower
    for m in range(2,6):
        ev=products(data,fr65,profiles[0],[I.iv(0),mag[m],I.iv(0)])
        minus,plus=ev['branches'];assert plus[0]>0
        if m==2:
            # Both exact normals are +/-old44: all three frame products and
            # all six white products vanish, and their positive complement
            # is one dimensional. The old strict center selects old44.
            proof=dict(kind='old44_only',exact_linear_complement_argument=True)
        else:
            lo,hi=(theta12[1],theta13_lower) if m==3 else (theta11_upper,theta12[0])
            assert -minus[1]>lo and -minus[0]<hi
            proof=dict(kind='native_gap',labels=[12,13] if m==3 else [11,12],
                       left_margin=-minus[1]-lo,right_margin=hi+minus[0])
        high65.append(dict(label10=m,evaluation=ev,minus_branch=proof,plus_positive_margin=plus[0]))
    debit=old['exact_deletion_debit'];assert debit['lost_finite_count']==18 and debit['compact_faces_by_dimension']=={'2':12,'3':1}
    paths={Path(__file__),Path(I.__file__),Path(B.__file__),Path(B.T.E.__file__),Path(B.T.P.__file__)}
    paths.update(ROOT/f'{d}/{s}.{e}' for s in (B.TAG,I.TAG,I.PRIMARY) for d,e in (('scripts','py'),('results','json'),('results','md')))
    assert all(p.exists() for p in paths)
    return I.enc(dict(status='PASS_THREE_CHAIN_FIRST_CAP_EXACT_GATE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=[30,33,44],retained_old_walls=60,exact_deletion_debit=debit,
        proposed_compact10_chain=[7,11,8,'J','R','S'],first_cap_duties=[[8,10,'J'],[8,65,'J'],[63,65,'J']],
        projection=data,frames=frames,gate4_63_residual_cross=gate_cross,low_recipes=rows,low_branch_decisions=counts,
        symbolic_high_families=high,continuous_high8=high8,high65_native_first_frame=fr65,high65=high65,
        theta11_rational_upper=theta11_upper,theta12_interval=theta12,theta13_rational_lower=theta13_lower,
        complete_low_recipes=46,complete_low_signs=92,high_families=5,changed_native_caps=0,old44_control=True,
        new_normals=0,new_cones=0,scope='Conditional on the frozen exact physical control and complete seven-profile spherical classifier: every first cap satisfying the three finite duties and all retained native walls is old44. All46 low recipes fail exact native gates; all five high families are covered by a continuous necessary domain or native-first right constraint. No later S/R family, arbitrary-incidence nonexistence, new chamber or ambient finite volume is claimed.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'low',out['complete_low_signs'],'high',out['high_families'],'changed',out['changed_native_caps'],flush=True)
