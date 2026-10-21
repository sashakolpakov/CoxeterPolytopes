#!/usr/bin/env python3
"""Exact mandatory-cap gates for all central-triangle1+2+2 pentagon rotations."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_release12_47_52_pentagon_20261017 as S

G=S.G;P,U,I,D=G.P,G.U,G.I,G.D
ROOT=S.ROOT;TAG='fresh_three_cap_pentagon_20261017'
FRAMES=[[30,7,8],[7,8,65],[64,30,7]]
TARGETS=[4,8,9,10,17,18,27,53,64,65]

def exclusion(value):
    if value.l>0:return dict(kind='strict_obtuse',product=value.enc())
    for m in range(2,65):
        lower,upper=-D.theta(m+1),-D.theta(m)
        if lower.h<value.l and value.h<upper.l:
            return dict(kind='strict_consecutive_angle_gap',product=value.enc(),labels=[m,m+1],
                lower_threshold=lower.enc(),upper_threshold=upper.enc())
    return None

def run():
    control=json.loads(G.CONTROL.read_text())['after64'];sp=G.Space(control)
    current=json.loads(S.H.PARENT.read_text());assert current['deleted_old_walls']==[49,62]
    assert not set(TARGETS)&set(S.DELETE+[49,62])
    records=[];total_low=total_tail=0
    for frame in FRAMES:
        domain=S.domains(frame)
        assert all(r['profile']==0 for r in domain['low']+domain['tails'])
        gram=[[sp.rp(i,j)for j in frame]for i in frame];inv,pivots=P.C53.inverse(gram)
        # Fix one literal seed in each complementary positive line.
        seeds=[]
        for seed in[4,8,30,53,64,65]:
            if seed in frame:continue
            cov=[sp.rp(seed,i)for i in frame];coeff=U.imv(inv,cov)
            kq=U.rnd(sp.nodes[seed]['q']-U.idot(cov,coeff))
            if kq.l>0:seeds.append((seed,cov,coeff,kq))
        assert seeds
        seed,cov,kcoeff,kq=min(seeds,key=lambda v:float((v[3].h-v[3].l)/v[3].l))
        target_data={}
        for target in TARGETS:
            if target in frame:continue
            tc=[sp.rp(i,target)for i in frame]
            kb=kq if target==seed else U.rnd(sp.rp(seed,target)-U.idot(kcoeff,tc))
            target_data[target]=(tc,kb)
        rows=[]
        for recipe in domain['low']+domain['tails']:
            high=recipe['high_positions'];assert len(high)<=1
            if high:
                j=high[0];assert all(m==2 for i,m in enumerate(recipe['labels'])if i!=j)
                c=I(D.theta(6).l,F(2));contacts=[-c if i==j else I()for i in range(3)]
                # Exact whole real interval; no discrete high-label sampling.
                square=U.rnd(c*c*inv[j][j]);remaining=U.rnd(2-square)
            else:
                contacts=[-D.theta(m)for m in recipe['labels']]
                square=U.rnd(U.idot(contacts,U.imv(inv,contacts)));remaining=U.rnd(2-square)
            assert remaining.l>0,(frame,recipe,remaining.enc())
            scale=U.sqrt(remaining/kq);branches=[]
            for sign in(-1,1):
                rejected=[]
                for target,(tc,kb)in target_data.items():
                    if high:
                        value=U.rnd(-c*U.idot(inv[j],tc)+sign*scale*kb)
                    else:
                        value=U.rnd(U.idot(U.imv(inv,contacts),tc)+sign*scale*kb)
                    reason=exclusion(value)
                    if reason is not None:rejected.append(dict(wall=target,**reason))
                assert rejected,(frame,recipe,sign,'no exact rejection')
                branches.append(dict(sign=sign,retained_wall_rejections=rejected))
            rows.append(dict(recipe=recipe,remaining_norm=remaining.enc(),normalization_scale=scale.enc(),
                whole_cosine_interval=c.enc()if high else None,branches=branches))
        total_low+=len(domain['low']);total_tail+=len(domain['tails'])
        records.append(dict(frame=frame,complete_finite_domain=domain,gram=G.encmat(gram),inverse_pivots=pivots,
            kernel_seed=seed,kernel_square=kq.enc(),recipes=rows))
        print('THREE_CAP_MANDATORY_FRAME',frame,len(domain['low']),len(domain['tails']),flush=True)
    assert(total_low,total_tail)==(9,2)
    rotations=[]
    for i in range(5):
        one=[S.CYCLE[i],S.CYCLE[(i+1)%5]]
        left=[S.CYCLE[(i+j)%5]for j in(1,2,3)]
        right=[S.CYCLE[(i+j)%5]for j in(3,4,5)]
        matches=[j for j,fr in enumerate(FRAMES)if fr in[left,right]]
        assert matches
        rotations.append(dict(single_edge_cap=one,two_edge_caps=[left,right],excluded_frame_indices=matches))
    parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text())
    paths=[Path(__file__),Path(S.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,
        S.H.PARENT,S.H.BOUNDARY,ROOT/'results'/f'{S.TAG}.json',
        ROOT/'results/fresh_nonarithmetic_pentagon_two_center_gate_20261017.json']
    return dict(status='PASS_EXACT_ALL_LABEL_THREE_CAP_PENTAGON_CENTRAL_TRIANGLE_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,retained_boundary_cycle=S.CYCLE,white_support=[0,1,2,3,5,6],
        mandatory_frame_gates=records,all_five_rotations=rotations,
        counts=dict(mandatory_frames=3,low_recipes=total_low,low_norm_branches=2*total_low,
                    complete_high_label_families=total_tail,whole_interval_branches=2*total_tail,rotations=5),
        exact_old_deletion_cost=parent['deletion_cost'],
        completeness='Each mandatory cap must make two prescribed old boundary edges finite. Finite H4+A2 extension classification forces the displayed profile0 domains. Labels m>=6 occur only as isolated I2(m), so symbolic6 covers all m>=6. Both real norm branches are rejected over the entire closed c∈[sqrt3,2] interval for each tail.',
        topology='Three new actual facets whose dual vertices form a central triangle, with consecutive boundary-arc sizes1,2,2 on the genuine old pentagon. The5 rotations list every assignment up to renaming the two2-edge caps. Each contains an already impossible mandatory cap; no mutual-label or third-wall solution can repair that failure.',
        scope='No such prescribed central-triangle compact disk closure while all other current63 native walls remain. Not an exclusion of other three-cap distributions, more caps, different boundary incidence, a further native-wall release/deformation, or general finite-volume completion. No new chamber or numerical geometry promotion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
