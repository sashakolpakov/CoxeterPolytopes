#!/usr/bin/env python3
"""Exact all-label obstruction for the joint-release pentagon two-center disk.

Shared frozen physical interval engine and finite Coxeter classification;
new complete four-contact norm and all-right branch checks.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_release12_47_52_pentagon_20261017 as S

G=S.G;P,U,I,D=G.P,G.U,G.I,G.D
ROOT=S.ROOT;TAG='fresh_nonarithmetic_pentagon_two_center_gate_20261017'

def run():
    old=json.loads(G.CONTROL.read_text())['after64'];sp=G.Space(old)
    current=json.loads(S.H.PARENT.read_text())
    assert current['deleted_old_walls']==[49,62]
    needed=[4,7,8,30,53,64,65]
    assert not set(needed)&set(S.DELETE+[49,62])
    records=[];low_count=tail_count=0
    for at in (1,2,3,4):
        frame=[S.CYCLE[(at+j)%5]for j in(1,2,3,4)]
        domain=S.domains(frame)
        gram=[[sp.rp(i,j)for j in frame]for i in frame];inv,pivots=P.C53.inverse(gram)
        low=[];tails=[]
        for recipe in domain['low']:
            assert recipe['profile']==0
            c=[-D.theta(m)for m in recipe['labels']]
            norm=U.rnd(U.idot(c,U.imv(inv,c)));res=U.rnd(norm-2)
            assert res.h<0 or res.l>0,(at,recipe,res.enc())
            low.append(dict(recipe=recipe,norm_square=norm.enc(),residual_from_two=res.enc(),strict_sign=1 if res.l>0 else -1))
        for recipe in domain['tails']:
            assert recipe['profile']==0 and len(recipe['high_positions'])==1
            j=recipe['high_positions'][0]
            assert all(m==2 for i,m in enumerate(recipe['labels'])if i!=j)
            a=inv[j][j];lower=U.rnd(3*a-2);upper=U.rnd(4*a-2)
            proof=dict(recipe=recipe,coefficient_square=a.enc(),value_at_sqrt3=lower.enc(),value_at_two=upper.enc())
            assert a.l>0
            if upper.h<0:proof['reason']='norm_below_two_on_entire_closed_cosine_interval'
            elif lower.l>0:proof['reason']='norm_above_two_on_entire_closed_cosine_interval'
            else:
                # For this unique tail, every c>0 gives an obtuse old8 contact.
                target=8;cov=[sp.rp(i,target)for i in frame]
                coefficient=U.rnd(-U.idot(inv[j],cov))
                assert coefficient.l>0,(at,recipe,coefficient.enc())
                proof.update(reason='strict_obtuse_retained_wall_for_all_positive_cosines',wall=target,product_per_positive_cosine=coefficient.enc())
            tails.append(proof)
        records.append(dict(split=at,frame=frame,gram=G.encmat(gram),inverse_pivots=pivots,
            low_recipes=low,all_high_label_tails=tails))
        low_count+=len(low);tail_count+=len(tails)
        print('FOUR_CONTACT_GATE',at,frame,len(low),len(tails),flush=True)
    assert(low_count,tail_count)==(32,6)
    frame=[30,7,8];domain=S.domains(frame)
    assert len(domain['low'])==1 and not domain['tails']
    assert domain['low'][0]['profile']==0 and domain['low'][0]['labels']==[2,2,2]
    gram=[[sp.rp(i,j)for j in frame]for i in frame];inv,pivots=P.C53.inverse(gram)
    seed=53;cov=[sp.rp(seed,i)for i in frame];coeff=U.imv(inv,cov)
    kq=U.rnd(sp.nodes[seed]['q']-U.idot(cov,coeff));assert kq.l>0
    scale=U.sqrt(2/kq);branches=[]
    for sign in(-1,1):
        blockers=[]
        for target in needed:
            if target in frame:continue
            targetcov=[sp.rp(i,target)for i in frame]
            kernel_product=kq if target==seed else U.rnd(sp.rp(seed,target)-U.idot(coeff,targetcov))
            product=U.rnd(sign*scale*kernel_product)
            if product.l>0:blockers.append(dict(wall=target,product=product.enc()))
        assert blockers,(sign,'all-right branch lacks strict blocker')
        branches.append(dict(sign=sign,strict_obtuse_retained_products=blockers))
    all_right=dict(frame=frame,unique_recipe=domain['low'][0],gram=G.encmat(gram),inverse_pivots=pivots,
        kernel_seed=seed,kernel_square=kq.enc(),normalizing_scale=scale.enc(),branches=branches)
    paths=[Path(__file__),Path(S.__file__),Path(G.__file__),Path(G.E.__file__),Path(G.E.X.__file__),
        G.CONTROL,G.PARENT,S.H.PARENT,S.H.BOUNDARY,ROOT/'results'/f'{S.TAG}.json']
    return dict(status='PASS_EXACT_ALL_LABEL_PENTAGON_AT_MOST_TWO_CENTER_GATE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,retained_boundary_cycle=S.CYCLE,white_support=[0,1,2,3,5,6],
        four_contact_gates=records,all_right_three_contact_gate=all_right,
        counts=dict(four_contact_low_recipes=low_count,complete_high_label_tails=tail_count,all_right_branches=2),
        topology_reduction='For a genuine five-edge boundary disk with one or two new actual facet vertices, no old-only internal triangles or old chords, and simple finite endpoints, the two centers must be joined and each owns a consecutive boundary arc. Arc sizes are1+4 or2+3. The five2+3 splits are excluded above; a four-edge owner or single center contains excluded frame[64,30,7,8]. Thus this entire specified disk class needs at least three new actual facets.',
        high_label_completeness='Finite connected rank>=3 diagrams have no edge m>=6. In each required endpoint corner every such edge is an isolated I2(m). Enumerating symbolic6 with the complete rank<=5 finite catalogue therefore covers all labels m>=6. All six surviving tails reduce to one contact c in[sqrt3,2], and are excluded on that larger closed real interval.',
        scope='Exact obstruction to the stated compact pentagon disk closure while all other current63 native walls are retained. It is not an obstruction to three or more new facets, further deletions or changed retained roots, additional old-side intersections that change the specified disk, or an arbitrary finite-volume completion. No numerical root or boundary score is used for acceptance.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
