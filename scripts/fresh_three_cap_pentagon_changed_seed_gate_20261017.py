#!/usr/bin/env python3
"""Exact all-label mandatory-frame gate for the changed pentagon seed."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_three_cap_pentagon_changed_seed_20261017 as S
import fresh_three_cap_pentagon_20261017 as T

G,P,U,I,D,ROOT=S.G,T.P,T.U,T.I,T.D,S.ROOT
TAG='fresh_three_cap_pentagon_changed_seed_gate_20261017'
THREE=[[30,7,8],[7,8,44],[63,30,7]]
TARGETS=[4,8,9,10,17,18,27,30,44,53,63]

def exact_gram(sp,z,frame):
    out=[]
    for i in frame:
        row=[]
        for j in frame:
            v=U.rnd(G.E.X.literal(z,i,j)-G.white(sp.nodes[i]['cp'],sp.nodes[j]['cp']))
            assert G.overlap(v,sp.rp(i,j));row.append(v)
        out.append(row)
    return out

def three_gate(sp,z,frame):
    dom=S.domain(frame);gram=exact_gram(sp,z,frame);inv,pivots=P.C53.inverse(gram)
    seeds=[]
    for seed in TARGETS:
        if seed in frame:continue
        cc=[sp.rp(seed,i)for i in frame];kc=U.imv(inv,cc);kq=U.rnd(sp.nodes[seed]['q']-U.idot(cc,kc))
        if kq.l>0:seeds.append((seed,kc,kq))
    seed,kc,kq=min(seeds,key=lambda v:float((v[2].h-v[2].l)/v[2].l));records=[]
    for rec in dom['low']+dom['tails']:
        assert rec['profile']==0 and len(rec['high_positions'])<=1
        high=rec['high_positions'];c=I(D.theta(6).l,F(2))
        co=[-c if j in high else-D.theta(m)for j,m in enumerate(rec['labels'])]
        if high:
            j=high[0];assert all(m==2 for k,m in enumerate(rec['labels'])if k!=j)
            square=U.rnd(c*c*inv[j][j])
        else:square=U.rnd(U.idot(co,U.imv(inv,co)))
        remainder=U.rnd(2-square);assert remainder.l>0;scale=U.sqrt(remainder/kq);branches=[]
        for sign in[-1,1]:
            bad=[]
            for target in TARGETS:
                if target in frame:continue
                tc=[sp.rp(i,target)for i in frame];kb=kq if target==seed else U.rnd(sp.rp(seed,target)-U.idot(kc,tc))
                value=U.rnd((-c*U.idot(inv[j],tc)if high else U.idot(U.imv(inv,co),tc))+sign*scale*kb)
                e=T.exclusion(value)
                if e:bad.append(dict(wall=target,**e))
            assert bad,(frame,rec,sign)
            branches.append(dict(sign=sign,rejections=bad))
        records.append(dict(recipe=rec,remaining_norm=remainder.enc(),scale=scale.enc(),branches=branches,
                            whole_cosine_interval=c.enc()if high else None))
    return dict(frame=frame,gram=G.encmat(gram),pivots=pivots,kernel_seed=seed,kernel_square=kq.enc(),records=records)

def four_gate(sp,z,frame):
    dom=S.domain(frame);gram=exact_gram(sp,z,frame);inv,pivots=P.C53.inverse(gram);rows=[]
    for rec in dom['low']:
        cp=G.profile(rec['profile']);q=U.rnd(2-G.white(cp,cp))
        co=[U.rnd(-D.theta(m)-G.white(cp,sp.nodes[i]['cp']))for i,m in zip(frame,rec['labels'])]
        value=U.rnd(U.idot(co,U.imv(inv,co))-q)
        assert value.l>0 or value.h<0,(frame,rec,value.enc())
        rows.append(dict(recipe=rec,reason='strict_wrong_norm',norm_residual=value.enc()))
    tails=[]
    for rec in dom['tails']:
        assert rec['profile']==0
        high=rec['high_positions'];base=[I()if j in high else-D.theta(m)for j,m in enumerate(rec['labels'])]
        if len(high)==1:
            j=high[0];a=inv[j][j];b=U.rnd(-2*U.idot(inv[j],base));c=U.rnd(U.idot(base,U.imv(inv,base))-2)
            assert a.l>0 or a.h<0
            disc=U.rnd(b*b-4*a*c);record=dict(recipe=rec,polynomial=[v.enc()for v in[a,b,c]],discriminant=disc.enc())
            if a.l==a.h==F(1,2)and b.l==b.h==0 and c.l==c.h==-2:
                record['reason']='only_positive_norm_root_is_parallel_limit_c2_not_finite_endpoint'
            elif disc.h<0:record['reason']='negative_discriminant'
            else:
                assert disc.l>0;root_disc=U.sqrt(disc);roots=[]
                for sign in[-1,1]:
                    v=U.rnd((-b+sign*root_disc)/(2*a))
                    if v.h<D.theta(6).l:reason=dict(kind='below_complete_high_cosine_interval')
                    elif v.l>2:reason=dict(kind='above_complete_high_cosine_interval')
                    else:
                        reason=T.exclusion(-v)
                        assert reason is not None,(frame,rec,v.enc())
                    roots.append(dict(sign=sign,cosine=v.enc(),rejection=reason))
                record['roots']=roots
            tails.append(record)
        else:
            assert frame==[8,44,63,30]and len(high)==2 and all(m==2 for j,m in enumerate(rec['labels'])if j not in high)
            j,k=high;box=I(D.theta(6).l,F(2))
            if high!=[1,3]:
                value=U.rnd(inv[j][j]*box*box+2*inv[j][k]*box*box+inv[k][k]*box*box-2)
                assert value.h<0
                tails.append(dict(recipe=rec,reason='strict_wrong_norm_on_complete_two_cosine_box',norm_residual=value.enc()))
            else:
                # Root44 is exactly norm2 and right to the other three frame roots.
                assert all(gram[1][i].l==gram[1][i].h==0 for i in[0,2,3])
                assert inv[1][1].l==inv[1][1].h==F(1,2)and inv[1][3].l==inv[1][3].h==0
                a=inv[3][3];assert a.l>0
                c44=I(U.sqrt(4-8*a).l,U.sqrt(4-6*a).h)
                assert D.theta(15).h<c44.l and c44.h<D.theta(18).l
                branches=[]
                for label in[16,17]:
                    v=U.sqrt((2-D.theta(label)*D.theta(label)/2)/a);reason=T.exclusion(-v);assert reason is not None
                    branches.append(dict(label44=label,forced_cosine30=v.enc(),rejection=reason))
                tails.append(dict(recipe=rec,reason='finite_high_label_reduction_then_angle_gaps',coefficient30=a.enc(),cosine44_interval=c44.enc(),
                                  strict_outer_labels=[15,18],only_possible_labels44=[16,17],branches=branches))
    return dict(frame=frame,gram=G.encmat(gram),pivots=pivots,low=rows,tails=tails)

def run():
    old=json.loads(G.CONTROL.read_text())['after64'];sp=G.Space(old);z=G.E.X.geometry()
    assert not set(TARGETS)&set(S.DELETE+[49,62])
    three=[three_gate(sp,z,f)for f in THREE];four=[]
    for at in range(5):
        frame=[S.CYCLE[(at+j)%5]for j in range(4)];row=four_gate(sp,z,frame);four.append(row)
        print('CHANGED_EXACT_FOUR',frame,len(row['low']),len(row['tails']),flush=True)
    rotations=[]
    for at in range(5):
        frames=[[S.CYCLE[(at+j)%5]for j in js]for js in[(1,2,3),(3,4,5)]]
        matches=[j for j,f in enumerate(THREE)if f in frames];assert matches
        rotations.append(dict(one_edge=[S.CYCLE[at],S.CYCLE[(at+1)%5]],two_edge_frames=frames,excluded_frames=matches))
    parent=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text())
    paths=[Path(__file__),Path(S.__file__),Path(T.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,S.H.PARENT,S.H.BOUNDARY,ROOT/'results'/f'{S.TAG}.json']
    return dict(status='PASS_EXACT_ALL_LABEL_CHANGED_PENTAGON_MANDATORY_FRAME_GATE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,exact_old_deletion_cost=parent['exact_old_deletion_cost'],
        mandatory_three_contact_gates=three,all_five_four_contact_gates=four,central_triangle_122_rotations=rotations,
        counts=dict(three_contact_frames=3,three_contact_recipes=sum(len(r['records'])for r in three),four_contact_frames=5,
                    four_contact_low_recipes=sum(len(r['low'])for r in four),four_contact_all_high_families=sum(len(r['tails'])for r in four)),
        topology='Each one/two-center closure has a cap covering at least three consecutive old boundary edges, hence one of the five excluded4-frames. A three-center central triangle with positive arc lengths has1+1+3 or1+2+2. The former has an excluded4-frame; every rotation of the latter has an excluded3-frame. The retained chord30-44 and finite triangle30-44-63 lie outside this exterior disk.',
        scope='Exact all-label mandatory-cap obstruction for these specified disk incidence templates. Not a theorem excluding all three-facet disks, zero-arc centers, additional old intersections, further deletions/deformations, or general finite-volume completion. No new chamber. Shared frozen physical interval arithmetic; source independently reconstructs every displayed exact norm/gap test.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
