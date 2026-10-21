#!/usr/bin/env python3
"""Exact selected D14 long-leaf transfers with one restored A16 half.

This packet proves its arbitrary-support selected cases and stated interfaces.
The separate whole-face output inventory and changed forest count are not
inferred merely from these finite checks.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import fresh_finite_rank14_existing_payment_transfer_20261018 as H
import fresh_d12_half_face_composition_20261018 as G
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C

ROOT=H.ROOT; A=H.A; T=H.T; Q=H.Q; D=H.D
TAG='fresh_d14_long_leaf_a16_half_selected_transfer_20261018'
CHANGES={(0,12),(0,13)}
RESTORE=(1,14)
EXPECTED=[628,636,673,777,1354,1362,1558,1571,1640,1646,1666,1847,1920,1933,2240,2313,2326]
AW='fresh_finite_d16_alternate_half_zero_transfer_20261018'
AD='fresh_affine_d12_leaf_pair_transfer_20261017'
OLD_B='fresh_order15_affine_b_pair_reduction_20261017'
SCOUT='fresh_finite_d12_long_leaf_support_scout_20261018'
EXTRA=('fresh_finite_a14_near_endpoint_zero_transfer_20261018',
 'fresh_finite_a16_endpoint_sparse_half_transfer_20261018',
 'fresh_finite_a15_endpoint_coupon_zero_20261018',
 'fresh_finite8_a17_d11_coupon_budget_20261017',
 'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018')
STEMS=tuple(dict.fromkeys((*H.STEMS,H.TAG,G.TAG,C.TAG,AW,AD,OLD_B,SCOUT,*G.FINITE,*EXTRA)))
enc=H.enc


def same_images(g,white,marks,rank,kind):
    ids=sorted(set(white)|set(marks))
    if Q.inertia(A.sub(g,ids))!=[len(ids),0,0]:return []
    shared=next((p for p in A.components(g,ids) if set(marks)<=set(p)),None)
    if shared is None or len(shared)!=rank:return []
    return T.mark_pairs(g,shared,marks,rank,kind)


def payments(c,g,white):
    """A conservative subset of the current marked payments plus restoration."""
    rows=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids)
        aa,bb=[next(b for b in blocks if m in b['nodes'])for m in marks]
        val=F();reason='unused';images=[];roles=[]
        if aa is bb:
            typ=aa['type'];images=T.mark_pairs(g,aa['nodes'],marks,len(aa['nodes']),typ[0])
            if typ=='D14'and set(images)<=CHANGES:
                val=F(1,2);reason='new_D14_long_leaf_half'
            elif typ=='A14'and images==[(0,13)]:
                val=F(1);reason='retained_A14_endpoint_unit'
            elif typ=='A16'and set(images)<={(0,13),(2,15)}:
                val=F(1,2);reason='retained_A16_half'
            elif typ=='A16'and images==[RESTORE]:
                val=F(1,2);reason='restored_A16_half'
        else:
            roles=[dict(type=b['type'],images=T.mark_pairs(g,b['nodes'],[m],len(b['nodes']),b['type'][0]))
                   for m,b in zip(marks,(aa,bb))]
            for sm,la in (roles,roles[::-1]):
                if sm['type']=='A1'and la==dict(type='A13',images=[(0,),(12,)]):
                    val=F(1,2);reason='retained_whole_A1_A13_endpoint_half'
                if sm['type']=='A2'and la==dict(type='A13',images=[(0,),(12,)]):
                    val=F(1);reason='retained_whole_A2_A13_endpoint_unit'
                if sm==dict(type='A4',images=[(1,),(2,)])and la==dict(type='D13',images=[(11,),(12,)]):
                    val=F(1,2);reason='retained_whole_A4_D13_short_leaf_half'
        rows.append(dict(marks=marks,components=blocks,images=images,roles=roles,payment=val,reason=reason))
    selected=next(r for r in rows if r['marks']==c['preferred_corner'])
    assert selected['reason']=='new_D14_long_leaf_half'and selected['payment']==F(1,2)
    assert any(r['payment']>=F(1,2)for r in rows if r is not selected)
    assert sum(r['payment']for r in rows)>=1
    return rows


def selected(c,core,maps,target):
    g,white,left,right=H.basic(c,core);n=len(g)
    if len(c['marked_component'])<14:
        assert c['index']==1640 and len(c['marked_component'])==13
        proof=H.proper(c,core,maps,target,'D')
        assert len(proof['maps'])==4 and all(r['projection_lower_bound']==12 for r in proof['maps'])
        return proof
    assert len(c['marked_component'])==14
    branches=[(None,g,white)]
    if c['method']=='radical':
        assert Q.inertia(g)==[n-2,1,1]and len(core['radical'])==1
        z=list(map(F,core['radical'][0]))
        assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i]for i in allowed)
        du,dv=Q.graph_dist(g,left),Q.graph_dist(g,right);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        assert short
        seal=dict(method='radical_annihilation',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short)
    else:
        assert c['index']in(1646,1666,1933,2326)and n==17
        gg,J,du,dv,end=H.X.constraints(c,core);assert gg==g
        rem=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert sorted(J[i][i]for i in rem)==[2,2,6]
        equal=[i for i in rem if J[i][i]==2];assert equal==[3,13]and J[3][13]==-2
        assert all(J[i][j]>0 for ee in c['end_sets']for i in ee for j in ee)
        for p in equal:
            h=[r+[F(-int(i==p))]for i,r in enumerate(g)]
            h.append([F(-int(i==p))for i in range(n)]+[F(2)])
            assert Q.inertia(h)==[n-1,1,1]
            ww=white+[n]
            assert Q.inertia(A.sub(h,ww))==[len(ww),0,0]
            assert all(Q.inertia(A.sub(h,ww+[b]))==[len(ww)+1,0,0]for b in left+right)
            branches.append((p,h,ww))
        seal=dict(method='complete_saturated_lorentz_extensions',end_projection=end,
            remaining_inverse_diagonals=[[i,J[i][i]]for i in rem],saturated_contact_positions=equal,
            same_position_product='2',opposite_position_product='-2',maximum_touching_extra_whites=1,
            core_orthogonal_whites_orthogonal_to_saturated_white=True)
    variants=[]
    for p,h,ww in branches:
        cc=payments(c,h,ww)
        variants.append(dict(extra_contact=p,gram=h,white=ww,corners=cc,total=sum(r['payment']for r in cc)))
    return dict(index=c['index'],core_key=c['core_key'],method='whole_D14_transfer',
        gram=g,white=white,sealing=seal,variants=variants)


def literal(docs):
    old=docs[SCOUT]['literal_interfaces'];rows=[];changed=0
    for saved,(stem,case,g,white,cc)in zip(old,T.literal_inputs(docs),strict=True):
        assert (saved['source'],saved['case'])==(stem,enc(case))and len(saved['corners'])==len(cc)
        corners=[]
        for prev,(marks,*_)in zip(saved['corners'],cc):
            assert prev['marks']==marks
            images=same_images(g,white,marks,14,'D');hit=bool(set(images)&CHANGES)
            if hit:assert set(images)<=CHANGES and F(prev['current_lower'])==1
            val=F(prev['current_lower'])-(F(1,2)if hit else 0);changed+=hit
            assert val>=0
            corners.append(dict(marks=marks,D14_images=images,previous_lower=prev['current_lower'],
                debit=F(1,2)if hit else F(),retained_lower=val))
        total=sum(r['retained_lower']for r in corners);assert total>=1
        rows.append(dict(source=stem,case=case,corners=corners,total=total))
    assert len(rows)==216 and changed==16
    return dict(records=rows,changed_corners=changed,minimum_total=min(r['total']for r in rows),
        credit_for_restored_A16_counted=False)


def later(docs):
    packets=[]
    for stem,keys in G.FINITE.items():
        out=[]
        for key in keys:
            for path,c,gg,white,blocks in G.corners(docs[stem][key],key):
                g=A.matrix(gg);marks=c['marks'];ids=sorted(set(white)|set(marks))
                assert sorted(A.components(g,ids))==sorted(b['nodes']for b in blocks)
                images=same_images(g,white,marks,14,'D')
                assert not set(images)&CHANGES,(stem,path)
                out.append(dict(path=path,marks=marks,D14_images=images,new_debit='0'))
        assert out;packets.append(dict(source=stem,sections=list(keys),corners=out))
    assert len(packets)==22 and sum(len(r['corners'])for r in packets)==564
    return packets


def affine_D12(docs):
    out=[]
    for row in docs[AD]['two_A1_spectators']['rows']:
        g=A.matrix(row['gram']);white=row['white'];corners=[]
        for old in (row['corners'][0],row['corners'][2]):
            marks=old['marks'];ids=sorted(set(white)|set(marks))
            assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
            parts=D.finite_blocks(g,ids);aa,bb=[next(b for b in parts if m in b['nodes'])for m in marks]
            if aa is bb:
                assert aa['type']=='D14'
                images=T.mark_pairs(g,aa['nodes'],marks,14,'D');assert images==[(1,12),(1,13)]
                reason='retained_whole_D14_unit_used_as_half'
            else:
                assert sorted(b['type']for b in(aa,bb))==['A2','D12']
                la=aa if aa['type']=='D12'else bb;lm=next(m for m in marks if m in la['nodes'])
                images=T.mark_pairs(g,la['nodes'],[lm],12,'D');assert images==[(10,),(11,)]
                reason='retained_whole_A2_D12_short_leaf_half'
            corners.append(dict(marks=marks,components=parts,images=images,payment='1/2',reason=reason))
        assert corners[0]['marks']!=corners[1]['marks']
        out.append(dict(marks=row['marks'],attachments=row['leaf_attachments'],corners=corners,total='1'))
    assert len(out)==16
    return out


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    aw={tuple(r['pair']):F(r['weight'])for r in docs[AW]['A16_weighted_pairs']}
    dw={tuple(r['pair']):F(r['weight'])for r in next(r for r in docs[H.TAG]['rules']if r['kind']=='D')['weighted_pairs']}
    assert sum(aw.values())==F(C.actual_base(('A',16)),32)==F(441,4)and aw[RESTORE]==0
    assert sum(dw.values())==F(C.actual_base(('D',14)),32)==F(181,2)
    assert all(dw[p]==1 for p in CHANGES)
    aw[RESTORE]=F(1,2)
    for p in CHANGES:dw[p]=F(1,2)
    assert sum(aw.values())==F(443,4)and sum(dw.values())==F(179,2)
    rec=next(r for r in docs[H.R.TAG]['records']if r['kind']=='D');over=docs[A.OVER]
    patterns=[c for c in over['isolation']['cases']if (c.get('component_order_bound')or 0)>=14 and len(c['marked_component'])<=14]
    assert len(patterns)==len(rec['selected_patterns'])==600
    assert all({tuple(sorted(mp[i]for i in p))for p in CHANGES}==CHANGES for mp in rec['automorphisms'])
    assert all(dw[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    proofs=[];domain=[];total=0
    for c,saved in zip(patterns,rec['selected_patterns']):
        assert c['index']==saved['index']
        core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==saved['maps'];total+=len(maps)
        hits=[m for m in maps if tuple(m['marked_images'])in CHANGES]
        domain.append(dict(index=c['index'],maps=maps,affected_maps=hits))
        if hits:
            pr=selected(c,core,hits,rec['target']);pr['affected_maps']=hits;proofs.append(pr)
    assert total==1104 and [p['index']for p in proofs]==EXPECTED
    assert sum(len(p.get('variants',[]))for p in proofs)==24
    assert Counter(p['sealing']['method']for p in proofs if 'sealing'in p)==dict(radical_annihilation=12,complete_saturated_lorentz_extensions=4)
    old_b=docs[OLD_B]['affected_cases'];assert old_b==[1284,1623,1641,1646,1666,1933,2326]
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in [('scripts','py'),('results','json'),('results','md')]}
    return enc(dict(status='PASS_D14_LONG_LEAF_A16_HALF_SELECTED_TRANSFER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},baseline=G.TAG,
        D14_changed_pairs=sorted(CHANGES),D14_new_weight='1/2',D14_new_cost='179/2',
        A16_restored_pair=RESTORE,A16_new_weight='1/2',A16_new_cost='443/4',
        D14_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(dw.items())],
        A16_weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(aw.items())],
        selected_patterns=domain,complete_patterns=600,complete_maps=total,
        selected_proofs=proofs,affected_contexts=17,affected_maps=sum(len(p['affected_maps'])for p in proofs),
        whole_variants=24,proper_projection_exclusions=4,
        literal_interfaces=literal(docs),later_finite_interfaces=later(docs),
        affine_D12_two_A1_interfaces=affine_D12(docs),
        old_B15_dependency_partition=dict(new_selected=[1646,1666,1933,2326],
            unchanged_current_D14_center_half=[1284,1623],unchanged_A12_A13=[1641]),
        no_new_cross_coupon=True,ideal_weights_changed=False,whole_face_composition_proved_here=False,
        finite_forest_count_performed=False,new_dimension_bound=None,
        formal_A16_three_D14_value='1517/4',formal_witness_excess_over377='9/4',
        scope='Complete selected native domain and arbitrary-support17-case replacements,216 literal and564 later finite interfaces,16 affineD12twoA1 repairs. Universal output inventory and changed complete forest certificate remain separate obligations; no new ambient theorem or geometric example.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_contexts'],out['whole_variants'],flush=True)
