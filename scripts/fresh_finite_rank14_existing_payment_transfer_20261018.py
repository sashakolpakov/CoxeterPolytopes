#!/usr/bin/env python3
"""Finite B/D14 reductions with exact arbitrary-support alternatives."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_rank14_long_pair_support_scout_20261018 as R
import fresh_finite_d15_three_singleton_half_transfer_20261018 as L
import fresh_depth6_sparse_ambient60_forest_frontier_20261018 as C
import fresh_depth6_finite7_adaptive_transfer_20261018 as DEEP

ROOT=R.ROOT;T=R.T;A=R.A;P=R.P;Q=R.Q;D=L.D;X=L.X
TAG='fresh_finite_rank14_existing_payment_transfer_20261018'
CHANGES={'A':set(),'B':{(0,13),(1,13),(2,13),(3,13)},'D':{(0,11)}}
DWHOLE={554,616,1240,1276,1284,1342,1543,1623}
DPROPER={624,1158,1555,1618,1641}
EXPECTED={'A':[],'B':[],'D':sorted(DWHOLE|DPROPER)}
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,*L.STEMS,L.TAG,C.TAG,DEEP.TAG)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def support(rec,over):
    patterns=[c for c in over['isolation']['cases'] if (c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])<=14]
    assert len(patterns)==rec['complete_eligible_patterns']==600
    assert [c['index'] for c in patterns]==[r['index'] for r in rec['selected_patterns']]
    total=0;hits=[]
    for c,row in zip(patterns,rec['selected_patterns']):
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==row['maps'];total+=len(maps)
        affected=[m for m in maps if tuple(m['marked_images']) in CHANGES[rec['kind']]]
        if affected:hits.append(dict(index=c['index'],maps=affected))
    assert total==dict(A=552,B=260,D=1104)[rec['kind']]
    assert [r['index'] for r in hits]==EXPECTED[rec['kind']]
    return dict(kind=rec['kind'],complete_patterns=600,complete_maps=total,
                fresh_native_search=True,affected_maps=sum(len(r['maps']) for r in hits),affected=hits)


def basic(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    left=c['left'][1];right=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in left+right)
    actual=Q.inspect_corner(g,[1]*n,left+right,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    return g,white,left,right


def corner_payments(c,g,white,kind):
    rows=[];alternatives=0
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
        payment=F();reason=None;images=[];roles=[]
        if marks==c['preferred_corner']:
            assert aa is bb and aa['type']==kind+'14'
            images=T.mark_pairs(g,aa['nodes'],marks,14,kind)
            assert set(images)<=CHANGES[kind]
            payment=F(1,2) if kind=='D' else F();reason='new_selected_weight'
        elif aa is bb and kind=='D' and aa['type'] in ('D15','D16'):
            rank=int(aa['type'][1:]);images=T.mark_pairs(g,aa['nodes'],marks,rank,'D')
            assert images==([(1,12)] if rank==15 else [(2,13)])
            payment=F(1,2);reason='existing_whole_'+aa['type']+'_half';alternatives+=1
        elif aa is not bb:
            for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if kind=='D' and sm['type']=='A4' and la['type']=='D13':
                    si=T.mark_pairs(g,sm['nodes'],[m],4,'A');li=T.mark_pairs(g,la['nodes'],[lm],13,'D')
                    assert si==[(1,),(2,)] and li==[(10,)] and not payment
                    payment=F(1,2);reason='existing_whole_A4_D13_center_half'
                    roles=[dict(type='A4',images=si),dict(type='D13',images=li)];alternatives+=1
        rows.append(dict(marks=marks,components=blocks,images=images,payment=payment,reason=reason,roles=roles))
    assert alternatives==1 and sum(r['payment'] for r in rows)==1
    return rows


def proper(c,core,maps,target,kind):
    g,white,left,right=basic(c,core);gg,J,du,dv,end=X.constraints(c,core);assert g==gg
    nodes=c['marked_component'];records=[]
    for m in maps:
        image={t:i for i,t in zip(nodes,m['images'])};missing=sorted(set(range(14))-set(image))
        contacts={t:[(old,image[old],target[t][old]) for old in sorted(image) if target[t][old]] for t in missing}
        assert all(label==3 and i in c['possible_contacts'] for entries in contacts.values() for old,i,label in entries)
        bad=next(((t,old,i,label) for t,entries in contacts.items() for old,i,label in entries if J[i][i]>2),None)
        if bad:
            t,old,i,label=bad
            records.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,
                required_contacts=contacts,forced_node=i,projection_lower_bound=J[i][i],
                outcome='arbitrary_contact_projection_exceeds_norm2'))
            continue
        assert kind=='D' and c['index']==1641 and len(g)==16 and len(nodes)==12
        assert missing==[12,13] and len(contacts[12])==len(contacts[13])==1 and contacts[12]==contacts[13]
        old,i,label=contacts[12][0];assert old==11 and J[i][i]==2
        assert end['both_ends_shortcut_bound']<=8
        # Same-end positivity permits no other contact. Both norm-two leaves
        # are entirely in the Lorentz core and have identical projections.
        records.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,
            required_contacts=contacts,forced_node=i,projection_square='2',residual_square='0',
            forced_leaf_product='2',outcome='two_distinct_short_leaves_would_coincide'))
    return dict(index=c['index'],kind=kind,core_key=c['core_key'],gram=core['gram'],
        white=white,method='proper_completion_exclusion',end_projection=end,maps=records)


def whole(c,core,kind):
    g,white,left,right=basic(c,core);n=len(g);assert len(c['marked_component'])==14
    assert kind=='D'
    branches=[]
    if c['method']=='radical':
        assert kind=='D' and Q.inertia(g)==[n-2,1,1]
        z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
        du=Q.graph_dist(g,left);dv=Q.graph_dist(g,right);short=[]
        for i,j in combinations(allowed,2):
            if z[i]*z[j]<0:
                length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
                short.append(dict(nodes=[i,j],length=length))
        assert short
        seal=dict(method='radical_annihilation',radical=z,allowed_contacts=allowed,mixed_sign_shortcuts=short)
    else:
        gg,J,du,dv,end=X.constraints(c,core);assert gg==g
        remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
        assert all(J[i][i]>=2 for i in remaining)
        equal=[i for i in remaining if J[i][i]==2]
        if kind=='D':
            assert c['index'] in (1284,1623) and n==17 and equal==[3,13]
            assert J[3][13]==-2 and all(J[i][j]>0 for ee in c['end_sets'] for i in ee for j in ee)
            for p in equal:
                h=[row+[F(-1 if i==p else 0)] for i,row in enumerate(g)]
                h.append([F(-1 if i==p else 0) for i in range(n)]+[F(2)])
                assert Q.inertia(h)==[n-1,1,1]
                ww=white+[n]
                assert Q.inertia(A.sub(h,ww))==[len(ww),0,0]
                assert all(Q.inertia(A.sub(h,sorted(set(ww)|{b})))==[len(ww)+1,0,0] for b in left+right)
                branches.append(dict(extra_contact=p,gram=enc(h),white=ww,
                    sole_contact='-1',projection_square='2',residual_square='0',
                    corners=corner_payments(c,h,ww,kind),total='1'))
        seal=dict(method='positive_end_projection_with_complete_saturation_cases',
            end_projection=end,remaining_contact_nodes=remaining,
            remaining_inverse_diagonals=[[i,J[i][i]] for i in remaining],
            saturated_contact_positions=equal,saturated_same_position_product='2' if equal else None,
            saturated_opposite_position_product='-2' if equal else None,
            maximum_extra_whites_touching_core=1 if equal else 0,
            orthogonal_extra_whites_cannot_attach_saturated_white=bool(equal))
    corners=corner_payments(c,g,white,kind)
    return dict(index=c['index'],kind=kind,core_key=c['core_key'],gram=core['gram'],white=white,
        method='whole_component_transfer',sealing=seal,corners=corners,total='1',extra_white_branches=branches)


def literal(scout):
    rows=[];hits=Counter()
    for old in scout['literal_interfaces']:
        corners=[]
        for c in old['corners']:
            debit=F();which=[]
            for kind,changes in CHANGES.items():
                images=set(map(tuple,c['whole_rank14_images'].get(kind,[])))
                if images&changes:
                    assert images<=changes and c['baseline_payment']=='1'
                    debit+=F(1,2) if kind=='D' else F(1);which.append(kind);hits[kind]+=1
            value=F(c['baseline_payment'])-debit;assert value>=0
            corners.append(dict(marks=c['marks'],previous_lower=c['baseline_payment'],changed_types=which,
                debit=debit,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        rows.append(dict(source=old['source'],case=old['case'],corners=corners,total_lower_bound=total))
    assert len(rows)==216 and hits=={'D':8}
    return dict(configurations=216,changed_corners=dict(hits),rows=rows,minimum_total=min(F(r['total_lower_bound']) for r in rows))


def rejected_A14_shortcut(over):
    rows=[]
    for i in (1157,1617):
        c=over['isolation']['cases'][i];g,white,*_=basic(c,over['isolation']['cores'][c['core_key']])
        for old in c['all_corners']:
            marks=old['preferred_corner'];bs=D.finite_blocks(g,sorted(set(white)|set(marks)))
            aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
            for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if sm['type']=='A2' and la['type']=='A13':
                    si=T.mark_pairs(g,sm['nodes'],[m],2,'A');li=T.mark_pairs(g,la['nodes'],[lm],13,'A')
                    assert si==[(0,),(1,)] and li==[(1,),(11,)]
                    rows.append(dict(index=i,marks=marks,A2_images=si,A13_images=li,
                        existing_endpoint_coupon_does_not_apply=True))
    assert len(rows)==2
    return dict(proposed_zero_pairs=[[0,12],[1,13]],role_checks=rows,retained_A14_cost='91',
        conclusion='The proposed free transfer fails the exact marked-role check. No A14 weight or coupon is changed.')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[L.TAG]['new_cost']=='201/2'
    d15={tuple(r['pair']):F(r['weight']) for r in docs[L.TAG]['weighted_pairs']}
    d16={tuple(r['pair']):F(r['weight']) for r in docs[L.B.TAG]['D16_weighted_pairs']}
    assert d15[1,12]==d16[2,13]==F(1,2)
    coupon=docs[L.S.A15]['upgraded_coupon'];assert coupon['new_weight']=='1' and coupon['small_marks']==[0,1] and coupon['large_marks']==[0,12]
    broad=docs[L.S.U.E.N.COUPON]
    assert [r['rank'] for r in broad['small_components']]==[1,2,3,4,5]
    assert next(r for r in broad['old_large_components'] if r['type']=='D13')['eligible_marks']==3
    assert docs[L.S.TAG]['upgraded_coupon']['large_marks']==[11,12]
    scout=docs[R.TAG];over=docs[A.OVER];proofs=[];supports=[];rules=[]
    for rec in scout['records']:
        kind=rec['kind'];sup=support(rec,over);supports.append(sup)
        new=F(1,2) if kind=='D' else F();weights={p:(new if p in CHANGES[kind] else F(1)) for p in combinations(range(14),2)}
        assert C.actual_base((kind,14))==728 and C.deletions((kind,14))==0
        assert sum(weights.values())=={'A':F(91),'B':F(87),'D':F(181,2)}[kind]
        assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
        assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in rec['automorphisms'])
        for hit in sup['affected']:
            c=over['isolation']['cases'][hit['index']];core=over['isolation']['cores'][c['core_key']]
            pr=whole(c,core,kind) if len(c['marked_component'])==14 else proper(c,core,hit['maps'],rec['target'],kind)
            pr['affected_maps']=hit['maps'];proofs.append(pr)
        rules.append(dict(kind=kind,rank=14,target=rec['target'],automorphisms=rec['automorphisms'],
            old_cost='91',new_cost=sum(weights.values()),changed_pairs=sorted(CHANGES[kind]),new_weight=new,
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())]))
    later=[]
    for stem in L.INTERFACES:
        for proof in docs[stem].get('selected_proofs',[]):
            for corner in proof.get('corners',[]):
                if not F(corner.get('payment',0)):continue
                for block in corner.get('components',[]):
                    typ=block['type']
                    if typ in ('A14','B14','D14') and set(corner['marks'])<=set(block['nodes']):
                        images=T.mark_pairs(A.matrix(proof['gram']),block['nodes'],corner['marks'],14,typ[0])
                        assert not set(images)&CHANGES[typ[0]],(stem,proof['index'],images)
                        later.append(dict(source=stem,index=proof['index'],type=typ,images=images,payment=corner['payment']))
    cases=over['isolation']['cases'];counts=[sum(0<(c.get('component_order_bound') or 0)<=13 for c in cases),
        sum(c.get('component_order_bound') is None for c in cases),
        sum((c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])>14 for c in cases),600]
    assert len(cases)==2614 and counts==[1378,83,553,600] and sum(counts)==2614
    paths={Path(__file__),Path(R.__file__),Path(L.__file__),Path(C.__file__),Path(DEEP.__file__),Path(X.__file__),
        Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_RANK14_EXISTING_PAYMENT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        rules=rules,new_costs=dict(A14='91',B14='87',D14='181/2'),complete_selected_support=supports,
        selected_proofs=proofs,selected_contexts=len(proofs),whole_contexts=8,proper_contexts=5,
        saturated_extra_white_branches=4,combined_literal_interface=literal(scout),
        inherited_positive_rank14_interfaces=later,complete_high_overlap_partition=counts,
        existing_D15_half=[1,12],existing_D16_half=[2,13],existing_A2_A13_coupon=coupon,
        existing_A4_D13_coupon='any_A4_mark_vs_D13_fork_center_half',
        rejected_A14_shortcut=rejected_A14_shortcut(over),rank14_representative_upper_remains='91',
        current_depth6_theorem=DEEP.TAG,triangle_minimum=2,quadrilateral_minimum=1,
        new_cross_component_coupons=[],all_other_weights_unchanged=True,
        finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Fresh complete native maps for all600 selected patterns per finite rank14 type; thirteen selected D14 contexts, arbitrary-contact radical and Lorentz seals, proper projection and coincident-leaf exclusions, all saturated extra-white branches, and simultaneous216 current literal interfaces. Written universal/current-priority arguments remain explicit proof inputs. The A14 free-transfer shortcut is rejected by exact marked roles. No new ambient theorem or polytope.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],'contexts',out['selected_contexts'],flush=True)
