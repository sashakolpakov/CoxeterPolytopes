#!/usr/bin/env python3
"""Remove D16(1,11) using a whole D5-long/A11-second unit.

Complete saved D16 map validation, arbitrary-contact completion proofs,
two radical seals and all216 current literal interfaces. No forest solve.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_d17_fractional_ambient60_frontier_20261018 as C
import fresh_finite_d16_fork_exchange_20261018 as J

ROOT=C.ROOT;T=J.T;A=J.A;Q=J.Q;D=J.D;P=J.P;X=J.X
TAG='fresh_finite_d16_second_pair_zero_transfer_20261018'
CHANGE={(1,11)}
PROPER=[504,1031,1537,1614];WHOLE=[1188,2527]
STEMS=tuple(dict.fromkeys((C.TAG,J.TAG,J.V.TAG,T.TAG,A.OVER,C.Z.DEPTH)))


def coupon(g,ids,marks):
    bs=D.finite_blocks(g,ids);a,b=[next(b for b in bs if m in b['nodes']) for m in marks]
    if a is b:return []
    rows=[]
    for sm,la,m,lm in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='D5' and la['type']=='A11':
            si=T.mark_pairs(g,sm['nodes'],[m],5,'D');li=T.mark_pairs(g,la['nodes'],[lm],11,'A')
            if si==[(0,)] and li==[(1,),(9,)]:
                rows.append(dict(small_type='whole_D5',small_nodes=sm['nodes'],small_mark=m,small_images=si,
                    large_type='whole_A11',large_nodes=la['nodes'],large_mark=lm,large_images=li,weight='1'))
    assert len(rows)<=1
    return rows


def proper(c,core,maps,target):
    g,inv,du,dv,end=X.constraints(c,core);base=c['marked_component'];proofs=[]
    for mp in maps:
        image=dict(zip(mp['images'],base));missing=sorted(set(range(16))-set(image));entries=[]
        for new in missing:
            for old,node in image.items():
                if target[new][old]:
                    assert target[new][old]==3 and node in c['possible_contacts']
                    entries.append(dict(new_node=new,old_node=old,core_node=node,inverse_diagonal=inv[node][node]))
        assert entries
        row=dict(images=mp['images'],missing=missing,forced_entries=entries)
        if c['index'] in (504,1537):
            assert sorted(e['inverse_diagonal'] for e in entries)==[4,8]
            row.update(method='forced_projection_above_two',lower_bound='8')
        else:
            assert c['index'] in (1031,1614) and missing==[0,14,15]
            leaves=[e for e in entries if e['new_node'] in (14,15)]
            assert len(leaves)==2 and len({e['core_node'] for e in leaves})==1
            k=leaves[0]['core_node'];assert inv[k][k]==F(4,3) and target[14][15]==0
            same=next(e for e in c['end_sets'] if k in e)
            assert all(inv[i][j]>0 for i in same for j in same)
            assert all(min(du[i]+dv[j],du[j]+dv[i])+2<=8
                for i in c['end_sets'][0] for j in c['end_sets'][1])
            # Each forced simple contact lies in this same end. Further
            # contacts there only increase both projections and their product.
            row.update(method='two_orthogonal_leaves_cauchy',forced_contact_node=k,
                same_end=same,projection_lower_bound='4/3',residual_square_upper='2/3',
                residual_cross_magnitude_lower='4/3',residual_determinant_upper='-4/3')
        proofs.append(row)
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],inverse=core['inverse'],
        base_nodes=base,preferred_corner=c['preferred_corner'],end_projection=end,maps=proofs)


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g);assert n==19 and core['gauge']==[1]*n
    assert Q.inertia(g)==[17,1,1] and all(g[i][i]==2 for i in range(n))
    z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]];white=sorted(set(range(n))-set(u+v))
    assert len(white)==15 and Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in u+v)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):assert actual[k]==c[k]
    assert len(c['marked_component'])==16
    assert T.mark_pairs(g,c['marked_component'],c['preferred_corner'],16,'D')==[(1,11)]
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);short=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],dv[i]+du[j])+2;assert length<=8
            short.append(dict(nodes=[i,j],length=length))
    assert short;corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks));assert Q.inertia(A.sub(g,ids))==[17,0,0]
        cs=coupon(g,ids,marks);assert not cs or marks!=c['preferred_corner']
        corners.append(dict(marks=marks,components=D.finite_blocks(g,ids),coupons=cs,payment=len(cs)))
    assert sum(r['payment'] for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=z,white=white,
        selected_pair=[1,11],allowed_contacts=allowed,mixed_sign_shortcuts=short,
        all_extra_whites_orthogonal_unless_short8=True,corners=corners,total=1)


def literal(current):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.L.STEMS}
    previous={(r['source'],C.L.N.normalize(r['case'])):r for r in current['combined_literal_interface']['rows']}
    rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=previous[stem,C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();inc=F();images=[];coupons=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                bs=D.finite_blocks(g,ids);shared=next((b for b in bs if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']=='D16':
                    images=T.mark_pairs(g,shared['nodes'],marks,16,'D')
                    if set(images)&CHANGE:assert set(images)<=CHANGE;debit=F(1)
                coupons=coupon(g,ids,marks);inc=F(len(coupons))
            value=F(old['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,canonical_images=images,debit=debit,coupon_increment=inc,
                coupons=coupons,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['debit'] or c['coupon_increment'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    current=docs[C.TAG];assert current['finite']['exact_maximum']=='411'
    assert current['current_ideal_order17_costs']==dict(A='85',B='111',C='108',D='227/2')
    assert current['finite']['D17_cost']=='235/2'
    old=docs[J.TAG];weights={tuple(r['pair']):F(r['weight']) for r in old['D16_weighted_pairs']}
    assert sum(weights.values())==111 and weights[1,11]==1;weights[1,11]=F()
    assert sum(weights.values())==110
    rec=next(r for r in docs[T.TAG]['records'] if r['rank']==16);target=rec['target'];sup=rec['selected_support'];over=docs[A.OVER]
    assert J.V.validate_maps(sup,target,over)==578 and sup['complete_eligible_patterns']==394
    _,dist,_=P.distance_data(target);assert dist[1][11]==10
    assert all({tuple(sorted(m[i] for i in p)) for p in CHANGE}==CHANGE for m in P.embeddings(target,target))
    assert all(weights[tuple(p)]==1 for p in rec['basic_pairs_retained'])
    assert weights[1,14]==weights[1,15]==F(1,2)
    hits=[];proper_rows=[];whole=[]
    for row in sup['rows']:
        maps=[m for m in row['maps'] if tuple(m['marked_images']) in CHANGE]
        if not maps:continue
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']]
        hits.append(dict(index=c['index'],maps=maps))
        if c['index'] in PROPER:proper_rows.append(proper(c,core,maps,target))
        else:assert c['index'] in WHOLE;whole.append(seal(c,core))
    assert [r['index'] for r in hits]==sorted(PROPER+WHOLE) and sum(len(r['maps']) for r in hits)==10
    # Actual new cross cost2 fits the already counted19/2 representative.
    assert C.actual_cross(('D',5),('A',11))==0
    assert C.upper_cross(('A',5),('D',11))==38 and 0+8<=38
    assert docs[C.Z.DEPTH]['feature_depth']==10
    paths={Path(__file__),Path(C.__file__),Path(J.__file__),Path(J.V.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__),Path(X.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return J.enc(dict(status='PASS_FINITE_D16_SECOND_PAIR_ZERO_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        old_cost='111',new_cost='110',changed_pairs=sorted(CHANGE),old_weight='1',new_weight='0',
        target=target,weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        complete_selected_patterns=394,validated_saved_maps=578,affected_contexts=hits,
        proper_completion_proofs=proper_rows,radical_proofs=whole,
        new_coupon=dict(small='whole D5',small_marks=[0],large='whole A11',large_marks=[1,9],
            weight='1',component_pair_cost_increase='2',old_actual_cost='0',new_actual_cost='2',
            representative_pair=['A5','D11'],unchanged_representative_allowance='19/2'),
        literal_interface=literal(current),all_ideal_weights_unchanged=True,finite_D17_cost_unchanged='235/2',
        retained_rank15_D16_pairs=[[1,14],[1,15]],compatible_finite7_feature_depth=10,
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_bound_proved=False,new_dimension_bound=None,
        scope='Complete saved394/578 D16 map validation; four arbitrary-contact proper exclusions and two whole-core radical seals; new actual cost2 and216 current literal reconstruction. Historical all-label census and universal interfaces are inherited. No finite forest solve, new ideal cost or dimension bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'D16',out['new_cost'],'literal changes',len(out['literal_interface']['changed']),flush=True)
