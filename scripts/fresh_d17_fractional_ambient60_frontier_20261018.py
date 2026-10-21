#!/usr/bin/env python3
"""Count three disjoint D5/D12 upgrades and finite/ideal D17 reductions.

Fresh finite tables include the proved deep7 deletions in their upper rule.
Replay checks complete Bellman upper/attainment certificates. Ideal scalar
rows are reconstructed from integer convolutions in either mode.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_rank15_order15_ambient60_frontier_20261018 as Z
import fresh_finite_d17_d5_d12_unit_transfer_20261018 as L
import fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018 as I
import fresh_order17_d5_d12_center_quarter_transfer_20261018 as H

ROOT=Z.ROOT;N=60
TAG='fresh_d17_fractional_ambient60_frontier_20261018'
STEMS=(Z.TAG,L.TAG,I.TAG,H.TAG,Z.DEPTH,Z.Z.COMPLETE)


def deep(t):
    k,r=t
    return max(0,r-{'A':27,'B':28,'D':29}[k]) if k in ('A','B','D') else 0


# All finite quantities here use quarter units. The predecessor uses halves.
def base(t):return 2*Z.base(t)-4*deep(t)-(6 if t==('D',17) else 0)
def actual_cross(a,b):return 2*Z.actual_cross(a,b)+(17 if {a,b}=={('D',5),('D',12)} else 0)
def upper_cross(a,b):return 2*Z.upper_cross(a,b)+(17 if {a,b}=={('A',5),('D',12)} else 0)
def item(r,key):return 2*Z.item(r,key)-4*max(0,r-29)-(6 if r==17 else 0)+(17*(key[1]-key[3]) if r==5 else 0)


def table(key):
    values=[0]+[-10**9]*N
    for r in Z.Z.ITEMS:
        w=item(r,key)
        for n in range(r,N+1):values[n]=max(values[n],values[n-r]+w)
    last=[0]+[next(r for r in Z.Z.ITEMS if r<=n and values[n]==values[n-r]+item(r,key)) for n in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    count=0
    for n in range(1,61):
        for r in Z.Z.ITEMS:
            if r<=n:assert v[n]>=v[n-r]+item(r,key);count+=1
        r=last[n];assert r in Z.Z.ITEMS and r<=n and v[n]==v[n-r]+item(r,key)
    return count


def finite(saved=None):
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=N]
    for a in types:
        assert base(rep(a))>=base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('D',5),('D',12))==47
    assert upper_cross(('A',5),('D',12))==47
    domain=Z.Z.large_domain();keys=sorted({key for cs,used,key in domain})
    assert len(domain)==1044 and len(keys)==792
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {key:table(key) for key in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tables[key],key) for key in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in domain:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used
        assert key[1]-key[3]==cs[3] # c12; key stores c11+c12 and c11 separately.
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=60-used;t=tables[key];value=sum(base(x) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==60
        upper=sum(base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        actual=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            quadrupled_upper=value,quadrupled_actual=actual,deep7_deletions_included=sum(deep(t) for t in forest)))
    assert max(r['quadrupled_upper'] for r in rows)==1644
    witness=[('D',16),('D',16),('D',28)]
    assert sum(base(t) for t in witness)+sum(actual_cross(a,b) for i,a in enumerate(witness) for b in witness[i+1:])==1644
    return dict(exact_maximum='411',previous_maximum='413',curvature='405',finite_deficit='6',
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=sum(r['quadrupled_upper']>1620 for r in rows),
        actual_attainer=['D16','D16','D28'],D17_cost='235/2',D5_D12_actual_cost='47/4',finite_cost_unit='1/4',
        deep7_deletions_counted_in_upper=True,new_tables_generated=True,replay_optimizes=False)


def literal_interface(finite_doc):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in L.STEMS}
    previous={(r['source'],L.N.normalize(r['case'])):r for r in finite_doc['literal_rows']}
    ideal=json.loads((ROOT/f'results/{I.TAG}.json').read_text());target=ideal['target']
    changes={tuple(p):F(1,2) for p in ideal['changed_pairs']}
    center=json.loads((ROOT/f'results/{H.TAG}.json').read_text())
    center_pairs=set(map(tuple,center['changed_pairs']))
    assert not center_pairs&set(changes)
    changes.update({p:F(1,4) for p in center_pairs});rows=[];changed=[]
    for stem,index,g,white,cs in L.T.literal_inputs(docs):
        prev=previous[stem,L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));debit=F();inc=F();images=[];coupon=[]
            sig=L.Q.inertia(L.A.sub(g,ids))
            if sig==[len(ids),0,0]:
                bs=L.D.finite_blocks(g,ids);aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
                if aa is not bb:
                    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                        if sm['type']=='D5' and la['type']=='D12':
                            si=L.T.mark_pairs(g,sm['nodes'],[m],5);li=L.T.mark_pairs(g,la['nodes'],[lm],12)
                            if si==[(0,)] and li in ([(9,)],[(10,),(11,)]):
                                inc=F(1,4) if li==[(9,)] else F(1,2)
                                coupon.append(dict(small_images=si,large_images=li))
            elif sig[1]==0:
                shared=next((p for p in L.A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared)==17:
                    graph=L.T.T.labels(g,shared)
                    if graph is not None:
                        maps=L.P.embeddings(graph,target)
                        images=sorted({tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in maps})
                        if set(images)&set(changes):
                            assert set(images)<=set(changes)
                            debits={changes[p] for p in images};assert len(debits)==1;debit=debits.pop()
            value=F(prev['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,ideal_images=images,ideal_debit=debit,coupon_increment=inc,
                coupon_roles=coupon,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['ideal_debit'] or c['coupon_increment'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed_configurations=len(changed),changed=changed,rows=rows)


def ideal(old):
    domain=set(Z.Z.high_parts(59))
    for low in range(1,9):domain.update((low,)+p for p in Z.Z.high_parts(59-low))
    before={tuple(r['ranks']):r for r in old['ideal']['rows']};assert set(before)==domain and len(domain)==2113
    rows=[]
    for ranks in sorted(domain):
        row=Z.Z.scalar(ranks);s15=F(2*ranks.count(14),15)
        row['normalized_cost']=str(F(row['normalized_cost'])-s15)
        row['normalized_margin']=str(F(row['normalized_margin'])+s15)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['order15_normalized_saving']=str(s15);assert row==before[ranks]
        saving=F(5*ranks.count(16),34)
        row['normalized_cost']=str(F(row['normalized_cost'])-saving)
        row['normalized_margin']=str(F(row['normalized_margin'])+saving)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['order17_normalized_saving']=str(saving);rows.append(row)
    summary=Z.Z.summarize(rows)
    assert summary['negative_count']==28 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1117/7680'
    assert all(16 not in row['ranks'] for row in rows if F(row['normalized_margin'])<0)
    conditional=Z.Z.conditional(rows,{11:F(1,2),12:F(1,2),13:F(1,4),15:F(5,4)})
    assert conditional['summary']['negative_count']==conditional['summary']['zero_count']==0
    return dict(complete_profiles=2113,previous_negative_profiles=49,resolved_negative_profiles=21,
        rows=rows,summary=summary,remaining_conditional_target=conditional,low9_concentration_only=True,
        negative_terminal_absolute_margin_is_universal_bound=False)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[L.TAG]['new_cost']=='235/2' and docs[L.TAG]['upgraded_coupon']['component_pair_cost_increase']=='3'
    assert docs[I.TAG]['new_order17_costs']==dict(A=85,B=111,C=108,D=114)
    assert docs[I.TAG]['upgraded_coupon']['component_pair_cost_increase']=='1'
    assert docs[H.TAG]['new_cost']=='227/2'
    assert docs[H.TAG]['upgraded_coupon']['component_pair_cost_increase']=='1/4'
    a,b=docs[L.TAG]['upgraded_coupon'],docs[I.TAG]['upgraded_coupon']
    marks_a={(i,j) for i in a['small_marks'] for j in a['large_marks']}
    marks_b={(i,j) for i in b['small_marks'] for j in b['large_marks']}
    c=docs[H.TAG]['upgraded_coupon'];marks_c={(i,j) for i in c['small_marks'] for j in c['large_marks']}
    assert len(marks_a)==6 and len(marks_b)==2 and len(marks_c)==1
    assert not marks_a&marks_b and not marks_a&marks_c and not marks_b&marks_c
    assert all(docs[s]['triangle_minimum']==2 and docs[s]['quadrilateral_minimum']==1 for s in (L.TAG,I.TAG,H.TAG))
    assert docs[Z.DEPTH]['feature_depth']==10 and docs[Z.Z.COMPLETE]['dimension_bound']==60
    assert docs[Z.DEPTH]['deleted_pair_counts']=={k:'max(0,n-'+str(c)+')' for k,c in dict(A=27,B=28,D=29).items()}
    assert F(5*sum(min(j,59-j) for j in range(60))-59*60,2)==405
    paths={Path(__file__),Path(Z.__file__),Path(Z.Z.__file__),Path(Z.Z.B.__file__),Path(Z.Z.B.B.__file__),Path(L.__file__),Path(I.__file__),Path(H.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    assert sum(2*j-59 for j in range(60) if 2*j>59)==900
    return L.enc(dict(status='PASS_FINITE411_FRACTIONAL_IDEAL17_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        finite=finite(saved['finite'] if saved else None),ideal=ideal(docs[Z.TAG]),
        combined_literal_interface=literal_interface(docs[L.TAG]),
        disjoint_coupon_marked_pairs=dict(finite=sorted(marks_a),ideal_leaf=sorted(marks_b),ideal_center=sorted(marks_c)),
        current_ideal_order17_costs=dict(A='85',B='111',C='108',D='227/2'),
        combined_D5_D12_increment='17/4',combined_D5_D12_cost='47/4',strict_global_constant='4500',
        constant_counted_once=True,new_dimension_bound=None,new_finite_volume_example=False,
        scope='All three disjoint finite coupon increments included. Complete typed domination,1044/792 fresh finite certificates with deep7 upper deletions, all2113 exact ideal convolutions and216 simultaneous literal interfaces. Universal local proofs are frozen inputs; no new census or dimension59 theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text()) if a.verify else None
    out=run(saved)
    if a.verify:assert out==saved
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'negative',out['ideal']['summary']['negative_count'],
          'literal_changed',out['combined_literal_interface']['changed_configurations'],flush=True)
