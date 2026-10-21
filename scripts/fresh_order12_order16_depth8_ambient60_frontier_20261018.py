#!/usr/bin/env python3
"""Combine order12/order16 ideal improvements and depth8 finite7 zeros.

Generate complete rank60 Bellman certificates; replay all upper and attaining
predecessor checks. Reconstruct every rank59 terminal convolution and all216
current literal payments. This is not a new unrestricted dimension bound.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_order13_sparse_ambient60_frontier_20261018 as R
import fresh_order12_bc_near_terminal_half_transfer_20261018 as I12
import fresh_order16_bd_sparse_unit_transfer_20261018 as I16
import fresh_depth8_finite7_profile_transfer_20261018 as G

ROOT=R.ROOT;Z=R.Z;DIM=60
TAG='fresh_order12_order16_depth8_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((R.TAG,I12.TAG,I16.TAG,G.TAG,*R.STEMS)))
FINAL_MD={R.TAG:'0766d99165d8b5774e52453272bd22f584d17e8947abee2a747f3385858e2116',
 I12.TAG:'08c042d7282dbbd6f43ee2ef4c675dedf1dfbbc0b1a13951250eec344eaacb9f',
 I16.TAG:'c8fb3e418599cc00f6724163522d612f1efe01b5b602a5caca3fa507187aae63',
 G.TAG:'b7b6d23d2f1e705ecb1a62c6af1c0f69fcdf0fdb44b2dec11de9d412b7fd32ff'}


def additional_deletions(t):
    k,n=t
    if k not in 'ABD':return 0
    return max(0,n-{'A':23,'B':24,'D':25}[k])-max(0,n-{'A':27,'B':28,'D':29}[k])


def base(t):return R.base(t)-8*additional_deletions(t)
def actual_cross(a,b):
    return R.actual_cross(a,b)+(8 if {a,b}=={('D',5),('D',11)} else 0)+(4 if {a,b}=={('A',1),('B',11)} else 0)
def upper_cross(a,b):return R.upper_cross(a,b)
def item(r,key):return R.item(r,key)-8*additional_deletions(('A' if r<=8 else 'D',r))


def table(key):
    values=[0]+[-10**9]*DIM
    for r in Z.Z.ITEMS:
        w=item(r,key)
        for n in range(r,DIM+1):values[n]=max(values[n],values[n-r]+w)
    last=[0]+[next(r for r in Z.Z.ITEMS if r<=n and values[n]==values[n-r]+item(r,key)) for n in range(1,DIM+1)]
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


def actual_value(forest):
    return sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def finite(saved=None):
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=DIM]
    for a in types:
        assert base(rep(a))>=base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('D',5),('D',11))==68 and upper_cross(('A',5),('D',11))==76
    assert actual_cross(('A',1),('B',11))==4 and upper_cross(('A',1),('D',11))==20
    assert actual_cross(('A',4),('B',9))==upper_cross(('A',4),('B',9))==8
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
        assert sum(r for k,r in fixed)==used and key[5]==cs[0]
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=DIM-used;t=tables[key];value=sum(base(x) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==DIM
        upper=sum(base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual_value(forest)<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=value,octupled_actual=actual_value(forest),
            additional_depth8_deletions=sum(additional_deletions(t) for t in forest)))
    maximum=max(r['octupled_upper'] for r in rows);assert maximum==3272
    witness=[('A',1)]*4+[('D',5)]*4+[('D',12)]*3
    assert sum(r for k,r in witness)==60 and actual_value(witness)==maximum
    secondary=[('D',17),('D',17),('D',26)];assert actual_value(secondary)==3268
    return dict(exact_maximum='409',previous_maximum='819/2',curvature='405',finite_deficit='4',
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=sum(r['octupled_upper']>3240 for r in rows),
        actual_attainer=[k+str(r) for k,r in witness],actual_attaining_cost='409',
        secondary_actual_forest=['D17','D17','D26'],secondary_actual_cost='817/2',
        actual_coupons=dict(D5_D11='17/2',A1_B11='1/2',A4_B9='1'),
        representative_allowances=dict(A5_D11='19/2',A1_D11='5/2',A4_B9='1'),
        new_actual_coupon_increments=dict(D5_D11='1',A1_B11='1/2'),
        finite_cost_unit='1/8',new_tables_generated=True,replay_optimizes=False,
        formal_forest_attainment_not_geometric_realization=True)


def ideal(docs):
    previous=R.ideal(docs);assert previous==docs[R.TAG]['ideal'];rows=[]
    for old in previous['rows']:
        row=dict(old);s12=F(row['ranks'].count(11),12);s16=F(row['ranks'].count(15),8)
        row['normalized_cost']=str(F(row['normalized_cost'])-s12-s16)
        row['normalized_margin']=str(F(row['normalized_margin'])+s12+s16)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['new_order12_normalized_saving']=str(s12);row['new_order16_normalized_saving']=str(s16);rows.append(row)
    summary=Z.Z.summarize(rows)
    assert len(rows)==2113 and summary['negative_count']==summary['zero_count']==0
    assert F(summary['normalized_minimum'])>0 and F(summary['absolute_minimum'])>0
    return dict(complete_profiles=2113,previous_negative_profiles=10,rows=rows,summary=summary,
        low9_concentration_only=True,positive_terminal_absolute_surplus_extends_to_all_actual_cusps=True,
        old_conditional_targets_now_superseded_by_stronger_proved_local_rules=True)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in R.B.C.L.STEMS}
    previous={(r['source'],R.B.C.L.N.normalize(r['case'])):r for r in docs[R.TAG]['combined_literal_interface']['rows']}
    rules={12:docs[I12.TAG]['rules'],16:docs[I16.TAG]['records']}
    rows=[];changed=[];counts=dict(order12_debits=0,order16_debits=0,D5_D11_credits=0,A1_B11_credits=0,finite7_debits=0)
    finite_max_order=0
    for stem,index,g,white,cs in I16.T.literal_inputs(inputs):
        before=previous[stem,R.B.C.L.N.normalize(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            images={};debit12=F();debit16=F();credit=F();coupons=[]
            sig=I16.Q.inertia(I16.A.sub(g,ids))
            if sig==[len(ids),0,0]:
                # This complete saved literal domain is simply laced. In
                # particular it contains no whole B11 for the new coupon.
                bs=I16.D.finite_blocks(g,ids);finite_max_order=max(finite_max_order,max(len(b['nodes']) for b in bs))
                assert all(b['type'][0] in 'ADE' for b in bs)
                coupons=I16.coupon(g,ids,marks);credit=F(len(coupons),2)
            elif sig[1]==0:
                shared=next((p for p in I16.A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared) in rules:
                    graph=I16.T.T.labels(g,shared)
                    if graph is not None:
                        for rec in rules[len(shared)]:
                            pairs={tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in I16.P.embeddings(graph,rec['target'])}
                            images[rec['kind']]=sorted(pairs);changes=set(map(tuple,rec['changed_pairs']))
                            if pairs&changes:
                                assert pairs<=changes
                                if len(shared)==12:debit12=F(1,2) if rec['kind']=='B' else F(1)
                                else:debit16=F(1) if rec['kind']=='B' else F(1,2)
            value=F(before['corners'][j]['retained_payment_lower_bound'])+credit-debit12-debit16;assert value>=0
            counts['order12_debits']+=debit12>0;counts['order16_debits']+=debit16>0;counts['D5_D11_credits']+=credit>0
            corners.append(dict(marks=marks,canonical_images=images,order12_debit=debit12,order16_debit=debit16,
                D5_D11_increment=credit,coupons=coupons,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['order12_debit'] or c['order16_debit'] or c['D5_D11_increment'] for c in corners):changed.append(row)
    assert len(rows)==216 and counts==dict(order12_debits=0,order16_debits=4,D5_D11_credits=4,A1_B11_credits=0,finite7_debits=0)
    assert finite_max_order<=23 # no new depth8 finite7 pair can occur here
    return dict(configurations=216,counts=counts,maximum_whole_finite_component_order=finite_max_order,
        changed=changed,rows=rows,all_totals_at_least_one=True,
        inherited_D16_eighth_credits=docs[R.TAG]['combined_literal_interface']['inherited_D16_eighth_credits'])


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    for stem,h in FINAL_MD.items():assert sha256((ROOT/f'results/{stem}.md').read_bytes()).hexdigest()==h,stem
    assert docs[I12.TAG]['new_order12_family_costs']==dict(A=60,B=62,C=62,D=62)
    assert docs[I16.TAG]['new_order16_family_costs']==dict(A=80,B=102,C=103,D=103)
    assert docs[I12.TAG]['new_coupon']['component_pair_cost_increase']=='1/2'
    assert docs[I16.TAG]['upgraded_coupon']['component_pair_cost_increase']=='1'
    assert docs[G.TAG]['feature_depth']==8 and docs[G.TAG]['finite_distance']==7 and docs[G.TAG]['new_weight']=='0'
    assert docs[G.TAG]['deleted_pair_counts']['formulas']==dict(A='max(0,n-23)',B='max(0,n-24)',D='max(0,n-25)')
    for stem in (I12.TAG,I16.TAG,G.TAG):
        assert docs[stem]['triangle_minimum']==2 and docs[stem]['quadrilateral_minimum']==1
    assert docs[R.TAG]['completed_unrestricted_dimension_bound']==60
    assert F(5*sum(min(j,59-j) for j in range(60))-59*60,2)==405
    assert 5*sum(2*j-59 for j in range(60) if 2*j>59)==4500
    paths={Path(__file__),Path(R.__file__),Path(I12.__file__),Path(I16.__file__),Path(G.__file__),
           Path(Z.__file__),Path(Z.Z.__file__),Path(Z.Z.B.__file__),Path(Z.Z.B.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    fin=finite(saved['finite'] if saved else None);ide=ideal(docs);lit=literal(docs)
    surplus=F(ide['summary']['absolute_minimum'])
    return I16.enc(dict(status='PASS_ORDER12_ORDER16_DEPTH8_AMBIENT60_IDEAL_CLOSURE_AND_FINITE_FRONTIER',
        sources=sources,completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        finite=fin,ideal=ide,combined_literal_interface=lit,
        new_order12_family_costs=docs[I12.TAG]['new_order12_family_costs'],
        new_order16_family_costs=docs[I16.TAG]['new_order16_family_costs'],
        strict_global_constant='4500',constant_counted_once=True,
        necessary_ambient60_inequality=dict(finite_vertex_coefficient='4',constant='4500',ideal_vertex_coefficient=str(surplus)),
        all_ideal_dimension_bound=59,new_unrestricted_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete combined local-payment composition,1044 fixed vectors/792 fresh finite tables,2113 ideal convolutions and216 literal interfaces. Formal finite maximum409 exceeds405; all ideal surpluses are positive. All-ideal ceiling59 follows using completed60 higher tail, but unrestricted ceiling remains60.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if a.verify else None
    out=run(saved)
    if a.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal',out['ideal']['summary'],flush=True)
