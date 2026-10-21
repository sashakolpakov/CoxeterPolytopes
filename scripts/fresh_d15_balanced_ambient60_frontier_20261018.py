#!/usr/bin/env python3
"""Count D15=100.5 with rank15 representative allowance101 after balancing.

Only fixed rank15 changes. The complete792 residual certificates are reused
and checked, not optimized again. Actual and representative bases differ.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_center_reversion_rank16_balance_ambient60_frontier_20261018 as B
import fresh_finite_d15_three_singleton_half_transfer_20261018 as L

ROOT=B.ROOT;Z=B.Z;DIM=60
TAG='fresh_d15_balanced_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((B.TAG,L.TAG,*B.STEMS)))
FINAL_MD={B.TAG:'741fc285012fe3626dbe61a5aa2bce6adb230379a23140aec00173b0f9a987bc',
 L.TAG:'7dcfbd0056b95d4b53a04c34eddce21d15dc0895e89b89d724e3e0d918ebe70b'}


def actual_base(t):return B.base(t,1)-(12 if t==('D',15) else 0)
def upper_base(t):return B.base(t,1)-(8 if t==('D',15) else 0)
def actual_cross(a,b):return B.actual_cross(a,b)
def upper_cross(a,b):return B.upper_cross(a,b)
def item(r,key):return B.item(r,key,1)


def actual_value(forest):
    return sum(actual_base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def finite(previous,saved=None):
    assert 15 not in Z.Z.ITEMS
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=DIM]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert [actual_base((k,15)) for k in 'ABD']==[808,776,804]
    assert upper_base(('D',15))==808
    domain=Z.Z.large_domain();keys=sorted({key for cs,used,key in domain})
    assert len(domain)==1044 and len(keys)==792
    tabs=saved['tables'] if saved else previous['tables'];assert tabs==previous['tables']
    tables={tuple(t['key']):t for t in tabs};assert set(tables)==set(keys) and len(tabs)==len(keys)
    checks=sum(B.check_table(tables[k],k,1) for k in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    before={tuple(r['counts']):r for r in previous['alternatives']};rows=[];items=0
    assert set(before)=={cs for cs,used,key in domain}
    for cs,used,key in domain:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key)==upper_base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=DIM-used;tab=tables[key];value=sum(upper_base(x) for x in fixed)+tab['values'][rem]
        assert value==before[cs]['octupled_upper']-8*cs[6]
        forest=fixed[:]
        while rem:
            r=tab['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==DIM
        upper=sum(upper_base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual_value(forest)<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=value,octupled_actual=actual_value(forest),
            octupled_fixed_rank15_upper_saving=8*cs[6]))
    maximum=max(r['octupled_upper'] for r in rows);assert maximum==3260
    witness=[('D',17),('D',17),('D',26)];assert actual_value(witness)==maximum
    assert actual_value([('D',15)]*4)==3216 and actual_value([('A',15)]*4)==3232
    return dict(exact_maximum='815/2',previous_maximum='408',curvature='405',finite_deficit='5/2',
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=tabs,alternatives=rows,above_curvature=sum(r['octupled_upper']>3240 for r in rows),
        actual_attainer=['D17','D17','D26'],actual_attaining_cost='815/2',
        actual_rank15_costs=dict(A15='101',B15='97',D15='201/2'),rank15_representative_upper='101',
        former_four_D15_actual_cost='402',four_A15_actual_cost='404',
        D5_D12_actual_cost='23/2',all_cross_component_costs_unchanged=True,
        finite_cost_unit='1/8',residual_tables_identical_to_frozen_predecessor=True,
        new_dynamic_program_run=False,replay_optimizes=False,formal_forest_attainment_not_geometric_realization=True)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in B.C.R.B.C.L.STEMS}
    prev={(r['source'],B.C.R.B.C.L.N.normalize(r['case'])):r for r in docs[B.TAG]['stages'][1]['literal_interface']['rows']}
    local={(r['source'],B.C.R.B.C.L.N.normalize(r['case'])):r for r in docs[L.TAG]['combined_literal_interface']['rows']}
    rows=[];changed=[];debits=0
    for stem,index,g,white,cs in L.T.literal_inputs(inputs):
        key=(stem,B.C.R.B.C.L.N.normalize(index));before=prev[key];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();images=[]
            if L.Q.inertia(L.A.sub(g,ids))==[len(ids),0,0]:
                block=next((b for b in L.D.finite_blocks(g,ids) if b['type']=='D15' and set(marks)<=set(b['nodes'])),None)
                if block:
                    images=L.T.mark_pairs(g,block['nodes'],marks,15,'D')
                    if set(images)&L.CHANGES:assert set(images)<=L.CHANGES;debit=F(1,2)
            assert debit==F(local[key]['corners'][j]['debit'])
            value=F(before['corners'][j]['retained_payment_lower_bound'])-debit;assert value>=0;debits+=debit>0
            corners.append(dict(marks=marks,D15_images=images,D15_debit=debit,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['D15_debit'] for c in corners):changed.append(row)
    assert len(rows)==216 and debits==8
    return dict(configurations=216,debited_corners=debits,changed=changed,rows=rows,all_totals_at_least_one=True,
        earlier_center_reversion_and_finite_balance_preserved=True)


def run(saved=None):
    assert all(h!='PENDING' for h in FINAL_MD.values()),'Await final frozen local and baseline reports.'
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    for stem,h in FINAL_MD.items():assert sha256((ROOT/f'results/{stem}.md').read_bytes()).hexdigest()==h,stem
    assert docs[B.TAG]['final_finite_maximum']=='408' and docs[L.TAG]['new_cost']=='201/2'
    assert docs[L.TAG]['no_new_coupons'] and docs[L.TAG]['required_rank15_upper_base']=='101'
    assert docs[L.TAG]['rank15_actual_costs']==dict(A15='101',B15='97',D15='201/2')
    assert sum(F(r['weight']) for r in docs[L.TAG]['weighted_pairs'])==F(actual_base(('D',15)),8)
    d16={tuple(r['pair']):F(r['weight']) for r in docs[B.TAG]['local_proofs']['finite_weighted_pairs']['D16']}
    assert d16[1,12]==d16[1,13]==F(1,2)
    assert not set(map(tuple,docs[L.TAG]['existing_D16_half_pairs']))&B.PAIRS['D16']
    assert docs[L.TAG]['existing_A2_D13_role']=='fork_center10'
    assert docs[L.TAG]['triangle_minimum']==2 and docs[L.TAG]['quadrilateral_minimum']==1
    ide=B.ideal(docs);assert ide==docs[B.TAG]['ideal']
    assert ide['summary']['normalized_minimum']=='31/6000' and ide['summary']['absolute_minimum']=='54'
    paths={Path(__file__),Path(B.__file__),Path(L.__file__),Path(L.T.__file__),Path(L.D.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return L.enc(dict(status='PASS_D15_BALANCED_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        finite=finite(docs[B.TAG]['stages'][1]['finite'],saved['finite'] if saved else None),
        ideal=ide,combined_literal_interface=literal(docs),strict_global_constant='4500',constant_counted_once=True,
        necessary_ambient60_inequality=dict(finite_vertex_coefficient='5/2',constant='4500',ideal_vertex_coefficient='54'),
        all_ideal_dimension_bound=59,new_unrestricted_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete fixed-vector upper update with distinct actualD15=100.5 and representative101; all unchanged residual certificates checked. Current216literal interfaces and2113idealconvolutions reconstructed. Exactformal407.5 remainsabove405; unrestricted60 unchanged.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if a.verify else None
    out=run(saved)
    if a.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal',out['ideal']['summary'],flush=True)
