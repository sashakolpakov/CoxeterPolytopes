#!/usr/bin/env python3
"""Twelve affine-D16 pair weights reduced to one half, using existing coupons."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_order13_affine_b_singleton_transfer_20261017 as L

ROOT=L.ROOT;P=L.P;Q=L.Q
TAG='fresh_d17_partial_weight_transfer_20261017'
OLD='fresh_d17_singleton_pair_coupon_20261017'
OVER='fresh_distance8_high_overlap_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
HALF_PAIRS={(0,10),(2,12),(0,11),(1,12),(1,15),(1,16),
            (11,13),(11,14),(2,15),(2,16),(10,13),(10,14)}
EXPECTED=[1311,1314,1974,1978,2367,2371,2494,2504,2511,2529,2535,2539]


def inspect(c,core):
    g=[[F(x) for x in row] for row in core['gram']];N=set(range(len(g)))
    assert len(g)==19 and core['gauge']==[1]*19 and Q.inertia(g)==[17,1,1]
    assert c['core_key']=='(25, 25, ((9, 9),))' and c['method']=='radical'
    z=list(map(F,core['radical'][0]));assert z[3]*z[13]<0
    assert all(sum(g[i][j]*z[j] for j in N)==0 for i in N)
    black=set(c['left'][1])|{c['right_copy_map'][i] for i in c['right'][1]}
    white=sorted(N-black);assert len(black)==4 and Q.inertia(L.sub(g,white))==[15,0,0]
    for b in black:
        assert Q.inertia(L.sub(g,sorted(set(white)|{b})))==[16,0,0]
    ids=c['marked_component'];assert ids==sorted(N-{3,13})
    assert set(ids)==set(white)|set(c['preferred_corner'])
    assert Q.inertia(L.sub(g,ids))==[16,0,1]
    k=list(map(F,c['blocks'][0]['positive_kernel']));assert all(x>0 for x in k)
    assert all(sum(g[i][j]*k[t] for t,j in enumerate(ids))==0 for i in ids)
    assert c['possible_contacts']==[3,13] and c['shortcut_bound']==2
    rows=[]
    for corner in c['all_corners']:
        marks=corner['preferred_corner'];ns=sorted(set(white)|set(marks));sig=Q.inertia(L.sub(g,ns))
        if sig==[16,0,1]:
            assert ns==ids and marks==c['preferred_corner']
            rows.append(dict(marks=marks,signature=sig,kind='affineD16_order17',
                             ideal_weight='1/2',finite_unit=0,finite_half='0'))
            continue
        assert sig==[17,0,0]
        cc=L.blocks(g,ns,core['gauge']);assert sorted(r['nodes'] for r in cc)==sorted(r['nodes'] for r in corner['blocks'])
        selected=[next(r for r in cc if m in r['nodes']) for m in marks]
        unit=int(selected[0] is selected[1]);coupon=F(0);roles=[]
        for i,j in ((0,1),(1,0)):
            a,b=selected[i],selected[j]
            role=b['roles'].get(str(marks[j]))
            if len(a['nodes'])<=5 and b['type'] in {'D'+str(t) for t in range(11,17)} and role in ('fork_center','short_fork_leaf'):
                coupon=F(1,2)
                roles.append(dict(small_type=a['type'],small_mark=marks[i],large_type=b['type'],
                                  large_mark=marks[j],large_role=role))
        rows.append(dict(marks=marks,signature=sig,kind='finite',components=cc,
                         marked_component_indices=[cc.index(r) for r in selected],ideal_weight='0',
                         finite_unit=unit,finite_half=str(coupon),coupon_roles=roles))
    assert len(rows)==4 and sum(F(r['ideal_weight']) for r in rows)==F(1,2)
    assert sum(F(r['finite_half']) for r in rows)>=F(1,2)
    total=sum(F(r['ideal_weight'])+r['finite_unit']+F(r['finite_half']) for r in rows)
    assert total>=1
    return dict(index=c['index'],left=c['left'],right=c['right'],match=c['match'],core_key=c['core_key'],
                gram=L.encode(g),signature=[17,1,1],radical=list(map(str,z)),white=white,
                affine_nodes=ids,positive_affine_kernel=list(map(str,k)),omitted=[3,13],
                omitted_radical_coefficients=[str(z[i]) for i in (3,13)],corners=rows,
                ideal_total='1/2',finite_coupon_total=str(sum(F(r['finite_half']) for r in rows)),
                finite_unit_total=sum(r['finite_unit'] for r in rows),payment_total=str(total))


def run():
    stems=[OLD,OVER,SMALL,'fresh_finite8_ideal5_combined_payments_20261017',
           'fresh_finite8_a17_d11_coupon_budget_20261017']
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,j in inputs.items():
        assert j['status'].startswith('PASS'),stem
        for name,h in j['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    old=inputs[OLD];over=inputs[OVER]
    rec=next(r for r in inputs[SMALL]['records'] if r['order']==17 and r['kind']=='D')
    target=rec['target'];autos=P.embeddings(target,target);_,ds,features=P.distance_data(target)
    paid=set(map(tuple,old['new_paid_pairs']));assert len(paid)==126
    assert len(HALF_PAIRS)==12 and HALF_PAIRS<=paid
    assert all(ds[a][b]>9 for a,b in HALF_PAIRS)
    assert all({tuple(sorted(a[i] for i in p)) for p in HALF_PAIRS}==HALF_PAIRS for a in autos)
    affected=[r['case'] for r in rec['case_embeddings'] if set(map(tuple,r['marked_images']))&HALF_PAIRS]
    assert affected==EXPECTED
    rows=[inspect(over['isolation']['cases'][i],over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in affected]
    weights={p:(F(1,2) if p in HALF_PAIRS else F(1) if p in paid else F(0))
             for p in combinations(range(17),2)}
    assert sum(weights.values())==120
    assert weights[(1,11)]==0
    triangles=set(map(tuple,old['triangle_pairs_retained']))
    assert all(weights[p]==1 for p in triangles)
    prior_paid=set(map(tuple,rec['retained_paid_pairs']))
    short8={p for p in combinations(range(17),2) if ds[p[0]][p[1]]<=8}
    assert short8<=prior_paid and all(weights[p]==1 for p in short8)
    orbits=[];todo=set(HALF_PAIRS)
    while todo:
        p=min(todo);o={tuple(sorted(a[i] for i in p)) for a in autos};assert o<=todo;todo-=o
        cases=[r['case'] for r in rec['case_embeddings'] if set(map(tuple,r['marked_images']))&o]
        orbits.append(dict(representative=list(p),distance=ds[p[0]][p[1]],pairs=[list(p) for p in sorted(o)],
                           old_weight='1',new_weight='1/2',cases=cases))
    assert [len(r['pairs']) for r in orbits]==[2,2,4,4]
    assert [len(r['cases']) for r in orbits]==[2,2,4,4]
    paths=[Path(__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_AFFINE_D16_TWELVE_HALF_WEIGHT_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='D16',diagram_order=17,target=target,automorphisms=autos,
                inherited_cost=126,new_cost='120',half_weight_pairs=[list(p) for p in sorted(HALF_PAIRS)],
                weighted_pairs=[dict(pair=list(p),weight=str(w)) for p,w in sorted(weights.items())],
                half_weight_orbits=orbits,affected_cases=affected,cases=rows,
                triangle_pairs_retained_at_weight_one=[list(p) for p in sorted(triangles)],
                short8_pairs_retained_at_weight_one=[list(p) for p in sorted(short8)],
                added_finite_coupons=[],
                scope='Complete saved-map support and twelve actual sealed-core replacements. Each selected ideal half combines with an existing finite half. This is a weighted local rule, not a binary paid-pair set. No new ambient bound or scalar comparison is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'12 half-weight pairs,12 sealed contexts,48 corners; total120',flush=True)
