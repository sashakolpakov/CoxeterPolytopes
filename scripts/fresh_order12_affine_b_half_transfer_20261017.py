#!/usr/bin/env python3
"""Two affine-B11 terminal/fork pairs have weight one half, total cost65.

Complete selected support, actual low-core variants, and the universal
immediate-D13 provenance are all retained. No new finite coupon is used.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as sp
import fresh_order13_affine_b_singleton_transfer_20261017 as L

ROOT=L.ROOT;Q=L.Q;P=L.P
TAG='fresh_order12_affine_b_half_transfer_20261017'
DIAG='fresh_order12_payment_support_20261017'
OVER='fresh_distance8_high_overlap_20261017'
D13='fresh_affine_d13_long_pair_transfer_20261017'
DROP={(9,10),(9,11)}
EXPECTED=[98,100,119,121,164,167,171,1452,1728,1734,1745,2121,2127,2138]


def sub(g,ids):
    return [[g[i][j] for j in ids] for i in ids]


def encode(g):
    return [[str(v) for v in row] for row in g]


def finite_corners(g,gauge,white,c):
    rows=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        signature=Q.inertia(sub(g,ids))
        if signature[2]:
            assert marks==c['preferred_corner'] and signature==[len(ids)-1,0,1]
            rows.append(dict(marks=marks,signature=signature,kind='affineB11_order12',
                finite_unit=0,finite_half='0',ideal_weight='1/2'))
            continue
        assert signature==[len(ids),0,0]
        cc=L.blocks(g,ids,gauge)
        a,b=[next(r for r in cc if m in r['nodes']) for m in marks]
        unit=int(a is b)
        assert not unit or len(a['nodes'])<=13
        half=F()
        for sm,lg,lm in ((a,b,marks[1]),(b,a,marks[0])):
            if len(sm['nodes'])<=5 and lg['type']=='D11' and lg['roles'].get(str(lm)) in ('fork_center','short_fork_leaf'):
                half=F(1,2)
        rows.append(dict(marks=marks,signature=signature,kind='finite',components=cc,
            finite_unit=unit,finite_half=str(half),ideal_weight='0'))
    finite=sum(F(r['finite_unit'])+F(r['finite_half']) for r in rows)
    ideal=sum(F(r['ideal_weight']) for r in rows)
    assert finite>=F(1,2) and finite+ideal>=1
    return dict(corners=rows,finite_total=str(finite),ideal_total=str(ideal),payment_total=str(finite+ideal))


def constraints(c,core):
    g=[[F(x) for x in row] for row in core['gram']];gauge=list(map(F,core['gauge']))
    n=len(g);assert all(g[i][i]==2*gauge[i] for i in range(n))
    left=c['left'][1];right=[c['right_copy_map'][i] for i in c['right'][1]]
    black=left+right;assert len(set(black))==4
    white=sorted(set(range(n))-set(black))
    assert Q.inertia(sub(g,white))==[len(white),0,0]
    for b in black:
        assert Q.inertia(sub(g,white+[b]))==[len(white)+1,0,0]
    checked=Q.inspect_corner(g,gauge,black,*c['preferred_corner'])
    for k,v in checked.items():assert c[k]==v,(c['index'],k)
    du=Q.graph_dist(g,left);dv=Q.graph_dist(g,right)
    return g,gauge,white,du,dv


def inspect_case(c,core):
    g,gauge,white,du,dv=constraints(c,core)
    n=len(g);allowed=c['possible_contacts'];base=c['marked_component']
    row=dict(index=c['index'],method=c['method'],core_key=c['core_key'],
        gram=encode(g),gauge_squared=list(map(str,gauge)),white=white,
        preferred_corner=c['preferred_corner'],base_nodes=base,allowed_contacts=allowed)
    if c['method']=='radical':
        assert n==15 and len(base)==12 and Q.inertia(g)==[13,1,1]
        radical=list(map(F,core['radical'][0]))
        assert len(core['radical'])==1
        assert all(sum(g[i][j]*radical[j] for j in range(n))==0 for i in range(n))
        assert all(radical[i] for i in allowed)
        pos=[i for i in allowed if radical[i]>0];neg=[i for i in allowed if radical[i]<0]
        shortcut=min(max(du[i] for i in pos)+max(dv[i] for i in neg)+2,
            max(dv[i] for i in pos)+max(du[i] for i in neg)+2)
        assert shortcut==c['shortcut_bound']<=3
        # Under the unpaid-short-path hypothesis every additional white is
        # orthogonal to the entire indefinite core. No indirect attachment
        # can enlarge a core component when every first contact vanishes.
        row.update(outcome='sealed_radical_core',radical=list(map(str,radical)),
            positive_contact_domain=pos,negative_contact_domain=neg,
            mixed_sign_shortcut=shortcut,sealed_payment=finite_corners(g,gauge,white,c))
        return row
    assert c['method']=='lorentz_residual' and n==14 and Q.inertia(g)==[13,1,0]
    inv=Q.inverse(g);assert encode(inv)==core['inverse']
    ends=c['end_sets'];assert set(ends[0]+ends[1])==set(allowed)
    assert set(ends[0]).isdisjoint(ends[1])
    shortcut=max(du[i] for i in ends[0])+max(dv[i] for i in ends[1])+2
    assert shortcut<=8
    assert all(inv[i][j]>0 for end in ends for i in end for j in end)
    projections={i:gauge[i]*inv[i][i] for i in allowed}
    pairs=[]
    for end in ends:
        for i,j in combinations(end,2):
            lower=projections[i]+projections[j]+2*min(gauge[i],gauge[j])*inv[i][j]
            assert lower>2
            pairs.append(dict(nodes=[i,j],projection_lower_bound=str(lower)))
    candidates=[i for i in allowed if projections[i]<=2]
    assert len(candidates)==1 and candidates==[r['contact_node'] for r in c['single_contact_candidates']]
    contact=candidates[0]
    assert projections[contact]==2 and gauge[contact]==2
    assert contact in c['omitted_black'] and contact not in base
    assert c['single_contact_candidates'][0]['residual_square']=='0'
    row.update(inverse=encode(inv),end_sets=ends,both_ends_shortcut=shortcut,
        normalized_singleton_projections={str(i):str(p) for i,p in projections.items()},
        same_end_pair_exclusions=pairs,unique_possible_contact=contact,
        contact_native_label=3,contact_normalized_projection='2',contact_residual='0')
    if len(base)==11:
        assert c['index'] in (164,1452)
        corner=sorted(set(white)|set(c['preferred_corner']))
        cc=L.blocks(g,corner,gauge)
        marked=next(r for r in cc if set(c['preferred_corner'])<=set(r['nodes']))
        assert marked['type']=='B11' and marked['nodes']==base
        assert not(set(base)&set(candidates))
        row.update(outcome='whole_finite_B11_before_any_affine_completion',
            whole_finite_component=marked,finite_unit=1,payment_total='1')
        return row
    assert len(base)==12 and c['index'] in (167,171,1745,2138)
    # At most one nonorthogonal new white exists: its norm2 residual is0
    # in the positive-definite complement of the Lorentz core. Identical
    # projection forces identical normal; every core-orthogonal white is
    # orthogonal to it. Label>=4 would project to square>=4>2.
    before=finite_corners(g,gauge,white,c)
    # Scale this new root by sqrt2. Its norm is4, and its one native unit
    # contact to the already sqrt2-scaled old root becomes Gram entry-2.
    extended=[r+[F(-2 if i==contact else 0)] for i,r in enumerate(g)]
    extended.append([F(-2 if i==contact else 0) for i in range(n)]+[F(4)])
    gauge2=gauge+[F(2)]
    assert Q.inertia(extended)==[13,1,1]
    after=finite_corners(extended,gauge2,white+[n],c)
    row.update(outcome='zero_or_one_sealed_unit_contact_root',
        no_extra_root_payment=before,one_extra_root_gram=encode(extended),
        one_extra_root_gauge_squared=list(map(str,gauge2)),one_extra_root_payment=after)
    return row


def immediate_gap(saved):
    g=2*sp.eye(14)
    for i,j in [(i,i+1) for i in range(9)]+[(0,10),(0,11),(9,12),(9,13)]:
        g[i,j]=g[j,i]=-1
    nu=sp.Matrix([1]*10+[sp.Rational(1,2)]*4)
    assert g*nu==sp.zeros(14,1)
    t=sp.symbols('t',nonnegative=True);rows=[]
    for b in (9,12,13):
        a=0;R=9+2*(b!=9);W=sorted(set(range(14))-{a,b});J=g.extract(W,W).inv()
        contact={0:t,1:sp.sqrt(2)};w=sp.Matrix([contact.get(i,0) for i in W])
        far=-nu[b]*g.extract([b],W).T
        height=sum(nu[i]*x for i,x in contact.items())
        moment=sp.simplify(R*(far.T*J*w)[0])
        c0=sp.simplify(2-(w.T*J*w)[0]-moment**2/R)
        assert height==t+sp.sqrt(2) and moment==sp.sqrt(2) and c0==0
        cost=sp.simplify(2*moment/height+c0/height**2)
        assert cost.subs(t,0)==2 and sp.diff(cost,t)<0
        farmax=6 if b==9 else 8
        gap=F(1)-F(2+farmax,R)
        assert gap==(F(1,9) if b==9 else F(1,11))
        inherited=next(r for r in saved['formulas']['end_rows']
            if r['marks']==[a,b] and r['side']=='left' and r['profile']=='inner4')
        assert inherited['c0']=='0' and inherited['moment']=='sqrt(2)'
        rows.append(dict(marks=[a,b],R=R,height=str(height),moment=str(moment),c0=str(c0),
            cost=str(cost),maximum_cost='2',far_maximum_cost=farmax,strict_gap=str(gap)))
    return dict(rows=rows,reflection='The right-end inner4 profiles are the labelled reflections.',
        arbitrary_spectators='The inherited correction (x/h-y/k)^T G_Z^-1 (x/h-y/k)/R is nonnegative.',
        scope='All native opposite magnitudes t>=0 and all opposite-end profiles, including the former immediate profiles.')


def run():
    stems=(DIAG,OVER,D13,'fresh_finite8_ideal5_combined_payments_20261017',
        'fresh_finite_a17_endpoint_payment_support_20261017')
    inputs={}
    for stem in stems:
        d=json.loads((ROOT/f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS')
        for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
        inputs[stem]=d
    diag=inputs[DIAG];rec=next(r for r in diag['records'] if r['kind']=='B')
    target=rec['target'];autos=P.embeddings(target,target)
    assert all({tuple(sorted(mp[i] for i in pair)) for pair in DROP}==DROP for mp in autos)
    basic=set(map(tuple,rec['basic_pairs']))
    assert len(basic)==64 and basic==set(combinations(range(12),2))-DROP
    affected=[r for r in rec['high_case_embeddings'] if DROP&set(map(tuple,r['marked_images']))]
    assert [r['case'] for r in affected]==EXPECTED
    rows=[];mapcount=0
    for r in affected:
        c=inputs[OVER]['isolation']['cases'][r['case']]
        pattern=next(p for p in diag['high_patterns'] if p['index']==r['case'])
        maps=P.embeddings(pattern['labels'],target)
        assert maps==r['embeddings']
        images={tuple(sorted(mp[i] for i in pattern['marks'])) for mp in maps}
        assert images<=DROP and images
        result=inspect_case(c,inputs[OVER]['isolation']['cores'][c['core_key']])
        result.update(embeddings=maps,marked_images=sorted(images))
        rows.append(result);mapcount+=len(maps)
    assert Counter(r['outcome'] for r in rows)=={
        'sealed_radical_core':8,'whole_finite_B11_before_any_affine_completion':2,
        'zero_or_one_sealed_unit_contact_root':4}
    assert all(not(DROP&set(map(tuple,r['marked_images']))) for r in rec['low_case_embeddings'])
    assert len(diag['low_corner_domain']['cases'])==66
    assert sum(len(r['variants']) for r in diag['low_corner_domain']['cases'])==69
    immediate=[r for r in diag['immediate_D13_alternatives']
        if r['target_kind']=='B' and DROP&set(map(tuple,r['marked_images']))]
    assert len(immediate)==2 and all(r['source_profile']=='inner4' for r in immediate)
    assert [r['literal_marks'] for r in immediate]==[[12,14],[13,14]]
    _,ds,features=P.distance_data(target)
    triangles={p for p in combinations(range(12),2) if ds[p[0]][p[1]]<=5 or
        ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features)}
    assert triangles<=basic and set(map(tuple,rec['short8_pairs']))<=basic
    weights=[dict(pair=list(p),weight='1/2' if p in DROP else '1') for p in combinations(range(12),2)]
    assert sum(F(r['weight']) for r in weights)==65
    paths=[Path(__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER12_AFFINE_B_TWO_HALF_WEIGHTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='B11',diagram_order=12,target=target,automorphisms=autos,
        old_cost=66,new_cost='65',half_weight_pairs=sorted(DROP),weighted_pairs=weights,
        affected_cases=EXPECTED,affected_maps=mapcount,cases=rows,
        low_actual_domain_cases=len(diag['low_corner_domain']['cases']),
        low_actual_variants=sum(len(r['variants']) for r in diag['low_corner_domain']['cases']),
        low_affected_cases=[],immediate_D13_support=immediate,immediate_D13_strict_gap=immediate_gap(inputs[D13]),
        triangle_pairs_retained_at_weight_one=sorted(triangles),short8_pairs_retained_at_weight_one=rec['short8_pairs'],
        basic_pairs_retained_at_weight_one=sorted(basic),added_finite_coupons=[],
        scope='All14 selected contexts, complete actual low domain, and universal immediate-D13 provenance. Whole finite A1/A2 against D11 uses only inherited halves. No new finite coupon, A/C order12 rule or ambient dimension conclusion.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_maps'],'maps; B11 cost66 to65',flush=True)
