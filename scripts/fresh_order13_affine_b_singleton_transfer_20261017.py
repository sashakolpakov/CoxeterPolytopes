#!/usr/bin/env python3
"""One order13 affine-B payment removed using two sealed cores and a strict gap.

The four distance11 contexts are recorded separately as conditional follow-up;
only the canonical singleton (0,10) is removed by this package.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as sp
import fresh_ideal5_small_factor_payments_20261017 as S
import fresh_selective_distance10_dimension116_20261017_review as Q

P=S.P
ROOT=P.ROOT
TAG='fresh_order13_affine_b_singleton_transfer_20261017'
DIAG='fresh_order13_payment_support_20261017'
OVER='fresh_distance8_high_overlap_20261017'
D14='fresh_affine_d14_long_pair_transfer_20261017'
DROP=(0,10)


def sub(g,ids):
    return [[g[i][j] for j in ids] for i in ids]


def encode(g):
    return [[str(x) for x in row] for row in g]


def blocks(g,ids,gauge):
    adj={i:{j for j in ids if i!=j and g[i][j]} for i in ids}
    todo=set(ids)
    out=[]
    while todo:
        component={min(todo)}
        stack=list(component)
        while stack:
            i=stack.pop()
            for j in adj[i]-component:
                component.add(j)
                stack.append(j)
        todo-=component
        ns=sorted(component)
        deg={i:len(adj[i]) for i in ns}
        bonds=[(i,j,g[i][j]**2/F(gauge[i]*gauge[j]))
               for i in ns for j in adj[i] if i<j]
        assert len(bonds)==len(ns)-1
        assert all(v in (1,2) for i,j,v in bonds)
        doubled=[(i,j) for i,j,v in bonds if v==2]
        forks=[i for i in ns if deg[i]==3]
        if doubled:
            assert len(doubled)==1 and not forks and max(deg.values())<=2
            ends=[i for i in doubled[0] if deg[i]==1]
            assert len(ends)==1
            kind='B'+str(len(ns));roles={str(ends[0]):'doubled_terminal'}
        elif forks:
            assert len(forks)==1 and max(deg.values())==3
            leaves=[i for i in adj[forks[0]] if deg[i]==1]
            assert len(leaves)==2
            kind='D'+str(len(ns))
            roles={str(forks[0]):'fork_center',**{str(i):'short_fork_leaf' for i in leaves}}
        else:
            assert max(deg.values())<=2
            kind='A'+str(len(ns))
            roles={str(i):'endpoint' for i in ns if deg[i]<=1}
        assert Q.inertia(sub(g,ns))==[len(ns),0,0]
        out.append(dict(nodes=ns,type=kind,roles=roles))
    return out


def inspect_case(c,core):
    g=[[F(v) for v in row] for row in core['gram']]
    n=len(g);gauge=core['gauge']
    assert n==15 and Q.inertia(g)==[13,1,1]
    assert all(g[i][i]==2*gauge[i] for i in range(n))
    radical=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*radical[j] for j in range(n))==0 for i in range(n))
    black=set(c['left'][1])|{c['right_copy_map'][i] for i in c['right'][1]}
    assert len(black)==4
    white=sorted(set(range(n))-black)
    assert Q.inertia(sub(g,white))==[11,0,0]
    for b in black:
        assert Q.inertia(sub(g,sorted(set(white)|{b})))==[12,0,0]
    affine=c['marked_component'];omitted=sorted(set(range(n))-set(affine))
    assert omitted==sorted(c['omitted_black']) and len(affine)==13
    assert set(affine)==set(white)|set(c['preferred_corner'])
    assert Q.inertia(sub(g,affine))==[12,0,1]
    kernel=list(map(F,c['blocks'][0]['positive_kernel']))
    assert all(v>0 for v in kernel)
    assert all(sum(g[i][j]*kernel[k] for k,j in enumerate(affine))==0 for i in affine)
    assert radical[omitted[0]]*radical[omitted[1]]<0
    assert c['possible_contacts']==c['omitted_black'] and c['shortcut_bound']==2
    rows=[]
    for corner in c['all_corners']:
        marks=corner['preferred_corner'];ids=sorted(set(white)|set(marks))
        signature=Q.inertia(sub(g,ids))
        if signature==[12,0,1]:
            assert ids==affine
            rows.append(dict(marks=marks,signature=signature,kind='affineB12_order13',
                             finite_unit=0,old_half='0',new_half='0'))
            continue
        assert signature==[13,0,0]
        cc=blocks(g,ids,gauge)
        assert sorted(r['nodes'] for r in cc)==sorted(r['nodes'] for r in corner['blocks'])
        marked=[next(r for r in cc if m in r['nodes']) for m in marks]
        unit=int(marked[0] is marked[1])
        old=new=F(0)
        for i,j in ((0,1),(1,0)):
            a,b=marked[i],marked[j];bm=marks[j]
            if len(a['nodes'])<=5 and b['type']=='D12' and b['roles'].get(str(bm)) in ('fork_center','short_fork_leaf'):
                old=F(1,2)
            if a['type']=='A2' and b['type']=='B10' and b['roles'].get(str(bm))=='doubled_terminal':
                new=F(1,2)
        rows.append(dict(marks=marks,signature=signature,kind='finite',components=cc,
                         marked_component_indices=[cc.index(r) for r in marked],
                         finite_unit=unit,old_half=str(old),new_half=str(new)))
    return dict(index=c['index'],left=c['left'],right=c['right'],match=c['match'],
                core_key=c['core_key'],gram=encode(g),gauge_squared=gauge,signature=[13,1,1],
                radical=list(map(str,radical)),white=white,affine_nodes=affine,
                positive_affine_kernel=list(map(str,kernel)),omitted=omitted,
                omitted_radical_coefficients=[str(radical[i]) for i in omitted],corners=rows,
                finite_unit_total=sum(r['finite_unit'] for r in rows),
                old_half_total=str(sum(F(r['old_half']) for r in rows)),
                new_half_total=str(sum(F(r['new_half']) for r in rows)))


def immediate_gap(saved):
    g=2*sp.eye(15)
    for i,j in [(i,i+1) for i in range(10)]+[(0,11),(0,12),(10,13),(10,14)]:
        g[i,j]=g[j,i]=-1
    nu=sp.Matrix([1]*11+[sp.Rational(1,2)]*4)
    assert g*nu==sp.zeros(15,1)
    x=sp.symbols('x',nonnegative=True)
    rows=[]
    for b in (10,13,14):
        a=0;R=10+2*(b!=10)
        W=sorted(set(range(15))-{a,b});J=g.extract(W,W).inv()
        contacts={0:x,1:sp.sqrt(2)}
        w=sp.Matrix([contacts.get(i,0) for i in W]);vb=-nu[b]*g.extract([b],W).T
        height=sum(nu[i]*v for i,v in contacts.items())
        moment=sp.simplify(R*(vb.T*J*w)[0])
        c0=sp.simplify(2-(w.T*J*w)[0]-moment**2/R)
        assert height==x+sp.sqrt(2) and moment==sp.sqrt(2) and c0==0
        cost=sp.simplify(2*moment/height+c0/height**2)
        assert cost.subs(x,0)==2 and sp.diff(cost,x)<0
        far_max=6 if b==10 else 8
        gap=F(1)-F(2+far_max,R)
        assert gap==(F(1,5) if b==10 else F(1,6))
        inherited=next(r for r in saved['profile_formulas']['end_rows']
                       if r['marks']==[a,b] and r['side']=='left' and r['profile']=='inner4')
        assert inherited['c0']=='0' and inherited['moment']=='sqrt(2)'
        rows.append(dict(marks=[a,b],resistance=R,height=str(height),moment=str(moment),c0=str(c0),
                         cost=str(cost),cost_maximum='2',far_cost_maximum=far_max,
                         normalized_determinant_lower_bound=str(gap)))
    assert saved['profile_formulas']['determinant']=='1-(2*l/h+c0/h^2+2*m/k+d0/k^2)/R'
    return dict(profile='inner4',rows=rows,
                reflected_cases='Identical after the labelled reflection i -> 10-i, leaves11/12 <->13/14.',
                all_label_scope='x ranges over all nonnegative reals; far-center bound6 includes both former immediate profiles of costs<=2 and<=4.',
                spectator_correction='(1/R)*(xZ/h-yZ/k)^T GZ^-1*(xZ/h-yZ/k) >= 0')


def run():
    stems=[DIAG,OVER,D14,'fresh_finite8_ideal5_combined_payments_20261017']
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,data in inputs.items():
        assert data['status'].startswith('PASS'),stem
        for name,h in data['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    diag=inputs[DIAG];over=inputs[OVER]
    assert diag['eligible_patterns']==801
    rec=next(r for r in diag['records'] if r['kind']=='B')
    target=rec['target'];autos=P.embeddings(target,target);_,ds,features=P.distance_data(target)
    paid=set(map(tuple,rec['required_union']));basic=set(map(tuple,rec['basic_pairs']))
    assert len(paid)==78 and len(basic)==73 and DROP in paid-basic and ds[0][10]==10
    assert all(tuple(sorted(a[i] for i in DROP))==DROP for a in autos)
    affected=[r['case'] for r in rec['case_embeddings'] if DROP in set(map(tuple,r['marked_images']))]
    assert affected==[160,1206]
    rows=[]
    for index in affected:
        c=over['isolation']['cases'][index]
        r=inspect_case(c,over['isolation']['cores'][c['core_key']])
        assert r['finite_unit_total']==0 and r['old_half_total']=='1/2' and r['new_half_total']=='1/2'
        rows.append(r)
    immediate=[r for r in diag['immediate_D14_alternatives'] if DROP in set(map(tuple,r['marked_images']))]
    assert len(immediate)==1 and immediate[0]['source_profile']=='inner4' and immediate[0]['literal_marks']==[10,15]
    gap=immediate_gap(inputs[D14])
    reduced=paid-{DROP};assert len(reduced)==77 and basic<=reduced
    triangle={p for p in combinations(range(13),2) if ds[p[0]][p[1]]<=5 or
              (ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
    assert triangle<=reduced and all(ds[a][b]>8 for a,b in paid-reduced)
    # Diagnostic only: these extra pairs are NOT removed here.
    long_pairs={(10,11),(10,12)}
    long_indices=[r['case'] for r in rec['case_embeddings'] if long_pairs&set(map(tuple,r['marked_images']))]
    assert long_indices==[166,170,1744,2137]
    follow=[]
    for index in long_indices:
        c=over['isolation']['cases'][index]
        r=inspect_case(c,over['isolation']['cores'][c['core_key']])
        if index in (166,1744):
            assert r['finite_unit_total']==1
            r['conditional_followup']='Already has a whole finite B13 unit after sealing.'
        else:
            assert r['finite_unit_total']==0 and r['old_half_total']=='1/2'
            half=next(q for q in r['corners'] if q['old_half']=='1/2')
            d=next(p for p in half['components'] if p['type']=='D12')
            dm=next(m for m in half['marks'] if m in d['nodes'])
            assert d['roles'][str(dm)]=='short_fork_leaf'
            r['conditional_followup']='A full whole-A1 / whole-D12 short-leaf coupon would pay; it is not supplied by this package.'
        follow.append(r)
    paths=[Path(__file__),Path(Q.__file__),Path(S.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER13_AFFINE_B_SINGLETON_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='B12',diagram_order=13,old_cost=78,new_cost=77,
                target=target,automorphisms=autos,removed_pairs=[list(DROP)],
                new_paid_pairs=[list(p) for p in sorted(reduced)],basic_pairs_retained=True,
                triangle_pairs_retained=[list(p) for p in sorted(triangle)],
                selected_support=dict(eligible_patterns=801,affected_cases=affected,cases=rows),
                immediate_D14_support=immediate,immediate_D14_strict_gap=gap,
                new_finite_coupon=dict(weight='1/2',small_type='whole A2',small_marks='either endpoint',
                                       large_type='whole B10',large_mark='unique doubled terminal'),
                conditional_long_orbit=dict(pairs=[list(p) for p in sorted(long_pairs)],
                                           removed=False,cases=follow),
                scope='Both selected-core and universal immediate-D14 provenances of the singleton are discharged. The new coupon needs a separate finite-vertex budget; no global dimension bound or removal of the distance11 orbit is claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'2 sealed contexts; 3 strict immediate gaps; 4 conditional contexts',flush=True)
