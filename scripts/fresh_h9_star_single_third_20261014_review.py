#!/usr/bin/env python3
"""Independent literal-domain and eleven-polynomial single-third review."""
import argparse
import ast
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from math import comb
from pathlib import Path
import fresh_h9_star_lastblock_20261014_review as R

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_third_20261014'
S=[3,4,6,7,9,10,12,13];CORE=sorted(S+[11,14,17])
Z,O=R.Z,R.O
def k(x):return R.scale(O,F(x))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(v):return tuple(map(F,v))
def enc(v):return list(map(str,v))
def clean(p):return{m:v for m,v in p.items()if v!=Z}
def pa(p,q):
    out=dict(p)
    for m,v in q.items():out[m]=R.add(out.get(m,Z),v)
    return clean(out)
def ps(p,x):return clean({m:R.mul(v,x)for m,v in p.items()})
def pm(p,q):
    out={}
    for m,v in p.items():
        for n,w in q.items():
            e=tuple(a+b for a,b in zip(m,n));out[e]=R.add(out.get(e,Z),R.mul(v,w))
    return clean(out)
def pp(p,n):
    if not p:assert n>0;return{}
    out={tuple(0 for _ in next(iter(p))):O}
    for _ in range(n):out=pm(out,p)
    return out
def fpow(x,n):
    out=O
    for _ in range(n):out=R.mul(out,x)
    return out
def shifted(p,shift):
    out={};zero=(0,)*len(shift)
    for ex,coef in p.items():
        term={zero:coef}
        for i,n in enumerate(ex):
            factor={}
            for j in range(n+1):
                e=list(zero);e[i]=j;factor[tuple(e)]=R.scale(fpow(shift[i],n-j),comb(n,j))
            term=pm(term,factor)
        out=pa(out,term)
    return out
def detpoly(M):
    out={}
    for p in permutations(range(3)):
        parity=(-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))
        out=pa(out,ps(pm(pm(M[0][p[0]],M[1][p[1]]),M[2][p[2]]),k(parity)))
    return out
def divide_den_squared(p):
    # Long division in c, with coefficients polynomial in beta and D.
    den2={(0,0,0):O,(0,2,0):k(-2),(0,4,0):O};out={};p=dict(p)
    while p:
        ex=max(p,key=lambda e:(e[1],e[0],e[2]));assert ex[1]>=4
        m=(ex[0],ex[1]-4,ex[2]);v=p[ex]
        out=pa(out,{m:v});p=pa(p,ps(pm({m:O},den2),R.neg(v)))
    return out
def parse_poly(text,one,variables):
    def go(t):
        if isinstance(t,ast.Name):return variables[t.id]
        if isinstance(t,ast.Constant):return ps(one,k(t.value))
        if isinstance(t,ast.UnaryOp)and isinstance(t.op,ast.USub):return ps(go(t.operand),k(-1))
        assert isinstance(t,ast.BinOp)
        if isinstance(t.op,ast.Add):return pa(go(t.left),go(t.right))
        if isinstance(t.op,ast.Sub):return pa(go(t.left),ps(go(t.right),k(-1)))
        if isinstance(t.op,ast.Mult):return pm(go(t.left),go(t.right))
        if isinstance(t.op,ast.Pow):return pp(go(t.left),t.right.value)
        assert isinstance(t.op,ast.Div)and isinstance(t.right,ast.Constant)
        return ps(go(t.left),k(F(1,t.right.value)))
    return go(ast.parse(text,mode='eval').body)

def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):read(f)
        return d
    d=read(f'results/{TAG}.json')
    assert d['status']=='PASS_ALL_LABEL_SINGLE_THIRD_SHARED_ANCHOR_CORE'
    g=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[set(v['incident'])for v in g['vertices']];MN=set(map(frozenset,g['minimal_nonfaces']))
    assert len(V)==378 and len(MN)==23 and d['rank_core']==CORE
    assert any(set(S)<=v for v in V)and all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
    faces={frozenset(v&set(CORE))for v in V};maximal={v for v in faces if not any(v<w for w in faces)}
    assert len(maximal)==9 and maximal==set(map(frozenset,d['complete_maximal_core_faces']))
    base={(3,4):4,(6,7):7,(9,10):7,(12,13):7,(3,14):2,(4,14):5,
          (11,12):2,(11,13):3,(12,17):2,(13,17):3}
    expected=[]
    def expect(kind,nodes,labels,cycle,condition):
        expected.append((kind,tuple(nodes),tuple(labels),cycle,json.dumps(condition,sort_keys=True)))
    for i,j in combinations(range(4),2):
        for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
            aa=next(v for v in S[2*i:2*i+2]if v!=a);bb=next(v for v in S[2*j:2*j+2]if v!=b)
            expect('cross_block_right',[aa,a,b,bb],[4,3,4],False,dict(trial_nonright=[a,b]))
    for anchor in(11,17):
        for block in((6,7),(9,10)):
            for contact in block:
                other=next(x for x in block if x!=contact)
                expect('high_middle_forces_anchor_right',[other,contact,anchor],[6,3],False,dict(trial_nonright=[contact,anchor]))
        for contact in(3,4):
            expect('owner_foreign_requires_internal_at_most_five',[7-contact,contact,anchor],[6,3],False,
                   dict(trial_nonright=[contact,anchor],trial_owner_internal_at_least=6))
        expect('third_anchor_label_at_most_three',[4,14,anchor],[4,4],False,dict(trial_label_at_least_four=[14,anchor]))
        expect('nonright_own_root_foreign_forces_third_right',[4,14,anchor],[3]*3,True,
               dict(owner_foreign=[4,anchor],trial_nonright=[14,anchor]))
    expect('shared_anchors_are_right',[11,13,17],[3]*3,True,dict(trial_nonright=[11,17]))
    for a,b in product((3,4),repeat=2):
        nodes=[a,11,13,17]+([]if a==b else[b])
        expect('at_most_one_owner_foreign_across_both_anchors',nodes,[3]*len(nodes),True,
               dict(two_foreign_contacts=[[a,11],[b,17]]))
    expect('third_cannot_contact_both_anchors',[14,11,13,17],[3]*4,True,dict(two_third_contacts=[[14,11],[14,17]]))
    for a,b in((11,17),(17,11)):
        expect('nonright_own_root_foreign_forces_other_third_right',[4,14,b,13,a],[3]*5,True,
               dict(owner_foreign=[4,a],trial_nonright=[14,b]))
    saved=[(r['kind'],tuple(r['ordered_support']),tuple(r['labels']),r['cycle'],json.dumps(r['conditional'],sort_keys=True))for r in d['positive_kernel_guards']]
    assert Counter(saved)==Counter(expected)and len(saved)==48
    guard_counts=Counter()
    for row in d['positive_kernel_guards']:
        nodes=row['ordered_support'];assert len(set(nodes))==len(nodes)
        assert set(nodes)==set(row['support'])<=V[row['actual_vertex']]and set(nodes)<=set(CORE)
        facts=dict(base)
        for name,value in row['conditional'].items():
            if name=='trial_owner_internal_at_least':facts[3,4]=value
            elif name in('two_foreign_contacts','two_third_contacts'):
                for p in value:facts[tuple(sorted(p))]=3
            else:
                assert name in('trial_nonright','trial_label_at_least_four','owner_foreign')
                facts[tuple(sorted(value))]=4 if name=='trial_label_at_least_four' else 3
        edges={tuple(sorted((nodes[i],nodes[(i+1)%len(nodes)]))):m for i,m in enumerate(row['labels'])}
        assert all(facts.get(edge,0)>=label for edge,label in edges.items())
        # Guard arithmetic is Qsqrt2 or Qsqrt3; this is separate from the
        # Qsqrt2sqrt5 determinant arithmetic below.
        rad=3 if 6 in row['labels'] else 2
        def qa(x,y):return x[0]+y[0],x[1]+y[1]
        def qm(x,y):return x[0]*y[0]+rad*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
        def qs(x,c):return x[0]*c,x[1]*c
        def qparse(text):
            def go(t):
                if isinstance(t,ast.Constant):return F(t.value),F(0)
                if isinstance(t,ast.UnaryOp):assert isinstance(t.op,ast.USub);return qs(go(t.operand),-1)
                if isinstance(t,ast.Call):assert t.func.id=='sqrt'and t.args[0].value==rad;return F(0),F(1)
                assert isinstance(t,ast.BinOp)
                if isinstance(t.op,ast.Add):return qa(go(t.left),go(t.right))
                if isinstance(t.op,ast.Sub):return qa(go(t.left),qs(go(t.right),-1))
                if isinstance(t.op,ast.Mult):return qm(go(t.left),go(t.right))
                assert isinstance(t.op,ast.Div)and isinstance(t.right,ast.Constant)
                return qs(go(t.left),F(1,t.right.value))
            return go(ast.parse(text,mode='eval').body)
        cos={3:(F(1,2),F(0)),4:(F(0),F(1,2)),6:(F(0),F(1,2))}
        G=[[(F(1),F(0))if i==j else qs(cos[edges[tuple(sorted((a,b)))]],-1)
            if tuple(sorted((a,b)))in edges else(F(0),F(0))for j,b in enumerate(nodes)]for i,a in enumerate(nodes)]
        assert G==[[qparse(x)for x in r]for r in row['gram']]
        w=list(map(qparse,row['positive_kernel']));assert all(R.sign_quad(*v,rad)>0 for v in w)
        for r in G:
            v=(F(0),F(0))
            for a,b in zip(r,w):v=qa(v,qm(a,b))
            assert v==(0,0)
        guard_counts[row['kind']]+=1
    budgets={11:(12,13),17:(12,13),14:(3,4)}
    for row in d['projection_budgets']:
        leaf=row['leaf'];drop,keep=budgets.pop(leaf)
        assert(row['removed_own_root'],row['retained_own_root'])==(drop,keep)
        assert set(row['support'])==(set(S)-{drop})|{leaf}and len(row['support'])==8
        assert set(row['support'])<=V[row['actual_vertex']]
        if leaf!=14:assert F(row['own_cost_lower'])==F(1,4)
    assert not budgets
    assert R.sign(R.sub(R.mul(R.COS[5],R.COS[5]),k(F(1,2))))>0
    assert F(1,4)+2*F(1,2)>1 and F(1,2)+2*F(1,4)==1

    # Independent Chebyshev recurrence, root identification and strict tail.
    def ra(a,b):
        c=[F(0)]*max(len(a),len(b))
        for i,x in enumerate(a):c[i]+=x
        for i,x in enumerate(b):c[i]+=x
        return c
    def rm(a,b):
        c=[F(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        return c
    t0,t1=[F(1)],[F(0),F(1)]
    for _ in range(2,8):t0,t1=t1,ra(rm([F(0),F(2)],t1),[-x for x in t0])
    cubic=list(map(F,[-1,-4,4,8]))
    assert ra(t1,[-F(1)])==rm([-F(1),F(1)],rm(cubic,cubic))
    at=F(31,50);value=sum(x*at**i for i,x in enumerate(cubic))
    assert value==F(d['tail']['value_at_test'])and value<0
    assert 24*F(1,2)**2+8*F(1,2)-4>0 and 48*F(1,2)+8>0
    assert(1+at)/2==F(81,100)and 1/(1-F(81,100))==F(100,19)

    one={(0,0,0):O};beta={(1,0,0):O};c={(0,1,0):O};D={(0,0,1):O}
    den=pa(one,ps(pp(c,2),k(-1)))
    patterns={(foreign,dx,dy)for foreign in(None,0,1)for dx,dy in product((F(0),F(1,2)),repeat=2)
              if not(dx and dy)and not(foreign==1 and(dx or dy))}
    assert len(patterns)==7 and d['residual_patterns']==7
    expectedcases={(f,dx,dy,m)for f,dx,dy in patterns for m in(['all finite >=4']if f is None else[4,5])}
    got=set();checks=[];coeffs=0
    for row in d['continuous_cases']:
        foreign=row['foreign_owner_root'];dx,dy=map(F,row['third_anchor_cosines']);internal=row['owner_internal']
        key=(foreign,dx,dy,internal);assert key in expectedcases and key not in got;got.add(key)
        aa=pm(den,pa(one,ps(D,k(F(-1,4)))))
        ab=ps(pm(den,D),k(F(-1,4)))
        at=ps(den,k(-dx));bt=ps(den,k(-dy))
        if foreign is not None:at=pa(at,ps(pm(beta,c)if foreign==0 else beta,k(F(-1,2))))
        G=[[pa(aa,ps(one,k(F(-1,4))))if foreign is not None else aa,ab,at],
           [ab,aa,bt],[at,bt,pa(den,ps(pp(beta,2),k(-1)))]]
        P=divide_den_squared(detpoly(G))
        assert P==parse_poly(row['scaled_determinant'],one,{'beta':beta,'c':c,'D':D})
        assert max(ex[2]for ex in P)==1
        if foreign is None:
            result=shifted(P,[R.COS[5],R.COS[4],k(F(100,19))]);variables=['x','y','z'];count=10
        else:
            reduced={}
            for(exb,exc,exD),v in P.items():
                reduced=pa(reduced,{(exb,exD):R.mul(v,fpow(R.COS[internal],exc))})
            result=shifted(reduced,[R.COS[5],k(F(100,19))]);variables=['x','z'];count=6
        saved={tuple(r['exponents']):dec(r['coefficient'])for r in row['shifted_terms']}
        assert len(saved)==len(row['shifted_terms'])==count and saved==result and row['variable_order']==variables
        assert all(R.sign(x)>0 for x in result.values())
        assert result[(0,)*len(variables)]==dec(row['shifted_constant'])
        coeffs+=count;checks.append(dict(foreign_owner_root=foreign,third_cosines=list(map(str,(dx,dy))),
          owner_internal=internal,independent_polynomial=[dict(exponents=list(ex),coefficient=enc(v))for ex,v in sorted(P.items())],
          shifted_terms=[dict(exponents=list(ex),coefficient=enc(v))for ex,v in sorted(result.items())]))
    assert got==expectedcases and len(got)==11 and coeffs==d['positive_coefficient_count']==78

    fixed=[[3,14,2],[11,12,2],[11,13,3],[12,17,2],[13,17,3]]
    bounds=[[3,4,4],[6,7,7],[9,10,7],[12,13,7],[4,14,5]]
    assert len(d['families'])==len(d['literal_orientations'])==4
    for i,(owner,last)in enumerate(product((False,True),repeat=2)):
        perm=list(range(18))
        if owner:perm[3],perm[4]=4,3
        if last:perm[12],perm[13]=13,12
        assert perm==d['literal_orientations'][i]
        assert{frozenset(perm[x]for x in v)for v in V}==set(map(frozenset,V))
        assert{frozenset(perm[x]for x in m)for m in MN}==MN
        for key,rows in(('fixed',fixed),('bounds',bounds)):
            assert sorted([sorted((perm[a],perm[b]))+[m]for a,b,m in rows])==d['families'][i][key]
    batch=read('results/fresh_h9_star_third34_20261014_search_batch_012.json')
    labels={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
    assert any(all(labels[a,b]==m for a,b,m in f['fixed'])and all(labels[a,b]>=m for a,b,m in f['bounds'])for f in d['families'])
    assert batch['cumulative_additional_nodes']==d['covered_saved_cursor']==26241
    paths=[f'results/{TAG}.json',f'scripts/{TAG}.py',f'scripts/{TAG}_review.py','scripts/fresh_h9_star_lastblock_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_SINGLE_THIRD_LITERAL_DOMAINS_AND_POLYNOMIALS',sources={p:sha(p)for p in paths},
      arithmetic='No producer/SymPy imports. Independent Fraction Qsqrt2sqrt5 operations and sparse polynomial permutation determinants; separate Qsqrt2/Qsqrt3 kernel arithmetic.',
      recursive_manifests=len(seen),vertices=378,core=CORE,core_maximal_faces=9,
      literal_guard_counts=dict(guard_counts),literal_guards=48,projection_budgets=3,
      rational_tail=dict(cos2pi7_test='31/50',cubic_value=str(value),D_strict_lower='100/19'),
      independent_cases=checks,polynomials=11,positive_coefficients=78,literal_families=4,covered_cursor=26241,
      scope='Owner finite internal>=4, otherthree>=7; third own(2,k),k>=5,eitherorientation; shared singlepi3 anchors. All noncore contacts free. Full new guard and polynomial arithmetic reviewed; no prior bank replay and no whole-incidence exclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['polynomials'],d['positive_coefficients'])
