#!/usr/bin/env python3
"""Independent Fraction/Qsqrt2 polynomial and literal-domain review."""
import argparse
import ast
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_third34_20261014'
S=[3,4,6,7,9,10,12,13];CORE=sorted(S+[11,14,17])
Z=(F(0),F(0));ONE=(F(1),F(0));RT=(F(0),F(1))
def k(a=0,b=0):return F(a),F(b)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],-a[1]
def sub(a,b):return add(a,neg(b))
def mul(a,b):return a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def sc(a,c):return a[0]*c,a[1]*c
def power(a,n):
    q=ONE
    for _ in range(n):q=mul(q,a)
    return q
def sign(v):
    a,b=v
    if not a:return(b>0)-(b<0)
    if not b:return(a>0)-(a<0)
    if a*b>0:return(a>0)-(a<0)
    d=a*a-2*b*b;assert d
    return((a>0)-(a<0))*((d>0)-(d<0))
def dec(v):return tuple(map(F,v))
def enc(v):return list(map(str,v))
def total(v):
    out=Z
    for q in v:out=add(out,q)
    return out
def clean(p):return{m:v for m,v in p.items()if v!=Z}
def pa(p,q):
    out=dict(p)
    for m,v in q.items():out[m]=add(out.get(m,Z),v)
    return clean(out)
def ps(p,c):return clean({m:mul(v,c)for m,v in p.items()})
def pm(p,q):
    out={}
    for m,v in p.items():
        for n,w in q.items():
            e=tuple(a+b for a,b in zip(m,n));out[e]=add(out.get(e,Z),mul(v,w))
    return clean(out)
def pp(p,n):
    if not p:
        assert n>0
        return {}
    out={tuple(0 for _ in next(iter(p))):ONE}
    for _ in range(n):out=pm(out,p)
    return out
def shifted(p,shift):
    out={};zero=tuple(0 for _ in shift)
    for ex,c in p.items():
        term={zero:c}
        for i,n in enumerate(ex):
            f={}
            for j in range(n+1):
                e=list(zero);e[i]=j;f[tuple(e)]=sc(power(shift[i],n-j),comb(n,j))
            term=pm(term,f)
        out=pa(out,term)
    return out
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

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
    assert d['status']=='PASS_ALL_LABEL_DENSE_THIRD_SINGLE_ANCHOR_CORE'
    g=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[set(r['incident'])for r in g['vertices']];missing=set(map(frozenset,g['minimal_nonfaces']))
    assert len(V)==378 and d['rank_core']==CORE and any(set(S)<=v for v in V)
    assert all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
    faces={frozenset(v&set(CORE))for v in V}
    maximal={f for f in faces if not any(f<t for t in faces)}
    assert len(maximal)==9 and maximal==set(map(frozenset,d['complete_maximal_core_faces']))
    bounds={(3,4):4,(6,7):4,(9,10):4,(12,13):7,(3,14):3,(4,14):4}
    cos={3:k(F(1,2)),4:k(0,F(1,2))}
    cross={tuple(sorted((a,b)))for i,j in combinations(range(4),2)
           for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2])}
    gotcross=set();uppers=set();conditionals=set();guards=[]
    for row in d['positive_kernel_guards']:
        nodes=row['ordered_support'];kind=row['kind'];trial=tuple(row['test_pair'])
        assert len(set(nodes))==len(nodes) and set(nodes)==set(row['support'])<=V[row['actual_vertex']]
        assert set(nodes)<=set(CORE)
        edges={tuple(sorted((nodes[i],nodes[i+1]))):m for i,m in enumerate(row['labels'])}
        condition=None
        if kind=='cross_block_right':
            assert trial in cross and row['labels']==[4,3,4]and not row['cycle'];gotcross.add(trial)
            facts=bounds|{trial:3}
        elif kind=='third_anchor_label_at_most_three':
            assert nodes in([4,14,11],[4,14,17])and row['labels']==[4,4]and not row['cycle']
            assert trial==tuple(sorted(nodes[1:]));uppers.add(trial);facts=bounds|{trial:4}
        else:
            assert kind=='owner_contact_forces_third_contact_right'
            root,third,leaf=nodes;assert root in(3,4)and third==14 and leaf in(11,17)
            assert row['labels']==[3,3]and row['cycle']and trial==tuple(sorted((third,leaf)))
            condition=tuple(sorted((root,leaf)));conditionals.add(condition)
            edges[condition]=3;facts=bounds|{trial:3,condition:3}
        assert all(facts.get(p,0)>=m for p,m in edges.items())
        A=[[ONE if a==b else neg(cos[edges[tuple(sorted((a,b)))]])
            if tuple(sorted((a,b)))in edges else Z for b in nodes]for a in nodes]
        assert A==[[dec(v)for v in r]for r in row['gram']]
        w=list(map(dec,row['positive_kernel']));assert all(sign(v)>0 for v in w)
        assert all(total(mul(a,b)for a,b in zip(r,w))==Z for r in A)
        guards.append(dict(kind=kind,support=nodes,actual_vertex=row['actual_vertex'],
                           conditional_nonright_pair=condition,kernel=list(map(enc,w))))
    assert gotcross==cross and uppers=={(11,14),(14,17)}
    assert conditionals=={(3,11),(4,11),(3,17),(4,17)}and len(guards)==30
    budgets={(11,12,13):F(1,4),(17,12,13):F(1,4),(17,13,12):F(1,4),(14,3,4):F(1,2)}
    for row in d['projection_budgets']:
        key=row['leaf'],row['removed_own_root'],row['retained_own_root']
        assert F(row['own_cost_lower'])==budgets.pop(key)
        support=(set(S)-{key[1]})|{key[0]}
        assert len(support)==8 and support==set(row['support'])<=V[row['actual_vertex']]
    assert not budgets
    assert F(1,4)+2*F(1,2)>1 # two foreign pi3 contacts or one foreign >=4
    assert F(1,2)+2*F(1,4)==1 # any third foreign contact fails strict PD
    assert 2*(F(1,4)+F(1,2))>1 # X>1 even without positive cross term
    assert F(3)>F(3,2)**2 # 2+sqrt3>7/2

    # Three-variable physical polynomial reconstruction, no symbolic engine.
    one={(0,0,0):ONE};a={(1,0,0):ONE};b={(0,1,0):ONE};c={(0,0,1):ONE}
    den=pa(one,ps(pp(c,2),k(-1)))
    NX=pa(pa(pp(a,2),pp(b,2)),ps(pm(pm(a,b),c),k(2)))
    multiplicities=Counter((v.count(0),v.count(1))for v in product((None,0,1),repeat=2))
    case_results=[];got=set()
    for row in d['continuous_cases']:
        n0,n1=row['owner_contact_counts'];assert(n0,n1)not in got;got.add((n0,n1))
        assert row['anchor_assignment_multiplicity']==multiplicities[n0,n1]
        delta=F(2-n0-n1,2);assert F(row['ambient_contact_sum_upper'])==delta
        NF=pa(ps(one,k(F(n0*n0+n1*n1,4))),ps(c,k(F(n0*n1,2))))
        NT=ps(pa(pa(ps(a,k(n0)),ps(pm(c,b),k(n0))),
                 pa(ps(b,k(n1)),ps(pm(c,a),k(n1)))),k(F(1,2)))
        wedge=ps(pp(pa(ps(a,k(n1)),ps(b,k(-n0))),2),k(F(1,4)))
        assert pa(pm(NX,NF),ps(pp(NT,2),k(-1)))==pm(den,wedge)
        P=pa(pa(ps(pa(NX,ps(den,k(-1))),k(F(3,2))),wedge),
             pa(ps(NF,k(-1)),pa(ps(NT,k(-2*delta)),ps(den,k(-delta*delta)))))
        # Independent AST parser checks the saved unshifted expression too.
        def parse(t):
            if isinstance(t,ast.Name):return {'alpha':a,'beta':b,'c':c}[t.id]
            if isinstance(t,ast.Constant):return ps(one,k(t.value))
            if isinstance(t,ast.UnaryOp)and isinstance(t.op,ast.USub):return ps(parse(t.operand),k(-1))
            assert isinstance(t,ast.BinOp)
            if isinstance(t.op,ast.Add):return pa(parse(t.left),parse(t.right))
            if isinstance(t.op,ast.Sub):return pa(parse(t.left),ps(parse(t.right),k(-1)))
            if isinstance(t.op,ast.Mult):return pm(parse(t.left),parse(t.right))
            if isinstance(t.op,ast.Pow):return pp(parse(t.left),t.right.value)
            assert isinstance(t.op,ast.Div)and isinstance(t.right,ast.Constant)
            return ps(parse(t.left),k(F(1,t.right.value)))
        assert parse(ast.parse(row['polynomial'],mode='eval').body)==P
        shiftedP=shifted(P,[k(F(1,2)),k(0,F(1,2)),k(0,F(1,2))])
        saved={tuple(t['exponents']):dec(t['coefficient'])for t in row['shifted_terms']}
        assert len(saved)==len(row['shifted_terms'])==11 and saved==shiftedP
        assert all(sign(v)>0 for v in shiftedP.values())
        assert shiftedP[0,0,0]==dec(row['shifted_constant'])and sign(shiftedP[0,0,0])>0
        case_results.append(dict(owner_counts=[n0,n1],multiplicity=multiplicities[n0,n1],
          unshifted=[dict(exponents=list(m),coefficient=enc(v))for m,v in sorted(P.items())],
          shifted=[dict(exponents=list(m),coefficient=enc(v))for m,v in sorted(shiftedP.items())]))
    assert got==set(multiplicities)and len(got)==6 and sum(multiplicities.values())==9
    assert d['positive_coefficient_count']==66

    # Eight formal variables X,F,T,Y,M,mu,e,delta: verify the determinant
    # difference without the producer's rational-expression simplifier.
    n=8;u={tuple([0]*n):ONE}
    var=[]
    for j in range(n):
        ex=[0]*n;ex[j]=1;var.append({tuple(ex):ONE})
    X,FF,T,Y,MM,MU,E,DEL=var
    xm=pa(X,ps(u,k(-1)))
    residual=pa(pm(xm,pa(pa(Y,FF),pa(MM,pa(ps(MU,k(2)),ps(u,k(-2)))))),
                ps(pp(pa(E,T),2),k(-1)))
    lower=pa(pm(xm,pa(FF,ps(u,k(F(3,2))))),ps(pp(pa(DEL,T),2),k(-1)))
    target=pa(pm(xm,pa(pa(Y,ps(u,k(F(-7,2)))),pa(MM,ps(MU,k(2))))),
              pm(pa(DEL,ps(E,k(-1))),pa(pa(DEL,E),ps(T,k(2)))))
    assert pa(residual,ps(lower,k(-1)))==target

    expected=[]
    for opposite,owner_reverse,last_reverse in product((False,True),repeat=3):
        p=list(range(18))
        if owner_reverse:p[3],p[4]=4,3
        if last_reverse:p[12],p[13]=13,12
        assert{frozenset(p[i]for i in v)for v in V}==set(map(frozenset,V))
        assert{frozenset(p[i]for i in v)for v in missing}==missing
        fixed={(11,12):2,(11,13):3,(12,17):3 if opposite else 2,(13,17):2 if opposite else 3}
        item={key:sorted([sorted((p[i],p[j]))+[m]for(i,j),m in vals.items()])
              for key,vals in [('fixed',fixed),('bounds',bounds)]}
        expected.append((p,item))
    assert len(d['families'])==len(d['literal_orientations'])==8
    for(p,item),perm,family in zip(expected,d['literal_orientations'],d['families']):
        assert p==perm and all(family[key]==value for key,value in item.items())
    batch=read('results/fresh_h9_star_core_transfers_20261014_search_batch_011.json')
    labels={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
    assert any(all(labels[i,j]==m for i,j,m in f['fixed'])and all(labels[i,j]>=m for i,j,m in f['bounds'])for f in d['families'])
    assert batch['cumulative_additional_nodes']==d['covered_saved_cursor']==20837
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_THIRD34_CONTINUOUS_CORE_REVIEW',sources={p:sha(p)for p in paths},
      arithmetic='Standalone Fraction Qsqrt2 arithmetic and sparse polynomials; no producer, SymPy or prior-bank arithmetic imports.',
      recursive_manifests=len(seen),actual_vertex_count=len(V),core=CORE,
      guard_checks=guards,literal_budgets=4,core_maximal_faces=9,
      continuous_cases=case_results,positive_coefficients=66,determinant_difference_identity=True,
      literal_families=8,covered_cursor=20837,
      scope='All finite internal labels with first3>=4,last>=7; third own>=3,>=4 in eitherorientation; shared/opposite single ownpi3 anchors. All noncore contacts free. Six polynomials and all new domains independently reconstructed; no incidence-wide exclusion or new label search.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['positive_coefficients'])
