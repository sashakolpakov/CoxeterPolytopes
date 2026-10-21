#!/usr/bin/env python3
"""Independent Golden-arithmetic census, reversed coordinate order and cone DD."""
import argparse,json,sys
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.golden import Golden as K
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone
TAG='fresh_h9_beta55_pilot_20261014_review'
PRIMARY='results/fresh_h9_beta55_pilot_20261014.json'
GEOMETRY='results/fresh_h9_beta55_pilot_20261014_geometry.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(v):
    if isinstance(v,K):return v.to_json()
    if isinstance(v,F):return str(v)
    if isinstance(v,(list,tuple)):return [enc(x)for x in v]
    if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
    return v
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def inverse(M):
    n=len(M);a=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(M)]
    for c in range(n):
        j=next(i for i in range(c,n)if a[i][c]);a[c],a[j]=a[j],a[c]
        a[c]=[x/a[c][c]for x in a[c]]
        for i in range(n):
            if i!=c:
                f=a[i][c];a[i]=[x-f*y for x,y in zip(a[i],a[c])]
    return [r[n:]for r in a]
def dot(x,y):return sum((a*b for a,b in zip(x,y)),K())
def mv(M,v):return [dot(r,v)for r in M]
def pair(x,y,G):return dot(x,mv(G,y))
def key(r):return tuple((x.a,x.b)for x in r)
def run():
    d,g=read(PRIMARY),read(GEOMETRY);assert F(d['height_bound'])==20
    phi=K(0,1);delta=1+8*phi;alpha=2+phi;beta=4-phi
    G=[[K.from_json(x)for x in r]for r in d['form']]
    C=[[x.a for x in r[1:9]]for r in G[1:9]];J=inverse(C)
    assert all(x.b==0 for r in G[1:9]for x in r[1:9])
    assert all(x>=0 and x.denominator==1 for r in J for x in r)
    expected=[[K()for _ in range(10)]for _ in range(10)]
    for i in range(8):
        for j in range(8):expected[i+1][j+1]=K(C[i][j])
    expected[0][0]=-2*phi;expected[0][9]=expected[9][0]=K(-1);expected[9][9]=K(4)
    assert G==expected and delta==phi*alpha*beta
    norms=[K(2),2*alpha,2*beta,2*alpha*beta]
    assert [K.from_json(q)for q in d['norm_classes']]==norms
    assert alpha.norm()==5 and beta.norm()==11 and delta.norm()==-55
    assert all(q.sign()>0 and q.sign(1)>0 and q<18 and q.conjugate()<13 for q in norms)
    assert delta>13 and -delta.conjugate()>3
    # Independently larger boxes: x<11, |x'|<5 => |a|,|b|<10;
    # |l|<39, |l'|<8 => |a|,|b|<26. Lambda bounds are smaller.
    assert F(4*18*20,13)<121 and 4*18*21<39*39 and F(4*13,3)<25 and 4*13<64
    assert K(0,2)-1>K(F(11,5)) and phi<K(F(17,10)) and phi-1<K(F(7,10))
    assert F(16,1)/F(11,5)<10 and (F(17,10)*5+F(7,10)*11)/F(11,5)<10
    assert F(47,1)/F(11,5)<26 and (F(17,10)*8+F(7,10)*39)/F(11,5)<26
    controls=[]
    for q in norms:
        for a in range(-10,11):
            for b in range(-10,11):
                x=K(a,b);h=delta*x*x/(4*q)
                if x>0 and h<=20 and (4*q+delta*x*x).sign(1)>=0:controls.append((h,q,x))
    controls.sort(key=lambda r:r[0]);assert len(controls)==58
    assert len({r[0]for r in controls})==len(controls)
    assert [(K.from_json(r['height']),K.from_json(r['norm']),K.from_json(r['control']))for r in d['levels']]==controls
    # Inverse of positive form diag(1,4J), natural lambda order unlike producer.
    Binv=[[F(1)]+[F(0)]*8]+[[F(0)]+[x/4 for x in r]for r in C]
    forms=[inverse([r[:n]for r in Binv[:n]])for n in range(1,10)]
    pool=[K(a,b)for a in range(-26,27)for b in range(-26,27)]
    squares={z:z*z for z in pool}
    N=[[K(-int(i==j))for i in range(10)]for j in range(1,9)]
    assert N==[[K.from_json(x)for x in r]for r in d['zero_simple_walls']]
    total=Counter();levels=[]
    for index,(height,q,x)in enumerate(controls):
        target=4*q+delta*x*x;b=q/2
        divisors=[z for z in pool if (z/b).is_integral]
        signed=[z for z in divisors if (target-squares[z]).sign()>=0 and (target-squares[z]).sign(1)>=0 and ((x-z)/4).is_integral]
        positive=[z for z in divisors if z>=0 and (target-2*squares[z]).sign()>=0 and (target-2*squares[z]).sign(1)>=0]
        assert set(signed)=={K.from_json(z)for z in d['levels'][index]['signed_l_alphabet']}
        assert set(positive)=={K.from_json(z)for z in d['levels'][index]['lambda_alphabet']}
        prior=N[8:];rows=[[r[0]-4*r[9]]+[4*r[i]for i in range(1,9)]for r in prior]
        bounds=[delta*x*r[0]for r in prior]
        assert all(z>=0 for r in rows for z in r[1:])
        accepted=[];count=Counter()
        def visit(v,partials):
            n=len(v);count['nodes']+=1
            if n:
                M=forms[n-1];minimum=sum((v[i]*v[j]*M[i][j]for i in range(n)for j in range(n)),K())
                if (target-minimum).sign()<0 or (target-minimum).sign(1)<0:count['schur']+=1;return
                partials=[s+v[-1]*r[n-1]for s,r in zip(partials,rows)]
                if any(s>bound for s,bound in zip(partials,bounds)):count['wall']+=1;return
            if n==9:
                if minimum!=target:return
                y=mv(J,v[1:]);r=[x,*y,(x-v[0])/4]
                assert all(z.is_integral for z in r)and pair(r,r,G)==q
                assert all((z/b).is_integral for z in mv(G,r))
                assert all(pair(r,t,G)<=0 for t in N)
                accepted.append(r);return
            for z in (signed if n==0 else positive):visit(v+[z],partials)
        visit([], [K()for _ in rows])
        assert {key(r)for r in accepted}=={key([K.from_json(x)for x in r])for r in d['levels'][index]['roots']}
        assert all(pair(r,s,G)<=0 for i,r in enumerate(accepted)for s in accepted[i+1:])
        N+=accepted;total.update(count);levels.append(dict(index=index,nodes=count['nodes'],roots=accepted))
    assert {key(r)for r in N}=={key([K.from_json(x)for x in r])for r in d['roots']}
    H=[K.from_json(x)for x in d['strict_center']]
    assert pair(H,H,G)<0 and all(pair(H,r,G)<0 for r in N)
    for r,p in zip([[K.from_json(x)for x in r]for r in d['roots']],d['actual_facet_points']):
        v=[K.from_json(x)for x in p];assert pair(v,v,G)<0 and pair(v,H,G)<0
        assert all(pair(v,n,G).sign()==(0 if n==r else -1)for n in N)
    # Match producer ordering only for incidence comparison, never for enumeration.
    N=[[K.from_json(x)for x in r]for r in d['roots']]
    cone=diagnose_chamber_cone(N,G,H,9)
    actual={'finite':set(cone.finite_vertices),'ideal':set(cone.cusps),'hyperideal':set(cone.hyperideal_rays)}
    for kind,incidences in actual.items():assert incidences=={tuple(r['incident'])for r in g['all_extreme_rays']if r['kind']==kind}
    assert not cone.opposite_sheet_rays
    E=[K.from_json(x)for x in g['strict_escape']['point']]
    assert pair(E,E,G)>0 and pair(E,H,G)<0 and all(pair(E,r,G)<0 for r in N)
    paths=[PRIMARY,GEOMETRY,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/golden.py','src/coxeter_polytopes/verify/chamber.py']
    return enc(dict(status='PASS_INDEPENDENT_BETA55_T20_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Project Golden arithmetic, no primary or primary-arithmetic helper imports.',
        coordinate_order=['l']+[f'lambda{i}'for i in range(1,9)],controls=len(controls),counts=dict(total),
        levels=levels,geometry={k:len(v)for k,v in actual.items()},
        scope='Independent complete bounded enumeration and complete cone replay. Ring/unit and zero-sector proofs are reviewed mathematically; no host reflectivity conclusion.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['controls'],d['counts'],d['geometry'])
