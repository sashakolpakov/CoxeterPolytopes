#!/usr/bin/env python3
"""Producer-free Fraction/Qsqrt11 replay of all32 bridge projections."""
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE='fresh_h22_coupled_cusp_20261013'
SOURCE=ROOT/'results'/f'{PRE}.json'
OUT=ROOT/'results'/f'{PRE}_independent.json'
ZERO=(F(0),F(0))
def lift(x):return x if isinstance(x,tuple)else(F(x),F(0))
def add(x,y):x,y=lift(x),lift(y);return x[0]+y[0],x[1]+y[1]
def neg(x):x=lift(x);return-x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    x,y=lift(x),lift(y);return x[0]*y[0]+11*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def sign(x):
    a,b=lift(x)
    if not b:return(a>0)-(a<0)
    if not a:return(b>0)-(b<0)
    if a*b>0:return(a>0)-(a<0)
    norm=a*a-11*b*b;assert norm
    return((a>0)-(a<0))*((norm>0)-(norm<0))
def dot(x,y):
    out=ZERO
    for a,b in zip(x,y,strict=True):out=add(out,mul(a,b))
    return out
def mv(M,v):return tuple(dot(row,v)for row in M)
def inv(M):
    n=len(M);A=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
    for i in range(n):
        j=next(j for j in range(i,n)if A[j][i]);A[i],A[j]=A[j],A[i]
        p=A[i][i];A[i]=[v/p for v in A[i]]
        for j in range(n):
            if i!=j:
                c=A[j][i];A[j]=[x-c*y for x,y in zip(A[j],A[i])]
    return tuple(tuple(row[n:])for row in A)
def ldl(M):
    A=[list(row)for row in M];out=[]
    for k in range(len(A)):
        p=A[k][k];assert p;out.append(p)
        for i in range(k+1,len(A)):
            for j in range(k+1,len(A)):A[i][j]-=A[i][k]*A[k][j]/p
    return out
def enc(x):return list(map(str,lift(x)))
def dec(x):return tuple(map(F,x))


def check():
    d=json.loads(SOURCE.read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    raw=json.loads((ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input']
    q=tuple(tuple(-F(x)for x in row)for row in raw['basis_gram'])
    pivots=ldl(q);assert all(p>0 for p in pivots[:-1])and pivots[-1]<0
    Q=tuple(tuple(row)+(F(0),)for row in q)+((F(0),)*22+(F(11),),)
    walls=[tuple(map(F,row))+(F(0),)for row in raw['wall_coordinates_as_columns_transposed']]
    assert len(walls)==210 and len(Q)==23
    pair=lambda x,y:dot(x,mv(Q,y))
    loc=json.loads((ROOT/'results/fresh_h22_lower_new_cusp_20261013.json').read_text())['representatives'][0]
    delta=tuple(map(F,loc['null_ray']));assert pair(delta,delta)==ZERO and delta[-1]==2
    dw=[pair(delta,w)for w in walls];assert all(x[1]==0 and sign(x)<=0 for x in dw)
    assert[i for i,x in enumerate(dw)if x==ZERO]==loc['incident_old_wall_indices']
    chains=[]
    for comp in loc['components']:
        if not comp['type'].startswith('finite'):continue
        ids=[int(s[3:])for s in comp['walls']]
        adj={i:[j for j in ids if j!=i and pair(walls[i],walls[j])==lift(-1)]for i in ids}
        assert all(pair(walls[i],walls[j])in(lift(2),lift(-1),ZERO)for i in ids for j in ids)
        if len(ids)==1:order=ids
        else:
            ends=sorted(i for i in ids if len(adj[i])==1);assert len(ends)==2
            order=[ends[0]]
            while len(order)<len(ids):
                nxt=[j for j in adj[order[-1]]if j not in order];assert len(nxt)==1;order+=nxt
        chains.append(order)
    chains.sort(key=lambda row:(len(row),row[0]));assert sorted(map(len,chains))==[1,3,3,4,5]
    basis=sum(chains,[])+[26];assert len(set(basis))==17 and all(dw[j]==ZERO for j in basis)
    B=tuple(tuple(pair(walls[i],walls[j])[0]for j in basis)for i in basis)
    assert all(p>0 for p in ldl(B));BI=inv(B)
    assert all(B[-1][i]==0 for i in range(16))and B[-1][-1]==2
    e=json.loads((ROOT/'results/fresh_h22_new_cusp_completion_20261013_lower_escape.json').read_text())
    X=tuple(lift(F(v,e['denominator']))for v in e['numerator'][:-1])+((F(0),F(e['numerator'][-1],11*e['denominator'])),)
    xnorm=pair(X,X);dx=pair(delta,X)
    assert xnorm==lift(F(e['square']))and sign(xnorm)>0 and sign(dx)<0
    assert all(sign(pair(w,X))<0 for w in walls)
    expected={tuple(sorted((a,b)))for i,j in combinations(range(5),2)
              for a in set((chains[i][0],chains[i][-1]))for b in set((chains[j][0],chains[j][-1]))}
    assert len(expected)==32 and{tuple(row['endpoints'])for row in d['certificates']}==expected
    rows=[]
    for row in d['certificates']:
        endpoints=set(row['endpoints']);p=tuple(F(-1 if j in endpoints else 0)for j in basis)
        rho=sub(2,dot(p,mv(BI,p)));assert rho[1]==0 and sign(rho)>0
        weights=[(int(j),F(w))for j,w in row['weights']]
        assert all(w>0 and sign(dw[j])<0 for j,w in weights)and sum(w for j,w in weights)==1
        lam=[mul(dx,w/dw[j][0])for j,w in weights];assert all(sign(v)>0 for v in lam)
        z=[]
        for k in range(23):
            v=X[k]
            for (j,_),ll in zip(weights,lam):v=sub(v,mul(ll,walls[j][k]))
            z.append(v)
        assert pair(z,delta)==ZERO
        zp=tuple(pair(walls[j],z)for j in basis);bz=mv(BI,zp)
        kappa=dot(p,bz);zeta=sub(pair(z,z),dot(zp,bz));margin=sub(mul(kappa,kappa),mul(rho,zeta))
        assert sign(kappa)<0 and sign(zeta)>=0 and sign(margin)>0
        assert lam==list(map(dec,row['positive_multipliers']))and rho==lift(F(row['residual_norm']))
        assert kappa==dec(row['kappa'])and zeta==dec(row['zeta'])and margin==dec(row['kappa_squared_minus_rho_zeta'])
        rows.append({'endpoints':row['endpoints'],'multipliers':list(map(enc,lam)),
                     'rho':enc(rho),'kappa':enc(kappa),'zeta':enc(zeta),'margin':enc(margin)})
    return {'status':'PASS_INDEPENDENT_ALL32_COUPLED_BRIDGE_PROJECTIONS',
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'arithmetic':'Standalone Fraction and Q(sqrt11) operations, independent rational Gaussian inverse and LDL. No producer imports or SymPy.',
        'pairs':32,'old_target_inequalities':210,'positive_projection_rank':17,'rebuilt_chains':chains,
        'target_norm':enc(xnorm),'delta_target_pairing':enc(dx),'certificates':rows,
        'scope':'For the fixed delta4 and prescribed two-chain endpoint pairings, every old-compatible real norm-two normal misses the exact target, including all null shifts and residual directions. Only finite families of these bridges are excluded from eliminating the strict escape. Other cusp orbits, different incidence profiles and removed old walls are not excluded.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=check()
    if args.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
