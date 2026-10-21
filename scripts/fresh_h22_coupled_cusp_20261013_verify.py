#!/usr/bin/env python3
"""Exact all-shift/all-direction obstruction for coupled affine-A20 bridges.

No numerical optimization enters verification. Each of the 32 endpoint pairs
has a dual certificate supported on one or two of the original 210 walls.
All algebraic inequalities are checked in Q(sqrt(11)) by rational squaring.
"""
import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import sympy as sp

P=Path(__file__).resolve().parents[1]
PRE='fresh_h22_coupled_cusp_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
LOCAL='results/fresh_h22_lower_new_cusp_20261013.json'
ESCAPE='results/fresh_h22_new_cusp_completion_20261013_lower_escape.json'
PROPOSALS='data/enumerations/'+PRE+'_dual_proposals.json'
OUT='results/'+PRE+'.json'

@dataclass(frozen=True)
class K:
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
    def __add__(self,y):
        y=lift(y);return K(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,y):return self+-lift(y)
    def __rsub__(self,y):return lift(y)+-self
    def __mul__(self,y):
        y=lift(y);return K(self.a*y.a+11*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=lift(y);d=y.a*y.a-11*y.b*y.b;assert d
        return self*K(y.a/d,-y.b/d)
    def sign(self):
        if self.b==0:return (self.a>0)-(self.a<0)
        if self.a==0:return (self.b>0)-(self.b<0)
        if self.a*self.b>0:return (self.a>0)-(self.a<0)
        d=self.a*self.a-11*self.b*self.b;assert d
        return ((self.a>0)-(self.a<0))*((d>0)-(d<0))
    def json(self):return [str(self.a),str(self.b)]

def lift(x):return x if isinstance(x,K)else K(x)
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(x,y):return sum((a*b for a,b in zip(x,y,strict=True)),K())
def mv(M,x):return [dot(row,x)for row in M]
def rational_matrix(M):return [[F(int(x.p),int(x.q))for x in row]for row in M.tolist()]

def run():
    raw=read(RAW)['coxeter_input'];q=sp.Matrix([[-int(x)for x in row]for row in raw['basis_gram']])
    minors=[q[:j,:j].det()for j in range(1,23)]
    assert all(x>0 for x in minors[:-1])and minors[-1]<0
    # The first 21 pivots are positive, the last negative; added coordinate is positive.
    Q=rational_matrix(sp.diag(q,11));W=[[F(x)for x in r]+[F(0)]for r in raw['wall_coordinates_as_columns_transposed']]
    assert len(W)==210 and len(Q)==23
    loc=read(LOCAL)['representatives'][0];delta=[F(x)for x in loc['null_ray']];dc=mv(Q,delta)
    assert dot(delta,dc)==K()and delta[-1]==2
    dw=[dot(w,dc)for w in W];assert all(x.b==0 and x.sign()<=0 for x in dw)
    assert [i for i,x in enumerate(dw)if x==K()]==loc['incident_old_wall_indices']
    chains=[]
    for comp in loc['components']:
        if not comp['type'].startswith('finite'):continue
        ids=[int(x[3:])for x in comp['walls']]
        G=[[dot(W[i],mv(Q,W[j])).a for j in ids]for i in ids]
        assert all(G[i][i]==2 for i in range(len(ids)))
        assert all(G[i][j]in(-1,0)for i in range(len(ids))for j in range(i))
        if len(ids)==1:ordered=ids
        else:
            ends=[i for i in range(len(ids))if sum(v==-1 for v in G[i])==1];assert len(ends)==2
            ordered=[ids[ends[0]]]
            while len(ordered)<len(ids):
                js=[j for j in ids if j not in ordered and dot(W[ordered[-1]],mv(Q,W[j]))==K(-1)]
                assert len(js)==1;ordered+=js
        chains.append(ordered)
    assert sorted(map(len,chains))==[1,3,3,4,5]
    finite=sum(chains,[]);basis=finite+[26];assert len(set(basis))==17
    B=sp.Matrix([[dot(W[i],mv(Q,W[j])).a for j in basis]for i in basis])
    assert all(B[:j,:j].det()>0 for j in range(1,18))
    Bi=rational_matrix(B.inv());assert all(B[-1,j]==0 for j in range(16))and B[-1,-1]==2
    ex=read(ESCAPE);X=[K(F(x,ex['denominator']))for x in ex['numerator'][:-1]]+[K(0,F(ex['numerator'][-1],11*ex['denominator']))]
    xc=mv(Q,X);xq=dot(X,xc);assert xq==K(F(ex['square']))and xq.sign()>0
    assert all(dot(w,xc).sign()<0 for w in W)
    dx=dot(delta,xc);assert dx.sign()<0
    expected=set()
    for i,j in combinations(range(5),2):
        for a in set((chains[i][0],chains[i][-1])):
            for b in set((chains[j][0],chains[j][-1])):expected.add(tuple(sorted((a,b))))
    assert len(expected)==32
    proposals=read(PROPOSALS)['certificates'];assert len(proposals)==32
    assert {tuple(sorted(r['endpoints']))for r in proposals}==expected
    rows=[]
    for row in sorted(proposals,key=lambda r:sorted(r['endpoints'])):
        ends=tuple(sorted(row['endpoints']));p=[F(-1 if j in ends else 0)for j in basis]
        bp=mv(Bi,p);rho=K(2)-dot(p,bp);assert rho.b==0 and rho.sign()>0
        weights=[(int(j),F(w))for j,w in row['weights']]
        assert all(w>0 and dw[j].sign()<0 for j,w in weights)and sum(w for _,w in weights)==1
        lambdas=[dx*(w/dw[j].a)for j,w in weights];assert all(v.sign()>0 for v in lambdas)
        z=[X[k]-sum((v*W[j][k]for (j,_),v in zip(weights,lambdas,strict=True)),K())for k in range(23)]
        zc=mv(Q,z);assert dot(delta,zc)==K()
        zp=[dot(W[j],zc)for j in basis];bz=mv(Bi,zp)
        kappa=dot(p,bz);zeta=dot(z,zc)-dot(zp,bz)
        margin=kappa*kappa-rho*zeta
        assert kappa.sign()<0 and zeta.sign()>=0 and margin.sign()>0
        rows.append({'endpoints':list(ends),'weights':[[j,str(w)]for j,w in weights],
                     'positive_multipliers':[v.json()for v in lambdas],
                     'residual_norm':str(rho.a),'kappa':kappa.json(),
                     'zeta':zeta.json(),'kappa_squared_minus_rho_zeta':margin.json(),
                     'exact_signs':[-1,zeta.sign(),1]})
    return {'status':'PASS_EXACT_ALL_32_COUPLED_AFFINE_A20_BRIDGES_MISS_ESCAPE',
            'sources':{s:sha(s)for s in(RAW,LOCAL,ESCAPE,PROPOSALS)},
            'field':'a+b*sqrt(11), positive sqrt(11); last coordinate is coefficient of sqrt(11)',
            'ambient_signature':[22,1],'old_form_leading_minors':list(map(str,minors)),
            'chains':chains,'orthogonal_affine_A1_old_wall':26,
            'target_square':xq.json(),'delta_target_pairing':dx.json(),
            'certificates':rows,
            'scope':'Every norm-two new wall perpendicular to delta4 and old26, pairing -1 with endpoints of two distinct stated finite chains and zero with the other finite-chain walls, and nonpositively with the original210 walls, has strictly negative pairing with the current exact escape. The result includes all real shifts and all four-dimensional tangent directions, with no field or dihedral-label restriction on the new normal. Thus any finite union of these bridges, in particular every coherent affine A20 completion of these five chains by five bridges, leaves this escape. It does not exclude other diagram extensions, nonincident caps or changes to retained walls.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();out=run()
    if args.verify:assert out==read(OUT)
    else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print('Exact certificates:',len(out['certificates']))
