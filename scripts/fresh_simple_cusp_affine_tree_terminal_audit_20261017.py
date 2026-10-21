#!/usr/bin/env python3
"""Independent physical-matrix audit of terminal and exceptional affine-tree cases.

No primary imports. All-rank coverage is a written classification/network proof;
finite matrix samples confirm its formulas without serving as rank cutoffs.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_tree_terminal_audit_20261017'
Q=S.Rational;T=S.sqrt(2)

def clean(x):return S.simplify(x)
def eq(x,y):assert clean(x-y)==0,(x,y)
def wire(x):return str(clean(x))
def mat(n,edges):
    G=2*S.eye(n)
    for i,j,w in edges:G[i,j]=G[j,i]=-w
    return G

def pd(G):
    # Sylvester's criterion, distinct from Schur pivots in the primary.
    return all(G[:i,:i].det()>0 for i in range(1,G.rows+1))

def affine(kind,n):
    nu=S.ones(n,1)
    if kind=='C':
        E=[(i,i+1,T if i in (0,n-2) else 1) for i in range(n-1)]
        nu[0]=nu[n-1]=T/2
    elif kind=='B':
        E=[(i,i+1,T if i==n-4 else 1) for i in range(n-3)]+[(0,n-2,1),(0,n-1,1)]
        nu[n-3]=T/2;nu[n-2]=nu[n-1]=Q(1,2)
    else:
        assert kind=='D'
        E=[(i,i+1,1) for i in range(n-5)]+[(0,n-4,1),(0,n-3,1),(n-5,n-2,1),(n-5,n-1,1)]
        for i in range(n-4,n):nu[i]=Q(1,2)
    G=mat(n,E);assert G*nu==S.zeros(n,1)
    D=S.diag(*list(nu));L=D*G*D
    assert all(clean(sum(L[i,j] for j in range(n)))==0 for i in range(n))
    return G,nu,E

def context(kind,n,a,b):
    G,nu,E=affine(kind,n);W=[i for i in range(n) if i not in (a,b)]
    H=G.extract(W,W);J=H.inv();assert H*J==S.eye(n-2)
    va=-G.extract([a],W).T;vb=-G.extract([b],W).T
    alpha=clean(nu[a]**2*(2-(va.T*J*va)[0]))
    eq(alpha,nu[b]**2*(2-(vb.T*J*vb)[0]));eq(alpha,nu[a]*nu[b]*(va.T*J*vb)[0])
    def profile(contacts):
        w=S.Matrix([contacts.get(i,0) for i in W])
        B=clean(nu[b]*(vb.T*J*w)[0]);C=clean(2-(w.T*J*w)[0]);l=clean(B/alpha)
        return clean(C-l*l*alpha),l,clean((nu.extract(W,[0]).T*w)[0]),nu[a]
    return G,nu,W,J,alpha,profile


def symbolic_det(alpha,A,B,C,D,E):
    h,k,u,v=S.symbols('h k u v',positive=True)
    G=S.Matrix([[alpha,-alpha,-h+B,-A],[-alpha,alpha,-B,-k+A],[-h+B,-B,C,-E],[-A,-k+A,-E,D]])
    f=S.factor(G.det()/(h*h*k*k)).subs({h:1/u,k:1/v}).expand()
    return S.factor(f),u,v

def run():
    records=[]
    for n in (18,19,26):
        for kind,a,b,port,nxt,fork in [('C',0,6,1,2,False),('B',n-3,n-9,n-4,n-5,False),
                  ('B',n-2,4,n-1,0,True),('D',n-4,4,n-3,0,True)]:
            G,nu,W,J,alpha,p=context(kind,n,a,b)
            eq(alpha,Q(1,6))
            offsets=[0,Q(1,2),T/2,1] if fork else [0,1,T,1]
            scale=Q(1,2) if fork else T/2
            profiles=[{}, {port:1},{port:T},{nxt:1}]
            for i,contacts in enumerate(profiles):
                data=p(contacts);expect=(S.Integer([2,1,0,0][i]),[0,1,T,2][i],offsets[i],scale)
                for x,y in zip(data,expect):eq(x,y)
                records.append(dict(kind=kind,order=n,marked=[a,b],profile=i,contacts=[[j,wire(v)] for j,v in contacts.items()],
                    alpha=wire(alpha),c0=wire(data[0]),moment=wire(data[1]),height_offset=wire(data[2]),opposite_contact_scale=wire(data[3])))
    # Literal exceptional white blocks, independent of any full-tree helper.
    D3=mat(3,[(0,1,1),(0,2,1)])
    A3=mat(3,[(0,1,1),(1,2,1)])
    W=S.diag(D3,A3);J=W.inv();assert pd(W)
    va=S.Matrix([1,0,0,1,0,0]);vb=S.Matrix([0,0,0,0,0,1]);wd=S.Matrix([0,0,1,0,0,0])
    eq(1-(vb.T*J*vb)[0],Q(1,4)) # far outer end has unit cost
    eq(2-(va.T*J*va)[0],Q(1,4));eq((va.T*J*vb)[0],Q(1,4))
    A=clean((va.T*J*wd)[0]);D=clean(2-(wd.T*J*wd)[0]);eq(A,Q(1,2));eq(D,Q(5,4))
    main=[]
    expected=[('isolated',S.Matrix([0,0,0,0,0,0]),Q(1,2),Q(1,18)),
              ('inner',S.Matrix([0,0,0,1,0,0]),Q(1,4),Q(11,144)),
              ('center_join',S.Matrix([1,0,0,1,0,0]),Q(1,4),Q(1,18)),
              ('same_leaf_join',S.Matrix([0,0,1,1,0,0]),Q(4,25),Q(14,225)),
              ('other_leaf_join',S.Matrix([0,1,0,1,0,0]),Q(4,25),None)]
    for name,w,umax2,minimum in expected:
        B=clean((vb.T*J*w)[0]);C=clean(2-(w.T*J*w)[0]);E=clean((w.T*J*wd)[0])
        f,u,v=symbolic_det(Q(1,4),A,B,C,D,E)
        target={'isolated':1-v-v*v/16-u*u/2,
          'inner':1-u/2-v-(2*u-v)**2/16,
          'center_join':1-u/2-v-v*v/16,
          'same_leaf_join':1-u/2-v-(u+v)**2/16,
          'other_leaf_join':1-u/2-v-(u-v)**2/16}[name]
        eq(f,target)
        umax=S.sqrt(umax2);values=[clean(f.subs({u:uu,v:vv})) for uu in (0,umax) for vv in (0,Q(2,3))]
        assert all(x>0 for x in values)
        if minimum is not None:assert min(values)==minimum
        main.append(dict(profile=name,B=wire(B),C=wire(C),E=wire(E),determinant=str(f),rectangle_values=list(map(wire,values))))
    # Wleft is A6: mate L -- center0 -- main1 -- main2 -- main3 -- main4.
    H=mat(6,[(i,i+1,1) for i in range(5)]);J=H.inv();assert pd(H)
    va=S.Matrix([0,1,0,0,0,0]);vb=S.Matrix([0,0,0,0,0,1]);wd=S.Matrix([1,0,0,0,0,0])
    alpha=clean((2-(va.T*J*va)[0])/4);eq(alpha,Q(1,7))
    eq(1-(vb.T*J*vb)[0],alpha);eq((va.T*J*vb)[0]/2,alpha)
    A=clean((va.T*J*wd)[0]/2);D=clean(2-(wd.T*J*wd)[0]);eq(A,Q(5,14));eq(D,Q(8,7))
    fork=[]
    for name,w in [('isolated',S.zeros(6,1)),('center3',va)]:
        B=clean((vb.T*J*w)[0]);C=clean(2-(w.T*J*w)[0]);E=clean((w.T*J*wd)[0])
        f,u,v=symbolic_det(alpha,A,B,C,D,E)
        target=1-5*v/7-v*v/28-(2*u*u/7 if name=='isolated' else 4*u/7)
        eq(f,target)
        if name=='isolated':eq(f.subs({u:Q(4,3),v:Q(2,3)}),0)
        else:eq(f.subs({u:Q(2,3),v:Q(2,3)}),Q(8,63))
        fork.append(dict(profile=name,B=wire(B),C=wire(C),E=wire(E),determinant=str(f),
            weak_rectangle_minimum='0, strict by u<4/3' if name=='isolated' else '8/63'))
    # Opposite component identifications that corrected the original draft.
    D8=mat(8,[(i,i+1,1) for i in range(5)]+[(0,6,1),(0,7,1)])
    E8=mat(8,[(i,i+1,1) for i in range(5)]+[(1,6,1),(0,7,1)])
    assert pd(D8) and pd(E8);eq(D8.det(),4);eq(E8.det(),1)
    witness=[]
    for name,n,edges,kernel in [
       ('affine_B4',5,[(0,1,1),(0,2,1),(0,3,1),(3,4,T)],[2,1,1,2,T]),
       ('affine_B3',4,[(0,1,1),(0,2,1),(0,3,T)],[2,1,1,T]),
       ('affine_D4',5,[(0,i,1) for i in range(1,5)],[2,1,1,1,1]),
       ('affine_A2',3,[(0,1,1),(1,2,1),(0,2,1)],[1,1,1])]:
        G=mat(n,edges);v=S.Matrix(kernel);assert G*v==S.zeros(n,1)
        assert all(x>0 for x in v)
        assert all(pd(G.extract([j for j in range(n) if j!=i],[j for j in range(n) if j!=i])) for i in range(n))
        witness.append(dict(name=name,edges=[[i,j,wire(w)] for i,j,w in edges],positive_kernel=list(map(wire,kernel))))
    # New full affine replacements really are full-order affine corners.
    replacements=[]
    for kind,n,a,b,port,nxt,fork in [('C',18,0,6,1,2,False),('B',18,15,9,14,13,False),('B',18,16,4,17,0,True),('D',18,14,4,15,0,True)]:
        G,nu,E=affine(kind,n);ids=[i for i in range(n) if i!=a]
        K=G.extract(ids,ids).row_join(S.zeros(n-1,1));K=K.col_join(S.zeros(1,n));K[-1,-1]=2
        contact=port if fork else nxt;w=T if fork else 1;ix=ids.index(contact);K[ix,-1]=K[-1,ix]=-w
        null=K.nullspace();assert len(null)==1
        z=null[0];z=z if z[0]>0 else -z;assert all(x>0 for x in z)
        assert pd(K[:n-1,:n-1])
        replacements.append(dict(old_kind=kind,new_contact=contact,weight=wire(w),order=n,kernel=list(map(wire,z))))
    paths=[f'scripts/{TAG}.py',
      'scripts/fresh_simple_cusp_affine_tree_completion_20261017.py',
      'results/fresh_simple_cusp_affine_tree_completion_20261017.json',
      'results/fresh_simple_cusp_affine_tree_short_ends_20261017.json',
      'results/fresh_simple_cusp_affine_a_quadrilateral_20261017.json',
      'results/fresh_simple_cusp_affine_exclusion_20261017.json']
    return dict(status='PASS_INDEPENDENT_TERMINAL_AND_TWO_E8_LITERAL_AUDIT',
      sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
      terminal_confirmations=records,main_E8_exception=main,fork_E8_exception=fork,
      fork_exception_white_A6=[[wire(x) for x in row] for row in H.tolist()],
      opposite_components={'mate3':'finite D8','center3':'finite E8','mate4':'contains proper affine B3'},
      affine_kernel_witnesses=witness,missing_full_affine_profiles=replacements,
      scope='Independent literal normalization/determinant/correction checks only. Classification and all-rank Green identities are reviewed in prose. No all-finite-corner transfer or ambient dimension conclusion.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['terminal_confirmations']),flush=True)
