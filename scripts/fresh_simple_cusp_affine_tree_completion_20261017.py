#!/usr/bin/env python3
"""Exact checks for the all-rank affine-tree quadrilateral argument, n>=18."""
import argparse
from collections import deque
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

import fresh_simple_cusp_affine_a_quadrilateral_20261017 as A
import fresh_simple_cusp_affine_tree_short_ends_20261017 as E

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_simple_cusp_affine_tree_completion_20261017'
Q,Z,SQ=A.Q,A.Z,A.SQ

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def inv(M):
    n=len(M);T=[r[:]+[Q(i==j) for j in range(n)] for i,r in enumerate(M)]
    for k in range(n):
        assert T[k][k].sign()>0
        d=T[k][k];T[k]=[x/d for x in T[k]]
        for i in range(n):
            if i!=k:
                d=T[i][k]
                if d!=0:T[i]=[x-d*y for x,y in zip(T[i],T[k])]
    return [r[n:] for r in T]

def affine(kind,n):
    nu=[Q(1)]*n;G=[[Q(2) if i==j else Z for j in range(n)] for i in range(n)]
    if kind=='C':
        edges=[(i,i+1,SQ if i in(0,n-2) else Q(1)) for i in range(n-1)]
        nu[0]=nu[-1]=SQ/2;ends=[set(range(6)),set(range(n-6,n))]
    elif kind=='B':
        edges=[(i,i+1,SQ if i==n-4 else Q(1)) for i in range(n-3)]+[(0,n-2,Q(1)),(0,n-1,Q(1))]
        nu[n-3]=SQ/2;nu[-2]=nu[-1]=Q(F(1,2))
        ends=[set(range(4))|{n-2,n-1},set(range(n-8,n-2))]
    else:
        assert kind=='D'
        edges=[(i,i+1,Q(1)) for i in range(n-5)]+[(0,n-4,Q(1)),(0,n-3,Q(1)),(n-5,n-2,Q(1)),(n-5,n-1,Q(1))]
        nu[-4:]=[Q(F(1,2))]*4
        ends=[set(range(4))|{n-4,n-3},set(range(n-8,n-4))|{n-2,n-1}]
    for i,j,w in edges:G[i][j]=G[j][i]=-w
    assert all(sum((G[i][j]*nu[j] for j in range(n)),Z)==0 for i in range(n))
    return G,nu,ends

def path(G,a,b):
    todo=deque([[a]]);seen={a}
    while todo:
        p=todo.popleft()
        if p[-1]==b:return p
        for j,w in enumerate(G[p[-1]]):
            if j!=p[-1] and w!=0 and j not in seen:seen.add(j);todo.append(p+[j])
    assert False

def context(kind,n,a,b):
    G,nu,ends=affine(kind,n);W=[i for i in range(n) if i not in(a,b)]
    J=inv([[G[i][j] for j in W]for i in W]);pos={i:j for j,i in enumerate(W)}
    def dot(w,z):return sum((x*J[pos[i]][pos[j]]*y for i,x in w.items() for j,y in z.items()),Z)
    va={i:-G[a][i] for i in W if G[a][i]!=0};vb={i:-G[b][i]for i in W if G[b][i]!=0}
    alpha=nu[a]**2*(2-dot(va,va))
    assert alpha==nu[b]**2*(2-dot(vb,vb))==nu[a]*nu[b]*dot(va,vb)
    p=path(G,a,b);distance=len(p)-1
    R=Q(distance+int(nu[a]==Q(F(1,2)))+int(nu[b]==Q(F(1,2))))
    assert alpha==Q(1)/R
    def profile(w,l,c0):
        B=nu[b]*dot(vb,w);C=2-dot(w,w)
        assert B==l/R and C==c0+l*l/R
        return dict(w=w,l=l,c0=c0,height=sum((nu[i]*v for i,v in w.items()),Z))
    return G,nu,W,dot,alpha,R,profile

def profiles(G,nu,a,b,R,profile):
    p=path(G,a,b);near=p[1];out=[]
    def add(name,w,l,c0,hoff=None):
        t=profile(w,l,c0);t['name']=name
        if hoff is not None:assert t['height']==hoff
        out.append(t)
    add('isolated',{},Z,Q(2),Z)
    if nu[a]==SQ/2:
        add('terminal_long3',{near:Q(1)},Q(1),Q(1),Q(1))
        add('terminal_long4',{near:SQ},SQ,Z,SQ)
        add('terminal_next3',{p[2]:Q(1)},Q(2),Z,Q(1))
    elif nu[a]==Q(F(1,2)):
        mates=[i for i in range(len(G)) if i!=a and G[near][i]!=0 and nu[i]==Q(F(1,2))]
        assert len(mates)==1;mate=mates[0]
        add('fork_mate3',{mate:Q(1)},Q(1),Q(1),Q(F(1,2)))
        add('fork_mate4',{mate:SQ},SQ,Z,SQ/2)
        add('fork_center3',{near:Q(1)},Q(2),Z,Q(1))
    else:
        assert nu[a]==1
        # Outer components: the part beyond a, away from b.
        outer=set();todo=[j for j in range(len(G)) if j not in(a,near) and G[a][j]!=0]
        while todo:
            i=todo.pop()
            if i in outer:continue
            outer.add(i)
            todo.extend(j for j in range(len(G)) if j!=i and j!=a and G[i][j]!=0 and j not in outer)
        outer=sorted(outer);r=len(outer);assert r>=1
        Jout=inv([[G[i][j] for j in outer]for i in outer]);idx={i:k for k,i in enumerate(outer)}
        v={i:-G[a][i] for i in outer if G[a][i]!=0}
        def odot(w,z):return sum((x*Jout[idx[i]][idx[j]]*y for i,x in w.items()for j,y in z.items()),Z)
        assert odot(v,v)==1
        for i in outer:
            q=odot({i:Q(1)},{i:Q(1)});t=odot(v,{i:Q(1)})
            assert t.sign()>0
            add('outer_basis_'+str(i),{i:Q(1)},Z,2-q,t)
        add('inner3',{near:Q(1)},Q(1),Q(1),Q(1))
        kind='B' if any(nu[i]==SQ/2 for i in outer) else 'D'
        ports=[i for i in outer if G[a][i]!=0]
        if kind=='B':
            if r==1:joins=[{outer[0]:Q(1)},{outer[0]:SQ}]
            elif r==2:joins=[{i:Q(1)}for i in outer]
            else:assert len(ports)==1;joins=[{ports[0]:Q(1)}]
        else:
            if r==2:joins=[{i:Q(1)}for i in outer]+[{i:Q(1)for i in outer}]
            elif r in(3,4):
                leaves=[i for i in outer if sum(G[i][j]!=0 for j in outer if j!=i)==1]
                joins=[{i:Q(1)}for i in sorted(set(leaves+ports))]
            else:assert len(ports)==1;joins=[{ports[0]:Q(1)}]
        for j,w in enumerate(joins):
            q,t=odot(w,w),odot(v,w)
            add('join_'+kind+str(r)+'_'+str(j),w|{near:Q(1)},Q(1),1-q,1+t)
    return out

def schur_checks():
    records=[];count=0;digest=hashlib.sha256()
    # Literal finite checks of formulas whose all-rank proof is the Laplacian identity.
    for n in (18,19,23):
        for kind in ('B','C','D'):
            G,nu,ends=affine(kind,n)
            for a,b in combinations(range(n),2):
                P=len(path(G,a,b))-1
                if P<=3 or any({a,b}<=e for e in ends):continue
                G,nu,W,dot,alpha,R,pc=context(kind,n,a,b)
                _,_,_,_,_,_,pd=context(kind,n,b,a)
                cs=profiles(G,nu,a,b,R,pc);ds=profiles(G,nu,b,a,R,pd)
                for c,d in product(cs,ds):
                    assert dot(c['w'],d['w'])==c['l']*d['l']/R
                    count+=1
                row=dict(type=kind,n=n,a=a,b=b,distance=P,R=R.wire(),
                    cusp_coefficients=[nu[a].wire(),nu[b].wire()],c_profiles=len(cs),d_profiles=len(ds))
                records.append(row);digest.update((json.dumps(row,sort_keys=True)+'\n').encode())
            print('TREE_SCHUR_PASS',kind,n,len(records),flush=True)
    return dict(contexts=records,profile_pairs=count,digest=digest.hexdigest())

def polynomial(alpha,A0,B,C,D,E0):
    # Coefficients of detS/(h²k²) in 1,u,v,u²,v²,uv, u=1/h,v=1/k.
    return [Q(1),-2*B,-2*A0,B*B-alpha*C,A0*A0-alpha*D,2*A0*B-2*alpha*E0]
def evaluate(co,u,v):return co[0]+co[1]*u+co[2]*v+co[3]*u*u+co[4]*v*v+co[5]*u*v

def exception_checks():
    out=[]
    for kind in ('B','D'):
        n=20
        leaf=n-1 if kind=='B' else n-3
        for name,a,b,cases in [
          ('central_D3_P4',1,5,[('isolated',{},[1,0,-1,F(-1,2),F(-1,16),0],SQ/2),
            ('inner',{2:Q(1)},[1,F(-1,2),-1,F(-1,4),F(-1,16),F(1,4)],Q(F(1,2))),
            ('join_center',{0:Q(1),2:Q(1)},[1,F(-1,2),-1,0,F(-1,16),0],Q(F(1,2))),
            ('join_same_leaf',{leaf:Q(1),2:Q(1)},[1,F(-1,2),-1,F(-1,16),F(-1,16),F(-1,8)],Q(F(2,5))),
            ('join_other_leaf',{leaf-1:Q(1),2:Q(1)},[1,F(-1,2),-1,F(-1,16),F(-1,16),F(1,8)],Q(F(2,5)))]),
          ('fork_P6',leaf-1,5,[('isolated',{},[1,0,F(-5,7),F(-2,7),F(-1,28),0],Q(F(4,3))),
            ('center',{0:Q(1)},[1,F(-4,7),F(-5,7),0,F(-1,28),0],2-SQ)])]:
            G,nu,W,dot,alpha,R,pc=context(kind,n,a,b)
            va={i:-G[a][i]for i in W if G[a][i]!=0};vb={i:-G[b][i]for i in W if G[b][i]!=0}
            w_d={leaf:Q(1)};aa=nu[a]*dot(va,w_d);dd=2-dot(w_d,w_d)
            for tag,w,expected,umax in cases:
                bb=nu[b]*dot(vb,w);cc=2-dot(w,w);ee=dot(w,w_d)
                coeff=polynomial(alpha,aa,bb,cc,dd,ee)
                assert coeff==list(map(Q,expected)),(kind,name,tag,[x.wire()for x in coeff])
                values=[evaluate(coeff,u,v)for u,v in product((Z,umax),(Z,Q(F(2,3))))]
                assert all(x.sign()>=0 for x in values)
                strict=(name,tag)!=('fork_P6','isolated')
                if strict:assert all(x.sign()>0 for x in values)
                else:assert values[-1]==0 # actual u<4/3, so the physical value is strictly positive.
                out.append(dict(type=kind,case=name,profile=tag,coefficients=[x.wire()for x in coeff],
                    umax=umax.wire(),vmax='2/3',rectangle_corner_values=[x.wire()for x in values],
                    strict_boundary=strict,strict_physical=True))
    return out

def scalar_checks():
    rows=[]
    def bound(name,c0,l,u,limit):
        value=2*l*u+c0*u*u;assert (Q(limit)-value).sign()>=0
        rows.append(dict(name=name,c0=c0.wire(),l=l.wire(),umax=u.wire(),value=value.wire(),bound=str(limit)))
    bound('D2_single_outer_conservative',Q(2),Z,Q(1),2) # separate 8/d²−2<14/9 at d>3/2.
    assert 8/F(3,2)**2-2==F(14,9)
    bound('D2_single_leaf_join',Q(F(1,2)),Q(1),Q(F(2,3)),F(14,9))
    bound('D2_double_leaf_join',Z,Q(1),Q(F(1,2)),1)
    bound('double_isolated',Q(2),Z,Q(1),2)
    bound('double_long3',Q(1),Q(1),2-SQ,2)
    bound('double_long4',Z,SQ,SQ/2,2)
    bound('double_next3',Z,Q(2),2-SQ,F(5,2))
    bound('fork_isolated_conservative',Q(2),Z,SQ,4)
    bound('fork_mate3',Q(1),Q(1),2*(SQ-1),4)
    bound('fork_mate4',Z,SQ,SQ,4)
    bound('fork_center3',Z,Q(2),Q(1),4)
    margins=dict(central=1-F(28,9)/4,double_central=1-(F(5,2)+F(14,9))/6,
                 fork_central=1-(4+F(14,9))/6,both_double=1-F(5,17),
                 double_fork=1-F(13,34),both_fork=1-F(8,17))
    assert all(x>0 for x in margins.values())
    assert ((36*SQ-40)/63).sign()>0
    return dict(bounds=rows,margins={k:str(v)for k,v in margins.items()})

def run():
    prev=json.loads((ROOT/f'results/{E.TAG}.json').read_text())
    assert prev['status']=='PASS_ALL_RANK_AFFINE_TREE_SHORT_OUTER_END_PAYMENT'
    for p,h in prev['sources'].items():assert sha(p)==h,p
    schur=schur_checks();exceptions=exception_checks();scalars=scalar_checks()
    b4=A.affine_witness('B4_missing_terminal_profiles',5,
        [(0,1,Q(1)),(0,2,Q(1)),(0,3,Q(1)),(3,4,SQ)],
        [Q(2),Q(1),Q(1),Q(2),SQ])
    d4=A.affine_witness('D4_center_contact_zero',5,
        [(0,i,Q(1))for i in range(1,5)],[Q(2)]+[Q(1)]*4)
    paths=[f'scripts/{TAG}.py',f'scripts/{E.TAG}.py',f'results/{E.TAG}.json',f'results/{E.TAG}.md',
        f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',f'results/{A.TAG}.md',
        'results/fresh_simple_cusp_affine_exclusion_20261017.json',
        'results/fresh_simple_cusp_affine_exclusion_20261017.md',
        'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.json',
        'results/fresh_simple_cusp_vinberg_weight_frontier_20261017.md']
    return dict(status='PASS_ALL_RANK_CLASSICAL_AFFINE_TREE_QUADRILATERAL_PAYMENT_N18',
        sources={p:sha(p)for p in paths},minimum_dimension=18,
        exact_schur_formula_checks=schur,exceptional_E8_checks=exceptions,scalar_bounds=scalars,
        extra_affine_witnesses=[b4,d4],formal_determinant_identity=A.determinant_identity(),
        all_rank_proof='Cusp-weighted Gram is a tree Laplacian; series resistance and branch Green formulas give the displayed Schur identities. Finite/affine one-node classification, not a rank cutoff, gives the complete profile lists. The report handles D2 centers, both terminal affine replacements and both E8 overlaps.',
        scope='All quadrilaterals having a classical affine-tree corner in an all-simple finite-volume Coxeter polytope n>=18 are paid. All-finite-corner transfer is separate. No ambient dimension bound asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],len(result['exact_schur_formula_checks']['contexts']),
        result['exact_schur_formula_checks']['profile_pairs'],len(result['exceptional_E8_checks']),flush=True)
