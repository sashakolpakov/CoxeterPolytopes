#!/usr/bin/env python3
"""Two exact all-rank families refuting a proposed affine-extension reduction.

The all-rank arguments are in the accompanying report; exact finite instances
check the displayed graph formulas and continuation obstructions.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_large_affine_white_contact_counterexamples_20261017'


def field(d):
    class Q:
        def __init__(self,a=0,b=0):
            if isinstance(a,Q):self.a,self.b=a.a,a.b
            else:self.a,self.b=F(a),F(b)
        def __add__(self,x):x=Q(x);return Q(self.a+x.a,self.b+x.b)
        __radd__=__add__
        def __neg__(self):return Q(-self.a,-self.b)
        def __sub__(self,x):return self+-Q(x)
        def __rsub__(self,x):return Q(x)+-self
        def __mul__(self,x):x=Q(x);return Q(self.a*x.a+d*self.b*x.b,self.a*x.b+self.b*x.a)
        __rmul__=__mul__
        def __truediv__(self,x):
            x=Q(x);norm=x.a*x.a-d*x.b*x.b;assert norm
            return self*Q(x.a/norm,-x.b/norm)
        def __bool__(self):return bool(self.a or self.b)
        def __eq__(self,x):x=Q(x);return self.a==x.a and self.b==x.b
        def sign(self):
            if not self.b:return (self.a>0)-(self.a<0)
            if not self.a:return (self.b>0)-(self.b<0)
            if self.a*self.b>0:return (self.a>0)-(self.a<0)
            delta=self.a*self.a-d*self.b*self.b
            return ((self.a>0)-(self.a<0))*((delta>0)-(delta<0))
        def out(self):return [str(self.a),str(self.b)]
    return Q


def sub(G,ids):return [[G[i][j]for j in ids]for i in ids]


def pd_pivots(G,Q):
    H=[[Q(x)for x in row]for row in G];piv=[]
    for k in range(len(H)):
        z=H[k][k];assert z.sign()>0;piv.append(z)
        for i in range(k+1,len(H)):
            for j in range(k+1,len(H)):H[i][j]-=H[i][k]*H[k][j]/z
    return piv


def determinant(G,Q):
    H=[[Q(x)for x in row]for row in G];out=Q(1)
    for j in range(len(H)):
        p=next((i for i in range(j,len(H))if H[i][j]),None)
        if p is None:return Q()
        if p!=j:H[p],H[j]=H[j],H[p];out=-out
        z=H[j][j];out*=z
        for i in range(j+1,len(H)):
            f=H[i][j]/z
            for k in range(j+1,len(H)):H[i][k]-=f*H[j][k]
    return out


def graph(n,edges):
    g={i:{}for i in range(n)}
    for i,j,m in edges:g[i][j]=m;g[j][i]=m
    return g


def induced(g,ids):
    ids=set(ids);return {i:{j:m for j,m in g[i].items()if j in ids}for i in ids}


def pieces(g):
    unseen=set(g);out=[]
    while unseen:
        todo=[min(unseen)];part=set()
        while todo:
            i=todo.pop()
            if i in part:continue
            part.add(i);todo.extend(g[i])
        unseen-=part;out.append(induced(g,part))
    return out


def tree_code(g):
    assert sum(map(len,g.values()))==2*(len(g)-1)
    def rec(i,p):return tuple(sorted((m,rec(j,i))for j,m in g[i].items()if j!=p))
    return min(rec(i,None)for i in g)


def finite_models(n):
    out={}
    def add(name,m,es):out[tree_code(graph(m,es))]=name
    add('A1',1,[]);add('I2(5)',2,[(0,1,5)])
    add('B3',3,[(0,1,4),(1,2,3)]);add('H4',4,[(0,1,5),(1,2,3),(2,3,3)])
    for m in range(4,n+1):add('D'+str(m),m,[(i,i+1,3)for i in range(m-2)]+[(1,m-1,3)])
    return out


def minimal_affine_certificate(G,nu,Q):
    n=len(G);assert len(nu)==n and all(x.sign()>0 for x in nu)
    assert all(sum((G[i][j]*nu[j]for j in range(n)),Q())==0 for i in range(n))
    edges=[(i,j)for i in range(n)for j in range(i)if G[i][j]]
    assert all(G[i][j].sign()<0 for i,j in edges)
    # D_nu G D_nu is a connected weighted Laplacian.  This proves PSD,
    # corank1, and positive definiteness of every proper principal block.
    seen={0};todo=[0]
    for i in todo:
        for j in range(n):
            if j!=i and G[i][j]and j not in seen:seen.add(j);todo.append(j)
    assert len(seen)==n
    assert all(G[i][i]*nu[i]*nu[i]==sum((-G[i][j]*nu[i]*nu[j]for j in range(n)if j!=i),Q())for i in range(n))
    return dict(positive_kernel=[x.out()for x in nu],connected=True,laplacian_identity=True,
        inertia=[n-1,0,1],every_proper_principal_block_PD=True)


def example(kind,n):
    assert n>=18
    Q=field(5 if kind=='white5_leaf'else 2);sqrt=Q(0,1)
    p,t=0,1;a,c=n-4,n-3;z=n;y=n+1
    oldedges=[(i,i+1,3)for i in range(n-5)]+[(0,a,3),(0,c,3),(n-5,n-2,3),(n-5,n-1,3)]
    extras=[(z,a,5)]if kind=='white5_leaf'else[(z,a,4),(z,c,3)]
    g=graph(n+1,oldedges+extras);G=[[Q(2)if i==j else Q()for j in range(n+1)]for i in range(n+1)]
    values={3:Q(1),5:Q(F(1,2),F(1,2))}if kind=='white5_leaf'else{3:Q(1),4:sqrt}
    for i,j,m in oldedges+extras:G[i][j]=G[j][i]=-values[m]
    oldnu=[Q(1)]*(n-4)+[Q(F(1,2))]*4
    original=minimal_affine_certificate(sub(G,range(n)),oldnu,Q)
    black=[p,t]if kind=='white5_leaf'else[a,p]
    assert any(i not in black for i in g[z])
    models=finite_models(n);deletions=[]
    for b in black:
        ids=[i for i in range(n+1)if i!=b]
        forest=sorted(models[tree_code(h)]for h in pieces(induced(g,ids)))
        expected=sorted(['I2(5)','A1','D'+str(n-3)]if b==p and kind=='white5_leaf'else
                        ['H4','D'+str(n-4)]if kind=='white5_leaf'else
                        ['D'+str(n)]if b==a else['B3','D'+str(n-3)])
        assert forest==expected
        piv=pd_pivots(sub(G,ids),Q)
        deletions.append(dict(deleted=b,forest=forest,pivots=[x.out()for x in piv]))
    assert determinant(G,Q).sign()<0
    # Interlacing with a checked PD n-root deletion gives exactly one negative
    # direction; determinant<0 excludes a radical.
    rec=dict(kind=kind,order=n,old_affine_type='D_tilde_order_'+str(n),
        old_edges=oldedges,added_edges=extras,black=black,white_contacts=sorted(set(g[z])-set(black)),
        roles=dict(p=p,t=t,a=a,c=c,z=z),original_affine=original,spherical_deletions=deletions,
        whole_determinant=determinant(G,Q).out(),whole_inertia=[n,1,0])
    H=[[Q(2)if i==j else Q()for j in range(n+2)]for i in range(n+2)]
    for i in range(n+1):
        for j in range(n+1):H[i][j]=G[i][j]
    H[y][z]=H[z][y]=Q(-1)
    if kind=='white5_leaf':
        tau=values[5]*values[5];assert tau==Q(F(3,2),F(1,2))
        ids=[y,z,a,p];det=determinant(sub(H,ids),Q)
        assert det==Q(3,-2)and det.sign()<0
        rec['no_small_affine_proof']='Any connected z-containing proper old block has Schur2-tau*r with r positive rational, hence nonzero; the full old affine block plusz is indefinite. Only the original ordern affine occurs.'
        rec['continuation_failure']=dict(deletion=t,principal_roots=ids,labels=[3,5,3],determinant=det.out(),reason='negative determinant')
    else:
        ids=[i for i in range(n+1)if i!=c]
        nu={i:Q(1)for i in ids};nu[z]=sqrt/2;nu[n-2]=nu[n-1]=Q(F(1,2))
        rec['second_large_affine']=dict(omitted=c,type='B_tilde_order_'+str(n),
            **minimal_affine_certificate(sub(G,ids),[nu[i]for i in ids],Q))
        cycle=[z,a,p,c];value=sum((G[i][j]for i in cycle for j in cycle),Q())
        assert value==Q(2,-2)and value.sign()<0
        rec['cycle_negative_vector']=dict(roots=cycle,coefficients=[1,1,1,1],value=value.out())
        rec['no_small_affine_proof']='A PSD subdiagram containingz must omit a,p,or c, since the four-cycle is indefinite. The first two full deletions are spherical; the c deletion is minimal affine of ordern. Only the two displayed ordern affines occur.'
        ids=[z,a,c,y];nu=[Q(1),sqrt/2,Q(F(1,2)),Q(F(1,2))]
        obstruction=minimal_affine_certificate(sub(H,ids),nu,Q)
        rec['continuation_failure']=dict(deletion=p,principal_roots=ids,labels=[4,3,3],reason='affine doubled star',**obstruction)
    return rec


def run():
    checks=[example(kind,n)for kind in ('white5_leaf','white3_black4_switch')for n in (18,19,23,31)]
    p=Path(__file__)
    return dict(status='PASS_EXACT_LARGE_AFFINE_WHITE_CONTACT_COUNTEREXAMPLES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()},checks=checks,
        all_rank_formulas=True,minimum_order=18,number_of_all_rank_families=2,
        continuation_tests='Both displayed families fail even the weakest nonright extra leaf y--3--z, with y right to every old affine root.',
        scope='Counterexamples to the proposed graph-theoretic contact-only-black implication from two spherical deletions and no affine subdiagram of order<=17. The formulas and written proof apply to every order>=18. No finite-volume ambient realization or counterexample to a dimension theorem is asserted. No complete arbitrary-contact or continuation classification.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=json.loads(json.dumps(run()));f=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(f.read_text())
    else:f.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['checks']),flush=True)
