#!/usr/bin/env python3
"""Exact rational H9 polytope with an actual all-ideal seven-simplex face.

Independent literal Gram and exhaustive 220-support cone census. No project
arithmetic imports. This is a scope example, not a dimension/novelty claim.
"""
import argparse
from functools import reduce
from itertools import combinations
import hashlib
import json
from math import gcd, lcm
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_all_ideal_seven_simplex_realization_20261017'

def enc(M):
    if isinstance(M,s.MatrixBase):return [[str(x) for x in M.row(i)] for i in range(M.rows)]
    return str(M)

def inertia(A):
    """Exact rational symmetric congruence, including zero-diagonal pivots."""
    A=s.Matrix(A);p=n=z=0
    while A.rows:
        k=next((i for i in range(A.rows) if A[i,i]),None)
        if k is not None:
            ids=[k]+[i for i in range(A.rows) if i!=k];A=A.extract(ids,ids)
            a=A[0,0];p+=int(bool(a>0));n+=int(bool(a<0))
            A=A[1:,1:]-A[1:,:1]*A[:1,1:]/a
        else:
            pair=next(((i,j) for i,j in combinations(range(A.rows),2) if A[i,j]),None)
            if pair is None:z+=A.rows;break
            ids=list(pair)+[i for i in range(A.rows) if i not in pair];A=A.extract(ids,ids)
            H=A[:2,:2];p+=1;n+=1
            A=A[2:,2:]-A[2:,:2]*H.inv()*A[:2,2:]
    return [p,n,z]

def primitive(v):
    den=lcm(*(int(x.q) for x in v));v=[int(x*den) for x in v]
    g=reduce(gcd,(abs(x) for x in v));assert g
    return tuple(x//g for x in v)

def components(A,ids):
    todo=set(ids);out=[]
    while todo:
        seed=min(todo);todo.remove(seed);part={seed};stack=[seed]
        while stack:
            i=stack.pop()
            for j in list(todo):
                if A[i,j]:todo.remove(j);part.add(j);stack.append(j)
        out.append(sorted(part))
    return out

def affine_name(H):
    n=H.rows
    if n==2:assert H==s.Matrix([[2,-2],[-2,2]]);return 'affine A1'
    adj={i:[j for j in range(n) if i!=j and H[i,j]] for i in range(n)}
    assert all(H[i,j]==-1 for i in range(n) for j in adj[i])
    degrees=sorted(map(len,adj.values()))
    if degrees==[2]*n:return f'affine A{n-1}'
    forks=[i for i in adj if len(adj[i])==3]
    if len(forks)==2:
        assert degrees==[1]*4+[2]*(n-6)+[3]*2
        assert all(sum(len(adj[j])==1 for j in adj[i])==2 for i in forks)
        return f'affine D{n-1}'
    assert len(forks)==1
    arms=[];f=forks[0]
    for j in adj[f]:
        prev=f;cur=j;length=1
        while len(adj[cur])==2:
            nxt=next(k for k in adj[cur] if k!=prev);prev,cur=cur,nxt;length+=1
        assert len(adj[cur])==1;arms.append(length)
    arms.sort();assert arms in ([1,2,5],[1,3,3])
    return 'affine E8' if arms==[1,2,5] else 'affine E7'

def sources():
    tag='fresh_finite_edge_cusp_tradeoff10_20261017'
    todo=[f'scripts/{TAG}.py',f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md']
    seen={}
    while todo:
        p=todo.pop()
        if p in seen:continue
        f=ROOT/p;seen[p]=hashlib.sha256(f.read_bytes()).hexdigest()
        if p.endswith('.json'):
            for child,h in json.loads(f.read_text()).get('sources',{}).items():
                assert hashlib.sha256((ROOT/child).read_bytes()).hexdigest()==h,child
                todo.append(child)
    return seen

def run():
    G=2*s.eye(10);edges=[(i,(i+1)%8) for i in range(8)]+[(8,7),(9,3)]
    for i,j in edges:G[i,j]=G[j,i]=-1
    J=G.inv();assert G.det()==-16 and inertia(G)==[9,1,0]
    assert [J[i,i] for i in range(10)]==[0]*8+[s.Rational(1,2)]*2
    assert J[8:,8:]==s.Matrix([[s.Rational(1,2),-s.Rational(1,2)],[-s.Rational(1,2),s.Rational(1,2)]])
    # Normal rows are coordinates in the ten original-root basis.
    N=s.eye(10).col_join((-2*J[:,8]).T).col_join((-2*J[:,9]).T)
    GG=N*G*N.T;assert all(GG[i,i]==2 for i in range(12))
    assert all(GG[i,j] in (0,-1,-2) for i,j in combinations(range(12),2))
    assert [(i,j) for i,j in combinations(range(12),2) if GG[i,j]==-2]==[(8,10),(9,11),(10,11)]
    H=-J*s.ones(10,1);assert (H.T*G*H)[0]==-70
    assert list(N*G*H)==[-1]*10+[-8,-8]
    # y=-Gx; A y >=0 is equivalent to all twelve root halfspaces.
    A=s.eye(10).col_join(-J[8,:]).col_join(-J[9,:])
    assert N*G*(-J)==-s.diag(*([1]*10+[2]*2))*A
    assert all(v>0 for v in A*s.ones(10,1))
    rays={};support_trials=[]
    for ids in combinations(range(12),9):
        B=A[list(ids),:];ker=B.nullspace()
        rec={'support':list(ids),'rank':B.rank()}
        if len(ker)!=1:rec['status']='rank_below9';support_trials.append(rec);continue
        v=ker[0];pos=all(a>=0 for a in A*v);neg=all(a<=0 for a in A*v)
        assert not (pos and neg)
        if not pos and not neg:rec['status']='infeasible';support_trials.append(rec);continue
        if neg:v=-v
        key=primitive(v);rec.update(status='ray',ray=list(key));support_trials.append(rec)
        rays.setdefault(key,[]).append(list(ids))
    assert len(support_trials)==220 and len(rays)==23
    expected={tuple(s.eye(10)[:,j]) for j in range(8)}
    for i in (8,9):
        for j in range(8):
            if J[i,j]<0:expected.add(primitive(s.eye(10)[:,j]-2*J[i,j]*s.eye(10)[:,i]))
    expected.add(tuple([0]*8+[1,1]));assert set(rays)==expected
    records=[]
    for key in sorted(rays):
        y=s.Matrix(key);x=-J*y;q=(x.T*G*x)[0];assert q<=0
        future=(x.T*G*H)[0];assert future<0
        active=[i for i,a in enumerate(A*y) if not a]
        assert N[active,:].rank()==9
        rec=dict(y=list(key),x=enc(x),norm=str(q),center_product=str(future),active=active,
                 type='ideal' if q==0 else 'finite',support_bases=rays[key])
        if q<0:
            assert len(active)==9 and inertia(GG.extract(active,active))==[9,0,0]
        else:
            parts=components(GG,active);cusps=[]
            for part in parts:
                B=GG.extract(part,part);assert inertia(B)==[len(part)-1,0,1]
                ker=B.nullspace();assert len(ker)==1
                h=ker[0]
                if all(t<0 for t in h):h=-h
                assert all(t>0 for t in h)
                null=N[part,:].T*h
                c=next(null[i]/x[i] for i in range(10) if x[i])
                assert c>0 and null==c*x
                cusps.append(dict(nodes=part,type=affine_name(B),positive_kernel=enc(h),
                                  rank=len(part)-1,null_ray_scale=str(c)))
            assert sum(c['rank'] for c in cusps)==8
            rec['cusp_components']=cusps
        records.append(rec)
    assert sum(r['type']=='finite' for r in records)==14
    assert sum(r['type']=='ideal' for r in records)==9
    # Every face is an intersection of facets; enumerate all 2^12 subsets,
    # closing them by their complete vertex incidence. The y>=0 section is
    # bounded, so this enumerates the entire projective face lattice.
    allmask=(1<<23)-1
    masks=[sum(1<<j for j,r in enumerate(records) if i in r['active']) for i in range(12)]
    closed={}
    for subset in range(1<<12):
        mask=allmask
        for i in range(12):
            if subset>>i&1:mask &= masks[i]
        if not mask or mask in closed:continue
        vids=[j for j in range(23) if mask>>j&1]
        active=sorted(set.intersection(*(set(records[j]['active']) for j in vids)))
        rank=N[active,:].rank() if active else 0;dimension=9-rank
        assert s.Matrix.hstack(*(s.Matrix(records[j]['x']) for j in vids)).rank()==dimension+1
        closed[mask]=dict(vertices=vids,active=active,dimension=dimension)
    fv=[sum(f['dimension']==k for f in closed.values()) for k in range(10)]
    assert fv[0]==23 and fv[8]==12 and fv[9]==1
    facets=[]
    for i in range(12):
        f=closed[masks[i]];assert f['dimension']==8
        hi=H-(N[i,:]*G*H)[0]*N[i,:].T/2
        dots=N*G*hi;assert dots[i]==0 and all(dots[j]<0 for j in range(12) if j!=i)
        assert (hi.T*G*hi)[0]<0
        vids=f['vertices'];subfaces=[f for mask,f in closed.items() if mask&masks[i]==mask]
        localfv=[sum(f['dimension']==k for f in subfaces) for k in range(9)]
        facet=dict(index=i,vertices=vids,finite_vertices=sum(records[j]['type']=='finite' for j in vids),
                   ideal_vertices=sum(records[j]['type']=='ideal' for j in vids),f_vector=localfv,
                   strict_relative_interior=enc(hi),interior_norm=str((hi.T*G*hi)[0]))
        if i in (10,11):
            assert len(vids)==9 and facet['finite_vertices']==7 and facet['ideal_vertices']==2
            assert localfv==[int(s.binomial(9,k+1)) for k in range(9)]
        facets.append(facet)
    targetmask=masks[8]&masks[9];target=closed[targetmask]
    assert target['dimension']==7 and len(target['vertices'])==8 and target['active']==[8,9]
    assert all(records[j]['type']=='ideal' for j in target['vertices'])
    assert {tuple(records[j]['y']) for j in target['vertices']}=={tuple(s.eye(10)[:,j]) for j in range(8)}
    targetfv=[sum(f['dimension']==k for mask,f in closed.items() if mask&targetmask==mask) for k in range(8)]
    assert targetfv==[int(s.binomial(8,k+1)) for k in range(8)]
    projected=G[:8,:8]-G[:8,8:]*G[8:,8:].inv()*G[8:,:8]
    assert inertia(projected)==[7,1,0]
    assert projected[2,3]==-1 and projected[2,2]==2 and projected[3,3]==s.Rational(3,2)
    # Its angle has cos^2=1/3, strictly between the consecutive m3/m4
    # values 1/4 and1/2; thus the induced simplex is not Coxeter.
    assert s.Rational(1,4)<s.Rational(1,3)<s.Rational(1,2)
    edge_types={}
    for f in closed.values():
        if f['dimension']!=1:continue
        assert len(f['vertices'])==2
        key='-'.join(sorted(records[j]['type'] for j in f['vertices']))
        edge_types[key]=edge_types.get(key,0)+1
    return dict(status='PASS_EXACT_FINITE_VOLUME_H9_WITH_ALL_IDEAL_SEVEN_SIMPLEX',
        scope='Explicit rational scope example; no new dimension, novelty, or classification identification claim.',
        sources=sources(),original_gram=enc(G),original_inverse=enc(J),original_inertia=inertia(G),
        root_coordinates=enc(N),full_gram=enc(GG),common_center=enc(H),common_center_norm='-70',
        y_inequality_matrix=enc(A),support_trials=support_trials,vertices=records,
        counts=dict(walls=12,finite_vertices=14,ideal_vertices=9,exterior_vertices=0,
                    support_subsets=220,face_subsets=4096),f_vector=fv,edge_types=edge_types,
        facets=facets,faces=sorted(closed.values(),key=lambda f:(f['dimension'],f['vertices'])),
        all_ideal_seven_simplex=dict(support=[8,9],vertices=target['vertices'],f_vector=targetfv,
            projected_black_gram=enc(projected),noncoxeter_angle_cosine_squared='1/3'),
        finite_volume_reason='Complete pointed cone has only future finite/ideal rays; every ideal link is a full-rank product of bounded affine simplex alcoves.',
        relation_to_tradeoff10='The all-eight-ideal alternative is realized. It has no finite edge and does not challenge coefficient10.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(out.read_text())==d
    else:out.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],'saved_replay=',a.verify,'f_vector=',d['f_vector'])
