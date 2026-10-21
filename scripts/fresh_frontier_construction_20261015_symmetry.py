#!/usr/bin/env python3
"""One exact finite chamber symmetry and its transported fourteenth wall."""
import argparse,json
from pathlib import Path
from hashlib import sha256
import fresh_frontier_construction_20261015 as A
K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_frontier_construction_20261015_symmetry'
SOURCE='results/fresh_frontier_construction_20261015.json'
GEO='results/fresh_frontier_construction_20261015_geometry.json'
HELPER='scripts/fresh_frontier_construction_20261015.py'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def transpose(M):return list(map(list,zip(*M)))
def mm(a,b):return [[A.dot(r,c)for c in zip(*b)]for r in a]
def inverse(M):
    n=len(M);a=[r[:]+[K(i==j)for j in range(n)]for i,r in enumerate(M)]
    for i in range(n):
        j=next(j for j in range(i,n)if a[j][i]);a[i],a[j]=a[j],a[i]
        q=a[i][i];a[i]=[x/q for x in a[i]]
        for j in range(n):
            if j!=i:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return [r[n:]for r in a]
def determinant(M):
    a=[r[:]for r in M];d=K(1)
    for i in range(len(a)):
        j=next(j for j in range(i,len(a))if a[j][i])
        if i!=j:a[i],a[j]=a[j],a[i];d=-d
        q=a[i][i];d*=q
        for j in range(i+1,len(a)):
            z=a[j][i]/q
            for k in range(i+1,len(a)):a[j][k]-=z*a[i][k]
    return d
def run():
    p=json.loads((ROOT/SOURCE).read_text());g=json.loads((ROOT/GEO).read_text())
    for d in (p,g):
        for f,h in d['sources'].items():assert sha(f)==h,f
    G=[[A.dec(x)for x in r]for r in p['form']]
    N=[list(map(A.dec,r))for r in p['roots']];H=list(map(A.dec,p['basepoint']))
    V=list(map(A.dec,g['finite_vertices'][2]['point']))
    src=list(range(9));dst=list(range(7))+[10,11]
    assert sorted(dst)==g['finite_vertices'][2]['incident']
    gs=[[A.pair(N[i],N[j],G)for j in src]for i in src]
    gt=[[A.pair(N[i],N[j],G)for j in dst]for i in dst]
    assert gs==gt
    J=inverse(gt);assert all(x.sign()>=0 for row in J for x in row)
    B=transpose([N[i]for i in src]+[H]);D=transpose([N[i]for i in dst]+[V])
    scale=determinant(B)/determinant(D)
    if scale.sign()<0:scale=-scale
    assert scale==7+5*A.S and A.pair(H,H,G)==scale*scale*A.pair(V,V,G)
    D=transpose([N[i]for i in dst]+[[scale*x for x in V]])
    T=mm(D,inverse(B));I=[[K(i==j)for j in range(10)]for i in range(10)]
    assert all(x.integral()for row in T for x in row)
    assert mm(T,T)==I and mm(mm(transpose(T),G),T)==G
    assert determinant(T)==1 and sum(T[i][i]for i in range(10))==6
    assert all(A.mv(T,N[i])==N[j]for i,j in zip(src,dst))
    assert A.mv(T,H)==[scale*x for x in V] and A.pair(V,H,G).sign()<0
    # Exact cone-persistence decomposition for the target corner.
    coefficient=A.pair(H,V,G)/A.pair(V,V,G)
    b=A.mv(J,[A.pair(H,N[i],G)for i in dst])
    assert coefficient.sign()>0 and all(x.sign()<=0 for x in b)
    assert H==[coefficient*V[k]+sum(b[i]*N[dst[i]][k]for i in range(9))for k in range(10)]
    R=A.mv(T,N[9]);assert R not in N
    N.append(R);permutation=[]
    for n in N:
        image=A.mv(T,n);assert image in N;permutation.append(N.index(image))
    gram=[[A.pair(a,b,G)for b in N]for a in N]
    assert all(gram[i][j].sign()<=0 for i in range(14)for j in range(i))
    C=list(map(A.dec,p['strict_center']));TC=A.mv(T,C)
    C=[a+b for a,b in zip(C,TC)]
    assert A.pair(C,C,G).sign()<0 and A.pair(C,H,G).sign()<0
    assert all(A.pair(C,n,G).sign()<0 for n in N)
    facets=[]
    for i,n in enumerate(N):
        v=[x-A.pair(C,n,G)/gram[i][i]*y for x,y in zip(C,n)]
        assert A.pair(v,v,G).sign()<0 and A.pair(v,H,G).sign()<0
        assert all(A.pair(v,m,G).sign()==(0 if i==j else-1)for j,m in enumerate(N))
        facets.append(v)
    E=list(map(A.dec,g['strict_escape']['point']))
    assert A.pair(E,E,G).sign()>0 and A.pair(E,C,G).sign()<0
    assert all(A.pair(E,n,G).sign()<0 for n in N)
    assert all((4*gram[-1][j]*gram[-1][j]/(gram[-1][-1]*gram[j][j])-4).sign()>0
        for j in range(13)if gram[-1][j])
    return A.enc(dict(status='PASS_FINITE_CHAMBER_INVOLUTION_AND_FOURTEEN_ACTUAL_WALLS',
        sources={f:sha(f)for f in [f'scripts/{TAG}.py',HELPER,SOURCE,GEO]},
        source_indices=src,target_indices=dst,target_vertex=2,matrix=T,inverse=T,
        determinant=K(1),trace=K(6),order=2,positive_corner_scale=scale,
        target_basepoint_decomposition=dict(vertex_coefficient=coefficient,normal_coefficients=b),
        normals=N,gram=gram,normal_permutation=permutation,new_wall=R,
        new_wall_height=A.DELTA*R[0]*R[0]/(2*A.pair(R,R,G)),
        strict_center=C,actual_facet_points=facets,
        strict_escape=dict(point=E,norm=A.pair(E,E,G),center_pairing=A.pair(E,C,G),
            wall_pairings=[A.pair(E,n,G)for n in N]),
        scope='One order-two full reflection-chamber symmetry and one actual transported wall. The fourteen-wall chamber has a strict spacelike escape. This is not an infinite stabilizer, a nonreflectivity proof, or a complete fourteen-wall vertex census.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS symmetry saved replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
