#!/usr/bin/env python3
"""One exact infinite-order symmetry of the full norm31 reflection chamber."""
import argparse,json
from pathlib import Path
from hashlib import sha256
import fresh_frontier_construction_20261015 as A
import fresh_frontier_construction_20261015_symmetry as B
K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_norm31_continuation_20261015_symmetry'
SOURCE='results/fresh_norm31_continuation_20261015.json'
CORNERS='results/fresh_norm31_continuation_20261015_corners.json'
HELPER='scripts/fresh_frontier_construction_20261015.py'
MATRIX='scripts/fresh_frontier_construction_20261015_symmetry.py'
SRC=list(range(9));DST=[5,1,4,3,2,0,21,20,7]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pd(M):
    M=[r[:]for r in M];p=[]
    for k in range(len(M)):
        q=M[k][k];assert q.sign()>0;p.append(q)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):M[i][j]-=M[i][k]*M[k][j]/q
    return p
def run():
    d=json.loads((ROOT/SOURCE).read_text());c=json.loads((ROOT/CORNERS).read_text())
    for x in (d,c):
        for f,h in x['sources'].items():assert sha(f)==h,f
    assert d['height_bound']=='100'and len(d['levels'])==32 and len(d['roots'])==28
    G=[[A.dec(x)for x in r]for r in d['form']]
    N=[list(map(A.dec,r))for r in d['roots']];H=list(map(A.dec,d['basepoint']))
    raw=next(v for v in c['finite_vertices']if v['incident']==sorted(DST))
    V=list(map(A.dec,raw['point']))
    SG=[[A.pair(N[i],N[j],G)for j in SRC]for i in SRC]
    assert SG==[[A.pair(N[i],N[j],G)for j in DST]for i in DST]
    pivots=pd(SG);J=B.inverse(SG);assert all(x.sign()>=0 for r in J for x in r)
    corners=[]
    for point,ids in ((H,SRC),(V,DST)):
        norm=A.pair(point,point,G);hp=A.pair(H,point,G)
        assert norm.sign()<0 and hp.sign()<0
        vals=[A.pair(point,n,G)for n in N]
        assert [i for i,x in enumerate(vals)if not x]==sorted(ids)
        assert all(x.sign()<=0 for x in vals)
        a=hp/norm;b=A.mv(J,[A.pair(H,N[i],G)for i in ids])
        assert a.sign()>0 and all(x.sign()<=0 for x in b)
        rebuilt=[a*point[k]+sum(b[i]*N[ids[i]][k]for i in range(9))for k in range(10)]
        assert rebuilt==H
        corners.append(dict(ordered_indices=ids,point=point,norm=norm,basepoint_pairing=hp,
            wall_pairings=vals,preservation_decomposition=dict(vertex_coefficient=a,normal_coefficients=b)))
    F=B.transpose([N[i]for i in SRC]+[H]);D=B.transpose([N[i]for i in DST]+[V])
    scale=B.determinant(F)/B.determinant(D)
    if scale.sign()<0:scale=-scale
    assert scale==18+14*A.S and scale.sign()>0
    for r in D:r[-1]*=scale
    T=B.mm(D,B.inverse(F));TI=B.inverse(T)
    assert all(x.integral()for r in T+TI for x in r)
    assert B.mm(B.mm(B.transpose(T),G),T)==G
    I=[[K(i==j)for j in range(10)]for i in range(10)]
    assert B.mm(T,TI)==B.mm(TI,T)==I
    assert B.determinant(T)==-1
    assert all(A.mv(T,N[i])==N[j]for i,j in zip(SRC,DST))
    assert A.mv(T,H)==[scale*x for x in V]
    trace=sum(T[i][i]for i in range(10))
    assert trace==9+5*A.S and (trace-10).sign()>0
    conjugate_pivots=pd([[x.conjugate()for x in r]for r in G])
    center=list(map(A.dec,d['strict_center']))
    assert A.pair(center,center,G).sign()<0 and A.pair(center,H,G).sign()<0
    assert all(A.pair(center,n,G).sign()<0 for n in N)
    return A.enc(dict(status='PASS_EXACT_INFINITE_ORDER_NORM31_CHAMBER_SYMMETRY',
        sources={f:sha(f)for f in [f'scripts/{TAG}.py',HELPER,MATRIX,SOURCE,CORNERS]},
        source_corner=corners[0],target_corner=corners[1],common_spherical_gram=SG,
        spherical_pivots=pivots,inverse_spherical_gram=J,positive_perpendicular_scale=scale,
        matrix=T,inverse_matrix=TI,determinant=K(-1),trace=trace,trace_minus_ten=trace-10,
        conjugate_form_positive_pivots=conjugate_pivots,strict_center=center,
        theorem='The literal discriminant31 lattice is nonreflective: its full integral reflection chamber has an infinite-order integral future stabilizer and infinite volume.',
        scope=dict(norm_alphabet='All primitive integral reflectors, using the frozen complete two-class norm alphabet and zero sector.',
            corner_persistence='Both complete outward spherical sectors are corners of the full reflection chamber by the saved exact basepoint decompositions and acuteness.',
            search='One exact witness, not an exhaustive chamber-automorphism census.',
            host_wide='No finite-volume Coxeter chamber can have all wall reflections preserving this literal O-lattice.',
            limitation='No exclusion of other lattices, rational similarity classes in general, or nonintegral geometric reflection walls.')))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS infinite symmetry saved replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
