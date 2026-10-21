#!/usr/bin/env python3
"""A literal infinite-order symmetry of the full discriminant131 chamber.

Reconstruct one map from two complete outward corner frames. No exhaustive
map search, higher-height enumeration, or classification is claimed.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A

K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta131_h300_20261016_symmetry'
PREFIX='results/fresh_h9_beta131_h300_20261016.json'
GEOMETRY='results/fresh_h9_beta131_h300_20261016_finite.json'
PILOT='results/fresh_h9_beta131_pilot_20261016.json'
POOL='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
SOURCE_VERTEX=1
TARGET_VERTEX=56
SOURCE_IDS=[0,1,2,3,4,5,6,7,10]
TARGET_IDS=[5,2,4,3,1,244,42,11,7]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def tr(M):return list(map(list,zip(*M)))
def mm(X,Y):return [[A.dot(r,c)for c in zip(*Y)]for r in X]
def identity(n):return [[K(i==j)for j in range(n)]for i in range(n)]
def inverse_det(M):
    n=len(M);rows=[r[:]+identity(n)[i]for i,r in enumerate(M)];det=K(1)
    for i in range(n):
        j=next(j for j in range(i,n)if rows[j][i])
        if j!=i:rows[i],rows[j]=rows[j],rows[i];det=-det
        pivot=rows[i][i];det*=pivot;rows[i]=[v/pivot for v in rows[i]]
        for j in range(n):
            if i!=j:
                c=rows[j][i];rows[j]=[a-c*b for a,b in zip(rows[j],rows[i])]
    return [r[n:]for r in rows],det
def pd_pivots(M,conjugate=False):
    a=[[x.conjugate()if conjugate else x for x in row]for row in M];out=[]
    for k in range(len(a)):
        v=a[k][k];assert v.sign()>0;out.append(v)
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/v
    return out

def run():
    d,g,pilot,pool=map(read,(PREFIX,GEOMETRY,PILOT,POOL))
    assert d['status']=='PASS_COMPLETE_BOUNDED_BETA131_PREFIX'
    assert d['height_bound']=='300' and len(d['levels'])==200
    assert g['status']=='PASS_COMPLETE_BETA131_H300_FINITE_CORNER_CENSUS'
    assert len(d['roots'])==282 and len(g['finite_vertices'])==70
    assert len(d['norm_classes'])==2 and len(d['ordinary_zero_roots'])==240
    G=[[A.dec(x)for x in row]for row in d['form']]
    N=[list(map(A.dec,r))for r in d['roots']]
    B=list(map(A.dec,d['basepoint']))
    assert B==[K(6)]+[K()]*8+[K(1)]
    assert G[0][0]==K(0,-2) and G[0][9]==G[9][0]==-1 and G[9][9]==6
    assert A.pair(B,B,G)==K(-6,-72)
    gram=[[A.dec(x)for x in row]for row in d['gram']]
    assert all(gram[i][j].sign()<=0 for i in range(len(N))for j in range(i))
    bpairs=[A.pair(B,n,G)for n in N]
    assert [i for i,x in enumerate(bpairs)if not x]==list(range(8))
    assert all(x.sign()<=0 for x in bpairs)
    S=[[gram[i][j]for j in SOURCE_IDS]for i in SOURCE_IDS]
    assert S==[[gram[i][j]for j in TARGET_IDS]for i in TARGET_IDS]
    assert S==[[A.pair(N[i],N[j],G)for j in SOURCE_IDS]for i in SOURCE_IDS]
    assert S==[[A.pair(N[i],N[j],G)for j in TARGET_IDS]for i in TARGET_IDS]
    assert all(S[i][j].sign()<=0 for i in range(9)for j in range(i))
    pivots=pd_pivots(S);J,detS=inverse_det(S)
    assert all(x.sign()>=0 for row in J for x in row)
    corners=[]
    for vertex,ids in ((SOURCE_VERTEX,SOURCE_IDS),(TARGET_VERTEX,TARGET_IDS)):
        raw=g['finite_vertices'][vertex];V=list(map(A.dec,raw['point']))
        norm=A.pair(V,V,G);bv=A.pair(B,V,G)
        assert norm==A.dec(raw['norm']) and norm.sign()<0 and bv.sign()<0
        vals=[A.pair(V,n,G)for n in N]
        assert [i for i,x in enumerate(vals)if not x]==raw['incident']==sorted(ids)
        assert all(x.sign()<=0 for x in vals)
        a=bv/norm;b=A.mv(J,[bpairs[i]for i in ids])
        assert a.sign()>0 and all(x.sign()<=0 for x in b)
        rebuilt=[a*V[k]+sum((b[j]*N[ids[j]][k]for j in range(9)),K())for k in range(10)]
        assert rebuilt==B
        corners.append(dict(vertex=vertex,incident=raw['incident'],ordered_frame=ids,
            point=V,norm=norm,basepoint_pairing=bv,pairings=vals,
            preservation_decomposition=dict(a=a,coefficients=b,reconstruction=rebuilt)))
    V,P=[r['point']for r in corners]
    F=tr([N[i]for i in SOURCE_IDS]+[V]);FI,detF=inverse_det(F)
    R=tr([N[i]for i in TARGET_IDS]+[P]);_,detR=inverse_det(R)
    scale=detF/detR
    if scale.sign()<0:scale=-scale
    assert scale==K(A.F(-6920,899),A.F(4300,899)) and scale.sign()>0
    for row in R:row[-1]*=scale
    T=mm(R,FI);TI,detT=inverse_det(T)
    assert all(x.integral()for row in T+TI for x in row)
    assert detT==-1 and mm(tr(T),mm(G,T))==G
    assert mm(T,TI)==mm(TI,T)==identity(10)
    for i,j in zip(SOURCE_IDS,TARGET_IDS):assert A.mv(T,N[i])==N[j]
    assert A.mv(T,V)==[scale*x for x in P]
    assert A.pair(V,P,G).sign()<0
    trace=sum((T[i][i]for i in range(10)),K())
    assert trace==K(32,53) and (trace-10).sign()>0
    conjugate_pivots=pd_pivots(G,True)
    C=list(map(A.dec,d['strict_center']))
    assert A.pair(C,C,G).sign()<0 and A.pair(C,B,G).sign()<0
    cpairs=[A.pair(C,n,G)for n in N];assert all(x.sign()<0 for x in cpairs)
    sources={p:sha(p)for p in (f'scripts/{TAG}.py',HELPER,PREFIX,GEOMETRY,PILOT,POOL)}
    return A.enc(dict(status='PASS_EXACT_INFINITE_ORDER_BETA131_CHAMBER_SYMMETRY',sources=sources,
        arithmetic='Fraction arithmetic in basis 1,phi; map reconstructed from literal root-plus-apex frames',
        source_corner=corners[0],target_corner=corners[1],basepoint=B,basepoint_pairings=bpairs,
        common_spherical_gram=S,pivots=pivots,inverse_gram=J,gram_determinant=detS,
        frames=dict(source=F,source_inverse=FI,source_determinant=detF,
            target_determinant_before_scaling=detR,positive_perpendicular_scale=scale),
        matrix=T,inverse_matrix=TI,determinant=detT,trace=trace,trace_minus_dimension=trace-10,
        conjugate_form_positive_pivots=conjugate_pivots,strict_center=C,strict_center_pairings=cpairs,
        exact_checks=dict(root_images=9,source_and_target_wall_signs=564,isometry_entries=100,
            integral_matrix_and_inverse_entries=200,finite_corner_decomposition_coordinates=20),
        theorem='The literal discriminant131 lattice is nonreflective: its full integral reflection chamber has an infinite-order integral future stabilizer and infinite volume.',
        scope=dict(norm_alphabet='Both primitive reflection norm classes, at arbitrary height; completeness inherited from the height20 norm/zero-sector proof.',
            prefix_dependency='Complete height300 Vinberg prefix. Future simple roots cannot cut either checked finite corner.',
            witness_search='One exact map. No exhaustive corner-map or full stabilizer classification claimed.',
            host_wide_consequence='No finite-volume Coxeter polytope can have all wall reflections preserving this literal lattice.',
            limitation='Other lattices and arbitrary nonintegral/geometric reflection walls are not excluded.')))

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved exact beta131 symmetry replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS beta131 infinite-order chamber symmetry: trace 32+53phi')
if __name__=='__main__':main()
