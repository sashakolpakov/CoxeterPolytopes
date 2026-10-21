#!/usr/bin/env python3
"""Exact infinite-order symmetry of the full beta11 reflection chamber.

The finite corner is checked directly; no exhaustive corner-map search is
claimed or required. Arithmetic helper is the frozen independent Fraction
implementation used in the complete height300 replay.
"""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A

K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta11_corner_symmetries_20261014'
PREFIX='results/fresh_h9_beta11_vinberg_h300_20261014.json'
GEOMETRY='results/fresh_h9_beta11_geometry_h300_20261014.json'
REPLAY='results/fresh_h9_beta11_vinberg_h300_20261014_review.json'
ZERO='results/fresh_h9_beta11_zero_sector_20261014.json'
POOL='results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
SOURCE_IDS=list(range(9))
TARGET_IDS=[7,4,6,5,20,29,30,11,2]
TARGET_VERTEX=735

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
                coeff=rows[j][i];rows[j]=[a-coeff*b for a,b in zip(rows[j],rows[i])]
    return [r[n:]for r in rows],det

def pd_pivots(M,conjugate=False):
    a=[[x.conjugate()if conjugate else x for x in row]for row in M];out=[]
    for k in range(len(a)):
        v=a[k][k];assert v.sign()>0;out.append(v)
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/v
    return out

def run():
    d,g,review,zero,pool=map(read,(PREFIX,GEOMETRY,REPLAY,ZERO,POOL))
    assert d['status']=='PASS_COMPLETE_BOTH_NORM_VINBERG_PREFIX'
    assert g['status']=='PASS_COMPLETE_BETA11_PREFIX_CONE_INFINITE_VOLUME'
    assert review['status']=='PASS_INDEPENDENT_COMPLETE_BOTH_NORM_HEIGHT300_REPLAY'
    assert len(d['roots'])==35 and len(g['finite_vertices'])==880
    G=[[A.dec(x)for x in row]for row in d['form']]
    N=[list(map(A.dec,r))for r in d['roots']]
    H=list(map(A.dec,d['basepoint']))
    P=list(map(A.dec,g['finite_vertices'][TARGET_VERTEX]['point']))
    target_incident=g['finite_vertices'][TARGET_VERTEX]['incident']
    assert set(target_incident)==set(TARGET_IDS) and len(target_incident)==9
    assert H==[K(2)]+[K()]*8+[K(1)]
    assert H==list(map(A.dec,g['finite_vertices'][0]['point']))
    assert g['finite_vertices'][0]['incident']==SOURCE_IDS
    gram=[[A.pair(r,s,G)for s in N]for r in N]
    assert gram==[[A.dec(x)for x in row]for row in d['gram']]
    assert all(gram[i][j].sign()<=0 for i in range(35)for j in range(i))
    hnorm=A.pair(H,H,G);pnorm=A.pair(P,P,G);hp=A.pair(H,P,G)
    assert hnorm.sign()<0 and pnorm.sign()<0 and hp.sign()<0
    assert pnorm==A.dec(g['finite_vertices'][TARGET_VERTEX]['norm'])
    hpairs=[A.pair(H,n,G)for n in N];ppairs=[A.pair(P,n,G)for n in N]
    assert [i for i,x in enumerate(hpairs)if not x]==SOURCE_IDS
    assert [i for i,x in enumerate(ppairs)if not x]==target_incident
    assert all(x.sign()<=0 for x in hpairs+ppairs)
    S=[[gram[i][j]for j in SOURCE_IDS]for i in SOURCE_IDS]
    Tgram=[[gram[i][j]for j in TARGET_IDS]for i in TARGET_IDS]
    assert S==Tgram
    pivots=pd_pivots(S);J,detS=inverse_det(S)
    assert all(x.sign()>=0 for row in J for x in row)
    # Check the finite-vertex preservation decomposition physically.
    a=hp/pnorm;b=A.mv(J,[hpairs[i]for i in TARGET_IDS])
    assert a.sign()>0 and all(x.sign()<=0 for x in b)
    decomposition=[a*x+sum((b[j]*N[TARGET_IDS[j]][i]for j in range(9)),K())for i,x in enumerate(P)]
    assert decomposition==H
    # Reconstruct, rather than import, the complete map from its root images.
    F=tr([N[i]for i in SOURCE_IDS]+[H]);FI,detF=inverse_det(F)
    R=tr([N[i]for i in TARGET_IDS]+[P]);_,detR=inverse_det(R)
    scale=detF/detR
    if scale.sign()<0:scale=-scale
    assert scale==K(8,14) and scale.sign()>0
    for row in R:row[-1]*=scale
    T=mm(R,FI);TI,detT=inverse_det(T)
    assert all(x.integral()for row in T+TI for x in row)
    assert detT==1 and mm(tr(T),mm(G,T))==G
    assert mm(T,TI)==identity(10)
    for i,j in zip(SOURCE_IDS,TARGET_IDS):assert A.mv(T,N[i])==N[j]
    assert A.mv(T,H)==[scale*x for x in P]
    trace=sum((T[i][i]for i in range(10)),K())
    assert trace==K(4,8) and (trace-10).sign()>0
    assert all((G[i][i]/2).integral()for i in range(10))
    conjugate_pivots=pd_pivots(G,True)
    # A concrete strict chamber point is retained as a check that our side
    # conventions agree with the source prefix.
    C=list(map(A.dec,d['strict_center']))
    assert A.pair(C,C,G).sign()<0 and A.pair(C,H,G).sign()<0
    cpairs=[A.pair(C,n,G)for n in N];assert all(x.sign()<0 for x in cpairs)
    sources={p:sha(p)for p in (f'scripts/{TAG}.py',HELPER,PREFIX,GEOMETRY,REPLAY,ZERO,POOL)}
    return A.enc(dict(status='PASS',sources=sources,
        arithmetic='Fraction arithmetic in the basis 1,phi; literal frame reconstruction',
        source_corner=dict(vertex=0,incident=SOURCE_IDS,point=H,norm=hnorm,pairings=hpairs),
        target_corner=dict(vertex=TARGET_VERTEX,incident=target_incident,ordered_images=TARGET_IDS,
            point=P,norm=pnorm,future_pairing=hp,pairings=ppairs),
        common_spherical_gram=S,pivots=pivots,inverse_gram=J,gram_determinant=detS,
        preservation_decomposition=dict(a=a,coefficients=b,reconstruction=decomposition),
        frames=dict(source=F,source_inverse=FI,source_determinant=detF,
            target_determinant_before_scaling=detR,positive_perpendicular_scale=scale),
        matrix=T,inverse_matrix=TI,determinant=detT,trace=trace,trace_minus_dimension=trace-10,
        conjugate_form_positive_pivots=conjugate_pivots,strict_center=C,strict_center_pairings=cpairs,
        exact_checks=dict(root_images=9,source_and_target_wall_signs=70,isometry_entries=100,
            integral_matrix_and_inverse_entries=200,finite_corner_decomposition_coordinates=10),
        theorem='An infinite-order future integral isometry stabilizes the full integral reflection chamber; the lattice is nonreflective.',
        scope=dict(full_reflection_group='Both complete primitive norm classes 2 and 2*beta, at all heights.',
            prefix_dependency='Complete Vinberg prefix through height300, separately replayed; later walls cannot cut either finite corner.',
            witness_search='One exact corner map; no exhaustive symmetry-group or corner-map classification claimed.',
            geometric_limit='This closes reflectivity of this literal arithmetic lattice, not existence of arbitrary compact H9 or nonarithmetic H20 Coxeter polytopes.')))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=ROOT/'results'/f'{TAG}.json';d=run()
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved exact chamber symmetry replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS exact infinite-order chamber symmetry; trace=4+8phi')
if __name__=='__main__':main()
