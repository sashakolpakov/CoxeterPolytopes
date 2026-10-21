#!/usr/bin/env python3
"""Two literal integral corner maps whose product stabilizes the full chamber."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A

K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_witt_twist_pilot_20261016_symmetry'
PREFIX='results/fresh_h9_witt_twist_pilot_20261016.json'
GEOMETRY='results/fresh_h9_witt_twist_pilot_20261016_geometry.json'
SETUP='results/fresh_h9_witt_twist_20261016.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
RECIPES=[
    (11,24,[0,1,2,4,6,7,8,9,11],[4,3,2,0,6,7,8,11,9],K(3,4),K(-1),K(4)),
    (2,22,[0,1,2,3,4,5,6,8,10],[10,14,2,3,4,5,6,8,0],K(2,2),K(1),K(6)),
]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def tr(M):return list(map(list,zip(*M)))
def mm(X,Y):return [[A.dot(r,c) for c in zip(*Y)] for r in X]
def identity(n):return [[K(i==j) for j in range(n)] for i in range(n)]
def inverse_det(M):
    n=len(M);rows=[r[:]+identity(n)[i] for i,r in enumerate(M)];det=K(1)
    for i in range(n):
        j=next(j for j in range(i,n) if rows[j][i])
        if j!=i:rows[i],rows[j]=rows[j],rows[i];det=-det
        pivot=rows[i][i];det*=pivot;rows[i]=[v/pivot for v in rows[i]]
        for j in range(n):
            if i!=j:
                c=rows[j][i];rows[j]=[a-c*b for a,b in zip(rows[j],rows[i])]
    return [r[n:] for r in rows],det
def pd_pivots(M):
    a=[r[:] for r in M];out=[]
    for k in range(len(a)):
        v=a[k][k];assert v.sign()>0;out.append(v)
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/v
    return out

def run():
    d,g,setup=map(read,(PREFIX,GEOMETRY,SETUP))
    assert d['status']=='PASS_COMPLETE_BOUNDED_WITT_TWIST_PREFIX'
    assert d['height_bound']=='20' and d['counts']['controls']==13 and len(d['roots'])==16
    assert g['status']=='PASS_COMPLETE_WITT_TWIST_PILOT_CONE'
    phi=K(0,1);beta=4-phi;diagonal=[-phi]+[K(1)]*7+[beta]*2
    G=[[diagonal[i] if i==j else K() for j in range(10)] for i in range(10)]
    assert G==[[A.dec(x) for x in row] for row in d['form']]
    assert G==[[A.dec(x) for x in row] for row in setup['form']]
    assert list(map(A.dec,d['norm_classes']))==[K(1),K(2),beta,2*beta]
    assert setup['norm_and_zero']['zero_root_count']==106
    N=[list(map(A.dec,r)) for r in d['roots']]
    H=list(map(A.dec,d['basepoint']));assert H==[K(1)]+[K()]*9
    gram=[[A.pair(a,b,G) for b in N] for a in N]
    assert gram==[[A.dec(x) for x in row] for row in d['gram']]
    assert all(gram[i][j].sign()<=0 for i in range(16) for j in range(i))
    hp=[A.pair(H,n,G) for n in N]
    assert [i for i,v in enumerate(hp) if not v]==list(range(9))
    assert all(v.sign()<0 for v in hp[9:])
    generators=[]
    for sv,tv,sids,tids,expected_scale,expected_det,expected_trace in RECIPES:
        S=[[gram[i][j] for j in sids] for i in sids]
        assert S==[[gram[i][j] for j in tids] for i in tids]
        pivots=pd_pivots(S);J,detS=inverse_det(S)
        assert all(v.sign()>=0 for row in J for v in row)
        corners=[]
        for vertex,ids in ((sv,sids),(tv,tids)):
            raw=g['finite_vertices'][vertex];V=list(map(A.dec,raw['point']))
            norm=A.pair(V,V,G);hv=A.pair(H,V,G)
            assert norm==A.dec(raw['norm']) and norm.sign()<0 and hv.sign()<0
            vals=[A.pair(V,n,G) for n in N]
            assert [i for i,v in enumerate(vals) if not v]==raw['incident']==sorted(ids)
            assert all(v.sign()<=0 for v in vals)
            a=hv/norm;b=A.mv(J,[hp[i] for i in ids])
            assert a.sign()>0 and all(v.sign()<=0 for v in b)
            rebuilt=[a*V[k]+sum((b[j]*N[ids[j]][k] for j in range(9)),K()) for k in range(10)]
            assert rebuilt==H
            corners.append(dict(vertex=vertex,incident=raw['incident'],ordered_frame=ids,
                point=V,norm=norm,basepoint_pairing=hv,pairings=vals,
                preservation_decomposition=dict(a=a,coefficients=b,reconstruction=rebuilt)))
        V,P=[r['point'] for r in corners]
        F=tr([N[i] for i in sids]+[V]);FI,detF=inverse_det(F)
        R=tr([N[i] for i in tids]+[P]);_,detR=inverse_det(R)
        scale=detF/detR
        if scale.sign()<0:scale=-scale
        assert scale==expected_scale and scale.sign()>0
        assert A.pair(V,V,G)==scale*scale*A.pair(P,P,G)
        for row in R:row[-1]*=scale
        T=mm(R,FI);TI,detT=inverse_det(T)
        assert all(v.integral() for row in T+TI for v in row)
        assert detT==expected_det and mm(tr(T),mm(G,T))==G
        assert mm(T,TI)==mm(TI,T)==identity(10)
        assert all(A.mv(T,N[i])==N[j] for i,j in zip(sids,tids))
        assert A.mv(T,V)==[scale*v for v in P]
        assert A.pair(V,P,G).sign()<0
        trace=sum((T[i][i] for i in range(10)),K());assert trace==expected_trace
        generators.append(dict(source_corner=corners[0],target_corner=corners[1],
            common_spherical_gram=S,positive_pivots=pivots,inverse_gram=J,gram_determinant=detS,
            positive_perpendicular_scale=scale,matrix=T,inverse_matrix=TI,determinant=detT,trace=trace))
    T=mm(generators[0]['matrix'],generators[1]['matrix']);TI,detT=inverse_det(T)
    assert all(v.integral() for row in T+TI for v in row)
    assert mm(tr(T),mm(G,T))==G and mm(T,TI)==mm(TI,T)==identity(10)
    trace=sum((T[i][i] for i in range(10)),K())
    assert trace==K(50,76) and (trace-10).sign()>0 and detT==-1
    assert A.pair(H,A.mv(T,H),G).sign()<0
    conjugate_pivots=pd_pivots([[v.conjugate() for v in row] for row in G])
    C=list(map(A.dec,d['strict_center']));cp=[A.pair(C,n,G) for n in N]
    assert A.pair(C,C,G).sign()<0 and A.pair(C,H,G).sign()<0 and all(v.sign()<0 for v in cp)
    sources={p:sha(p) for p in (f'scripts/{TAG}.py',HELPER,PREFIX,GEOMETRY,SETUP)}
    return A.enc(dict(status='PASS_EXACT_INFINITE_ORDER_WITT_TWIST_CHAMBER_SYMMETRY',sources=sources,
        arithmetic='Fraction arithmetic in 1,phi; each map reconstructed from literal nine-root-plus-apex frames',
        form=G,basepoint=H,basepoint_pairings=hp,generators=generators,
        product_order=[0,1],matrix=T,inverse_matrix=TI,determinant=detT,trace=trace,
        trace_minus_dimension=trace-10,conjugate_form_positive_pivots=conjugate_pivots,
        strict_center=C,strict_center_pairings=cp,
        exact_checks=dict(outward_root_images=18,four_corner_wall_signs=64,
            integral_matrix_and_inverse_entries=600,isometry_entries=300,basepoint_decomposition_coordinates=40),
        theorem='The literal diagonal lattice <-phi> + I7 + (4-phi) I2 is nonreflective: its full integral reflection chamber has an infinite-order integral future stabilizer and infinite volume.',
        scope=dict(norm_alphabet='All four primitive reflection norm classes at every height, from the bound setup.',
            prefix_dependency='Complete height20 Vinberg prefix. Every later simple wall is strict at the four certified finite corners.',
            witness='Two exact generators and their product. Completeness of the finite-corner census or map screen is not needed.',
            host_wide_consequence='No finite-volume Coxeter polytope has all wall reflections preserving this literal lattice.',
            limitation='Other lattices and arbitrary nonintegral geometric reflections are not excluded.')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('PASS saved exact Witt-twist symmetry replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS Witt-twist infinite-order symmetry: trace 50+76phi')
