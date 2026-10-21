#!/usr/bin/env python3
"""Independent Golden reconstruction of the beta55 infinite chamber symmetry."""
import argparse,json
from fractions import Fraction as F
import fresh_h9_beta55_vinberg_h100_20261014_review as R
K,ROOT,sha,enc,read,pair,mv,dot=R.K,R.ROOT,R.sha,R.enc,R.read,R.pair,R.mv,R.dot
TAG='fresh_h9_beta55_corner_symmetries_20261014_review'
PREFIX='results/fresh_h9_beta55_vinberg_h100_20261014.json'
REVIEW='results/fresh_h9_beta55_vinberg_h100_20261014_review.json'
GEOM='results/fresh_h9_beta55_geometry_h100_20261014.json'
PRIMARY='results/fresh_h9_beta55_corner_symmetries_20261014.json'
def dm(m):return [[K.from_json(x)for x in r]for r in m]
def dv(v):return [K.from_json(x)for x in v]
def transpose(m):return [list(x)for x in zip(*m)]
def mm(a,b):return [[dot(r,c)for c in zip(*b)]for r in a]
def inverse(m):
    n=len(m);a=[r[:]+[K(i==j)for j in range(n)]for i,r in enumerate(m)]
    det=K(1)
    for j in range(n):
        i=next(i for i in range(j,n)if a[i][j])
        if i!=j:a[j],a[i]=a[i],a[j];det=-det
        q=a[j][j];det*=q;a[j]=[z/q for z in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return [r[n:]for r in a],det

def run():
    d,g,review,primary=map(read,(PREFIX,GEOM,REVIEW,PRIMARY))
    assert review['status']=='PASS_INDEPENDENT_BETA55_T100_REVIEW'
    G,N,B=dm(d['form']),dm(d['roots']),dv(d['basepoint']);H=dv(d['strict_center'])
    assert len(N)==44 and len(g['finite_vertices'])==129 and B==[K(4)]+[K()]*8+[K(1)]
    source,target=3,42;images=[25,0,35,2,3,4,5,24,7];scale=K(2,2)
    left,right=g['finite_vertices'][source],g['finite_vertices'][target]
    ids=left['incident'];v,w=dv(left['point']),dv(right['point'])
    assert ids==[0,1,2,3,4,5,6,9,10]and set(images)==set(right['incident'])
    U=transpose([v]+[N[i]for i in ids]);V=transpose([[scale*x for x in w]]+[N[i]for i in images])
    Uinv,_=inverse(U);Vinv,_=inverse(V);T=mm(V,Uinv);Ti=mm(U,Vinv)
    I=[[K(i==j)for j in range(10)]for i in range(10)]
    assert mm(T,Ti)==I==mm(Ti,T)and all(z.is_integral for row in T+Ti for z in row)
    assert mm(mm(transpose(T),G),T)==G
    assert mv(T,v)==[scale*x for x in w]and scale>0
    assert all(mv(T,N[i])==N[j]for i,j in zip(ids,images))
    assert pair(mv(T,H),H,G)<0
    trace=sum((T[i][i]for i in range(10)),K());_,det=inverse(T)
    assert trace==K(5,8)and trace>10 and det==-1
    # The matrix is reconstructed from the literal complete corner sectors.
    # Equality with the producer's serialized witness is checked below after
    # locating its one explicitly named matrix field.
    matrices=[]
    def scan(x):
        if isinstance(x,dict):
            for k,z in x.items():
                if k in ('matrix','isometry','T')and isinstance(z,list)and len(z)==10:
                    try:matrices.append(dm(z))
                    except (TypeError,KeyError):pass
                scan(z)
        elif isinstance(x,list):
            for z in x:scan(z)
    scan(primary);assert T in matrices
    assert dm(primary['inverse_matrix'])==Ti and K.from_json(primary['trace'])==trace
    assert K.from_json(primary['determinant'])==det
    inertias=[]
    for embedding in (0,1):
        M=[[z if embedding==0 else z.conjugate()for z in r]for r in G];signs=[]
        for k in range(10):
            q=M[k][k];assert q;signs.append(q.sign())
            for i in range(k+1,10):
                for j in range(k+1,10):M[i][j]-=M[i][k]*M[k][j]/q
        inertias.append([signs.count(1),signs.count(-1)])
    assert inertias==[[9,1],[10,0]]
    corners=[]
    for row in (left,right):
        ids=row['incident'];p=dv(row['point']);norm=pair(p,p,G)
        assert norm<0 and pair(p,B,G)<0
        vals=[pair(p,n,G)for n in N]
        assert all(z<=0 for z in vals)and [i for i,z in enumerate(vals)if not z]==ids
        S=[[pair(N[i],N[j],G)for j in ids]for i in ids]
        Si,_=inverse(S);assert all(z>=0 for r in Si for z in r)
        M=[r[:]for r in S]
        for k in range(9):
            q=M[k][k];assert q>0
            for i in range(k+1,9):
                for j in range(k+1,9):M[i][j]-=M[i][k]*M[k][j]/q
        a=pair(B,p,G)/norm;b=mv(Si,[pair(B,N[i],G)for i in ids])
        assert a>0 and all(z<=0 for z in b)
        reconstructed=[a*x+sum((c*N[i][k]for c,i in zip(b,ids)),K())for k,x in enumerate(p)]
        assert reconstructed==B
        corners.append(dict(incident=ids,point=p,norm=norm,basepoint_vertex_coefficient=a,
                            basepoint_normal_coefficients=b,inverse_incident_gram=Si))
    paths=[PREFIX,GEOM,REVIEW,PRIMARY,f'scripts/{TAG}.py',
           'scripts/fresh_h9_beta55_vinberg_h100_20261014_review.py','src/coxeter_polytopes/exact/golden.py']
    return enc(dict(status='PASS_INDEPENDENT_BETA55_FULL_CHAMBER_INFINITE_SYMMETRY',sources={p:sha(p)for p in paths},
        source_vertex=source,target_vertex=target,outward_normal_images=images,positive_scale=scale,
        matrix=T,inverse=Ti,determinant=det,trace=trace,embedding_inertias=inertias,complete_corners=corners,
        complete_sector_lemma='B=a*v+sum(b_i*n_i), a>0,b_i<=0. Every later simple root has B.r<0 and n_i.r<=0, so a*v.r=B.r-sum(b_i*n_i.r)<0. Thus neither corner is cut or touched by a later wall. Local finiteness gives the complete chamber sector.',
        chamber_consequence='Integral future T permutes the full reflection arrangement and maps one complete outward corner sector to the other, hence stabilizes the full chamber. Trace>10 proves infinite order. Discreteness supplies infinitely many disjoint equal-volume balls within that chamber.',
        scope='The full integral reflection group of this literal lattice has infinite covolume; no finite-volume Coxeter polytope has all wall reflections preserving it. Other lattices or nonintegral walls are outside this conclusion.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['trace'],d['determinant'])
