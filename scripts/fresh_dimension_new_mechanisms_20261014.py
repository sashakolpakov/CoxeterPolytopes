#!/usr/bin/env python3
"""Exact compact H7 source and its coupled B2-supported H9 lift.

This is a geometric construction/census, not a compactness claim.
All field operations use the separately frozen Q(sqrt5) implementation.
"""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from fresh_h9_makarov_face_extension_20261013 import K,Z,O,PHI,ident,inv,null1,inertia,mv,pair,dot,matrix_enc
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_dimension_new_mechanisms_20261014'
EDGES=[(0,1,5),(1,2,3),(2,3,3),(3,4,3),(4,5,3),(5,6,5),(3,7,3),(1,9,3),(5,10,3),(8,9,4),(8,10,4)]
def enc(x):return [a.enc() for a in x]
def smv(V):return [sum((v[i] for v in V),Z) for i in range(len(V[0]))]
def scale(a,v):return [a*x for x in v]
def add(a,b):return [x+y for x,y in zip(a,b)]
def classify(p,a,b):
    if not p:return '2'
    assert p.sign()<0
    c=p*p/(a*b)
    if c==K(F(1,4)):return '3'
    if c==K(F(1,2)):return '4'
    if c==K(F(3,8),F(1,8)):return '5'
    if c==K(F(3,4)):return '6'
    if c==K(F(5,8),F(1,8)):return '10'
    if c==O:return 'parallel'
    assert (c-O).sign()>0
    return 'dotted'
def build():
    B=[[2*x for x in row]for row in ident(8)]
    for a,b,m in EDGES:
        if a<8 and b<8:B[a][b]=B[b][a]=-({3:O,5:PHI}[m])
    assert inertia(B)==[7,1]
    J=inv(B)
    V=ident(8)+[[-PHI*J[i][7] for i in range(8)],[-J[i][1] for i in range(8)],[-J[i][5] for i in range(8)]]
    q=[pair(v,v,B)for v in V]
    assert q==[K(2)]*8+[O,K(2),K(2)]
    G=[[pair(v,w,B)for w in V]for v in V]
    assert G[7][8]==-PHI and G[9][10]==-2*PHI
    expected={tuple(sorted((a,b))):str(m)for a,b,m in EDGES}
    expected[(7,8)]=expected[(9,10)]='dotted'
    assert all(classify(G[i][j],q[i],q[j])==expected.get((i,j),'2')for i,j in combinations(range(11),2))
    # A timelike positive sum makes its -1 projective slice bounded.
    dual=smv(V);assert pair(dual,dual,B).sign()<0
    verts=[]
    for ids in combinations(range(11),7):
        v=null1([mv(B,V[i])for i in ids])
        if v is None:continue
        ps=[pair(v,w,B)for w in V]
        if all(x.sign()>=0 for x in ps):v=scale(K(-1),v);ps=[-x for x in ps]
        if any(x.sign()>0 for x in ps):continue
        den=-pair(dual,v,B);assert den.sign()>0
        v=scale(1/den,v)
        assert pair(v,v,B).sign()<0
        active=[i for i,x in enumerate(ps)if not x]
        assert active==list(ids)
        assert inertia([[G[i][j]for j in ids]for i in ids])==[7,0]
        verts.append(dict(point=v,incident=active))
    h=smv([v['point']for v in verts])
    assert pair(h,h,B).sign()<0 and all(pair(h,v,B).sign()<0 for v in V)
    # Ambient B2 basis s0^2=2,s1^2=1,s0.s1=-1; d0=s0+s1, d1=s0+2s1.
    Q=[row+[Z,Z]for row in B]+[[Z]*8+[K(2),K(-1)],[Z]*8+[K(-1),O]]
    W=[v+[Z,Z]for v in V]
    W[8][-2:]=[K(-1),K(-1)]
    W += [[Z]*8+[O,Z],[Z]*8+[Z,O]]
    nq=[pair(w,w,Q)for w in W]
    assert nq==[K(2)]*12+[O] and inertia(Q)==[9,1]
    GG=[[pair(v,w,Q)for w in W]for v in W]
    labels={(i,j):classify(GG[i][j],nq[i],nq[j])for i,j in combinations(range(13),2)}
    assert labels[(7,8)]=='5' and labels[(8,9)]==labels[(8,10)]==labels[(8,11)]=='3' and labels[(11,12)]=='4'
    eps=-pair(h,V[8],B)/4
    H=h+[-2*eps,-3*eps]
    # This perturbation can be shrunk until it is timelike.
    while pair(H,H,Q).sign()>=0:
        eps/=2;H=h+[-2*eps,-3*eps]
    assert all(pair(H,w,Q).sign()<0 for w in W)
    rays=[]
    for vi,vrec in enumerate(verts):
        v=vrec['point'];L=-pair(v,V[8],B);assert L.sign()>=0
        options=[(0,[Z,Z])] if not L else [(0,[Z,Z]),(1,[-L,-L]),(2,[-L,-2*L])]
        for fiber,u in options:
            x=v+u;ps=[pair(x,w,Q)for w in W];assert all(p.sign()<=0 for p in ps)
            active=[i for i,p in enumerate(ps)if not p]
            assert len(active)==9
            assert null1([mv(Q,W[i])for i in active]) is not None
            norm=pair(x,x,Q)
            kind={-1:'finite',0:'ideal',1:'spacelike'}[norm.sign()]
            if kind=='finite':assert inertia([[GG[i][j]for j in active]for i in active])==[9,0]
            rays.append(dict(source_vertex=vi,fiber=fiber,point=x,incident=active,norm=norm,kind=kind))
    assert len({tuple(r['incident'])for r in rays})==len(rays)
    outside=next(r['point']for r in rays if r['kind']=='spacelike')
    a=K(F(1,1024));E=add(outside,scale(a,H))
    while pair(E,E,Q).sign()<=0:a/=2;E=add(outside,scale(a,H))
    assert all(pair(E,w,Q).sign()<0 for w in W) and pair(E,H,Q).sign()<0
    return B,V,G,verts,Q,W,GG,labels,H,rays,E

def rank(A):
    M=[list(r)for r in A];k=0
    for j in range(len(M[0]) if M else 0):
        z=next((i for i in range(k,len(M))if M[i][j]),None)
        if z is None:continue
        M[k],M[z]=M[z],M[k];d=M[k][j]
        for i in range(k+1,len(M)):
            if M[i][j]:
                t=M[i][j]/d
                M[i]=[x-t*y for x,y in zip(M[i],M[k])]
        k+=1
        if k==len(M):break
    return k

def clip_cone(Q,N,rays,caps):
    dual=smv(N);assert pair(dual,dual,Q).sign()<0
    def normal(x):
        d=-pair(x,dual,Q);assert d.sign()>0
        return scale(1/d,x)
    current=[dict(point=normal(r['point']),incident=r['incident'])for r in rays]
    trace=[];walls=list(N)
    for cap in caps:
        prod=[pair(r['point'],cap,Q)for r in current]
        pos=[i for i,p in enumerate(prod)if p.sign()>0]
        neg=[i for i,p in enumerate(prod)if p.sign()<0]
        nxt={};crossed=0;rank_cache={}
        def insert(v):
            v=normal(v);key=tuple((x.a,x.b)for x in v)
            ps=[pair(v,n,Q)for n in walls+[cap]]
            assert all(x.sign()<=0 for x in ps)
            active=[i for i,p in enumerate(ps)if not p]
            assert rank([(walls+[cap])[i]for i in active])==9
            nxt[key]=dict(point=v,incident=active)
        for i,p in enumerate(prod):
            if p.sign()<=0:insert(current[i]['point'])
        for i in pos:
            si=set(current[i]['incident'])
            for j in neg:
                common=tuple(sorted(si.intersection(current[j]['incident'])))
                if len(common)<8:continue
                if common not in rank_cache:rank_cache[common]=rank([walls[k]for k in common])
                if rank_cache[common]!=8:continue
                crossed+=1
                v=add(scale(-prod[j],current[i]['point']),scale(prod[i],current[j]['point']))
                insert(v)
        trace.append(dict(before=len(current),positive=len(pos),negative=len(neg),zero=len(current)-len(pos)-len(neg),crossing_edges=crossed,after=len(nxt)))
        walls.append(cap);current=list(nxt.values())
        print('clip',len(walls),trace[-1],flush=True)
    for r in current:
        r['norm']=pair(r['point'],r['point'],Q)
        r['kind']={-1:'finite',0:'ideal',1:'spacelike'}[r['norm'].sign()]
    return walls,current,trace

def polar_extension(B,V,G,Q,W,H,rays):
    # This exact finite alphabet is exhaustive because all squared cosines lie in Q(sqrt5).
    # cos(2pi/m) has degree phi(m)/2<=2; m8/12 have the wrong quadratic fields.
    trials=[];passed=[]
    for ri,row in enumerate(rays):
        if row['kind']!='spacelike':continue
        v=row['point'];q=row['norm'];relations=[];bad=[]
        for j,n in enumerate(W):
            p=pair(v,n,Q)
            try:r=classify(p,q,pair(n,n,Q))
            except AssertionError:r='forbidden';bad.append(j)
            relations.append(r)
        trials.append(dict(ray=ri,relations=relations,bad=bad))
        if not bad:passed.append(ri)
    mutual=[]
    for i,j in combinations(passed,2):
        a,b=rays[i],rays[j]
        lab=classify(pair(a['point'],b['point'],Q),a['norm'],b['norm'])
        assert lab=='dotted';mutual.append([i,j,lab])
    caps=[rays[i]['point']for i in passed]
    assert all(pair(H,p,Q).sign()<0 for p in caps)
    NN,rr,trace=clip_cone(Q,W,rays,caps)
    counts=Counter(r['kind']for r in rr)
    outside=next(r['point']for r in rr if r['kind']=='spacelike')
    e=K(F(1,1024));E=add(outside,scale(e,H))
    while pair(E,E,Q).sign()<=0:e/=2;E=add(outside,scale(e,H))
    assert pair(E,H,Q).sign()<0 and all(pair(E,n,Q).sign()<0 for n in NN)
    # Actual facets follow explicitly by orthogonal projection of the common center.
    facets=[]
    for i,n in enumerate(NN):
        x=add(H,scale(-pair(H,n,Q)/pair(n,n,Q),n))
        assert pair(x,x,Q).sign()<0
        assert all(pair(x,m,Q).sign()==(0 if i==j else -1)for j,m in enumerate(NN))
        facets.append(x)
    # Single permitted source double; every surviving or reflected facet has an odd ridge.
    D=[v for i,v in enumerate(V)if i!=8]+[add(V[i],scale(-2*G[i][8],V[8]))for i in [7,9,10]]
    dq=[pair(v,v,B)for v in D];DG=[[pair(a,b,B)for b in D]for a in D]
    dl={(i,j):classify(DG[i][j],dq[i],dq[j])for i,j in combinations(range(13),2)}
    odd={str(i):next([a,b,int(m)]for(a,b),m in dl.items()if i in(a,b)and m in ['3','5'])for i in range(13)}
    return dict(polar_trials=trials,accepted_rays=passed,mutual_relations=mutual,normals=matrix_enc(NN),
        clipping_trace=trace,complete_ray_census=[dict(point=enc(r['point']),incident=r['incident'],norm=r['norm'].enc(),kind=r['kind'])for r in rr],
        counts=dict(counts),actual_facet_points=matrix_enc(facets),strict_escape=enc(E),escape_norm=pair(E,E,Q).enc(),escape_pairings=enc([pair(E,n,Q)for n in NN]),
        source_double=dict(normals=matrix_enc(D),gram=matrix_enc(DG),odd_ridge_witnesses=odd,conclusion='The unique permitted first source double has no doubling wall and no B2-special facet.'))

def determinant(A):
    M=[list(r)for r in A];d=O
    for j in range(len(M)):
        k=next(i for i in range(j,len(M))if M[i][j])
        if k!=j:M[j],M[k]=M[k],M[j];d=-d
        pivot=M[j][j];d*=pivot
        for i in range(j+1,len(M)):
            t=M[i][j]/pivot
            M[i]=[a-t*b for a,b in zip(M[i],M[j])]
    return d

def arithmetic_host(Q,W,extension):
    def dec(M):return [[K(F(x[0]),F(x[1]))for x in row]for row in M]
    N=dec(extension['normals'])
    basis=ident(10);basis[0]=W[8]
    A=[[basis[j][i]for j in range(10)]for i in range(10)]
    Ai=inv(A);G=[[pair(a,b,Q)for b in basis]for a in basis]
    def integral(x):return (x.a-x.b).denominator==1 and (2*x.b).denominator==1
    assert determinant(A)==PHI/2
    assert determinant(Q)==-4/PHI and determinant(G)==-PHI
    assert all(integral(x)for row in G for x in row)
    assert inertia([[K(x.a,-x.b)for x in row]for row in G])==[10,0]
    assert G[9][9]==O # odd, not an even unimodular lattice
    refs=[]
    for i,n in enumerate(N):
        c=mv(Ai,n);row=mv(G,c);q=pair(n,n,Q)
        R=[[K(int(a==b))-2*c[a]*row[b]/q for b in range(10)]for a in range(10)]
        assert all(integral(x)for r in R for x in r)
        refs.append(dict(wall=i,normal_coordinates=enc(c),norm=q.enc(),reflection=matrix_enc(R)))
    return dict(basis=matrix_enc(basis),basis_determinant=(PHI/2).enc(),gram=matrix_enc(G),determinant=(-PHI).enc(),conjugate_inertia=[10,0],odd_unit_vector=9,integral_reflections=refs,
      primary_theorem_url='https://www.mathnet.ru/php/getFT.phtml?jrnid=vmumm&option_lang=eng&paperid=3408&what=fullt',
      primary_theorem_location='Bugaenko1984, printedp6 (definition, uniqueness of odd unimodular admissible form, stated dimension split); pp10-12 (infinite chamber stabilizer for n>=8).',
      conclusion='The displayed19 reflections preserve an odd unimodular admissible Z[phi] lattice in signature(9,1). Bugaenko1984 identifies this lattice type and proves its full integral reflection group has infinite covolume. Hence no finite-volume completion using only reflections preserving this lattice exists.',
      limitation='This is not a no-go for arbitrary Q(sqrt5) normals, another commensurable lattice, or field-extension normals; such reflections need not preserve the specified lattice.')

def run():
    B,V,G,vs,Q,W,GG,labels,H,rays,E=build()
    # Every source facet except8 has a literal odd nonright edge. Such facets cannot be special.
    odd={str(i):next([a,b,m]for a,b,m in EDGES if i in(a,b)and m%2)for i in range(11)if i!=8}
    assert not any(8 in(a,b)and m!=4 for a,b,m in EDGES)
    extension=polar_extension(B,V,G,Q,W,H,rays)
    sources={f'scripts/{TAG}.py', 'scripts/fresh_h9_makarov_face_extension_20261013.py',f'references/{TAG}_h7_diagram.pdf'}
    return dict(status='PASS_COUPLED_B2_H7_H9_GEOMETRY_AND_UNIMODULAR_HOST_GATE',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(sources)},
      field='Q(sqrt5); coefficient pairs [a,b] mean a+b sqrt5',
      source_urls=['https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/7-11.pdf','https://arxiv.org/pdf/math/0510238','https://arxiv.org/pdf/0903.0138'],
      source_numbering={'0':'left outer','1':'left junction','2':'upper left','3':'upper center','4':'upper right','5':'right junction','6':'right outer','7':'upper hanging','8':'lower center','9':'lower left','10':'lower right'},
      source=dict(form=matrix_enc(B),normals=matrix_enc(V),gram=matrix_enc(G),finite_nonright_edges=EDGES,dotted_cosh_squares=[['7,8',K(F(3,4),F(1,4)).enc()],['9,10',(PHI*PHI).enc()]],vertices=[dict(point=enc(v['point']),incident=v['incident'])for v in vs],compact=True,positive_timelike_dual=enc(smv(V))),
      B2_profile_gate=dict(inverse_normalized=[['1','sqrt2/2'],['sqrt2/2','1']],face_profiles=['right/right','pi3/right','right/pi3'],special_finite_pair_rules={'special/ordinary':'source right or pi4; becomes right or pi3','same_color_special/special':'no finite meeting','opposite_color_special/special':'only source pi4; becomes right'},forced_ordinary_odd_witnesses=odd,special_source_facets=[8]),
      ambient=dict(form=matrix_enc(Q),normals=matrix_enc(W),gram=matrix_enc(GG),pair_counts=dict(Counter(labels.values())),center=enc(H),center_norm=pair(H,H,Q).enc(),rays=[dict(source_vertex=r['source_vertex'],fiber=r['fiber'],point=enc(r['point']),incident=r['incident'],norm=r['norm'].enc(),kind=r['kind'])for r in rays],ray_counts=dict(Counter(r['kind']for r in rays)),strict_escape=enc(E),escape_norm=pair(E,E,Q).enc(),escape_pairings=enc([pair(E,w,Q)for w in W])),
      extension=extension,arithmetic_host=arithmetic_host(Q,W,extension),
      full_ray_coverage='Projection to the compact source cone; transverse coordinates t0,t1>=0, t0+t1<=-v8.x. This simplex fiber decomposes over source extreme rays, so its extreme rays are exactly the three listed lifts, collapsed when v8.x=0.',
      scope='Exact coupled13-wall H9 chamber, all48 polar tests, coherent accepted-polar extension and complete projective cone clipping, with strict spacelike escapes. No compact/finite-volume H9 example; outer caps or changed sources remain open.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status']);print('sourcevertices',len(d['source']['vertices']),'ambient',d['ambient']['ray_counts'],'pairs',d['ambient']['pair_counts'])
