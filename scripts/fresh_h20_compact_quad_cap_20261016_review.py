#!/usr/bin/env python3
"""Separate new-cap and direct Lorentz-three face reconstruction.
Shares the frozen independent157-wall physical/interval engine, never imports
the new producer. New kernel seed, inverse, face coordinates and cone tests.
"""
import argparse,ast,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_compact_quad_cap_20261016_review'
PREVIOUS='scripts/fresh_h20_release144_cap_20261016_review.py'
PRIMARY='results/fresh_h20_compact_quad_cap_20261016.json'


def predecessor():
    path=ROOT/PREVIOUS;tree=ast.parse(path.read_text())
    run=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='run')
    stop=max(i for i,n in enumerate(run.body)if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='paths'for t in n.targets))
    run.body=run.body[:stop]+ast.parse('def f():\n return locals()').body[0].body
    ns={'__name__':'_reviewed157_geometry','__file__':str(path)}
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),ns)
    return ns,ns['run']()


def run():
    primary=json.loads((ROOT/PRIMARY).read_text());assert primary['status']=='PASS_EXACT_158_WALL_COMPACT_QUADRILATERAL_CAP'
    ns,L=predecessor()
    pt,add,sub,neg,mul,div,sqrt=[ns[k]for k in('pt','plus','minus','neg','times','div','root')]
    R,saved,vec=ns['R'],ns['saved'],ns['vec']
    ipair,lin,vadd,vscale=[L[k]for k in('ipair','ilinear','ivadd','ivscale')]
    N,qs,H,Hres,S16=[L[k]for k in('full157','norms157','H','Hres','S16')]
    def total(xs):
        z=pt(0)
        for x in xs:z=add(z,x)
        return z
    def dot(a,b):return total(mul(x,y)for x,y in zip(a,b))
    def mv(M,v):return[dot(r,v)for r in M]
    def inverse(M):
        n=len(M);A=[r[:]+[pt(i==j)for j in range(n)]for i,r in enumerate(M)];piv=[]
        for j in range(n):
            k=next(i for i in range(j,n)if A[i][j][0]>0 or A[i][j][1]<0)
            A[j],A[k]=A[k],A[j];d=A[j][j];piv.append(saved(d));A[j]=[div(x,d)for x in A[j]]
            for i in range(n):
                if i==j:continue
                d=A[i][j];A[i]=[sub(x,mul(d,y))for x,y in zip(A[i],A[j])]
        assert all(R.overlap(A[i][j],pt(i==j))for i in range(n)for j in range(n))
        return[r[n:]for r in A],piv
    def determinant3(M):
        return total(mul(pt(s),mul(M[0][i],sub(mul(M[1][j],M[2][k]),mul(M[1][k],M[2][j]))))for s,i,j,k in[(1,0,1,2),(-1,1,0,2),(1,2,0,1)])
    ids=[156,151,145,95];a=ipair(N[95],N[145]);assert a[1]<-2
    G=[[pt(2)if i==j else a if{i,j}=={95,145}else pt(-1)if{i,j}=={145,156}else pt(0)for j in ids]for i in ids]
    for i,row in enumerate(G):
        for j,g in enumerate(row):assert R.overlap(ipair(N[ids[i]],N[ids[j]]),g)
    inv,piv=inverse(G);target=[pt(0),pt(0),pt(-1),neg(sqrt(pt(3)))]
    p=lin([N[i]for i in ids],mv(inv,target))
    # A different, spacelike old-wall seed. It is exactly perpendicular toS16.
    seed=N[142];assert all(R.overlap(ipair(seed,N[i]),pt(0))for i in S16)
    k=vadd(seed,vscale(lin([N[i]for i in ids],mv(inv,[ipair(seed,N[i])for i in ids])),pt(-1)))
    kq=ipair(k,k);pq=dot(target,mv(inv,target));assert kq[0]>0
    rad=div(sub(pt(2),pq),kq);assert rad[0]>0
    rights=S16|{151,156};branches=[];accepted=[]
    for sign in(-1,1):
        r=vadd(p,vscale(k,mul(pt(sign),sqrt(rad))));pairs=[ipair(r,n)for n in N]
        bad=[j for j,g in enumerate(pairs)if j not in rights|{95,145}and not(g[1]<0 and g[1]**2>2*qs[j])]
        hc=ipair(r,H);branches.append(dict(sign=sign,failed_dotted=bad,center=saved(hc)))
        if not bad and hc[1]<0:accepted.append((sign,r,pairs))
    assert len(accepted)==1;sign,r,pairs=accepted[0]
    assert all(R.overlap(x,tuple(map(F,y)))for x,y in zip(r,primary['definition']['physical_normal_bounds']))
    assert R.overlap(ipair(r,r),pt(2))
    contacts=[]
    for j,g in enumerate(pairs):
        if j in rights:relation='right';assert R.overlap(g,pt(0))
        elif j==95:relation='pi/6';assert R.overlap(g,neg(sqrt(pt(3))))
        elif j==145:relation='pi/3';assert R.overlap(g,pt(-1))
        else:relation='dotted'
        assert relation==primary['all157_contacts'][j]['relation']
        contacts.append(dict(wall=j,relation=relation,pair=saved(g)))
    assert Counter(c['relation']for c in contacts)=={'right':18,'pi/6':1,'pi/3':1,'dotted':137}
    full=N+[r];norms=qs+[F(2)];hp=[ipair(H,n)for n in full];assert all(x[1]<0 for x in hp)
    facetq=[sub(ipair(H,H),div(mul(h,h),pt(q)))for h,q in zip(hp,norms)];assert all(q[1]<0 for q in facetq)
    # Direct projection to the three-dimensional complement ofS16+151+156.
    def project(v):
        return vadd(v,vscale(lin([N[151],N[156]],[div(ipair(v,N[j]),pt(2))for j in(151,156)]),pt(-1)))
    hf=project(Hres);sides=[95,142,145,157];P=[project(full[i])for i in sides]
    assert ipair(hf,hf)[1]<0
    for v in P+[hf]:
        assert all(R.overlap(ipair(v,full[i]),pt(0))for i in S16|{151,156})
    pg=[[ipair(v,w)for w in P]for v in P]
    expected=[(0,1,pt(0)),(1,2,neg(div(sqrt(pt(qs[142]/2)),pt(2)))),(2,3,pt(-1)),(3,0,neg(sqrt(pt(3))))]
    for i,j,v in expected:assert R.overlap(pg[i][j],v)
    for i,v in enumerate([pt(2),pt(3*qs[142]/4),pt(F(3,2)),pt(2)]):assert R.overlap(pg[i][i],v)
    F3=P[:3];Q3=[r[:3]for r in pg[:3]];det3=determinant3(Q3);assert det3[1]<0
    coords=[];Qi,_=inverse(Q3)
    for v in P:coords.append(mv(Qi,[ipair(v,b)for b in F3]))
    for v,c in zip(P,coords):assert all(R.overlap(x,y)for x,y in zip(v,lin(F3,c)))
    corners=[]
    for i,j in[(0,1),(1,2),(2,3),(3,0)]:
        M=[[pg[i][i],pg[i][j]],[pg[j][i],pg[j][j]]];dd=sub(mul(M[0][0],M[1][1]),mul(M[0][1],M[1][0]))
        assert M[0][0][0]>0 and dd[0]>0;J,_=inverse(M)
        v=vadd(hf,vscale(lin([P[i],P[j]],mv(J,[ipair(hf,P[i]),ipair(hf,P[j])])),pt(-1)))
        vq,vh=ipair(v,v),ipair(v,H);assert vq[1]<0 and vh[1]<0
        incidence=sorted(S16|{151,156,sides[i],sides[j]});assert len(incidence)==20
        products=[ipair(v,n)for n in full]
        for q,g in enumerate(products):
            if q in incidence:assert R.overlap(g,pt(0))
            else:assert g[1]<0,(i,j,q,saved(g))
        corners.append(dict(sides=[sides[i],sides[j]],active=incidence,point=vec(v),norm=saved(vq),all158_pairs=vec(products)))
    # Enumerate all six intersections directly in the Lorentz3 basis.
    cone=[];feasible=[]
    for i,j in combinations(range(4),2):
        u=[ipair(P[i],b)for b in F3];v=[ipair(P[j],b)for b in F3]
        c=[sub(mul(u[1],v[2]),mul(u[2],v[1])),sub(mul(u[2],v[0]),mul(u[0],v[2])),sub(mul(u[0],v[1]),mul(u[1],v[0]))]
        z=lin(F3,c);results=[]
        for orientation in(-1,1):
            w=vscale(z,pt(orientation));other=[ipair(w,P[t])for t in range(4)if t not in(i,j)]
            assert all(q[0]>0 or q[1]<0 for q in other)
            valid=all(q[1]<0 for q in other)
            if valid:
                assert ipair(w,w)[1]<0 and ipair(w,H)[1]<0;feasible.append([sides[i],sides[j]])
            results.append(dict(orientation=orientation,other_side_products=vec(other),feasible=valid))
        cone.append(dict(sides=[sides[i],sides[j]],tests=results))
    assert sorted(map(sorted,feasible))==sorted(map(sorted,[x['sides']for x in corners]))and len(feasible)==4
    cusps=[ipair(r,v)for v in L['cusp_vectors']+[L['nu2']]];assert len(cusps)==5 and all(g[1]<0 for g in cusps)
    ep=[ipair(L['E7'],n)for n in full];assert all(g[1]<0 for g in ep)and L['q7']>0
    triangles=[ipair([tuple(map(F,x))for x in rec['point']],r)for rec in L['alltrianglevertices']];assert all(g[1]<0 for g in triangles)
    neighbors=sorted(rights|{95,145});assert len(neighbors)==20
    rejected=[dict(omitted=j,dotted_witness=[95,145])for j in neighbors if j not in(95,145)];assert len(rejected)==18
    paths=[f'scripts/{TAG}.py',PREVIOUS,'results/fresh_h20_release144_cap_20261016_review.json',
           'scripts/fresh_h20_compact_quad_cap_20261016.py',PRIMARY]
    return dict(status='PASS_INDEPENDENT_158_WALL_CAP_AND_COMPLETE_COMPACT_QUADRILATERAL',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                shared_arithmetic='Frozen independent157-wall physical reconstruction and outward rational interval engine; no new producer import.',
                independent_definition=dict(frame=ids,kernel_seed=142,inverse_pivots=piv,particular_square=saved(pq),kernel_square=saved(kq),radical_square=saved(rad),branches=branches,selected_sign=sign,normal=vec(r)),
                all157_contacts=contacts,actual_facets=158,facet_squares=vec(facetq),
                compact_quad=dict(support=sorted(S16|{151,156}),sides=sides,projected_gram=[vec(row)for row in pg],
                                  basis=sides[:3],basis_determinant=saved(det3),vertices=corners,complete_six_pair_cone=cone,
                                  adjacent_angle_cosine_squares=['0','1/9','1/3','3/4']),
                new_finite_vertices=2,rejected_finite_supports=rejected,five_cusp_pairs=vec(cusps),old_triangle_pairs=vec(triangles),
                strict_escape=dict(norm=str(L['q7']),all158_pairs=vec(ep)),
                scope='Independent new-root reconstruction, all158 actual facets, exhaustive new-facet finite supports and full projectedH2 cone. Old157 geometry/parameter isolation is inherited and replayed through its independent assertions. The quadrilateral is not a Coxeter polygon. StrictE7 remains; no finite-volume completion is claimed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved independent H20 quad replay')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
