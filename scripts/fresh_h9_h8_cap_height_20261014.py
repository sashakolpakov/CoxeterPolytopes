#!/usr/bin/env python3
"""Maximal exceptional-cap heights and their uniform infinite-volume gate."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_cap_height_20261014'
BASE='results/fresh_h9_h8_disjoint_cap_20261014.json'
REVIEW='results/fresh_h9_h8_disjoint_cap_20261014_review.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d.get('sources',{}).items():assert sha(f)==h,f
    return d
def run():
    base,review,src=[read(p)for p in(BASE,REVIEW,SOURCE)]
    assert review['status'].startswith('PASS_')and review['sources'][BASE]==sha(BASE)
    a=A.E(A.Q(F(-1,4),F(1,4)));D=a+2;phi=A.PHI
    Q=A.dm(base['form']);Q[-1][-1]=D;N=A.dm(base['normals']);H=A.dv(base['center'])
    assert A.inertia(Q)==[9,1]and len(N)==15
    rows=[A.mv(Q,n)for n in N];G=[[A.dot(n,row)for row in rows]for n in N];norms=[G[i][i]for i in range(15)]
    assert norms[:11]==[A.E(2)]*11 and norms[11]==D and norms[12:]==[A.E(2)]*3
    assert A.pair(H,H,Q).sign()<0 and all(A.dot(H,row).sign()<0 for row in rows)
    pairs=[];counts=Counter()
    for i,j in combinations(range(15),2):
        kind,m=A.relation(G[i][j],norms[i],norms[j]);assert kind in('right','finite','dotted')
        key='pi'+str(m)if kind=='finite'else kind;counts[key]+=1
        pairs.append(dict(pair=[i,j],kind=kind,label=m,product=G[i][j]))
    assert dict(counts)==dict(right=79,pi3=12,pi5=5,dotted=9)
    facets=[]
    for i,n in enumerate(N):
        p=[x-A.dot(H,rows[i])/norms[i]*y for x,y in zip(H,n)]
        vals=[A.dot(p,row)for row in rows];q=A.pair(p,p,Q)
        assert q.sign()<0 and A.pair(p,H,Q).sign()<0 and not vals[i]
        assert all(v.sign()<0 for j,v in enumerate(vals)if j!=i)
        facets.append(dict(facet=i,point=p,norm=q,pairings=vals))
    source_vertices=[]
    for vi,row in enumerate(base['complete_finite_vertices'][:41]):
        assert row['source_vertex']==vi
        p=A.dv(row['point']);assert not p[-1]
        ids=src['geometry']['vertices'][vi];assert len(ids)==8
        vals=[A.dot(p,r)for r in rows]
        assert A.pair(p,p,Q).sign()<0 and [i for i,x in enumerate(vals)if not x]==ids+[11]
        assert all(x.sign()<=0 for x in vals)
        source_vertices.append(p)
    vertices=[];candidates=[]
    for vi,p in enumerate(source_vertices):
        vertices.append(dict(source_vertex=vi,cap=11,point=p))
    for cap in(12,13,14):
        for vi,v in enumerate(source_vertices):
            p=v[:];p[-1]=A.dot(v,rows[cap])/D
            assert not A.dot(p,rows[cap])
            vals=[A.dot(p,row)for row in rows];q=A.pair(p,p,Q)
            feasible=all(x.sign()<=0 for x in vals)
            assert q
            if q.sign()<0:
                assert feasible and A.pair(p,H,Q).sign()<0
                vertices.append(dict(source_vertex=vi,cap=cap,point=p))
            candidates.append(dict(cap=cap,source_vertex=vi,norm=q,feasible=feasible,
                kind='finite'if q.sign()<0 else'spacelike'))
    for row in vertices:
        p=row['point'];vals=[A.dot(p,r)for r in rows];ids=[i for i,x in enumerate(vals)if not x]
        expected=src['geometry']['vertices'][row['source_vertex']]+[row['cap']]
        assert ids==expected and len(ids)==9 and all(x.sign()<=0 for x in vals)
        assert A.inertia([[G[i][j]for j in ids]for i in ids])==[9,0]
        row.update(incident=ids,norm=A.pair(p,p,Q),pairings=vals)
    assert len(vertices)==62 and len({A.normkey(r['point'])for r in vertices})==62
    assert Counter(r['cap']for r in vertices)=={11:41,12:11,13:5,14:5}
    edges=[]
    for i,j in combinations(range(len(vertices)),2):
        support=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(support)==8:
            assert A.inertia([[G[x][y]for y in support]for x in support])==[8,0]
            edges.append(dict(vertices=[i,j],support=support))
    # Every finite vertex contains R or a cap: nine source-only normals
    # cannot have positive Gram in the source Lorentz9-space. No two of R
    # and the three caps can meet finitely. Projection of a finite vertex
    # to the source is timelike and has eight source incidences, hence is
    # one of its41vertices. This proves the finite census complete.
    exceptions=[A.dv(r['point'])for r in base['exceptional_caps']]
    E=[sum((v[i]for v in exceptions),A.Z)+H[i]/16 for i in range(9)]+[A.E(-5)]
    eq=A.pair(E,E,Q);ev=[A.dot(E,row)for row in rows]
    assert eq==A.E(A.Q(F(4363,256),F(-3057,512)))and eq.sign()>0
    assert all(x.sign()<0 for x in ev)and A.pair(E,H,Q).sign()<0
    # Exact all-label source contact alphabet, for any cap norm q>0.
    # Min contact1 gives q<=2; at q>1/2 it equals1/(2cos²(pi/m)).
    # At m3, phi/2=cos(pi5); for allm>=4, phi*cos(pi/m)>=phi/sqrt2>1.
    contacts=[]
    for row in base['exceptional_caps']:
        h=A.dv(row['contacts']);assert any(x==1 for x in h)
        assert sorted([i for i,x in enumerate(h)if x])in([4,9,10],[0,8,10],[0,8,9])
        assert all(x in(A.Z,A.ONE,phi,phi*phi*phi*phi)for x in h)
        contacts.append(dict(source_vertex=row['source_vertex'],values=h))
    assert(phi*phi-2).sign()>0 and(phi*phi*phi*phi-2).sign()>0
    distances=[]
    for i,j in combinations(range(3),2):
        b=-A.pair(exceptions[i]+[A.Z],exceptions[j]+[A.Z],Q)
        assert(b-a-4).sign()>0
        distances.append(dict(pair=[i,j],negative_source_product=b,strict_all_height_dotted_margin=b-a-4))
    # If0<qi<=2, eta_i=sqrt((a+qi)/(a+2))<=1. With E_last=-5,
    # (E,(vi,-eta_i))=(E,(vi,-1))+5D(eta_i-1)<0.
    # This isometry absorbs arbitraryd>0 and positive individual heights.
    sigmaQ=[[A.E(A.Q(x.a.a,-x.a.b))for x in row]for row in Q]
    assert all(not x.b for row in Q for x in row)and A.inertia(sigmaQ)==[10,0]
    fieldnorm=D.a.a*D.a.a-5*D.a.b*D.a.b;assert fieldnorm==F(11,4)
    paths=[BASE,REVIEW,SOURCE,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return A.enc(dict(status='PASS_MAXIMAL_HEIGHT_GEOMETRY_AND_ALL_HEIGHT_INFINITE_VOLUME',sources={p:sha(p)for p in paths},
        form=Q,vertical_coefficient=D,norms=norms,normals=N,gram=G,pairs=pairs,pair_counts=dict(counts),
        strict_center=dict(point=H,norm=A.pair(H,H,Q),pairings=[A.dot(H,r)for r in rows]),facet_points=facets,
        complete_finite_vertices=vertices,complete_finite_vertex_count=len(vertices),
        source_cap_vertex_candidates=candidates,complete_finite_edges=edges,complete_finite_edge_count=len(edges),
        open_finite_corner_germs=9*len(vertices)-2*len(edges),
        strict_escape=dict(point=E,norm=eq,pairings=ev,center_pairing=A.pair(E,H,Q)),
        general_height_gate=dict(source_contacts=contacts,cap_norm_domain='0<q<=1/2, orq=1/(2cos²(pi/m))forinteger m>=3',
            max_cap_norm=2,all_height_cap_pair_margins=distances,
            standard_form='Q8 plusD=a+2; arbitraryd>0 andhi>0 map eachcap to(vi,-eta_i), eta_i=sqrt((a+qi)/D), qi=dhi²-a.',
            escape_difference_formula='newpairing-maximalpairing=5D(eta_i-1)<=0',
            conclusion='Every Coxeter chamber retaining all11sourcewalls and the bottom, with onecap on each of the same three exceptional source directions and arbitrary positiveindividualheights, has infinitevolume. The same escape covers every subfamily.'),
        host_triage=dict(conjugate_signature=[10,0],coefficient_field_norm=str(fieldnorm),
            scope='This maximalhost is admissible at the other embedding. No reflectivity, arithmeticity or nonreflectivity classification is inferred.'),
        scope='Exact changed15wall maximalheight seed, complete62finitevertices and finiteedges, strictescape. Uniform all-real independent-height exclusion for these fixedthreeprojectedcapdirections only; newprojecteddirections, coupled sourcewalls or additional facets are outside. No finite-volume or compactH9example.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();d=run()
    out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text());print('SAVED_REPLAY_PASS')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['complete_finite_vertex_count'],d['complete_finite_edge_count'],d['open_finite_corner_germs'],flush=True)
