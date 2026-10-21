#!/usr/bin/env python3
"""Exact finite supports for a written universal equality proof, not a census."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cmp_to_key
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_ideal_ra4_equality_structure_20261018.json'
PINS={
 'results/ideal_ra4_second_volume_20261018.md':'9dff11c375bddcf62f16f3f01ddabd7f5b5c137756cf02ebff275076e6869c25',
 'results/fresh_ideal_ra4_double24_volume_20261018.json':'386d9650b8f669e3edcb3cdee56ed53d1c9ea4a90af2b4ba5c8c1b83f9d7bfa3',
}

def rank(rows):
    a=[[Q(x) for x in row] for row in rows]
    k=0
    if not a:return 0
    for j in range(len(a[0])):
        p=next((i for i in range(k,len(a)) if a[i][j]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k]
        d=a[k][j];a[k]=[x/d for x in a[k]]
        for i in range(k+1,len(a)):
            c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k

def solve(rows,rhs):
    a=[[Q(x) for x in row]+[Q(b)] for row,b in zip(rows,rhs)]
    n=len(a)
    for j in range(n):
        p=next((i for i in range(j,n) if a[i][j]),None)
        if p is None:return None
        a[j],a[p]=a[p],a[j];d=a[j][j];a[j]=[x/d for x in a[j]]
        for i in range(n):
            if i!=j:
                c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)

def dot(a,b,metric):return sum(x*y*m for x,y,m in zip(a,b,metric))

def canonical_cycle(c):
    return min(tuple(c[i:]+c[:i]) for c in (c,list(reversed(c))) for i in range(len(c)))

def ten_vertex_incidence():
    # Direct transcription of Kolpakov Figure7; completeness is published input.
    xy={'A':(0,0),'B':(6,0),'C':(6,6),'D':(0,6),'E':(2,2),
        'F':(4,2),'G':(5,3),'H':(3,3),'I':(2,4),'J':(1,3)}
    edges=sorted(tuple(sorted(e)) for e in ['AB','BC','CD','DA','AE','AJ','BF','BG',
              'CG','CI','DI','DJ','EF','FG','GH','HF','HE','EJ','JI','IH'])
    adj={v:[] for v in xy}
    for a,b in edges:adj[a].append(b);adj[b].append(a)
    assert set(map(len,adj.values()))=={4} and len(edges)==20
    def order(v):
        def compare(a,b):
            x,y=(xy[a][i]-xy[v][i] for i in range(2))
            X,Y=(xy[b][i]-xy[v][i] for i in range(2))
            h=int(y<0 or (y==0 and x<0));H=int(Y<0 or(Y==0 and X<0))
            if h!=H:return h-H
            z=x*Y-y*X
            assert z!=0
            return -1 if z>0 else 1
        return sorted(adj[v],key=cmp_to_key(compare))
    rotation={v:order(v) for v in xy};unseen={(a,b) for a,b in edges}|{(b,a) for a,b in edges}
    faces=[]
    while unseen:
        start=min(unseen);a,b=start;cycle=[]
        while True:
            unseen.remove((a,b));cycle.append(a)
            ns=rotation[b];c=ns[(ns.index(a)-1)%len(ns)];a,b=b,c
            if (a,b)==start:break
        faces.append(canonical_cycle(cycle))
    faces=sorted(faces);quads=[f for f in faces if len(f)==4]
    assert Counter(map(len,faces))=={3:8,4:4}
    assert set(quads)=={tuple('ABCD'),tuple('ABFE'),tuple('EHIJ'),tuple('CGHI')}
    face_edges=lambda f:{tuple(sorted((f[i],f[(i+1)%len(f)]))) for i in range(len(f))}
    matching=[(i,j) for i,j in combinations(range(4),2) if face_edges(quads[i])&face_edges(quads[j])]
    assert len(matching)==2 and sorted(i for p in matching for i in p)==list(range(4))
    return {'figure':'Kolpakov2012 Figure7','vertices':xy,'edges':edges,'faces':faces,
            'quadrilaterals':quads,'quadrilateral_edge_adjacency':matching,
            'historical_classification_replayed':False}

def double_oct_data(saved):
    d=saved['double'];walls=[[Q(x) for x in a] for a in d['walls']]
    vs=[[Q(x) for x in a] for a in d['vertices']];metric=(-1,1,1,1,1)
    facets=[{i for i,v in enumerate(vs) if dot(v,n,metric)==0} for n in walls]
    ridges=[set(a) for a in d['faces']['2']]
    ans=[]
    for row in d['facet_cusp_shapes']:
        if row['vertices']!=9:continue
        fi=row['facet'];f=facets[fi]
        ns={c['vertex'] for c in row['cusps'] if Q(c['aspect_ratio'])!=1}
        qs=[r for r in ridges if r<=f and len(r)==4];ts=[r for r in ridges if r<=f and len(r)==3]
        assert len(ns)==len(qs)==3 and len(ts)==8
        assert all(len(q&ns)==2 for q in qs) and all(len(t&ns)<=1 for t in ts)
        assert all(len(a&b)==1 and a&b<=ns for a,b in combinations(qs,2))
        widths=[]
        for vi in sorted(f):
            qn=[j for j,g in enumerate(facets) if j!=fi and f&g in qs and vi in g]
            tn=[j for j,g in enumerate(facets) if j!=fi and f&g in ts and vi in g]
            assert len(qn)==(2 if vi in ns else 1)
            if vi in ns:
                assert len(tn)==2
                assert dot(walls[qn[0]],walls[qn[1]],metric)==-1
                assert dot(walls[tn[0]],walls[tn[1]],metric)==-1
                kq=sum(walls[j][0] for j in qn);kt=sum(walls[j][0] for j in tn)
                assert kt==2*kq
                widths.append({'vertex':vi,'quad_pair_width':str(kq),'triangle_pair_width':str(kt)})
        ans.append({'facet':fi,'nonsquare_vertices':sorted(ns),'quads':[sorted(q) for q in qs],
                    'triangle_nonsquare_intersections':sorted(len(t&ns) for t in ts),'widths':widths})
    assert len(ans)==8
    return ans

def cuboctahedron():
    metric=(Q(-1),Q(1,2),Q(1,2),Q(1,2))
    vs=set()
    for ij in combinations(range(3),2):
        for signs in product((-1,1),repeat=2):
            x=[Q(0)]*3
            for i,s in zip(ij,signs):x[i]=Q(s)
            vs.add((Q(1),*x))
    walls=[]
    for i in range(3):
        for s in (-1,1):
            x=[Q(0)]*3;x[i]=2*s;walls.append((Q(1),*x))
    walls.extend((Q(1),*map(Q,s)) for s in product((-1,1),repeat=3))
    assert all(dot(v,v,metric)==0 for v in vs)
    assert all(dot(v,n,metric)<=0 for v in vs for n in walls)
    assert all(dot(a,b,metric)<=0 for a,b in combinations(walls,2))
    found=set()
    for ids in combinations(range(14),3):
        x=solve([[walls[i][j]*metric[j] for j in range(1,4)] for i in ids],[1]*3)
        if x is not None and all(dot((Q(1),*x),n,metric)<=0 for n in walls):found.add((Q(1),*x))
    assert found==vs
    vs=sorted(vs);facet_sizes=[sum(dot(v,n,metric)==0 for v in vs) for n in walls]
    assert facet_sizes==[4]*6+[3]*8
    cusps=[]
    for vi,v in enumerate(vs):
        ids=[i for i,n in enumerate(walls) if dot(v,n,metric)==0]
        assert len(ids)==4
        pairs=[(a,b) for a,b in combinations(ids,2) if dot(walls[a],walls[b],metric)<0]
        assert len(pairs)==2 and sorted(i for p in pairs for i in p)==ids
        kappas=[]
        for a,b in pairs:
            na=dot(walls[a],walls[a],metric);nb=dot(walls[b],walls[b],metric)
            assert na==nb and dot(walls[a],walls[b],metric)==-na
            sm=tuple(x+y for x,y in zip(walls[a],walls[b]))
            assert sm==tuple(sm[0]*x for x in v)
            kappas.append(sm[0]**2/na)
        assert max(kappas)/min(kappas)==2
        cusps.append({'vertex':vi,'opposite_pairs':pairs,'width_squares':list(map(str,kappas))})
    return {'scaled_metric':list(map(str,metric)),'vertices':[[str(x) for x in v] for v in vs],
      'walls':[[str(x) for x in n] for n in walls],'all_active_triples':364,
      'facet_sizes':facet_sizes,'cusps':cusps,'all_cusp_aspect_squared':'2'}

def six_cusp_graphs():
    all_edges=set(combinations(range(6),2));ans=[]
    for missing in combinations(sorted(all_edges),3):
        if sorted(v for e in missing for v in e)!=list(range(6)):continue
        edges=all_edges-set(missing)
        triangles=[t for t in combinations(range(6),3) if set(combinations(t,2))<=edges]
        assert len(triangles)==8
        signatures=[tuple(next(i for i,v in enumerate(p) if v in t) for p in missing) for t in triangles]
        assert set(signatures)==set(product((0,1),repeat=3))
        dual=[(i,j) for i,j in combinations(range(8),2) if len(set(triangles[i])&set(triangles[j]))==2]
        assert len(dual)==12
        assert all(sum(x!=y for x,y in zip(signatures[i],signatures[j]))==1 for i,j in dual)
        ans.append({'missing_matching':missing,'triangles':triangles,'cube_signatures':signatures,'dual_edges':dual})
    assert len(ans)==15
    return ans

def cube_algebra():
    ss=list(product((-1,1),repeat=3));equations=[]
    for i in range(3):
        other=[j for j in range(3) if j!=i]
        for sign in (-1,1):
            equations.append([s[other[0]]*s[other[1]] if s[i]==sign else 0 for s in ss])
    modes=[[1]*8]+[[s[i] for s in ss] for i in range(3)]
    assert rank(equations)==4 and rank(modes)==4
    assert all(sum(x*y for x,y in zip(e,m))==0 for e in equations for m in modes)
    gram=[[Q(-1+sum(x*y for x,y in zip(s,t)),2) for t in ss] for s in ss]
    assert rank(gram)==4 and all(gram[i][i]==1 for i in range(8))
    assert all(gram[i][j]==1-sum(x!=y for x,y in zip(ss[i],ss[j])) for i,j in combinations(range(8),2))
    return {'signs':ss,'square_equations':equations,'equation_rank':4,'affine_modes':modes,
            'normal_gram':[[str(x) for x in r] for r in gram],
            'orthogonal_affine_mode_squares':['-1/2','1/2','1/2','1/2']}

def build():
    for p,h in PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    saved=json.loads((ROOT/'results/fresh_ideal_ra4_double24_volume_20261018.json').read_text())
    partitions=[]
    for n in (7,8):
        def visit(prefix,lower,left):
            if len(prefix)==n:
                if left==0:partitions.append(prefix)
                return
            for x in range(lower,left//(n-len(prefix))+1):visit(prefix+[x],x,left-x)
        visit([],3,24)
    assert partitions==[[3,3,3,3,3,3,6],[3,3,3,3,3,4,5],[3,3,3,3,4,4,4],[3]*8]
    return {'status':'EXACT_FINITE_SUPPORTS_FOR_WRITTEN_EQUALITY_PROOF',
      'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'sources':PINS,
      'ten_vertex':ten_vertex_incidence(),'double_oct_local':double_oct_data(saved),
      'remaining_defect_partitions':partitions,'cuboctahedron':cuboctahedron(),
      'six_cusp_graphs':six_cusp_graphs(),'cube_algebra':cube_algebra(),
      'scope':{'historical_classification_replayed':False,'four_polytope_census':False,
               'universal_geometry_proved_by_this_computation':False}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=build()
    if a.verify:assert json.loads(OUT.read_text())==json.loads(json.dumps(d))
    else:
        assert not OUT.exists(), 'Refusing to overwrite an existing certificate; use --verify.'
        OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':'PASS','verify':a.verify,'ten_vertex_quads':4,'double_oct_facets':8,
      'cuboctahedron_cusps':12,'labeled_six_cusp_graphs':15,'cube_equation_rank':4}))

if __name__=='__main__':main()
