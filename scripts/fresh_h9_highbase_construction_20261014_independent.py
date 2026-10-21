#!/usr/bin/env python3
"""Independent stdlib convex-base and incidence audit; no producer imports.

For the 125-facet multiwedges, completeness of nine-faces is checked through
the full-block identity |face| <= 5 + |full base blocks| <= 9, rather than
attempting C(125,9) subsets. Equality reconstructs every archived vertex.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product, permutations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_highbase_construction_20261014'
FILES=[f'results/{TAG}_{s}.json' for s in ['cyclic','cell','refined','pentagon']]
OUT=ROOT/f'results/{TAG}_independent.json'


def ph(x):return tuple(map(F,x))
def plus(a,b):return(a[0]+b[0],a[1]+b[1])
def times(a,b):return(a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def inner(v,w):
    a=b=F(0)
    for (x,y),(z,t) in zip(v,w):a+=x*z+y*t;b+=x*t+y*z+y*t
    return(a,b)
def sign(v):
    a,b=2*v[0]+v[1],v[1]
    if b==0:return (a>0)-(a<0)
    if a==0:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mask(a):return sum(1<<i for i in a)


def run():
    cyclic,cell,refined,pentagon=map(read,FILES)
    for data in [cyclic,cell,refined,pentagon]:
        for p,h in data['sources'].items():assert sha(p)==h,p
    cyrows=[]
    for row in cyclic['cases']:
        b=row['base'];n=b['vertex_count'];d=b['dimension']
        facets=set(map(tuple,b['dual_facets']))
        # Gale evenness is independent of the producer's sign-polynomial
        # enumeration. Supporting polynomials are also checked directly.
        expected={f for f in combinations(range(n),d)
                  if all(sum(a<i<c for i in f)%2==0
                         for a,c in combinations([i for i in range(n) if i not in f],2))}
        assert facets==expected
        assert b['moment_curve_vertices']==[[t**k for k in range(1,d+1)] for t in range(n)]
        for f,p in zip(b['dual_facets'],b['exact_support_polynomials']):
            val=[sum(c*t**k for k,c in enumerate(p)) for t in range(n)]
            assert [i for i,x in enumerate(val) if x==0]==f
            assert all(x>0 for i,x in enumerate(val) if i not in f)
        faces={()}
        for f in facets:
            for k in range(1,d+1):faces.update(combinations(f,k))
        missing=[a for k in range(2,d+2) for a in combinations(range(n),k)
                 if a not in faces and all(a[:j]+a[j+1:] in faces for j in range(k))]
        assert missing==list(map(tuple,b['minimal_nonfaces']))
        total=row['facets_after_wedges'];hist=Counter();valid=[];count=0
        for cuts in combinations(range(1,total),n-1):
            bounds=(0,)+cuts+(total,);w=tuple(bounds[i+1]-bounds[i] for i in range(n))
            maximum=max(sum(w[i] for i in a) for a in missing);hist[maximum]+=1;count+=1
            if maximum<=5:valid.append(w)
        assert not valid and not row['all_admissible_weights']
        assert count==row['positive_weight_compositions']
        assert {str(k):v for k,v in hist.items()}==row['maximum_missing_order_histogram']
        cyrows.append(dict(base_dimension=d,vertices=n,facets=len(facets),all_weights=count,admissible=0))

    b=cell['base'];V=[tuple(map(ph,v)) for v in b['dual_vertices']];n=len(V);assert n==120
    G=[[inner(v,w) for w in V] for v in V]
    assert all(G[i][i]==(4,0) for i in range(n))
    assert all(tuple((-a,-c) for a,c in v) in V for v in V)
    adjacency=[{j for j in range(n) if i!=j and G[i][j]==(0,2)} for i in range(n)]
    missing={(i,j) for i,j in combinations(range(n),2) if j not in adjacency[i]}
    assert missing==set(map(tuple,b['minimal_nonfaces']))
    assert len(missing)==6420 and set(map(len,adjacency))=={12}
    # Different recursive clique enumeration from the producer.
    found=[];counts=Counter()
    def expand(prefix,available):
        counts[len(prefix)]+=1
        if not available:
            found.append(tuple(prefix));return
        while available:
            i=min(available);available.remove(i)
            expand(prefix+[i],available & adjacency[i])
    expand([],set(range(n)))
    tetra={f for f in found if len(f)==4}
    assert counts=={0:1,1:120,2:720,3:1200,4:600}
    assert tetra==set(map(tuple,b['dual_facets']))
    level=ph(b['support_level']);assert level==(4,6)
    for f,normal in zip(b['dual_facets'],b['exact_facet_normals']):
        u=tuple(map(ph,normal))
        assert u==tuple(tuple(sum(V[i][k][a] for i in f) for a in range(2)) for k in range(4))
        vals=[inner(u,v) for v in V]
        assert all(vals[i]==level for i in f)
        assert all(sign((level[0]-vals[i][0],level[1]-vals[i][1]))>0 for i in range(n) if i not in f)
    assert sign((4,-2))>0 and sign((4,6))>0 # tetra position Gram eigenvalues
    ridges=Counter(t for f in tetra for t in combinations(f,3));assert set(ridges.values())=={2}
    neighbors={f:set() for f in tetra};rd={t:[] for t in ridges}
    for f in tetra:
        for t in combinations(f,3):rd[t].append(f)
    for aa in rd.values():
        a,c=aa;neighbors[a].add(c);neighbors[c].add(a)
    reached={min(tetra)};todo=list(reached)
    for f in todo:
        for g in neighbors[f]-reached:reached.add(g);todo.append(g)
    assert reached==tetra
    group=list(map(tuple,b['tetrahedron_stabilizer']));f=tuple(b['chosen_tetrahedron'])
    assert len(set(group))==24
    assert {tuple(p[i] for i in f) for p in group}==set(permutations(f))
    for p in group:
        assert set(p)==set(range(n))
        assert all(G[p[i]][p[j]]==G[i][j] for i in range(n) for j in range(n))
    GS=set(group);assert all(tuple(p[q[i]] for i in range(n)) in GS for p in group for q in group)
    partition=[];left=set(range(n))-set(f)
    while left:
        orb={p[min(left)] for p in group};assert orb<=left
        partition.append(tuple(sorted(orb)));left-=orb
    assert partition==list(map(tuple,b['fifth_vertex_orbits'])) and len(partition)==14
    assert len(cell['cases'])==14 and cell['complete_selected_recipe_family']
    rows=[]
    a,c=pentagon['edge']
    assert c in adjacency[a]
    assert set(pentagon['doubled_vertices'])==adjacency[a]&adjacency[c]
    assert len(adjacency[a]&adjacency[c])==5
    for r in [*cell['cases'],dict(pentagon,case=14)]:
        chosen=set(r['doubled_vertices']);assert len(chosen)==5
        if r['case']<14:
            assert set(f)<=chosen
        else:
            assert all(len(adjacency[i]&chosen)==2 for i in chosen)
            assert not any(set(t)<=chosen for t in tetra)
        weights=r['weights'];assert weights==[2 if i in chosen else 1 for i in range(n)]
        blocks=list(map(tuple,r['blocks']));assert [i for b0 in blocks for i in b0]==list(range(125))
        expected_M={tuple(sorted(blocks[i]+blocks[j])) for i,j in missing}
        assert set(map(tuple,r['minimal_nonfaces']))==expected_M and max(map(len,expected_M))==4
        expected_V=set()
        for basefacet in tetra:
            core=[i for j in basefacet for i in blocks[j]]
            free=[blocks[j] for j in sorted(chosen-set(basefacet))]
            for representatives in product(*free):
                expected_V.add(tuple(sorted(core+list(representatives))))
        actual=set(map(tuple,r['actual_vertices']))
        assert actual==expected_V and len(actual)==len(r['actual_vertices'])
        assert all(len(a)==9 for a in actual)
        # Any missing-face-free set has full blocks forming a base clique.
        # Its size is <= sum(w-1)+|full blocks| <=5+4. Equality forces a
        # tetrahedron and exactly one omitted element in every other block.
        assert sum(w-1 for w in weights)==5 and max(counts)==4
        assert set(Counter(t for a in actual for t in combinations(a,8)).values())=={2}
        adj=[set() for _ in range(125)]
        for a in expected_M:
            for i in a:adj[i].update(set(a)-{i})
        reached={0};todo=[0]
        for i in todo:
            for j in adj[i]-reached:reached.add(j);todo.append(j)
        assert len(reached)==125
        proof=refined if r['case']==0 else r
        assert proof['verdict']=='EXACT_LANNER_VERTEX_OBSTRUCTION'
        assert not proof['state_limit_reached'] and not any(t['truncated'] for t in proof['stages'])
        assert proof['stages'][-1]['passing_states']==0
        rows.append(dict(case=r['case'],facets=125,dimension=9,vertices=len(actual),
                         missing_face_counts=dict(Counter(map(len,expected_M))),
                         complete_incidence=True,nonproduct=True,
                         label_certificate='refined' if r['case']==0 else ('pentagon' if r['case']==14 else 'primary')))
    files=FILES+[f'scripts/{TAG}_independent.py']
    return json.loads(json.dumps(dict(status='PASS_INDEPENDENT_HIGH_BASE_POLYTOPAL_INCIDENCE_AUDIT',
        sources={p:sha(p) for p in files},cyclic_cases=cyrows,cell_base_clique_counts=dict(counts),
        cell_cases=rows,total_cell_vertices=sum(r['vertices'] for r in rows),
        scope='Independent exact supporting facets, clique completeness, stabilizer/orbit structure, all multiwedge vertex incidences and ridge counts. '
              'Label-join truth is replayed separately; this checker does not implement a second Lanner catalogue or join.')))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['total_cell_vertices'])
