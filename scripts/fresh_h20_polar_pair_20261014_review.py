#!/usr/bin/env python3
"""Independent polar-pair geometry and complete new-vertex replay.

Uses the earlier independent sparse field, never the primary polar arithmetic.
The census is reconstructed by fixed-size clique recursion and exact tree
Schur elimination, independently of maximal cliques and dense leading pivots.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
import fresh_h20_bridge_faces_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_polar_pair_20261014_review'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
BRIDGE='results/fresh_h20_e7_bridge_face_20261014_exact.json'
PRE='results/fresh_h20_bridge_edge_closure_20261014.json'
SUM='results/fresh_h20_bridge_faces_20261014_review.json'
PRIMARY='results/fresh_h20_polar_pair_20261014.json'
CENSUS='results/fresh_h20_polar_pair_20261014_census.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rank(matrix):
    a=[list(map(F,r))for r in matrix];k=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(k,len(a))if a[i][j]),None)
        if pivot is None:continue
        a[k],a[pivot]=a[pivot],a[k];q=a[k][j];a[k]=[x/q for x in a[k]]
        for i in range(k+1,len(a)):
            q=a[i][j]
            if q:a[i]=[x-q*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k

def run():
    old,bridge,pre,summed,d,cn=[json.loads((ROOT/p).read_text())for p in[BASE,BRIDGE,PRE,SUM,PRIMARY,CENSUS]]
    for doc in[old,bridge,pre,summed,d,cn]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    ff=Field([parse(old['field']['discriminant'],0),
              parse(bridge['field']['new_positive_sqrt_radical'],1),
              parse(pre['field']['new_positive_sqrt_radical'],2)])
    zero=ff.zeros[3];c=ff.const
    def vp(v,k):return tuple(ff.lift(parse(x,k),k)for x in v)
    N=[vp(v,1)for v in old['normals']]+[vp(bridge['new_normal'],2)]+[vp(v,3)for v in pre['new_normals']]
    assert len(N)==98
    Q=[list(map(F,r))for r in old['form']]
    assert Q==[list(map(F,r))for r in d['form']]
    assert A.inertia(Q)==[20,1,0]
    caps=[list(map(F,r))for r in d['new_normals']]
    assert caps==[list(map(F,r))for r in cn['new_normals']]
    def dual(v):return[sum(x*q for x,q in zip(v,row))for row in Q]
    def rd(v,w):return sum(x*y for x,y in zip(v,dual(w)))
    def linear(v,w):return ff.total(ff.scale(x,q)for x,q in zip(w,dual(v))if q)
    qnew=[rd(v,v)for v in caps]
    assert qnew==[F(2),F(4)]and rd(*caps)==0
    norms=[c(F(q))for q in old['norms']]+[c(2)]+[parse(pre['new_normal_norm'],3)]*2
    newpairs=[[linear(r,n)for n in N]for r in caps]
    assert newpairs==[[parse(x,3)for x in row]for row in d['new_pairings']]
    def classify(g,p):
        if g==zero:return'right'
        assert ff.sign(g)<0
        # A strict interval inequality can certify dotted pairs without
        # multiplying very large boost-field coefficients unnecessarily.
        for bits in[192,384]:
            gl,gh=ff.interval(g,3,bits);pl,ph=ff.interval(p,3,bits)
            if gh<0 and gh*gh>ph:return'dotted'
        sq=ff.mul(g,g)
        for q,name in[(F(1,4),'pi3'),(F(1,2),'pi4'),(F(3,4),'pi6'),(F(1),'parallel')]:
            if sq==ff.scale(p,q):return name
        assert ff.sign(ff.add(sq,ff.scale(p,-1)))>0
        return'dotted'
    rels=[[classify(g,ff.scale(q,qnew[k]))for g,q in zip(row,norms)]for k,row in enumerate(newpairs)]
    assert rels==d['new_relations']
    counts=Counter(x for row in rels for x in row);counts['right']+=1
    assert dict(counts)==d['new_pair_counts']and sum(counts.values())==197
    finite={98+k:sorted([i for i,s in enumerate(row)if s in['right','pi3','pi4','pi6']]+[99-k])for k,row in enumerate(rels)}
    assert [finite[98],finite[99]]==d['complete_finite_neighbors']
    T=tuple(ff.total(n[k]for n in N[:95])for k in range(21))
    assert T==vp(d['strict_center']['point'],1)
    qT=ff.lift(parse(summed['positive_wall_sum_norm'],1),1)
    assert qT==ff.lift(parse(d['strict_center']['norm'],1),1)and ff.sign(qT)<0
    tp=[ff.lift(parse(x,1),1)for x in summed['positive_wall_sum_all_pairings']]
    assert [row['wall']for row in bridge['new_pairs']]==list(range(95))
    tp.append(ff.total(ff.lift(parse(row['pairing'],2),2)for row in bridge['new_pairs']))
    tp.extend(ff.total(parse(g,3)for g in row[:95])for row in pre['new_pairings'])
    newtp=[linear(r,T)for r in caps];tp+=newtp
    assert tp==[parse(x,3)for x in d['strict_center']['all100_pairings']]
    assert all(ff.sign(x)<0 for x in tp)
    for k,row in enumerate(d['new_facet_points']):
        alpha=ff.scale(newtp[k],1/qnew[k])
        P=tuple(ff.add(t,ff.scale(alpha,-x))for t,x in zip(T,caps[k]))
        assert P==vp(row['point'],1)
        qP=ff.add(qT,ff.scale(ff.mul(newtp[k],newtp[k]),-1/qnew[k]))
        assert qP==ff.lift(parse(row['norm'],1),1)and ff.sign(qP)<0
        # The strict signs follow directly from tp<0, alpha<0 and g<=0.
        # Independently check the saved pairings by the projection identity.
        signs=[parse(x,3)for x in row['all100_pairings']]
        full=newpairs[k]+[c(rd(caps[k],r))for r in caps]
        assert signs==[ff.add(t,ff.scale(ff.mul(alpha,g),-1))for t,g in zip(tp,full)]
        assert all(ff.sign(x)==(0 if i==98+k else -1)for i,x in enumerate(signs))
    allN=N+[tuple(c(x)for x in r)for r in caps]
    for row in d['preserved_bounded_cusps']:
        nu=vp(row['point'],0)
        gs=[linear(r,nu)for r in caps]
        assert gs==[ff.lift(parse(x,0),0)for x in row['new_pairings']]
        assert all(ff.sign(x)<0 for x in gs)
    E=list(map(F,d['strict_escape']['rational_point']))
    assert tuple(c(x)for x in E)==vp(pre['strict_escape']['point'],2)
    assert rd(E,E)==F(d['strict_escape']['norm'])>0
    assert [rd(E,r)for r in caps]==list(map(F,d['strict_escape']['new_pairings']))
    assert all(rd(E,r)<0 for r in caps)and ff.sign(linear(E,T))<0
    # Build the complete finite-neighbor pool from the verified full100 graph.
    pool=sorted((set(finite[98])|set(finite[99]))-{98,99})
    assert pool==cn['rational_old_wall_pool']
    rn={}
    for i in pool:
        v=[]
        for x in N[i]:
            for level in[3,2,1]:assert x[1]==ff.zeros[level-1];x=x[0]
            assert not any(x[1:]);v.append(x[0])
        rn[i]=v
    rn.update({98:caps[0],99:caps[1]})
    G={(i,j):rd(v,w)for i,v in rn.items()for j,w in rn.items()}
    graph={i:{j for j in rn if i!=j and G[i,j]**2<G[i,i]*G[j,j]}for i in rn}
    assert all(G[i,j]<=0 for i in rn for j in rn if i!=j)
    assert all(4*G[i,j]**2>=G[i,i]*G[j,j]for i in rn for j in rn if i!=j and G[i,j])
    assert {i:sorted(graph[i])for i in[98,99]}==finite
    assert {str(i):finite[i]for i in[98,99]}==cn['complete_finite_neighbors']
    zerosets=[[i for i,g in enumerate(row)if g==zero]for row in newpairs]
    assert zerosets==d['old98_zero_sets']
    zero_ranks=[rank([dual(rn[i])for i in row])for row in zerosets]
    assert zero_ranks==[20,20]
    # Direct fixed-size clique recursion: no maximal-clique data are needed.
    candidates=set()
    def extend(chosen,remaining,need):
        if need==0:candidates.add(tuple(sorted(chosen)));return
        while len(remaining)>=need:
            v=remaining[0];remaining=remaining[1:]
            extend(chosen+[v],[w for w in remaining if w in graph[v]],need-1)
    for i in[98,99]:extend([i],sorted(graph[i]),19)
    assert len(candidates)==cn['distinct_twenty_subsets']==20810
    # A connected graph with a cycle cannot be positive definite: after
    # normalizing diagonals to1, the all-ones vector on that cycle has square
    # <=0, since every nonzero finite edge has cosine at least1/2.
    # For a tree, exact leaf Schur elimination has no fill-in.
    cache={};kinds=Counter()
    def component_pd(ids):
        key=tuple(sorted(ids))
        if key in cache:return cache[key]
        adj={i:{j for j in ids if i!=j and G[i,j]}for i in ids}
        if sum(map(len,adj.values()))!=2*(len(ids)-1):
            kinds['cycle_obstruction']+=1;cache[key]=False;return False
        diagonal={i:G[i,i]for i in ids}
        while len(adj)>1:
            leaf=min(i for i,n in adj.items()if len(n)==1)
            if diagonal[leaf]<=0:
                kinds['nonpositive_leaf']+=1;cache[key]=False;return False
            parent=next(iter(adj[leaf]))
            diagonal[parent]-=G[leaf,parent]**2/diagonal[leaf]
            adj[parent].remove(leaf);del adj[leaf]
        ok=diagonal[next(iter(adj))]>0
        kinds['positive_tree'if ok else'nonpositive_root']+=1;cache[key]=ok;return ok
    def pd(ids):
        todo=set(ids)
        while todo:
            first=min(todo);todo.remove(first);part={first};queue=[first]
            while queue:
                i=queue.pop();extra={j for j in todo if G[i,j]}
                todo-=extra;part|=extra;queue.extend(extra)
            if not component_pd(part):return False
        return True
    positive={ids for ids in candidates if pd(ids)}
    assert positive=={tuple(v['incident'])for v in cn['vertices']}
    for row in cn['exact_component_checks']:
        assert component_pd(row['roots'])==row['positive']
    vertices=[];vcounts=Counter();points=set()
    for row in cn['vertices']:
        ids=tuple(row['incident']);x=list(map(F,row['primitive_rational_point']))
        assert tuple(x)not in points;points.add(tuple(x))
        qx=rd(x,x);assert qx==F(row['norm'])<0 and ff.sign(linear(x,T))<0
        signs=[ff.sign(linear(x,n))for n in allN]
        assert [i for i,s in enumerate(signs)if s==0]==list(ids)
        assert all(s<=0 for s in signs)
        for i in[98,99]:
            if i in ids:vcounts['on_'+str(i)]+=1
        if 98 in ids and 99 in ids:vcounts['on_both']+=1
        vertices.append(dict(incident=ids,norm=qx,zero_count=signs.count(0),strict_count=signs.count(-1)))
    assert len(vertices)==16 and dict(vcounts)==cn['vertex_counts']
    paths=[BASE,BRIDGE,PRE,SUM,PRIMARY,CENSUS,f'scripts/{TAG}.py',
           'scripts/fresh_h20_polar_pair_20261014.py','scripts/fresh_h20_polar_pair_20261014_census.py',
           'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
           'scripts/fresh_h20_bridge_faces_20261014_review.py',
           'results/fresh_h20_bridge_edge_closure_20261014_review.json']
    return A.enc(dict(status='PASS_INDEPENDENT_POLAR_PAIR_AND_COMPLETE_NEW_VERTEX_CENSUS',sources={p:sha(p)for p in paths},
        new_pair_counts=dict(counts),complete_finite_neighbors=finite,polar_old_zero_ranks=zero_ranks,
        candidate_twenty_cliques=len(candidates),independent_component_results=dict(kinds),
        checked_cached_component_count=len(cn['exact_component_checks']),positive_supports=vertices,
        vertex_counts=dict(vcounts),direct_vertex_inequalities=1600,preserved_old_cusps=3,
        inherited_old_finite_vertex_lower_bound=118,total_finite_vertex_lower_bound=134,
        escape_norm=rd(E,E),escape_new_pairings=[rd(E,r)for r in caps],maximum_interval_bits=ff.bits,
        scope='Independent197 new pairs, complete finite-neighbor lists, exact center/facet formulas, cusp and escape signs; direct fixed-size clique census, separate tree/cycle PD tests, and all1600 physical vertex signs. Frozen old98 full pair/corner census is not repeated.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['vertex_counts'],flush=True)
