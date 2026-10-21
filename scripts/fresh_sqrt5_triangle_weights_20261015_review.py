#!/usr/bin/env python3
"""Independent literal triangle/weight review using the separate Qsqrt5 helper."""
import argparse,importlib.util,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_triangle_weights_20261015_review'
PRIMARY='fresh_sqrt5_triangle_weights_20261015'
HELPER='fresh_sqrt5_quad_contacts_20261015_review'
spec=importlib.util.spec_from_file_location('independent_sqrt5_arithmetic',ROOT/f'scripts/{HELPER}.py')
Q=importlib.util.module_from_spec(spec);spec.loader.exec_module(Q)


def weights(nodes,g):
    adj=Q.graph(nodes,g);out={p:int(Q.distance(adj,*p)in(1,2))for p in combinations(sorted(nodes),2)}
    h4={p:out[p]for p in out};types=[]
    for c in Q.components(adj):
        assert Q.finite(c,g)
        five=[(i,j)for i,j in combinations(c,2)if g.get((i,j),2)==5]
        if not five:continue
        assert len(five)==1 and len(c)in(2,3,4);out[five[0]]+=1
        if len(c)==4:
            ends=tuple(i for i in c if len(adj[i])==1);assert len(ends)==2
            out[ends]+=1;h4[five[0]]+=1;h4[ends]+=1;types.append(('H4',c))
        else:types.append(('I2(5)'if len(c)==2 else'H3',c))
    return out,h4,types


def path(n,first=3):return {(i,i+1):first if i==0 else 3 for i in range(n-1)}
def arms(a):
    g={};k=1
    for length in a:
        prev=0
        for _ in range(length):g[(prev,k)]=3;prev,k=k,k+1
    return g


def model(name,n):
    if name.startswith('I2('):return path(2,int(name[3:-1]))
    if name[0]=='A':return path(n)
    if name[0]=='B':return path(n,4)
    if name[0]=='D':return arms((1,1,n-3))
    if name[0]=='E':return arms((1,2,n-4))
    if name=='F4':return {(0,1):3,(1,2):4,(2,3):3}
    if name[0]=='H':return path(n,5)
    raise AssertionError(name)


def run():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    cat=json.loads((ROOT/'results/fresh_sqrt5_lanner_catalogue_20261015.json').read_text())
    diagrams=[dict(order=3,labels=r['labels'])for r in cat['order_three']['accepted']]
    diagrams += [dict(order=b['order'],labels=r['labels'])for b in cat['higher_orders']for r in b['accepted_classes']]
    assert Counter(r['order']for r in diagrams)=={3:12,4:6,5:4}
    saved={(r['lanner_index'],tuple(r['black'])):r for r in primary['markings']};records=[];deficient=[]
    for index,row in enumerate(diagrams):
        n=row['order'];g=dict(zip(combinations(range(n),2),row['labels']))
        for black in combinations(range(n),3):
            record=saved[index,black];assert record['labels']==row['labels']and record['order']==n
            corners=[]
            for omit in black:
                nodes=[i for i in range(n)if i!=omit];marked=tuple(i for i in black if i!=omit)
                labels={Q.pair(i,j):g[Q.pair(nodes[i],nodes[j])]for i,j in combinations(range(n-1),2)}
                M=Q.tree_scaled(n-1,labels);assert M is not None and Q.inertia(M)==[n-1,0,0]
                w,h4,types=weights(nodes,g);dist=Q.distance(Q.graph(nodes,g),*marked)
                prior=next(c for c in record['corners']if c['omitted']==omit)
                assert w[marked]==prior['weight']and h4[marked]==prior['h4_only_weight']and dist==prior['distance']
                corners.append(dict(omitted=omit,marked=list(marked),distance=dist,weight=w[marked],
                                    h4_only_weight=h4[marked],independent_full_inertia=[n-1,0,0]))
            base=sum(c['distance']in(1,2)for c in corners);total=sum(c['weight']for c in corners)
            assert total>=2 and sum(c['h4_only_weight']for c in corners)>=2
            assert base==record['base_total']and total==record['total']
            if base<2:deficient.append(len(records))
            records.append(dict(lanner_index=index,order=n,labels=row['labels'],black=list(black),corners=corners,base=base,total=total))
    assert len(records)==76 and len(deficient)==7 and deficient==primary['base_deficient_marking_indices']
    finite=[]
    for row in primary['finite_component_totals']:
        n=row['rank'];name=row['type'];g=model(name,n);adj=Q.graph(range(n),g)
        assert Q.finite(list(range(n)),g)
        # A tree has one edge path per distance1 pair and one two-edge path
        # for every unordered pair of neighbors of a vertex.
        base=n-1+sum(len(v)*(len(v)-1)//2 for v in adj.values())
        bonus=1 if name in('I2(5)','H3')else 2 if name=='H4'else 0
        total=base+bonus;deficit=2*n-total
        assert total==row['total']and deficit==row['deficit']and deficit>=1
        assert deficit>=2 or name=='H4'and n==4
        if n<=8:
            M=Q.tree_scaled(n,g);assert M is not None and Q.inertia(M)==[n,0,0]
        finite.append(dict(type=name,rank=n,base=base,bonus=bonus,total=total,deficit=deficit))
    best=[0]
    for n in range(1,101):
        best.append(max(best[n-r['rank']]+r['total']for r in finite if r['rank']<=n))
        if n>=5:assert best[n]==2*n-2
    assert best==primary['rank_partition_maxima']
    five=[('I2(5)',2),('H3',3),('H4',4)];extensions=[]
    for a,n in five:
        g=model(a,n);w,_,_=weights(list(range(n)),g)
        for b,m in five:
            if n>m:continue
            h=model(b,m);z,_,_=weights(list(range(m)),h);count=0
            for image in permutations(range(m),n):
                if any(g.get((i,j),2)!=h.get(Q.pair(image[i],image[j]),2)for i,j in combinations(range(n),2)):continue
                assert all(z[Q.pair(image[i],image[j])]>=v for(i,j),v in w.items());count+=1
            extensions.append(dict(source=a,target=b,induced_embeddings=count))
        for perm in permutations(range(n)):
            relabel={Q.pair(perm[i],perm[j]):lab for(i,j),lab in g.items()};z,_,_=weights(list(range(n)),relabel)
            assert all(z[Q.pair(perm[i],perm[j])]==v for(i,j),v in w.items())
    assert extensions==primary['five_component_extensions']
    paths=[f'scripts/{TAG}.py',f'scripts/{HELPER}.py']
    for stem in(PRIMARY,'fresh_sqrt5_lanner_catalogue_20261015','fresh_sqrt5_lanner_catalogue_20261015_review'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json']
    manifests=Q.manifest_check(paths)
    return dict(status='PASS_INDEPENDENT_QSQRT5_TRIANGLE_WEIGHTS_AND_EXTENSION_STABILITY_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},recursive_json_manifests=manifests,
                markings=records,deficient_markings=deficient,component_totals=finite,partition_maxima=best,
                induced_five_component_embeddings=extensions,
                all_rank_proof='Only H4 has deficit1; all other irreducible finite components have deficit>=2. Rank>=5 implies total deficit>=2.',
                extension_proof='Connected old paths remain unique in every spherical forest. A finite5-edge component is I2(5),H3,H4; H4 cannot grow, so all bonuses persist.',
                scope='All compact arithmeticQsqrt5 triangular2faces, arbitrary white extensions. Vertexbound2n-2 for n>=5. No quadrilateral or dimension conclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text());print('PASS saved triangle review replay')
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
