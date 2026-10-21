#!/usr/bin/env python3
"""Integral local chamber automorphisms at the connected certified corners."""
from fresh_h9_quartic_edge_walk_20261013_next import *
from itertools import permutations
import hashlib

OUT=ROOT/'results/fresh_h9_quartic_next_frontier_20261013_symmetries.json'


def run():
    paths=['results/fresh_h9_quartic_first_endpoint_20261013'+s+'.json'for s in('','_second','_third')]
    paths += [f'results/fresh_h9_quartic_next_frontier_20261013_census_edge{i}.json'for i in(2,3)]
    g=Host(1).g;identity=tuple(tuple(K(i==j)for j in range(6))for i in range(6));corners=[];autos=[];seen=set()
    for index,path in enumerate(paths):
        data=json.loads((ROOT/path).read_text());edge=[tuple(map(K.from_json,r))for r in data['edge_roots']]
        kinds=['nearest']if index else['start','nearest']
        for kind in kinds:
            if kind=='start':
                base=json.loads((ROOT/data['base_source']).read_text())
                if data['base_source_case']is not None:base=base['cases'][data['base_source_case']]
                root=tuple(map(K.from_json,base['forced_H5_face']['all_walls'][data['start_corner_wall']]))
                point=tuple(map(K.from_json,data['start_corner']))
            else:root=tuple(map(K.from_json,data['nearest_root']));point=tuple(map(K.from_json,data['nearest_point']))
            roots=edge+[root];gram=tuple(tuple(pair(r,g,s)for s in roots)for r in roots)
            frame=tuple(zip(*roots,point));fi=inverse(frame);count=0;integral_count=0;corner_index=len(corners)
            for perm in permutations(range(5)):
                if any(gram[i][j]!=gram[perm[i]][perm[j]]for i in range(5)for j in range(5)):continue
                count+=1;target=tuple(zip(*(roots[i]for i in perm),point));F=matmul(target,fi);F_inv=inverse(F)
                assert matmul(matmul(tuple(zip(*F)),g),F)==g and mv(F,point)==point
                if not all(z.is_integral for row in F+F_inv for z in row):continue
                integral_count+=1
                if F in seen:continue
                seen.add(F);autos.append({'corner_index':corner_index,'permutation':list(perm),'matrix':menc(F),'inverse':menc(F_inv),
                                         'involution':matmul(F,F)==identity,'trace':sum((F[i][i]for i in range(6)),Z).to_json()})
            corners.append({'source':path,'kind':kind,'roots':menc(roots),'point':enc(point),'gram':menc(gram),
                            'gram_automorphisms':count,'integral_automorphisms':integral_count})
    products=[]
    for i,row in enumerate(autos):
        A=tuple(tuple(map(K.from_json,r))for r in row['matrix'])
        for j,second in enumerate(autos[:i]):
            B=tuple(tuple(map(K.from_json,r))for r in second['matrix']);F=matmul(A,B)
            tr=sum((F[k][k]for k in range(6)),Z)
            if tr>6 or tr<-6:products.append({'factors':[i,j],'matrix':menc(F),'trace':tr.to_json()})
    out={'status':'EXACT_LOCAL_AUTOMORPHISM_CENSUS','sources_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
         'corners':corners,'distinct_integral_automorphisms':autos,'products_with_trace_above_dimension':products,
         'scope':'Exact primitive-Gram automorphisms fixing one certified local corner and its full finite chamber sector. These corners are connected by the certified nearest-edge chain. Each integral automorphism therefore stabilizes the same full reflection chamber. A trace magnitude above six would certify infinite order; absence of such a two-factor product alone is not a reflectivity theorem.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print({'corner_counts':[(r['gram_automorphisms'],r['integral_automorphisms'])for r in corners],
           'distinct':len(autos),'infinite_products':len(products)},flush=True)


if __name__=='__main__':run()
