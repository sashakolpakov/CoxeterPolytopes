#!/usr/bin/env python3
"""Search exact sector maps in the transported cubic prefix, including fixed-point gates."""
import argparse,json
from hashlib import sha256
import networkx as nx
import numpy as np
import fresh_h9_cubic_corner_20261014 as P
K=P.K;Z,O=P.Z,P.O;ROOT=P.ROOT;TAG='fresh_h9_cubic_corner_20261014_symmetries'

def run(src='results/fresh_h9_cubic_corner_20261014_extension.json',esrc='results/fresh_h9_cubic_corner_20261014_edges.json'):
    d=json.loads((ROOT/src).read_text());ed=json.loads((ROOT/esrc).read_text())
    N=[[K.decode(x)for x in r]for r in d['normals']];a=K.decode(d['timelike_coefficient']);bound=d.get('complete_distance_bound',100)
    G=[[K.pairing(r,s,a)for s in N]for r in N];I=[[O if i==j else Z for j in range(10)]for i in range(10)]
    Q=[K.gram_vector(e,a)for e in I]
    oldT=[[K.decode(x)for x in r]for r in d['isometry']];centers=[I[0],P.mv(oldT,I[0])]
    keys={tuple(r):i for i,r in enumerate(N)};corners=[]
    for v in ed['finite_vertices']:
        ids=v['incident'];p=[K.decode(x)for x in v['point']];B=P.transpose([N[i]for i in ids]+[p]);inv,det=P.inverse(B)
        graph=nx.Graph()
        for j,i in enumerate(ids):graph.add_node(j,norm=G[i][i])
        graph.add_edges_from((j,k,{'pair':G[ids[j]][ids[k]]})for j in range(9)for k in range(j)if G[ids[j]][ids[k]]!=Z)
        corners.append(dict(ids=ids,p=p,inv=inv,det=det,graph=graph,full=v['full_chamber_corner_certified']))
    found=[];seen=set()
    for i,c in enumerate(corners):
        for j,e in enumerate(corners):
            if j<i:continue
            matcher=nx.algorithms.isomorphism.GraphMatcher(c['graph'],e['graph'],node_match=lambda x,y:x['norm']==y['norm'],edge_match=lambda x,y:x['pair']==y['pair'])
            for perm in matcher.isomorphisms_iter():
                B=P.transpose([N[e['ids'][perm[k]]]for k in range(9)]+[e['p']]);_,det=P.inverse(B);s=K.divide(c['det'],det)
                if K.sign(s)<0:s=K.neg(s)
                for k in range(10):B[k][9]=K.mul(B[k][9],s)
                T=P.mm(B,c['inv']);key=tuple(x for r in T for x in r)
                if key in seen or not P.integral(T):continue
                seen.add(key);inv,dt=P.inverse(T);assert P.integral(inv)and dt in(O,K.neg(O))and P.mm(P.mm(P.transpose(T),Q),T)==Q
                assert P.mv(T,c['p'])==[K.mul(s,x)for x in e['p']]and K.sign(s)>0
                accepted=c['full']and e['full'];fixed=None
                if not accepted and P.mm(T,T)==I:
                    H=[K.add(x,y)for x,y in zip(centers[0],P.mv(T,centers[0]))];q=K.pairing(H,H,a)
                    assert K.sign(q)<0 and P.mv(T,H)==H
                    heights=[K.sub(K.divide(K.mul(K.pairing(H,v,a),K.pairing(H,v,a)),K.mul(K.neg(q),a)),O)for v in centers]
                    ips=[K.pairing(H,r,a)for r in N]
                    if all(K.sign(x)<=0 for x in ips)and any(K.sign(K.sub(K.scale(O,bound),h))>0 for h in heights):
                        incident=[k for k,x in enumerate(ips)if x==Z]
                        targets=[keys.get(tuple(P.mv(T,N[k])))for k in incident]
                        if None not in targets and sorted(targets)==incident:
                            accepted=True;fixed=dict(point=H,norm=q,heights=heights,incident=incident,permutation=targets)
                found.append(dict(source_corner=i,target_corner=j,permutation=[perm[k]for k in range(9)],matrix=T,inverse=inv,
                                  certified_chamber_symmetry=accepted,fixed_point_certificate=fixed))
    good=[i for i,v in enumerate(found)if v['certified_chamber_symmetry']];witness=None
    fs={i:np.array([[P.floating(x)for x in r]for r in found[i]['matrix']])for i in good}
    for i in good:
        if witness:break
        for j in good:
            if abs(np.trace(fs[i]@fs[j]))<10.0001:continue
            prod=P.mm(found[i]['matrix'],found[j]['matrix']);tr=P.total(prod[k][k]for k in range(10))
            if K.sign(K.sub(tr,K.scale(O,10)))>0 or K.sign(K.add(tr,K.scale(O,10)))<0:
                witness=dict(map_indices=[i,j],trace=tr,product=prod);break
    paths=[src,esrc,'scripts/fresh_h9_cubic_corner_20261014.py',f'scripts/{TAG}.py']
    return P.encode(dict(status='EXACT_EXTENDED_PREFIX_SECTOR_MAP_SEARCH',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                         integral_maps=found,certified_map_indices=good,infinite_order_witness=witness,
                         scope='Exact maps between all finite prefix corners with literal norm-preserving diagram isomorphisms. Chamber certification uses known full corners or a near fixed-point sector. No conclusion from unvalidated maps.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--source',default='results/fresh_h9_cubic_corner_20261014_extension.json');ap.add_argument('--edges',default='results/fresh_h9_cubic_corner_20261014_edges.json');ap.add_argument('--output',default=f'results/{TAG}.json');args=ap.parse_args();out=run(args.source,args.edges);p=ROOT/args.output
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'maps',len(out['integral_maps']),'certified',len(out['certified_map_indices']),'infinite',bool(out['infinite_order_witness']))
