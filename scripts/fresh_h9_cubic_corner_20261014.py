#!/usr/bin/env python3
"""Exact bounded-distance corner and lattice-symmetry search on saved cubic hosts."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import networkx as nx
import audit_fresh_cubic7_h9_no_import_20261013 as K

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_cubic_corner_20261014'
Z,O=K.ZERO,K.ONE
def total(xs):
    z=Z
    for x in xs:z=K.add(z,x)
    return z
def transpose(A):return list(map(list,zip(*A)))
def mm(A,B):return [[total(K.mul(x,y)for x,y in zip(r,c))for c in zip(*B)]for r in A]
def mv(A,v):return [total(K.mul(x,y)for x,y in zip(r,v))for r in A]
def inverse(A):
    n=len(A);B=[list(r)+[O if i==j else Z for j in range(n)]for i,r in enumerate(A)];det=O
    for j in range(n):
        p=next(i for i in range(j,n)if B[i][j]!=Z)
        if p!=j:B[p],B[j]=B[j],B[p];det=K.neg(det)
        d=B[j][j];det=K.mul(det,d);inv=K.inverse(d);B[j]=[K.mul(x,inv)for x in B[j]]
        for i in range(n):
            if i==j or B[i][j]==Z:continue
            c=B[i][j];B[i]=[K.sub(x,K.mul(c,y))for x,y in zip(B[i],B[j])]
    return [r[n:]for r in B],det
def pd(A):
    A=[list(r)for r in A];piv=[]
    for j in range(len(A)):
        p=A[j][j];assert K.sign(p)>0;piv.append(p)
        for i in range(j+1,len(A)):
            for k in range(j+1,len(A)):A[i][k]=K.sub(A[i][k],K.divide(K.mul(A[i][j],A[j][k]),p))
    return piv
def encode(x):
    if isinstance(x,tuple)and len(x)==3 and all(hasattr(v,'denominator')for v in x):return K.encoded(x)
    if isinstance(x,(tuple,list)):return [encode(v)for v in x]
    if isinstance(x,dict):return {str(k):encode(v)for k,v in x.items()}
    return x
def floating(x):return sum(float(c)*(2*np.cos(2*np.pi/7))**i for i,c in enumerate(x))
def integral(A):return all(c.denominator==1 for r in A for x in r for c in x)

def run():
    scout=json.loads((ROOT/f'results/{TAG}_scout.json').read_text());rows=[]
    for proposal in scout['seeds']:
        source=proposal['source'];d=json.loads((ROOT/source).read_text());K.verify(d)
        a=K.decode(d['timelike_coefficient']);N=[[K.decode(x)for x in r]for r in d['all_walls']]
        Q=[[Z for _ in range(10)]for _ in range(10)]
        Q[0][0]=K.neg(a);Q[9][9]=O
        for i in range(1,9):Q[i][i]=K.scale(O,2)
        for i,j in K.EDGES:Q[i][j]=Q[j][i]=K.neg(O)
        cov=[K.gram_vector(r,a)for r in N];G=[[K.pairing(r,s,a)for s in N]for r in N]
        corners=[]
        for item in proposal['proposed_corners']+proposal['beyond_bound_proposals']:
            ids=item['incident'];B=[[G[i][j]for j in ids]for i in ids];piv=pd(B)
            inv,_=inverse([cov[i][1:]for i in ids]);p=[O]+mv(inv,[K.neg(cov[i][0])for i in ids]);q=K.pairing(p,p,a)
            assert K.sign(q)<0
            products=[K.pairing(p,r,a)for r in N]
            assert [i for i,x in enumerate(products)if x==Z]==ids and all(K.sign(x)<=0 for x in products)
            h=K.sub(K.divide(a,K.neg(q)),O);bounded=K.sign(K.sub(K.scale(O,100),h))>0
            frame=transpose([N[i]for i in ids]+[p]);fi,det=inverse(frame)
            graph=nx.Graph()
            for j,i in enumerate(ids):graph.add_node(j,norm=G[i][i])
            graph.add_edges_from((j,k,{'pair':G[ids[j]][ids[k]]})for j in range(9)for k in range(j)if G[ids[j]][ids[k]]!=Z)
            corners.append(dict(ids=ids,p=p,norm=q,height=h,pivots=piv,within_complete_height=bounded,frame=frame,inverse=fi,det=det,graph=graph))
        maps=[];keys=set()
        for i,c in enumerate(corners):
            for j,e in enumerate(corners):
                if j<i:continue
                matcher=nx.algorithms.isomorphism.GraphMatcher(c['graph'],e['graph'],
                    node_match=lambda x,y:x['norm']==y['norm'],edge_match=lambda x,y:x['pair']==y['pair'])
                for perm in matcher.isomorphisms_iter():
                    target=transpose([N[e['ids'][perm[k]]]for k in range(9)]+[e['p']]);_,tdet=inverse(target)
                    scale=K.divide(c['det'],tdet)
                    if K.sign(scale)<0:scale=K.neg(scale)
                    for k in range(10):target[k][9]=K.mul(target[k][9],scale)
                    T=mm(target,c['inverse'])
                    if not integral(T):continue
                    key=tuple(x for r in T for x in r)
                    if key in keys:continue
                    keys.add(key);Ti,dt=inverse(T)
                    assert integral(Ti)and dt in(O,K.neg(O))and mm(mm(transpose(T),Q),T)==Q
                    assert mv(T,c['p'])==[K.mul(scale,x)for x in e['p']]and K.sign(scale)>0
                    for k in range(9):assert mv(T,N[c['ids'][k]])==N[e['ids'][perm[k]]]
                    maps.append(dict(source_corner=i,target_corner=j,permutation=[perm[k]for k in range(9)],
                                     time_scale=scale,matrix=T,inverse=Ti,
                                     full_sectors_certified=c['within_complete_height']and e['within_complete_height']))
        witness=None
        floats=[np.array([[floating(x)for x in r]for r in t['matrix']])for t in maps]
        for i,Ti in enumerate(floats):
            if witness:break
            for j,Tj in enumerate(floats):
                if not maps[i]['full_sectors_certified']or not maps[j]['full_sectors_certified']:continue
                T=Ti@Tj
                if abs(np.trace(T))<10.0001:continue
                prod=mm(maps[i]['matrix'],maps[j]['matrix']);tr=total(prod[k][k]for k in range(10))
                if K.sign(K.sub(tr,K.scale(O,10)))>0 or K.sign(K.add(tr,K.scale(O,10)))<0:
                    witness=dict(map_indices=[i,j],trace=tr,product=prod);break
        rows.append(dict(seed=d['seed_index'],source=source,
                         exact_prefix_corners=[{k:v for k,v in c.items()if k in('ids','p','norm','height','pivots','within_complete_height')}for c in corners],
                         integral_sector_maps=maps,infinite_order_witness=witness,
                         status='EXACT_INFINITE_FULL_CHAMBER_STABILIZER'if witness else'NO_INFINITE_STABILIZER_FOUND_IN_PROPOSED_MAPS'))
        print('seed',d['seed_index'],'corners',len(corners),'integral maps',len(maps),'infinite witness',bool(witness),'pending',[(m['source_corner'],m['target_corner'])for m in maps if not m['full_sectors_certified']],flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{TAG}_scout.py',f'results/{TAG}_scout.json',
           'scripts/audit_fresh_cubic7_h9_no_import_20261013.py','scripts/search_fresh_cubic7_h9_20261013.py']+[x['source']for x in rows]
    return encode(dict(status='EXACT_BOUNDED_DISTANCE_CORNERS_AND_SECTOR_MAPS',
                       sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},seeds=rows,
                       scope='Corner proposals are numerical and not an exhaustive census. Accepted coordinates, heights and maps are exact. Full-chamber sector recognition additionally uses the inherited exhaustive Vinberg height100 prefixes and norm alphabets. Infinite trace gives whole-host nonreflectivity only for rows with witnesses.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(p.read_text())==out
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
