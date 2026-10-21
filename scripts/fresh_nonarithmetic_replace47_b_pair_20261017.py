#!/usr/bin/env python3
"""Minimal paired endpoint caps on the remaining two-ray B component."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release54_region47_20261017 as S
B,R,ROOT=S.B,S.R,S.ROOT
TAG='fresh_nonarithmetic_replace47_b_pair_20261017'
PARENT=ROOT/'results/fresh_nonarithmetic_replace47_close12_boundary_20261017.json'
ALPHABET=list(range(2,129))+[192,256,512,1024]

def setup():
    z=json.loads(PARENT.read_text())['after64']
    g=np.array([[B.mid(x)for x in row]for row in z['residual_gram']]);g=(g+g.T)/2
    vals,vecs=np.linalg.eigh(g);assert vals[0]<0<vals[1]
    transform=vecs@np.diag(1/np.sqrt(abs(vals)))
    labels=[r['wall']for r in z['projections']]
    cov=np.array([[B.mid(x)for x in r['covector']]for r in z['projections']]);normals=cov@transform@B.METRIC
    h=np.linalg.solve(transform,np.array([B.mid(x)for x in z['projected_common_center']]));h/=np.sqrt(-h@B.METRIC@h)
    return z,labels,dict(zip(labels,normals)),h

def compute():
    z,labels,ns,h=setup();old=np.array([ns[i]for i in labels]);groups=[];records=[];counts=Counter()
    profiles=[[(2,2,m)for m in [2,3,4]],[(m30,m47,2)for m30,m47 in product([2,3],ALPHABET)]+[(2,m47,3)for m47 in [2,3,4,5]]]
    frames=[[7,30,47],[30,47,52]]
    for group,(ids,prs)in enumerate(zip(frames,profiles)):
        frame=np.array([ns[i]for i in ids]);accepted=[]
        for ms in prs:
            for sign,n in R.roots(frame,[-B.theta(m)for m in ms]):
                status,rel=R.gate(n,old,labels,h,dict(zip(ids,ms)));counts[f'{group}_{status}']+=1
                if status!='accepted':continue
                accepted.append(len(records));records.append(dict(index=len(records),group=group,frame=ids,labels=list(ms),branch=sign,normal=n.tolist(),finite_contacts=rel))
        groups.append(accepted)
    pairs=[];relation_counts=Counter();ranges=[]
    for ia in groups[0]:
        a=np.array(records[ia]['normal']);products=[]
        for ib in groups[1]:
            b=np.array(records[ib]['normal']);x=float(a@B.METRIC@b);r=B.relation(x);relation_counts[str(r)]+=1;products.append(x)
            if not isinstance(r,int):continue
            actual=z['actual_sides']+[66,67];normals=np.array([ns[i]for i in z['actual_sides']]+[a,b])
            boundary=B.ledger(normals,actual)
            pairs.append(dict(indices=[ia,ib],product=x,label=r,boundary=boundary))
            print('FINITE_PAIR',ia,ib,r,boundary['open_ends'],flush=True)
        ranges.append(dict(left_index=ia,minimum=min(products),maximum=max(products)))
    return dict(status='BOUNDED_NUMERICAL_TWO_RAY_B_ENDPOINT_PAIR_SCREEN',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[PARENT,Path(S.__file__),Path(R.__file__),Path(B.__file__)]},
        left_profiles_complete=True,right_profiles=profiles[1],right_tail_sampled_labels=ALPHABET,
        frames=frames,records=records,groups=groups,counts=dict(counts),pair_relation_counts=dict(relation_counts),product_ranges=ranges,pairs=pairs,
        scope='Left cap handles both7/30 and7/47 endpoints, giving complete3 finite profiles. Right cap handles30/52 and47/52; its two isolated-label tails are sampled only. All64 native gates and full S6 boundary are tested; no all-label pair obstruction from the bounded sample.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(compute()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['pair_relation_counts'],out['product_ranges'],flush=True)

if __name__=='__main__':main()
