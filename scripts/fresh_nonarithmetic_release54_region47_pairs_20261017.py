#!/usr/bin/env python3
"""Two end-cap discovery on the changed12 three-ray exterior component."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release54_region47_20261017 as S
B,R,ROOT=S.B,S.R,S.ROOT
TAG='fresh_nonarithmetic_release54_region47_pairs_20261017'
ALPHABET=list(range(2,49))+[64,96,128,192,256]
FRAMES=[[7,12,47],[30,47,52]]


def compute():
    z,labels,ns,h=S.setup();nn=np.array([ns[i]for i in labels]);counts=Counter();records=[];groups=[]
    for group,ids in enumerate(FRAMES):
        frame=np.array([ns[i]for i in ids]);gram=frame@B.METRIC@frame.T
        inv=np.linalg.inv(gram);kernel=np.linalg.svd(frame@B.METRIC)[2][-1]
        kq=kernel@B.METRIC@kernel;assert kq>0
        kernel/=np.sqrt(kq);indices=[labels.index(i)for i in ids];accepted=[]
        for ms in product(ALPHABET,repeat=3):
            c=-np.array([B.theta(m)for m in ms])
            if any(np.linalg.eigvalsh(np.array([[gram[0,0],gram[0,j],c[0]],
                [gram[0,j],gram[j,j],c[j]],[c[0],c[j],2.]]))[0]<=1e-9 for j in [1,2]):continue
            counts['finite_profile']+=1;p=c@inv@frame;rem=2-p@B.METRIC@p
            if rem<=1e-9:counts['no_positive_norm_remainder']+=1;continue
            for sign in [-1,1]:
                n=p+sign*np.sqrt(rem)*kernel
                prods=nn@B.METRIC@n;prods[indices]=c
                if n@B.METRIC@h>=-1e-8:counts['wrong_center']+=1;continue
                if prods.max()>2e-7:counts['obtuse']+=1;continue
                rel=[B.relation(x)for x in prods]
                if any(x is None or x=='parallel_or_unresolved'for x in rel):counts['gap']+=1;continue
                counts['accepted']+=1
                row=dict(index=len(records),group=group,frame=ids,labels=list(ms),branch=sign,normal=n.tolist(),
                    finite_contacts=[[i,r]for i,r in zip(labels,rel)if isinstance(r,int)])
                records.append(row);accepted.append(row['index'])
        groups.append(accepted);print('group',group,len(accepted),dict(counts),flush=True)
    pairs=[];relations=Counter();smallest=[]
    if all(groups):
        aa=np.array([records[i]['normal']for i in groups[0]]);bb=np.array([records[i]['normal']for i in groups[1]])
        products=aa@B.METRIC@bb.T
        for ii,jj in product(range(len(aa)),range(len(bb))):
            x=float(products[ii,jj]);r=B.relation(x)
            relations[str(r)]+=1
            if r is None or r=='parallel_or_unresolved':continue
            if not isinstance(r,int):continue
            ia,ib=groups[0][ii],groups[1][jj]
            actual=z['actual_sides']+[66,67]
            roots=np.array([ns[i]for i in z['actual_sides']]+[records[ia]['normal'],records[ib]['normal']])
            ledger=B.ledger(roots,actual)
            pairs.append(dict(indices=[ia,ib],product=x,label=r,boundary=ledger))
            print('finite pair',ia,ib,r,ledger['open_ends'],flush=True)
        flat=np.argsort(products.ravel())[:12]
        smallest=[dict(indices=[groups[0][k//len(bb)],groups[1][k%len(bb)]],product=float(products.ravel()[k]))for k in flat]
    return dict(status='BOUNDED_NUMERICAL_CHANGED12_ENDPOINT_PAIR_SCREEN',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in [S.PARENT,Path(S.__file__),Path(B.__file__)]},
        alphabet=ALPHABET,frames=FRAMES,counts=dict(counts),records=records,groups=groups,
        pair_relation_counts=dict(relations),finite_pairs=pairs,smallest_products=smallest,
        scope='Both end caps terminate the two geodesic edges at their indicated old exterior endpoint. All64 native gates checked; finite label list bounded; finite mutual pairs scored on entire S6 boundary.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=compute()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['pair_relation_counts'],out['smallest_products'],flush=True)

if __name__=='__main__':main()
