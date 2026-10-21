#!/usr/bin/env python3
"""Finite-event replacement of the two walls exposed by deleting34."""
import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_release34_cavity_20261017 as B

ROOT=B.ROOT
TAG='fresh_nonarithmetic_release34_replace12_47_20261017'
ALPHABET=list(range(2,129))+[192,256,512]


def roots(frame, contacts):
    gram=frame@B.METRIC@frame.T
    p=np.array(contacts)@np.linalg.inv(gram)@frame
    k=np.linalg.svd(frame@B.METRIC)[2][-1]
    q=k@B.METRIC@k
    assert abs(q)>1e-10
    rem=(2-p@B.METRIC@p)/q
    if rem<=1e-10:return []
    out=[]
    for sign in [-1,1]:
        n=p+sign*np.sqrt(rem)*k
        assert abs(n@B.METRIC@n-2)<1e-6
        assert np.max(abs(frame@B.METRIC@n-contacts))<1e-7
        out.append((sign,n))
    return out


def gate(n, normals, labels, h, forced):
    if n@B.METRIC@h>=-1e-7:return 'center',None
    products=normals@B.METRIC@n
    for i,m in forced.items():
        j=labels.index(i)
        assert abs(products[j]+B.theta(m))<1e-6
        products[j]=-B.theta(m)
    if products.max()>2e-7:return 'obtuse',None
    relations=[B.relation(x)for x in products]
    if any(r is None or r=='parallel_or_unresolved'for r in relations):return 'gap',None
    return 'accepted',[[i,r]for i,r in zip(labels,relations)if isinstance(r,int)]


def compute(drop54=False):
    z,g,t,labels,normals,h,sides=B.setup()
    ns=dict(zip(labels,normals))
    actual=[labels[i]for i in sides]
    if drop54:
        assert 54 not in actual
        labels=[i for i in labels if i!=54]
    counts=Counter();first=[];second=[]
    ls=[i for i in labels if i!=12];nn=np.array([ns[i]for i in ls])
    # A replacement12 fixes its two surviving finite incidences. Either
    # finite event with52 or47 can terminate an exterior old edge.
    for anchor in [52,47]:
        frame=np.array([ns[i]for i in [7,8,anchor]])
        for m in ALPHABET:
            for sign,n in roots(frame,[0.,0.,-B.theta(m)]):
                status,rel=gate(n,nn,ls,h,{7:2,8:2,anchor:m})
                counts['first_'+status]+=1
                if status!='accepted':continue
                row=dict(index=len(first),anchor=anchor,label=m,branch=sign,
                         normal=n.tolist(),finite_contacts=rel)
                nmap={**ns,12:n}
                row['boundary']=B.ledger(np.array([nmap[i]for i in actual]),actual)
                first.append(row)
    print('first replacements',len(first),dict(counts),flush=True)
    # Independently move47 while fixing7/30 and giving it one finite contact
    # to new12. All retained walls and replacement12 are then checked again.
    for row in first:
        if row['anchor']!=52:continue
        nmap={**ns,12:np.array(row['normal'])}
        ls=[i for i in labels if i!=47];nn=np.array([nmap[i]for i in ls])
        frame=np.array([nmap[i]for i in [7,30,12]])
        for m in ALPHABET:
            for sign,n in roots(frame,[0.,0.,-B.theta(m)]):
                status,rel=gate(n,nn,ls,h,{7:2,30:2,12:m})
                counts['second_'+status]+=1
                if status!='accepted':continue
                new={**nmap,47:n}
                boundary=B.ledger(np.array([new[i]for i in actual]),actual)
                second.append(dict(first_index=row['index'],label=m,branch=sign,
                    normal=n.tolist(),finite_contacts=rel,boundary=boundary))
    first.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],r['label']))
    second.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],r['label']))
    return dict(status='BOUNDED_NUMERICAL_TWO_WALL_REPLACEMENT_EVENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in [B.PARENT,Path(B.__file__)]},
        alphabet=ALPHABET,counts=dict(counts),first_replacements=first,second_replacements=second,
        deleted_additional_old_walls=[54]if drop54 else [],
        baseline=B.ledger(normals[sides],actual),
        scope='Preserves all65 labels, replaces12 and optionally47. Numerical discovery only, complete boundary scored including deletion costs; finite event labels are bounded.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--drop54',action='store_true');args=p.parse_args()
    out=compute(args.drop54);path=ROOT/'results'/(TAG+('_drop54'if args.drop54 else '')+'.json')
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
    for key in ['first_replacements','second_replacements']:
        if out[key]:print(key,out[key][0],flush=True)


if __name__=='__main__':main()
