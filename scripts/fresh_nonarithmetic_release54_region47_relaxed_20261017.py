#!/usr/bin/env python3
"""Replace47 using any three boundary contacts, rather than preserving7/30."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release54_region47_20261017 as S
B,R,ROOT=S.B,S.R,S.ROOT
TAG='fresh_nonarithmetic_release54_region47_relaxed_20261017'
ALPHABET=list(range(2,17))+[24,32,64,128]


def compute():
    z,labels,ns,h=S.setup();remaining=[i for i in labels if i!=47];nn=np.array([ns[i]for i in remaining])
    frames=list(combinations([7,12,30,52],3));counts=Counter();blockers=Counter();records=[];seen=[]
    for ids in frames:
        frame=np.array([ns[i]for i in ids]);gram=frame@B.METRIC@frame.T
        inv=np.linalg.inv(gram);k=np.linalg.svd(frame@B.METRIC)[2][-1];kq=k@B.METRIC@k
        assert abs(kq)>1e-9
        for ms in product(ALPHABET,repeat=3):
            c=-np.array([B.theta(m)for m in ms]);p=c@inv@frame;rem=(2-p@B.METRIC@p)/kq
            if rem<=1e-9:continue
            if not any(np.linalg.eigvalsh(np.array([[gram[i,i],gram[i,j],c[i]],[gram[j,i],gram[j,j],c[j]],[c[i],c[j],2.]]))[0]>1e-8 for i,j in combinations(range(3),2)):continue
            for sign in [-1,1]:
                n=p+sign*np.sqrt(rem)*k;counts['branches']+=1
                if n@B.METRIC@h>=-1e-8:counts['center']+=1;continue
                products=nn@B.METRIC@n
                for i,v in zip(ids,c):products[remaining.index(i)]=v
                rel=[B.relation(x)for x in products]
                bad=[i for i,r in zip(remaining,rel)if r is None or r=='parallel_or_unresolved']
                blockers[tuple(bad)]+=1
                if len(bad)>1:continue
                # Permit deletion of one intrinsically redundant oldwall, but
                # never silently remove a compact-heptagon side.
                if any(i in z['actual_sides']for i in bad):continue
                if any(np.linalg.norm(n-q)<1e-7 for q in seen):continue
                seen.append(n);counts['accepted_up_to_one_causal_deletion']+=1
                new={**ns,47:n};actual=z['actual_sides'];boundary=B.ledger(np.array([new[i]for i in actual]),actual)
                row=dict(index=len(records),frame=list(ids),labels=list(ms),branch=sign,normal=n.tolist(),
                    additional_deleted_walls=bad,finite_contacts=[[i,r]for i,r in zip(remaining,rel)if isinstance(r,int)],boundary=boundary)
                records.append(row)
                if boundary['open_ends']<16:print('IMPROVEMENT',row,flush=True)
        print('frame',ids,dict(counts),'accepted',len(records),flush=True)
    return dict(status='BOUNDED_NUMERICAL_RELAXED47_REPLACEMENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in [S.PARENT,Path(S.__file__),Path(B.__file__)]},
        frames=[list(f)for f in frames],alphabet=ALPHABET,counts=dict(counts),records=records,
        blocker_sets=[dict(walls=list(k),count=v)for k,v in blockers.most_common()],
        scope='Any three contacts among7/12/30/52, bounded labels, actual finitecorner necessary.47removed; at mostone additional causally projected wall deletion explicitly recorded. Full S6 boundary scored numerically.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=compute()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],len(out['records']),flush=True)
    print(out['blocker_sets'][:12],flush=True)

if __name__=='__main__':main()
