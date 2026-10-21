#!/usr/bin/env python3
"""Bounded third-side search to close the whole47 facet after two endcaps."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release54_region47_20261017 as S
B,R,ROOT=S.B,S.R,S.ROOT
TAG='fresh_nonarithmetic_release54_region47_bridge_20261017'
BANK=ROOT/'results/fresh_nonarithmetic_release54_region47_pairs_20261017.json'

def compute():
    z,labels,ns,h=S.setup();x=json.loads(BANK.read_text())
    groups=[[i for i in g if max(x['records'][i]['labels'])<=8]for g in x['groups']]
    counts=Counter();records=[];best=[]
    old=np.array([ns[i]for i in labels]);actual=z['actual_sides']
    for ia,ib in product(*groups):
        a,b=[np.array(x['records'][i]['normal'])for i in [ia,ib]]
        frame=np.array([ns[47],a,b]);nn=np.vstack([old,a,b]);ls=labels+[66,67]
        for ma,mb in product(range(2,9),repeat=2):
            for sign,n in R.roots(frame,[0.,-B.theta(ma),-B.theta(mb)]):
                counts['branches']+=1
                status,rel=R.gate(n,nn,ls,h,{47:2,66:ma,67:mb})
                counts[status]+=1
                if status!='accepted':continue
                normals=np.array([ns[i]for i in actual]+[a,b,n]);sides=actual+[66,67,68]
                boundary=B.ledger(normals,sides)
                row=dict(index=len(records),endpoint_indices=[ia,ib],endpoint_labels=[x['records'][i]['labels']for i in [ia,ib]],
                    bridge_labels=[2,ma,mb],branch=sign,new_normals=[a.tolist(),b.tolist(),n.tolist()],
                    finite_contacts=rel,boundary=boundary)
                records.append(row)
                if boundary['open_ends']<16:
                    best.append(row['index']);print('IMPROVEMENT',row['index'],row['endpoint_labels'],[ma,mb],boundary['open_ends'],flush=True)
        if (ia==groups[0][0]):print('first A / B',ib,dict(counts),'accepted',len(records),flush=True)
    return dict(status='BOUNDED_NUMERICAL_THREE_CAP_47_FACE_CLOSURE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in [S.PARENT,BANK,Path(S.__file__),Path(R.__file__)]},
        endpoint_label_bound=8,bridge_labels=list(range(2,9)),endpoint_groups=groups,
        counts=dict(counts),records=records,improving_indices=best,
        scope='Third norm2 wall right47 and finite to the two selected end caps; all64 old plus2new native gates. Only bounded numerical discovery; entire S6 boundary scored.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=compute()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'accepted',len(out['records']),'improvements',len(out['improving_indices']),flush=True)

if __name__=='__main__':main()
