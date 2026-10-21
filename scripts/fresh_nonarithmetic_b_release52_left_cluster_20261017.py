#!/usr/bin/env python3
"""Release52: close the entire left-cap face by coordinated adjacent caps."""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_b_release52_patch_20261017 as P
ROOT=P.ROOT;B=P.B;J=P.J;M=B.METRIC
TAG='fresh_nonarithmetic_b_release52_left_cluster_20261017'
ALPHABET=[2,3,4,5,6,8,10,12]

def candidates(ids,frame,contacts_old,L,data):
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=52];old=np.array([ns[i]for i in retained]);out=[];counts=Counter()
    gram=frame@M@frame.T;gi=np.linalg.inv(gram);k=np.linalg.svd(frame@M)[2][-1];kq=k@M@k
    for pi,pr in enumerate(profiles):
        off=np.array(pr['offsets']);q=pr['square'];oo=np.array([off[labels.index(i)]if i>=0 else 0. for i in ids])
        for ms in product(ALPHABET,repeat=3):
            cc=-np.array([B.theta(m)for m in ms])-oo
            # Both specified old-end duties must be actual finite corners.
            if any(np.linalg.eigvalsh(np.block([[gram[np.ix_(ab,ab)],cc[list(ab),None]],[cc[None,list(ab)],np.array([[q]])]]))[0]<1e-8 for ab in[(0,1),(1,2)]):continue
            p=cc@gi@frame;rem=(q-p@M@p)/kq
            if rem<=1e-9:continue
            for sign in[-1,1]:
                y=p+sign*math.sqrt(rem)*k;counts['branches']+=1
                if y@M@h+pr['center_offset']>=-1e-8:continue
                products=old@M@y+off[[labels.index(i)for i in retained]]
                for i,m in zip(ids,ms):
                    if i>=0:products[retained.index(i)]=-B.theta(m)
                rel=[B.relation(t)for t in products]
                if any(t is None or t=='parallel_or_unresolved'for t in rel):continue
                out.append(dict(profile=pi,labels=list(ms),branch=sign,projection=y.tolist(),finite_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)]))
    return out,counts

def compute():
    data=P.setup();z,labels,ns,h,profiles=data;actual=[i for i in z['actual_sides']if i!=52]
    seeds=json.loads((ROOT/'results'/(P.TAG+'.json')).read_text())['candidates']
    seeds=[r for r in seeds if r['profile']==0 and r['frame']==[7,30,47]and r['labels']in[[2,2,2],[2,2,3],[2,2,4]]]
    assert len(seeds)==3
    records=[];batches=[]
    for s in seeds:
        L=np.array(s['projection'])
        left,cl=candidates([12,47,-1],np.array([ns[12],ns[47],L]),None,L,data)
        # Reorder C's frame so old30 is the shared endpoint side.
        right,cr=candidates([64,30,-1],np.array([ns[64],ns[30],L]),None,L,data)
        print('SEED',s['labels'],'left',len(left),'right',len(right),flush=True)
        relations=Counter();pairs=[]
        for ai,a in enumerate(left):
            av=np.array(a['projection']);ap=profiles[a['profile']];ac=np.array(ap['contacts'])
            for bi,b in enumerate(right):
                bv=np.array(b['projection']);bp=profiles[b['profile']];bc=np.array(bp['contacts'])
                product=float(av@M@bv+ac@J@bc);r=B.relation(product);relations[str(r)]+=1
                if not isinstance(r,int):continue
                corner=np.array([[2,-B.theta(a['labels'][2]),-B.theta(b['labels'][2])],[-B.theta(a['labels'][2]),ap['square'],product-ac@J@bc],[-B.theta(b['labels'][2]),product-ac@J@bc,bp['square']]])
                if np.linalg.eigvalsh(corner)[0]<1e-8:continue
                boundary=B.ledger(np.array([ns[i]for i in actual]+[L,av,bv]),actual+[66,67,68])
                row=dict(seed=s['labels'],left=a,right=b,mutual_label=r,projections=[L.tolist(),av.tolist(),bv.tolist()],boundary=boundary)
                pairs.append(row);batches.append(row)
        records.append(dict(seed=s,left=left,right=right,counts=[dict(cl),dict(cr)],pair_relations=dict(relations),finite_pairs=pairs))
    batches.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],-r['boundary']['finite']))
    return dict(status='BOUNDED_NUMERICAL_RELEASE52_LEFT_FACE_CLUSTER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(P.__file__),ROOT/'results'/(P.TAG+'.json'),P.C.PARENT]},
        deleted_wall=52,alphabet=ALPHABET,records=records,batches=batches,
        scope='Three complete left-cap options terminate7/30 and7/47. Two further caps have all seven A2 profiles, terminate12/47 and64/30 respectively, and meet the left cap at finite endpoints. Their finite mutual contact and PD triple would close the entire left-cap pentagon. Labels sampled as listed; all native old gates and complete numerical S6 boundary tested. No all-label or continuous exclusion.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(compute()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['batches']),flush=True)
    if out['batches']:print('BEST',out['batches'][0],flush=True)
if __name__=='__main__':main()
