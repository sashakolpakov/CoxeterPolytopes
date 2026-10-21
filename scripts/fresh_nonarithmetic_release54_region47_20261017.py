#!/usr/bin/env python3
"""Bounded numerical replacement47 screen on the frozen64 chamber.

All native wall gates are evaluated, including causal S6 projections. A
suggested deletion is explicitly recorded and its complete numerical S6
boundary is recomputed. This is discovery, not an exact acceptance proof.
"""
import argparse
import json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
import fresh_nonarithmetic_release34_cavity_20261017 as B
import fresh_nonarithmetic_release34_replace12_47_20261017 as R

ROOT=B.ROOT
TAG='fresh_nonarithmetic_release54_region47_20261017'
PARENT=ROOT/'results/fresh_nonarithmetic_release54_replace12_boundary_20261017.json'
GEOMETRY=ROOT/'results/fresh_nonarithmetic_release54_replace12_20261017.json'
ALPHABET=list(range(2,65))+[96,128,192,256,512]


def setup():
    z=json.loads(PARENT.read_text())['after64']
    g=np.array([[B.mid(x)for x in row]for row in z['residual_gram']])
    g=(g+g.T)/2
    vals,vecs=np.linalg.eigh(g)
    assert vals[0]<0<vals[1]
    transform=vecs@np.diag(1/np.sqrt(abs(vals)))
    assert np.max(abs(transform.T@g@transform-B.METRIC))<1e-12
    labels=[r['wall']for r in z['projections']]
    cov=np.array([[B.mid(x)for x in r['covector']]for r in z['projections']])
    normals=cov@transform@B.METRIC
    h=np.linalg.solve(transform,np.array([B.mid(x)for x in z['projected_common_center']]))
    h/=np.sqrt(-h@B.METRIC@h)
    return z,labels,dict(zip(labels,normals)),h


def compute():
    z,labels,ns,h=setup()
    actual=z['actual_sides']
    old=json.loads(GEOMETRY.read_text())
    finite=old['retained_old_finite_supports']+[r['support']for r in old['new_finite_vertices']]
    assert len(finite)==197
    baseline=B.ledger(np.array([ns[i]for i in actual]),actual)
    assert baseline['open_ends']==16
    counts=Counter();badsets=Counter();records=[];improvements=[]
    remaining=[i for i in labels if i!=47]
    nn=np.array([ns[i]for i in remaining])
    for anchor in [12,52,8,44,65]:
        frame=np.array([ns[i]for i in [7,30,anchor]])
        for m in ALPHABET:
            for branch,n in R.roots(frame,[0.,0.,-B.theta(m)]):
                counts['branches']+=1
                products=nn@B.METRIC@n
                for i,x in [(7,0.),(30,0.),(anchor,-B.theta(m))]:
                    products[remaining.index(i)]=x
                blockers=[];finite_contacts=[]
                for i,x in zip(remaining,products):
                    rel=B.relation(x)
                    if rel is None or rel=='parallel_or_unresolved':
                        blockers.append(dict(wall=i,product=float(x),kind='obtuse'if x>0 else'gap_or_parallel'))
                    elif isinstance(rel,int):finite_contacts.append([i,rel])
                badsets[tuple(r['wall']for r in blockers)]+=1
                if n@B.METRIC@h>=-1e-8:
                    counts['wrong_fixed_center']+=1
                    continue
                counts['common_center_branches']+=1
                if len(blockers)>3:continue
                deleted=[r['wall']for r in blockers]
                if any(i in [7,30,anchor]for i in deleted):raise AssertionError('forced blocker')
                sides=[i for i in actual if i not in deleted]
                nmap={**ns,47:n}
                boundary=B.ledger(np.array([nmap[i]for i in sides]),sides)
                affected=[s for s in finite if 47 in s or any(i in s for i in deleted)]
                row=dict(index=len(records),anchor=anchor,label=m,branch=branch,
                    normal=n.tolist(),center_product=float(n@B.METRIC@h),
                    blockers=blockers,additional_deleted_walls=deleted,
                    finite_contacts=finite_contacts,boundary=boundary,
                    old_finite_supports_incident_to_changed_or_deleted=affected,
                    preserves_heptagon_side_labels=not any(i in deleted for i in [7,8,10,11,12,44,52,65]))
                records.append(row)
                if boundary['open_ends']<16:
                    improvements.append(row['index'])
                    print('IMPROVEMENT',row['index'],anchor,m,branch,deleted,boundary['open_ends'],boundary['finite'],flush=True)
    return dict(status='BOUNDED_NUMERICAL_REPLACEMENT47_WITH_EXPLICIT_DELETION_COSTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in [PARENT,GEOMETRY,Path(B.__file__),Path(R.__file__)]},
        alphabet=ALPHABET,anchors=[12,52,8,44,65],baseline=baseline,counts=dict(counts),
        blocker_sets=[dict(walls=list(k),count=v)for k,v in badsets.most_common()],
        records=records,improving_indices=improvements,
        scope='Two fixed right contacts7/30 and one bounded finite event; all64 native contacts tested. Up to3 extra deletions scored. Numerical only; no inference outside this recipe domain.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=compute();path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improving_indices']),flush=True)
    print('top blocker sets',out['blocker_sets'][:12],flush=True)


if __name__=='__main__':main()
