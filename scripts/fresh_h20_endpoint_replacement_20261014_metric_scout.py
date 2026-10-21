#!/usr/bin/env python3
"""Numerical full-orbit test of the pi/5-tuned changed endpoint hierarchy."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse,json
import numpy as np
from fresh_h20_balanced_20261014 import data

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_endpoint_replacement_20261014_metric_scout'

def classify(cs):
    labels=np.full(len(cs),-1,dtype=int)
    labels[cs>1+1e-8]=0
    labels[np.abs(cs)<1e-8]=2
    labels[np.abs(cs-1)<1e-8]=1
    sel=(cs>1e-8)&(cs<1-1e-8)
    m=np.pi/np.arccos(cs[sel]);rounded=np.rint(m).astype(int)
    good=np.abs(m-rounded)<1e-6
    labels[np.where(sel)[0][good]]=rounded[good]
    return labels

def run():
    d=data();Z=np.array([v['point'] for v in d['vertices']]);V=d['vectors']
    idx=lambda a:a[0]*56+a[1]*8+a[2]
    c=(1+np.sqrt(5))/4
    r=float(max(np.roots([2-4*c,49*c-147,2058])))
    assert 16<r<49/3
    C=[np.zeros(19)]+list(Z);Q=[4.]+[r]*392
    rows=[dict(kind='central')]+[dict(kind='vertex',vertices=[i]) for i in range(392)]
    for kind in ('edge','triangle'):
        for f,n in enumerate((7,7,8)):
            for a in product(range(7),range(7),range(8)):
                b=list(a);b[f]=(b[f]+1)%n
                ii,jj=idx(a),idx(b)
                if kind=='edge':
                    C.append((Z[ii]+Z[jj])/2);Q.append(n*(n-1)/2-r)
                    rows.append(dict(kind=kind,factor=f,vertices=[ii,jj]))
                else:
                    cc=list(a);cc[f]=(cc[f]+2)%n;kk=idx(cc);t=r/(n*n)
                    C.append(t*Z[ii]+(1-2*t)*Z[jj]+t*Z[kk]);Q.append(4*r*r/(n*n)-r)
                    rows.append(dict(kind=kind,factor=f,vertices=[ii,jj,kk]))
    C=np.array(C);Q=np.array(Q);N=np.sum(C*C,axis=1)
    assert len(C)==2745 and min(Q)>0
    counts=Counter();bad=[]
    for i in range(len(C)):
        cos=(N[i]+N[:i]-2*C[:i]@C[i]-Q[i]-Q[:i])/(2*np.sqrt(Q[i]*Q[:i]))
        labels=classify(cos);counts.update(labels.tolist())
        bad += [dict(i=i,j=int(j),cosine=float(cos[j])) for j in np.where(labels<0)[0]]
    seed_bad=[]
    for i,row in enumerate((1-C@V.T)/np.sqrt(Q[:,None])):
        seed_bad += [dict(i=i,j=int(j),cosine=float(row[j])) for j in np.where(classify(row)<0)[0]]
    oldpath='results/fresh_h20_balanced_integer_exact_20261014.json'
    old=json.loads((ROOT/oldpath).read_text())
    W=np.array([[float(Fraction(x)) for x in row] for row in old['walls'][415:]])
    OC=W[:,:19]@d['cholesky'];ON=np.sum(OC*OC,axis=1);OQ=ON-2*W[:,19]
    old_cross_counts=Counter();old_bad=[]
    for i in range(len(C)):
        cos=(N[i]+ON-2*OC@C[i]-Q[i]-OQ)/(2*np.sqrt(Q[i]*OQ))
        labels=classify(cos);old_cross_counts.update(labels.tolist())
        old_bad += [dict(i=i,old_wall=int(j)+415,cosine=float(cos[j])) for j in np.where(labels<0)[0]]
    probes=[]
    for lengths in [(3,1,1),(4,1,1),(1,1,3),(1,1,4),(2,2,1),(2,1,2),(3,3,3),(4,4,4)]:
        ids=[idx(a) for a in product(*[range(n) for n in lengths])]
        x=.99*np.mean(Z[ids],axis=0)
        power=np.sum((C-x)**2,axis=1)-Q
        oldpower=np.sum((OC-x)**2,axis=1)-OQ
        radial=float(x@x-4-np.max((Z@x-4)**2)/192)
        probes.append(dict(face_lengths=list(lengths),point=x.tolist(),new_minimum_power=float(min(power)),
                           old_middle_minimum_power=float(min(oldpower)),
                           all_real_radial_lower_bound=radial,closest_new=rows[int(np.argmin(power))]))
    assert not bad and not seed_bad and not old_bad
    paths=[f'scripts/{TAG}.py','scripts/fresh_h20_balanced_20261014.py',oldpath]
    return dict(status='NUMERICAL_COHERENT_METRIC_HIERARCHY_WITH_GLOBAL_GAPS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        endpoint_radius_squared=r,caps=2745,walls_with_seed=2767,
        numerical_pair_count=sum(counts.values()),numerical_pair_labels={str(k):v for k,v in counts.items()},
        numerical_old_cross_count=sum(old_cross_counts.values()),numerical_old_cross_labels={str(k):v for k,v in old_cross_counts.items()},
        seed_failures=seed_bad,pair_failures=bad,old_cross_failures=old_bad,
        label_encoding='0=dotted, 1=parallel, integers>=2=pi/integer; -1=unrecognized',
        recipes=rows,coverage_probes=probes,
        scope='Floating proposal only. The exact companion proves a common escape for the entire parameter interval, so no finite-volume or new exact chamber checkpoint is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    f=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(f.read_text())==out
    else:f.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['numerical_pair_labels'],out['numerical_old_cross_labels'])
