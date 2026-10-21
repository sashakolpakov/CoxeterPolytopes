#!/usr/bin/env python3
"""Complete incident contact domain; optional numerical quotient relaxation."""
import argparse,json,math
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_null_link_profiles_20261016'
BASE='results/fresh_h20_affine_a17_completion_20261016.json'
RAD=(1,2,3,5,6,10,15,30)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def sq(d):
    den=1<<256;a=math.isqrt(d*den*den)
    return F(a,den),F(a+int(a*a!=d*den*den),den)
def sign(poly):
    poly={d:q for d,q in poly.items()if q}
    if not poly:return 0
    lo=hi=F(0)
    for d,q in poly.items():
        a,b=sq(d);lo+=q*(a if q>0 else b);hi+=q*(b if q>0 else a)
    assert lo>0 or hi<0,poly
    return 1 if lo>0 else -1
def mul(a,b):
    out={}
    for d,x in a.items():
        for e,y in b.items():
            g=math.gcd(d,e);r=d*e//g**2;out[r]=out.get(r,F(0))+g*x*y
    return out
def run():
    a=json.loads((ROOT/BASE).read_text());paths=a['definition']['old_paths'];white=sum(paths,[])
    sizes=list(map(len,paths));assert sizes==[2,5,4,1,1]
    J=[[F(0)for _ in white]for _ in white];off=0
    for n in sizes:
        for i in range(n):
            for j in range(n):J[off+i][off+j]=F((min(i,j)+1)*(n-max(i,j)),n+1)
        off+=n
    assert sorted(J[i][i]for i in range(13))[:4]==[F(1,2),F(1,2),F(2,3),F(2,3)]
    values={3:{1:F(1)},4:{2:F(1)},5:{1:F(1,2),5:F(1,2)},6:{3:F(1)}}
    prods={(m,n):mul(values[m],values[n])for m in values for n in values}
    records=[dict(walls=[],labels=[],cost=['0']*8)];rejected=0;trials=0
    for count in (1,2,3):
        for ids in combinations(range(13),count):
            for labels in product(range(3,7),repeat=count):
                trials+=1;q={}
                for i,m in zip(ids,labels):
                    for j,n in zip(ids,labels):
                        for d,c in prods[m,n].items():q[d]=q.get(d,F(0))+J[i][j]*c
                cmp=dict(q);cmp[1]=cmp.get(1,F(0))-2
                if sign(cmp)>0:rejected+=1;continue
                records.append(dict(walls=[white[i]for i in ids],labels=list(labels),cost=[str(q.get(d,F(0)))for d in RAD]))
    assert all(len(r['walls'])<=3 for r in records)
    return dict(status='PASS_COMPLETE_INCIDENT_CONTACT_DOMAIN_WITH_TWO_CONTINUOUS_TAILS',
        sources={BASE:sha(BASE),f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py')},
        white=white,components=paths,inverse_Gram=[[str(x)for x in row]for row in J],
        squarefree_basis=list(RAD),finite_trials=trials,rejected_over_budget=rejected,
        finite_profiles=len(records),profiles=records,
        continuous_tails=[dict(wall=i,other_incident_contacts='right',magnitude_interval=['sqrt3','2'],
                               overlaps_finite_label=6,includes_parallel=True)for i in (17,31)],
        proof=dict(cost_bound='d^t J d<=2 by Euclidean quotient projection',
                   nonnegative_inverse=True,minimum_diagonal='1/2',
                   four_smallest_diagonal_sum='7/3>2, so at most3 nonzero contacts',
                   nonisolated_diagonal='>=2/3, hence d<=sqrt3 and label<=6',
                   multiple_contact_isolated_bound='d^2/2+1/2<=2, hence label<=6',
                   dotted_incident_exclusion='d>2 implies diagonal cost>2'),
        scope='Complete all-label incident contact domain for norm2 acute walls through nu. Finite profiles plus two continuous overapproximating tails; not yet a Coxeter completion or E1 exclusion theorem.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['finite_profiles'],d['finite_trials'],flush=True)
