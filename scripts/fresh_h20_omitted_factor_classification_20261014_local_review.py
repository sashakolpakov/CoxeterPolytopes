#!/usr/bin/env python3
"""Independent Fraction/Q(sqrt(6)) replay of accepted local directions.

The source-level completeness proof and a primary saved replay are recorded
in the review note. This file replays every accepted record, not the whole
support/branch enumeration, and imports no producer arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_omitted_factor_classification_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
LOCAL=f'results/{TAG}.json'
Z=(F(0),F(0));O=(F(1),F(0))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+6*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,c):return(x[0]*c,x[1]*c)
def sub(x,y):return add(x,scale(y,-1))
def sgn(x):return int(x>0)-int(x<0)
def sign(x):
    a,b=x
    if not b:return sgn(a)
    if not a:return sgn(b)
    if sgn(a)==sgn(b):return sgn(a)
    return sgn(a)*sgn(a*a-6*b*b)
def dec(v):
    a=list(map(F,v));assert len(a)==8 and all(not a[i]for i in[1,2,4,5,6,7])
    return(a[0],a[3])
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def linear(a,b):
    z=Z
    for x,y in zip(a,b,strict=True):z=add(z,scale(y,x))
    return z
def kdot(a,b):
    z=Z
    for x,y in zip(a,b,strict=True):z=add(z,mul(x,y))
    return z
def pd(A):
    A=[[F(x)for x in row]for row in A];piv=[]
    while A:
        p=A[0][0];assert p>0;piv.append(p)
        A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,len(A))]for i in range(1,len(A))]
    return piv
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    saved=json.loads((ROOT/LOCAL).read_text())
    for p,h in saved['sources'].items():assert sha(p)==h,p
    raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q=[[-int(x)for x in r]for r in raw['basis_gram']]
    W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
    C=[[dot(q,w)for q in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
    walk=json.loads((ROOT/WALK).read_text())['corner_orbits'];byid={r['id']:r for r in walk}
    rows=[];groups={};support_counts=Counter();angle_counts=Counter();irrational=0
    for ix,r in enumerate(saved['local_directions']):
        ids=r['retained_walls'];a=list(map(dec,r['retained_pairings']));k=list(map(dec,r['positive_kernel']))
        e=list(map(dec,r['unit_normal_constant']));D=list(map(F,r['null_shift']));L=dec(r['parameter_lower_magnitude'])
        fac=next(f for f in byid[r['cusp']]['affine_components']if f['walls']==r['factor'])
        marks=dict(zip(fac['walls'],map(F,fac['positive_kernel']),strict=True));omit=r['omitted_wall']
        assert min(marks.values())>0 and F(fac['ray_multiple'])>0
        assert ids==[i for i in fac['walls']if i!=omit]
        B=[[G[i][j]for j in ids]for i in ids];piv=pd(B)
        assert min(map(sign,k))>0 and [linear(row,k)for row in B]==a and kdot(a,k)==O
        assert e==[linear([W[i][d]for i in ids],k)for d in range(22)]
        Qe=[linear(q,e)for q in Q]
        assert kdot(e,Qe)==O and linear(D,Qe)==Z
        assert D==[sum(marks[i]*W[i][d]for i in fac['walls'])for d in range(22)]
        assert D==[F(fac['ray_multiple'])*F(x)for x in byid[r['cusp']]['primitive_ray']]
        assert dot(D,[dot(q,D)for q in Q])==0
        assert [linear(C[i],e)for i in ids]==a
        support=[i for i,x in zip(ids,a,strict=True)if x!=Z];assert support==r['support']
        support_counts[len(support)]+=1;irrational+=int(any(x[1]for x in e))
        H=[[sub((F(B[i][j]),F(0)),mul(a[i],a[j]))for j in range(len(ids))]for i in range(len(ids))]
        assert all(sign(H[i][i])>0 and kdot(H[i],k)==Z for i in range(len(ids)))
        edges=[]
        for i in range(len(ids)):
            for j in range(i):
                p=H[i][j];nn=mul(H[i][i],H[j][j]);assert sign(p)<=0
                values=[c for c in[F(0),F(1,4),F(1,2)]if mul(p,p)==scale(nn,c)]
                assert len(values)==1;angle_counts[str(values[0])]+=1
                if values[0]:edges.append([ids[j],ids[i],str(values[0])])
        assert sorted(edges)==sorted(r['projected_edges_squared_cosines'])
        minima=[scale(x,1/marks[i])for i,x in zip(ids,k,strict=True)]
        assert sign(L)>0 and all(sign(sub(x,L))>=0 for x in minima)and L in minima
        assert minima[ids.index(r['lower_bound_witness_wall'])]==L
        coeff=[sub(x,scale(L,marks[i]))for i,x in zip(ids,k,strict=True)]
        assert all(sign(c)>=0 for c in coeff)
        e1=[sub(x,scale(L,d))for x,d in zip(e,D,strict=True)]
        assert e1==[sub(linear([W[i][d]for i in ids],coeff),scale(L,marks[omit]*W[omit][d]))for d in range(22)]
        pairs=[linear(C[i],e)for i in fac['walls']];o=fac['walls'].index(omit)
        assert sign(pairs[o])<0
        # Independently recover the first representative of this spatial
        # plane using proportional covectors on the whole affine factor.
        key=(r['cusp'],tuple(fac['walls']));match=None
        for oldid,oldpairs in groups.get(key,[]):
            if all(mul(x,oldpairs[o])==mul(y,pairs[o])for x,y in zip(pairs,oldpairs,strict=True)):
                match=oldid;break
        assert match==r['same_spatial_plane_as_row']
        if match is None:groups.setdefault(key,[]).append((ix,pairs))
        rows.append({'row':ix,'retained_LDL':piv,'positive_kernel':k,'unit_normal':e,
            'shift_lower_magnitude':L,'lower_endpoint_coefficients':coeff,'same_plane_as_row':match})
    assert len(rows)==175 and sum(map(len,groups.values()))==119 and irrational==48
    sources=[RAW,WALK,LOCAL,f'scripts/{TAG}.py',f'scripts/{TAG}_local_review.py']
    return enc({'status':'PASS_INDEPENDENT_ACCEPTED_LOCAL_DIRECTION_AND_SHIFT_REPLAY',
        'sources':{p:sha(p)for p in sources},'accepted_records':len(rows),'spatial_planes':119,
        'irrational_records':irrational,'support_sizes':dict(support_counts),'pair_counts':dict(angle_counts),
        'records':rows,'scope':'Separate Fraction/Q(sqrt(6)) replay of all accepted records, projected Gram identities, positivity, duplicate-plane relations and shift interval minima. Completeness of support/branch enumeration is source/proof reviewed and checked by the primary saved replay, not re-enumerated here.'})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();z=run()
    p=ROOT/f'results/{TAG}_local_review.json'
    if a.verify:assert z==json.loads(p.read_text())
    else:p.write_text(json.dumps(z,indent=2)+'\n')
    print(z['status']);print(z['accepted_records'],z['spatial_planes'],z['irrational_records'])
