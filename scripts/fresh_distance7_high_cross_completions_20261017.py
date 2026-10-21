#!/usr/bin/env python3
"""All labels and simultaneous subsets of the48 residual cross positions.

This closes extra cross contacts, not arbitrary supporting-root extensions.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_distance7_high_cross_domain_20261017 as C

H=C.H;ROOT=C.ROOT
TAG='fresh_distance7_high_cross_completions_20261017'


class Q2:
    __slots__=('a','b')
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __bool__(self):return bool(self.a or self.b)
    def __add__(self,y):
        if not isinstance(y,Q2):y=Q2(y)
        return Q2(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Q2(-self.a,-self.b)
    def __sub__(self,y):return self+-y
    def __mul__(self,y):
        if not isinstance(y,Q2):y=Q2(y)
        return Q2(self.a*y.a+2*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def __truediv__(self,y):
        if not isinstance(y,Q2):y=Q2(y)
        d=y.a*y.a-2*y.b*y.b;assert d
        return self*Q2(y.a/d,-y.b/d)
    def sign(self):
        s=lambda x:int(x>0)-int(x<0)
        if not self.b:return s(self.a)
        if not self.a:return s(self.b)
        if s(self.a)==s(self.b):return s(self.a)
        n=self.a*self.a-2*self.b*self.b;assert n
        return s(self.a)*s(n)
    def out(self):return [str(self.a),str(self.b)]


def surd_inertia(native):
    a=[[Q2(0,-1)if i!=j and x==-2 else Q2(x)for j,x in enumerate(row)]for i,row in enumerate(native)]
    pos=neg=zero=0;steps=[]
    while a:
        k=next((i for i in range(len(a))if a[i][i]),None)
        if k is not None:
            ids=[k]+[i for i in range(len(a))if i!=k];a=H.extract(a,ids);p=a[0][0];s=p.sign();assert s
            pos+=s>0;neg+=s<0;steps.append(dict(permutation=ids,pivot=p.out(),sign=s))
            a=[[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,len(a))]for i in range(1,len(a))]
        else:
            pair=next(((i,j)for i in range(len(a))for j in range(i+1,len(a))if a[i][j]),None)
            if pair is None:zero+=len(a);break
            i,j=pair;ids=[i,j]+[k for k in range(len(a))if k not in pair];a=H.extract(a,ids);p=a[0][1]
            pos+=1;neg+=1;steps.append(dict(permutation=ids,hyperbolic_pivot=p.out(),signature=[1,1]))
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/p for j in range(2,len(a))]for i in range(2,len(a))]
    return dict(method='exact_Q_sqrt2_congruence',signature=[pos,neg,zero],steps=steps)


@lru_cache(None)
def native_inertia(native):
    try:g,sc=H.native_to_rational(native)
    except AssertionError:return surd_inertia(native)
    return dict(method='positive_diagonal_rational_congruence',gram=H.encode(g),gauge_squared=sc,signature=H.inertia(g))


def bound_above_four(native,u,v,i,j):
    # This boundary occurs only in the all-simple base. For label>=5 the
    # norm-two contact is >=(1+sqrt5)/2 >3/2. A rational -3/2 comparison
    # with a nonnegative negative-square corner witness excludes all of it.
    assert {x for row in native for x in row}<={F(-1),F(0),F(2)}
    g=[r.copy()for r in native];g[i][j]=g[j][i]=-F(3,2)
    white=[k for k in range(len(g))if k not in u+v]
    for p,q in product(u,v):
        ids=sorted(white+[p,q]);o=C.D.strict_obstruction(H.extract(g,ids))
        if o:return dict(comparison_contact='3/2',corner=[p,q],nodes=ids,**o)
    raise AssertionError('missing all-label gate')


def run():
    cross=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    for rel,h in cross['sources'].items():assert sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
    ds={d['index']:d for d in json.loads((ROOT/f'results/{H.CENSUS}.json').read_text())['diagrams']}
    gg={i:H.bond_graph(d)for i,d in ds.items()};gates=[];rows=[];counts=Counter()
    for c in cross['cases']:
        if not c.get('remaining_positions'):continue
        b=c['base_context'];di,u=b['left'];dj,v0=b['right'];native,mp=H.unite(gg[di],gg[dj],b['match']);v=[mp[i]for i in v0]
        options=[]
        for i,j,dist in c['remaining_positions']:
            h=[r.copy()for r in native];h[i][j]=h[j][i]=F(-2)
            ob=H.edge_bad_native(h,u+v)
            if ob:
                allow=[3];gate=dict(context=c['index'],edge=[i,j],allowed_labels=allow,first_excluded_label=4,edge_witness=ob)
            else:
                ob=C.corner_bad(h,u,v)
                if ob:
                    allow=[3];gate=dict(context=c['index'],edge=[i,j],allowed_labels=allow,first_excluded_label=4,corner_witness=ob)
                else:
                    allow=[3,4];gate=dict(context=c['index'],edge=[i,j],allowed_labels=allow,first_excluded_label=5,
                        larger_label_corner_witness=bound_above_four(native,u,v,i,j))
            options.append((i,j,allow));gates.append(gate);counts['positions']+=1;counts['positions_allowing4']+=4 in allow
        for size in range(1,len(options)+1):
            for chosen in combinations(options,size):
                for labels in product(*(z[2]for z in chosen)):
                    h=[r.copy()for r in native]
                    for (i,j,_),label in zip(chosen,labels):h[i][j]=h[j][i]=F(-(label-2))
                    row=dict(base_context=c['index'],edges=[[i,j,m]for(i,j,_),m in zip(chosen,labels)])
                    dist=min(H.graph_dist(h,u)[i]for i in v)
                    if dist<=7:row.update(status='marked_short7',distance=dist)
                    else:
                        white=[i for i in range(len(h))if i not in u+v];checks=[];reason=None
                        for p in u+v:
                            ids=sorted(white+[p]);z=native_inertia(tuple(map(tuple,H.extract(h,ids))));checks.append(dict(kind='edge',marks=[p],nodes=ids,proof=z))
                            if z['signature']!=[len(ids),0,0]:reason='actual_edge_not_PD';break
                        if reason is None:
                            for p,q in product(u,v):
                                ids=sorted(white+[p,q]);z=native_inertia(tuple(map(tuple,H.extract(h,ids))));checks.append(dict(kind='corner',marks=[p,q],nodes=ids,proof=z))
                                if z['signature'][1]:reason='actual_corner_not_PSD';break
                        if reason is None:
                            z=native_inertia(tuple(map(tuple,h)));row['full_core_inertia']=z
                            assert z['signature'][1]>=1
                            if z['signature'][1]>1:reason='multiple_negative_directions'
                        row.update(status=reason or'survives_native_core_tests',minimum_marked_distance=dist,actual_principal_checks=checks)
                        if reason is None:row.update(native_graph=H.encode(h),left_marks=u,right_marks=v,right_copy_map=mp)
                    rows.append(row);counts[row['status']]+=1
    assert len(gates)==48 and counts['positions_allowing4']==2 and len(rows)==54
    paths=[Path(__file__),Path(C.__file__),Path(H.__file__),Path(C.D.__file__),Path(C.D.R.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths.extend(ROOT/f'{folder}/{tag}.{ext}'for tag in(C.TAG,H.CENSUS)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths=list(dict.fromkeys(paths))
    out=dict(status='PASS_DISTANCE7_ALL_LABEL_EXTRA_CROSS_COMPLETIONS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),all_label_gates=gates,nonempty_cross_assignments=rows,
        unchanged_zero_cross_contexts=6691,
        scope='Complete native labels and nonempty simultaneous subsets on all48 residual positions in44 high-overlap contexts. Higher labels including parallel/dotted contacts excluded by monotone edge/corner witnesses. Every remaining labelled Gram checked exactly, including Q(sqrt2) when no rational gauge exists. Zero-extra-edge base cases and arbitrary supporting-root extensions retain their separate obligations. No universal finite7 payment or dimension theorem.')
    return json.loads(json.dumps(out))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
