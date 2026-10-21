#!/usr/bin/env python3
"""Independent quarter-D10 finite arithmetic and fresh rank56 convolutions."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from functools import lru_cache
import json
from math import prod
import signal
import check_a16_d16_spectator_extension_independent_20261018 as Z

ROOT=Z.ROOT;W=Z.W;C=Z.C;K=Z.K;OLD=Z.OLD
TAG='fresh_order10_d_quarter_ideal57_scout_20261018'
LOCAL='fresh_order10_d_opposite_fork_fractional_transfer_20261018'


def object_at(path,marker):
    with path.open('rb')as stream:
        buf=b''
        while True:
            b=stream.read(65536)
            if not b:raise ValueError('Absent envelope object')
            buf+=b;i=buf.find(marker)
            if i<0:buf=buf[-len(marker):];continue
            buf=buf[i+len(marker)-1:]
            while True:
                try:return json.JSONDecoder().raw_decode(buf.decode())[0]
                except json.JSONDecodeError:buf+=stream.read(65536);assert len(buf)<1024**2


def increment(fs):
    co=Counter(map(tuple,fs));return co['D',9]*(12*co['A',1]+24*sum(co['A',r]for r in range(2,9)))


@lru_cache(None)
def parts(n,least=8):
    if n==0:return((),)
    return tuple((r,)+s for r in range(least,n+1)for s in parts(n-r,r))


def profiles(doc):
    env={int(k):F(v)for k,v in object_at(ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json',b'"current_envelopes": {').items()}
    assert env[9]==F(71,16)and env[16]==F(227,34)
    env[16]=F(113,17);oldenv=dict(env);env[9]=F(22,5)
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in parts(56-lo))
    ii=doc['ideal'];assert len(domain)==2270 and domain==[tuple(r['ranks'])for r in ii['rows']]
    bad=[];zeros=[];checks=0
    for rs,row in zip(domain,ii['rows']):
        h=[1]
        for r in rs:
            out=[0]*(len(h)+r)
            for i,a in enumerate(h):
                for j in range(r+1):out[i+j]+=a;checks+=1
            h=out
        v=prod(r+1 for r in rs);assert sum(rs)==56 and len(h)==57 and h==h[::-1]and sum(h)==v
        g=sum(min(i,56-i)*a for i,a in enumerate(h));curv=F(5*g-56*v,2*v)
        margin=curv-sum(env[r]for r in rs);previous=curv-sum(oldenv[r]for r in rs)
        assert margin==F(row['normalized_margin'])and previous==F(row['old_margin'])
        assert margin-previous==F(3*rs.count(9),80)==F(row['normalized_saving'])
        if margin<0:bad.append(row)
        if not margin:zeros.append(row)
    assert len(bad)==18 and not zeros and bad==ii['negative_profiles']
    assert min(F(r['normalized_margin'])for r in ii['rows'])==-F(17719,255552)
    assert ii['worst_profiles']==[[10,10,10,11,15]]
    return checks,len(bad)


def finite(doc):
    ff=doc['finite'];rep=C.representative;extra={}
    for a in C.TYPES:
        for b in C.TYPES:
            need=max(0,C.actual_cross(a,b)+increment([a,b])-C.C.upper_cross(rep(a),rep(b)))
            if need:key=tuple(sorted((rep(a),rep(b))));extra[key]=max(extra.get(key,0),need)
    assert extra=={tuple(sorted((t,('B',9)))):24 for t in(('A',5),('E',6),('A',7),('A',8))}
    domain={(cs,k,key[9]):key for cs,u,key in C.fixed_domain()if cs[0]for k in range((57-u)//5+1)}
    assert len(domain)==1181==ff['classes']and{(tuple(r['counts']),r['k5'],r['a12'])for r in ff['rows']}==set(domain)
    checks=preds=ordinary_options=quarter_options=0
    for r in ff['rows']:
        cs=tuple(r['counts']);k=r['k5'];key=domain[cs,k,r['a12']];rem=57-sum(a*b for a,b in zip(range(9,18),cs))-5*k
        assert rem==r['remaining_rank']and list(key)==r['key']
        ws={j:C.item(j,key,k)+(24*cs[0]if j in(6,7,8)else 0)for j in C.ITEMS};tabs={}
        for t in r['tables']:
            no2=t['exclude2'];js=[j for j in C.ITEMS if j<=rem and(not no2 or j!=2)]
            assert list(map(tuple,t['types']))==[rep(('A'if j<=8 else'D',j))for j in js]and t['weights']==[ws[j]for j in js]and len(t['values'])==rem+1
            a,b=K.bellman(t);checks+=a;preds+=b;tabs[no2]=t
        assert set(tabs)=={False,True}
        fixed=C.fixed_components(cs,key)+[('A',5)]*k;fc=C.upper_value(fixed)+24*cs[0]*k
        opts=[tabs[False]['values'][rem]]
        for a1 in(1,2):
            og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k)
            if rem<a1 or not(og or ng)or(a1==2 and not og):continue
            for a2 in range((rem-a1)//2+1):
                v=tabs[True]['values'][rem-a1-2*a2]
                if v is not None:opts.append(a1*ws[1]+a2*ws[2]+v+(36*a2*(a2-1)if og else 0)+(64*a2*k if ng else 0))
        ordinary_options+=len(opts);order17=16*k*key[9]+(8*k*rem if cs[2]else 0)
        ordinary=fc+max(v for v in opts if v is not None)+order17;quarter=None
        if cs[1]and rem>=7:
            for a1 in range(1,rem-5):
                for a4 in range(1,(rem-a1-2)//4+1):
                    for a2 in range(1,(rem-a1-4*a4)//2+1):
                        v=tabs[True]['values'][rem-a1-4*a4-2*a2]
                        if v is None:continue
                        q=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(36*a2*(a2-1)if a1>=2 else 0)+64*a2*k+64*a2*a4
                        quarter=q if quarter is None else max(quarter,q);quarter_options+=1
        assert F(ordinary,64)==F(r['ordinary_upper'])
        assert (None if quarter is None else F(quarter,64))==(None if r['quarter_upper']is None else F(r['quarter_upper']))
        assert F(max(ordinary,quarter if quarter is not None else -1),64)==F(r['upper'])
    assert(checks,preds,ordinary_options,quarter_options)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['ordinary_context_options'],ff['quarter_context_options'])
    assert max(F(r['upper'])for r in ff['rows'])==F(729,2)
    assert F(ff['coupled_A16_later_overlap_upper'])==F(741,2)+F(3,16)==F(5931,16)
    assert F(ff['guarded_D17_later_overlap_extra'])==max(F(3*x*(26-9*x),16)for x in(1,2))==F(51,16)
    assert F(ff['guarded_D17_later_overlap_upper'])==365+F(51,16)<F(1483,4)
    assert F(ff['complete_new_positive_family_upper'])==F(5931,16)<F(1483,4)
    assert Z.value(ff['attaining_forest'])+increment(ff['attaining_forest'])==64*F(1483,4)
    return checks,preds,ordinary_options,quarter_options


def local(doc):
    old=W.load(LOCAL);weights={tuple(r['pair']):F(r['weight'])for r in old['weighted_pairs']};changed={(6,8),(6,9),(7,8),(7,9)}
    assert all(weights[p]==F(27,32)for p in changed)and sum(F(3,4)if p in changed else w for p,w in weights.items())==44
    assert F(C.actual_base(('A',17)),64)==F(old['A17_new_cost'])-F(5,32)+F(1,2)
    for r in old['physical']['immediate']:assert F(3,4)+(F(1,2)if r['target']=='B'else F(3,4))>=1
    for r in old['whole_finite_coupons']:assert F(3,4)+(F(1,4)if '/wholeD9' in r['role']else 1)>=1
    assert old['triangle']['allowed_a']==[0,7]and old['triangle']['allowed_b']==[5,9]
    assert doc['local']['extra_A17_cost']=='0'


def main():
    doc=W.load(TAG);hashes=0
    for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    local(doc);a=profiles(doc);b=finite(doc);OLD.U.memory_guard()
    print('PASS_INDEPENDENT_ORDER10_D_QUARTER_IDEAL57',hashes,*a,*b,'1483/4','18',flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
