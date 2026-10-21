#!/usr/bin/env python3
"""Conservative exact-interval compatibility graph of 3040 archived cap IDs.

Includes all finite Coxeter labels: labels2..120 are enclosed explicitly,
and the whole interval from -2 to the label120 upper bound is retained.
Unknown pairs remain edges. Verification uses independent exact field signs
for coefficient enclosures and general four-product interval arithmetic.
"""
import argparse
from bisect import bisect_left
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_3034_safe_graph_20261013.json'
AUDIT=ROOT/'results/fresh_h9_split_leaf_3034_safe_graph_20261013_verify.json'
SCALE=10**24
PAIR_SCALE=SCALE*SCALE
FINE=10**80
TRIG_SCALE=10**70
ZERO=(F(0),F(0))


def ceildiv(a,b):return-((-a)//b)
def parse(x):return(F(x['a']),F(x['b']))
def digest(x):return hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()


class Enclosures:
    def __init__(self):
        q=isqrt(5*FINE*FINE)
        self.phi=((FINE+q)//2,(FINE+q+2)//2)
        assert self.phi[0]**2-self.phi[0]*FINE-FINE**2<0
        assert self.phi[1]**2-self.phi[1]*FINE-FINE**2>0
        self.kcache={};self.sqrtcache={};self.ecache={}
    def golden(self,x):
        if x in self.kcache:return self.kcache[x]
        a,b=x;p,q=self.phi if b>=0 else self.phi[::-1]
        lo=a*FINE+b*p;hi=a*FINE+b*q
        out=(lo.numerator//lo.denominator,ceildiv(hi.numerator,hi.denominator))
        self.kcache[x]=out;return out
    def sqrt(self,d):
        if d in self.sqrtcache:return self.sqrtcache[d]
        if d==ZERO:out=(0,0)
        else:
            lo,hi=self.golden(d);assert lo>=0
            out=(isqrt(lo*FINE),isqrt(hi*FINE)+1)
            assert out[0]**2<=lo*FINE and out[1]**2>=hi*FINE
        self.sqrtcache[d]=out;return out
    def element(self,x,d):
        key=(x,d)
        if key in self.ecache:return self.ecache[key]
        a,b=x;al,ah=self.golden(a);bl,bh=self.golden(b);dl,dh=self.sqrt(d)
        products=(bl*dl,bl*dh,bh*dl,bh*dh)
        lo=al+min(products)//FINE;hi=ah+ceildiv(max(products),FINE)
        step=FINE//SCALE
        out=(lo//step,ceildiv(hi,step));assert out[0]<=out[1]
        self.ecache[key]=out;return out


def atan_bounds(den,terms):
    x=F(1,den);square=x*x;term=x;s=F(0)
    for k in range(terms):
        s+=term/(2*k+1)*(1 if k%2==0 else-1)
        term*=square
    nxt=term/(2*terms+1)*(1 if terms%2==0 else-1)
    return min(s,s+nxt),max(s,s+nxt)


def cos_interval(xlo,xhi,terms=48):
    """Fixed-integer interval Taylor series; alternating tail included."""
    assert 0<=xlo<=xhi<2*TRIG_SCALE and terms>=2
    lo=hi=TRIG_SCALE;sl=sh=TRIG_SCALE
    for k in range(1,terms+1):
        denominator=TRIG_SCALE**2*(2*k-1)*(2*k)
        lo=lo*xlo*xlo//denominator
        hi=ceildiv(hi*xhi*xhi,denominator)
        if k==terms:
            # The next term bounds the entire remaining alternating tail.
            # Subsequent absolute terms decrease since x<2 and k>=2.
            if k%2==0:sh+=hi
            else:sl-=hi
        elif k%2==0:sl+=lo;sh+=hi
        else:sl-=hi;sh-=lo
    return sl,sh


def angle_enclosures():
    al,ah=atan_bounds(5,60);bl,bh=atan_bounds(239,18)
    pl,ph=16*al-4*bh,16*ah-4*bl
    assert 3<pl<ph<4
    # tan(4atan(1/5)-atan(1/239))=1; this angle lies in(0,pi/2).
    t2=2*F(1,5)/(1-F(1,25));t4=2*t2/(1-t2*t2)
    assert(t4-F(1,239))/(1+t4*F(1,239))==1
    rows=[]
    for m in range(2,121):
        xlo=(pl*TRIG_SCALE/m).numerator//(pl*TRIG_SCALE/m).denominator
        xhi=ceildiv((ph*TRIG_SCALE/m).numerator,(ph*TRIG_SCALE/m).denominator)
        cl,ch=cos_interval(xlo,xhi)
        low=(-2*ch*PAIR_SCALE)//TRIG_SCALE
        high=ceildiv(-2*cl*PAIR_SCALE,TRIG_SCALE)
        rows.append({'label':m,'lower':str(low),'upper':str(high)})
    assert int(rows[0]['lower'])<=0<=int(rows[0]['upper'])
    assert int(rows[1]['lower'])<=-PAIR_SCALE<=int(rows[1]['upper'])
    assert all(int(rows[i]['upper'])<int(rows[i-1]['lower'])for i in range(1,len(rows)))
    return {'pi_interval':[str(pl),str(ph)],'atan_terms':{'5':60,'239':18},
            'trig_scale':str(TRIG_SCALE),'cos_terms':48,'pair_scale':str(PAIR_SCALE),
            'finite_angles':rows,'tail_upper':rows[-1]['upper']}


class AngleGate:
    def __init__(self,proof):
        self.tail=int(proof['tail_upper'])
        intervals=list(reversed(proof['finite_angles'][:-1]))
        self.lower=[int(x['lower'])for x in intervals]
        self.upper=[int(x['upper'])for x in intervals]
    def rejection(self,low,high):
        if low>0:return'positive'
        if low<=self.tail:return None
        i=bisect_left(self.upper,low)
        if i<len(self.upper)and self.lower[i]<=high:return None
        return'finite_angle_gap'


def coefficient_intervals(source):
    enc=Enclosures();normals=[];covectors=[]
    for i,row in enumerate(source['caps']):
        assert row['index']==i;d=parse(row['radicand'])
        r=[enc.element(tuple(map(parse,x)),d)for x in row['normal']]
        g=[enc.element(tuple(map(parse,x)),d)for x in row['covector']]
        assert len(r)==len(g)==16 and all(hi<=0 for lo,hi in g)
        normals.append(r);covectors.append(g)
    return enc,normals,covectors


def interval_dot_fast(a,b):
    # All b intervals are nonpositive, checked during construction.
    low=high=0
    for(al,ah),(bl,bh)in zip(a,b):
        if al>=0:low+=ah*bl;high+=al*bh
        elif ah<=0:low+=ah*bh;high+=al*bl
        else:low+=ah*bl;high+=al*bl
    return low,high


def interval_dot_general(a,b):
    low=high=0
    for(al,ah),(bl,bh)in zip(a,b):
        p=(al*bl,al*bh,ah*bl,ah*bh)
        low+=min(p);high+=max(p)
    return low,high


def source_and_bounds():
    source=json.loads(SOURCE.read_text());assert len(source['caps'])==3040
    enc,normal,covector=coefficient_intervals(source)
    angles=angle_enclosures()
    return source,enc,normal,covector,angles


def build():
    start=time.time();source,enc,normal,covector,angles=source_and_bounds()
    n=len(normal);adj=[0]*n;bits=[1<<i for i in range(n)];counts=Counter();gate=AngleGate(angles)
    print('Exact coefficient and all-label angle enclosures ready',flush=True)
    for i in range(n):
        ai=adj[i];a=normal[i];bit=bits[i]
        for j in range(i+1,n):
            low,high=interval_dot_fast(a,covector[j]);reason=gate.rejection(low,high)
            if reason is None:ai|=bits[j];adj[j]|=bit;counts['possible_edges']+=1
            else:counts[reason]+=1
        adj[i]=ai
        if i%250==0:print('graph rows',i,'seconds',round(time.time()-start,1),dict(counts),flush=True)
    hexadecimal=[format(x,'x')for x in adj]
    assert sum(counts.values())==n*(n-1)//2
    out={'status':'COMPLETE_EXACT_SAFE_COMPATIBILITY_SUPERGRAPH','entry_count':n,
         'exact_distinct_mirrors':3034,'source':str(SOURCE.relative_to(ROOT)),
         'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         'coefficient_scale':str(SCALE),'fine_scale':str(FINE),
         'phi_fine_interval':list(map(str,enc.phi)),
         'angle_proof':angles,'pair_counts':dict(counts),
         'adjacency_hex':hexadecimal,'adjacency_sha256':digest(hexadecimal),
         'scope':'Every genuine compatible pair of DISTINCT cap mirrors is retained. Rejected intervals are provably positive or disjoint from all finite Coxeter values; dotted pairs and the entire label>=120 tail are retained. Unknown pairs may be edges. Six duplicate profile IDs are retained, but a geometric facet selection can use one representative per mirror. No clique, SAT, or volume conclusion follows from this graph alone.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print('COMPLETE',dict(counts),'seconds',round(time.time()-start,1),flush=True)


def verify():
    # Independent exact sign implementation from the completed normal audit.
    import fresh_h9_split_leaf_3040_independent_20261013 as Exact
    start=time.time();source,enc,normal,covector,angles=source_and_bounds()
    saved=json.loads(OUT.read_text());n=len(normal)
    assert saved['status']=='COMPLETE_EXACT_SAFE_COMPATIBILITY_SUPERGRAPH'
    assert saved['source']==str(SOURCE.relative_to(ROOT))
    assert saved['source_sha256']==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert saved['entry_count']==n==3040
    assert saved['coefficient_scale']==str(SCALE)and saved['fine_scale']==str(FINE)
    assert saved['phi_fine_interval']==list(map(str,enc.phi))
    assert saved['angle_proof']==angles
    # Certify every derived coefficient enclosure with the independent field
    # sign routine, including degenerate radicands and exact endpoint equality.
    seen=set();boundchecks=0
    for raw,nr,cg in zip(source['caps'],normal,covector):
        d=Exact.parse(raw['radicand']);field=Exact.Quadratic(d)
        for vals,bounds in((raw['normal'],nr),(raw['covector'],cg)):
            for value,(lo,hi)in zip(vals,bounds):
                x=tuple(map(Exact.parse,value));key=(x,d,lo,hi)
                if key in seen:continue
                seen.add(key)
                assert field.sign(field.sub(x,((F(lo,SCALE),F(0)),Exact.Z)))>=0
                assert field.sign(field.sub(((F(hi,SCALE),F(0)),Exact.Z),x))>=0
                boundchecks+=1
    print('Independent exact field check of all coefficient bounds PASS',boundchecks,flush=True)
    hexadecimal=saved['adjacency_hex'];assert len(hexadecimal)==n and digest(hexadecimal)==saved['adjacency_sha256']
    adj=list(map(lambda s:int(s,16),hexadecimal));mask=(1<<n)-1
    assert all(a>=0 and not a&~mask and not(a>>i&1)for i,a in enumerate(adj))
    counts=Counter();gate=AngleGate(angles);first_nonedge=None
    for i in range(n):
        for j in range(i+1,n):
            edge=bool(adj[i]>>j&1);assert edge==bool(adj[j]>>i&1)
            if edge:counts['possible_edges']+=1;continue
            # Rejected pairs are replayed with a separate, general interval
            # product implementation. Retained pairs need no exact gate.
            low,high=interval_dot_general(normal[i],covector[j])
            reason=gate.rejection(low,high);assert reason is not None,(i,j)
            counts[reason]+=1
            if first_nonedge is None:first_nonedge=(i,j,low,high)
        if i%500==0:print('verified nonedge rows',i,'seconds',round(time.time()-start,1),flush=True)
    assert dict(counts)==saved['pair_counts']
    # Forged rejection: zero must always be retained as a possible right angle.
    assert gate.rejection(0,0)is None
    # All high labels are deliberately retained without a degree bound.
    assert gate.rejection(-2*PAIR_SCALE,gate.tail)is None
    result={'status':'PASS_SAFE_SUPERGRAPH_REPLAY','standard_library_only':True,
            'entry_count':n,'exact_distinct_mirrors':3034,
            'distinct_coefficient_enclosures_checked_by_independent_exact_field_signs':boundchecks,
            'pair_counts':dict(counts),'adjacency_sha256':saved['adjacency_sha256'],
            'nonedges_checked_with_general_interval_products':counts['positive']+counts['finite_angle_gap'],
            'zero_rejection_disallowed':True,'unbounded_high_label_tail_retained':True,
            'scope':saved['scope']}
    AUDIT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    verify()if args.verify else build()
