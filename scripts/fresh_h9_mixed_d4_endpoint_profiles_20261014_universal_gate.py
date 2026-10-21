#!/usr/bin/env python3
"""Exact dyadic interval exclusion of old-edge endpoint caps cutting one gap.

The input finite profiles and isolated dihedral tails are already complete.
No floating-point value is used for acceptance or for a stopping condition.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate'
BASE = 'results/fresh_h9_mixed_d4_coupling_20261014.json'
PROFILES = 'results/fresh_h9_mixed_d4_endpoint_profiles_20261014.json'
BITS = 128
SCALE = 1 << BITS


class I:
    __slots__ = ('lo', 'hi')
    def __init__(self, x=0):
        if isinstance(x, I):
            self.lo, self.hi = x.lo, x.hi
        else:
            x = F(x)
            self.lo = x.numerator * SCALE // x.denominator
            self.hi = -((-x.numerator * SCALE) // x.denominator)
    @classmethod
    def raw(cls, lo, hi):
        assert lo <= hi
        r = object.__new__(cls); r.lo, r.hi = lo, hi
        return r
    def __add__(self, y):
        y = I(y); return I.raw(self.lo+y.lo, self.hi+y.hi)
    __radd__ = __add__
    def __neg__(self): return I.raw(-self.hi, -self.lo)
    def __sub__(self, y): return self+-I(y)
    def __rsub__(self, y): return I(y)+-self
    def __mul__(self, y):
        y = I(y)
        p = [self.lo*y.lo, self.lo*y.hi, self.hi*y.lo, self.hi*y.hi]
        return I.raw(min(p)//SCALE, -((-max(p))//SCALE))
    __rmul__ = __mul__
    def __truediv__(self, y):
        y = I(y)
        assert y.lo > 0 or y.hi < 0
        p = [F(x*SCALE, z) for x in (self.lo,self.hi) for z in (y.lo,y.hi)]
        l, h = min(p), max(p)
        return I.raw(l.numerator//l.denominator, -((-h.numerator)//h.denominator))
    def __rtruediv__(self, y): return I(y)/self
    def square(self):
        a = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        b = max(self.lo*self.lo, self.hi*self.hi)
        return I.raw(a//SCALE, -((-b)//SCALE))
    def sqrt(self):
        assert self.lo >= 0
        l, h = isqrt(self.lo*SCALE), isqrt(self.hi*SCALE)
        if h*h < self.hi*SCALE: h += 1
        return I.raw(l, h)
    def enc(self): return [str(self.lo), str(self.hi)]
    def overlap(self, y):
        y=I(y); return self.lo <= y.hi and y.lo <= self.hi


SQ5=I(5).sqrt(); SQ2=I(2).sqrt(); SQ3=I(3).sqrt()
PHI=(1+SQ5)/2
RR=((16+4*SQ5)/11).sqrt()
def qi(q): return I(q.a)+I(q.b)*SQ5
def ei(e): return qi(e.a)+qi(e.b)*RR
def jq(raw): return A.Q(*raw)
def isdot(x,y): return sum((a*b for a,b in zip(x,y)),I(0))


def atan_rational(q):
    """Alternating arctangent series, with an explicit next-term enclosure."""
    q=F(q);assert 0<q<1
    total=F(0);j=0
    while True:
        term=q**(2*j+1)/(2*j+1)
        total+=term if j%2==0 else -term
        j+=1
        following=q**(2*j+1)/(2*j+1)
        if following < F(1,1<<200):
            end=total+(following if j%2==0 else -following)
            return I.raw(I(min(total,end)).lo,I(max(total,end)).hi)


# Machin's identity: pi=16 atan(1/5)-4 atan(1/239).
PI=16*atan_rational(F(1,5))-4*atan_rational(F(1,239))
assert I(3).hi<PI.lo and PI.hi<I(F(22,7)).lo
def twice_cos(m):
    assert m>=2
    x=PI/m;assert x.hi<2*SCALE
    xx=x.square();term=I(1);s=term
    for j in range(1,65):
        term=term*(-xx)/((2*j-1)*(2*j));s+=term
    # For x<2 the next term is bounded by 2^130/130!.
    fact=1
    for j in range(1,131):fact*=j
    rem=F(1<<130,fact);assert rem<F(1,1<<140)
    err=I(rem)
    return 2*I.raw(s.lo-err.hi,s.hi+err.hi)

# A rational point selected from a numerical proposal. Every use below is exact.
ESCAPE = [
    '3.445009490103295','4.274481843159041','3.0011044891310377',
    '1.7539516778449502','.5329849937376314','.15617993073424302',
    '-1.048155530053127','-.7473950044921694','-.656929489224263',
    '-.665829155523975']


def coefficients():
    data=json.loads((ROOT/BASE).read_text())
    pp=json.loads((ROOT/PROFILES).read_text())
    for obj in (data,pp):
        assert obj['sources']
        for path,digest in obj['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest, path
    assert len(pp['profiles'])==pp['open_edge_count']==52
    assert sum(len(r['finite_profiles']) for r in pp['profiles'].values())==pp['finite_profile_count']==1825
    assert sum(len(r['dihedral_tail_profiles']) for r in pp['profiles'].values())==pp['unbounded_tail_count']==72
    M,N=A.dm(data['form']),A.dm(data['final_normals'])
    E=[A.E(F(s)) for s in ESCAPE]
    H=A.dv(data['center'])
    en=A.pair(E,E,M);eh=A.pair(E,H,M)
    ep=[A.pair(E,n,M) for n in N]
    assert en.sign()>0 and eh.sign()<0 and all(x.sign()<0 for x in ep)
    norm=[A.pair(n,n,M) for n in N]
    scales=[(2/ei(q)).sqrt() for q in norm]
    cross=[[ei(A.pair(a,b,M))*scales[i]*scales[j]
            for j,b in enumerate(N)] for i,a in enumerate(N)]
    crossE=[ei(x)*s for x,s in zip(ep,scales)]
    edges=[]
    for es,row in pp['profiles'].items():
        support=row['support']; outside=[i for i in range(len(N)) if i not in support]
        B=[[jq(v) for v in rr] for rr in row['normalized_Gram']]
        J=[[jq(v) for v in rr] for rr in row['normalized_inverse']]
        assert all(sum((B[i][k]*J[k][j] for k in range(8)),A.Q())==A.Q(i==j)
                   for i in range(8) for j in range(8))
        # Verify the saved normalized Gram by signs and squared exact values.
        for i,ii in enumerate(support):
            for j,jj in enumerate(support):
                raw=A.pair(N[ii],N[jj],M); want=A.E(B[i][j])
                assert raw.sign()==want.sign()
                assert 4*raw*raw==want*want*norm[ii]*norm[jj]
        assert A.inertia([[A.E(v) for v in rr] for rr in B])==[8,0]
        frame=row['oriented_residual_frame']
        X,Y=A.dv(frame['X']),A.dv(frame['Y'])
        xx=-A.pair(X,X,M); yy=A.pair(Y,Y,M)
        assert xx.sign()>0 and yy.sign()>0 and not A.pair(X,Y,M)
        assert A.pair(X,H,M).sign()<0
        assert all(A.pair(X,n,M).sign()<=0 for n in N)
        assert all(not A.pair(Z,N[i],M) for i in support for Z in (X,Y))
        assert A.pair(X,N[row['released_wall']],M)==0
        assert A.pair(Y,N[row['released_wall']],M).sign()<0
        dx,dy=ei(xx).sqrt(),ei(yy).sqrt()
        avec=[ei(A.pair(X,N[j],M))*scales[j]/dx for j in outside]
        cvec=[ei(A.pair(Y,N[j],M))*scales[j]/dy for j in outside]
        avec.append(ei(A.pair(X,E,M))/dx)
        cvec.append(ei(A.pair(Y,E,M))/dy)
        assert avec[-1].hi<0<cvec[-1].lo and (avec[-1]+cvec[-1]).hi<0
        # Every scalar pairing strictly decreases with h>=0.
        for a,c in zip(avec,cvec):
            assert a.hi<=0 and (a+c).hi<0
        cols=[[cross[i][j] for i in support] for j in outside]
        cols.append([crossE[i] for i in support])
        profiledata=[]
        for p in row['finite_profiles']:
            cc=[I(0)]*8
            active=[];cost=A.Q()
            for wall,m in p['contacts']:
                assert m in (3,4,5)
                j=support.index(wall);active.append(j)
                cc[j]={3:I(1),4:SQ2,5:PHI}[m]
                cost+=J[j][j]*{3:A.Q(1),4:A.Q(2),5:A.Q(F(3,2),F(1,2))}[m]
            assert all(not J[i][j] for i in active for j in active if i!=j)
            assert jq(p['residual_norm'])==2-cost
            weights=[sum((qi(J[i][j])*cc[j] for j in range(8)),I(0)) for i in range(8)]
            b=[-isdot(weights,col) for col in cols]
            rho=qi(jq(p['residual_norm']));assert rho.lo>0
            profiledata.append(dict(kind='finite',contacts=p['contacts'],b=b,rho=rho))
        for tail in row['dihedral_tail_profiles']:
            j=support.index(tail['wall'])
            assert J[j][j]==A.Q(F(1,2)) and all(not J[i][j] for i in range(8) if i!=j)
            weights=[qi(J[i][j]) for i in range(8)]
            b=[-isdot(weights,col) for col in cols]
            profiledata.append(dict(kind='tail',wall=tail['wall'],b=b))
        edges.append(dict(edge_id=int(es),support=support,outside=outside,a=avec,c=cvec,profiles=profiledata))
    return edges,dict(coordinates=ESCAPE,norm=A.enc(en),future_center_pairing=A.enc(eh),old_pairings=A.enc(ep))


ANGLE_TAIL=16
ALLOWED_POINTS=[I(0),I(-1),-SQ2,-PHI,-SQ3]+[-twice_cos(m) for m in range(7,ANGLE_TAIL)]
HIGH_THRESHOLD=-twice_cos(ANGLE_TAIL)
def forbidden(v):
    # An overapproximation of ALL normalized acute Coxeter pairings:
    # Labels 2..15 are isolated, while every label >=16, parallel and
    # dotted contact is retained in a continuous interval.
    return v.lo > HIGH_THRESHOLD.hi and not any(v.overlap(s) for s in ALLOWED_POINTS)


def profile_box(p,k):
    if p['kind']=='finite': return p['b'],p['rho']
    rho=2-k.square()/2
    assert rho.hi>=0
    rho=I.raw(max(0,rho.lo),rho.hi)
    return [v*k for v in p['b']],rho


def evaluate(a,b,c,rho,h):
    """Use monotonicity in h, so its two occurrences do not wrap."""
    lo,hi=I.raw(h.lo,h.lo),I.raw(h.hi,h.hi)
    vl=b+a*hi+c*(hi.square()+rho).sqrt()
    vh=b+a*lo+c*(lo.square()+rho).sqrt()
    # This is valid because a<0 and a+c<0, even when c<0.
    return I.raw(vl.lo,vh.hi)


def certify(edge,p,node_limit=2000000):
    is_tail=p['kind']=='tail'
    k=I.raw(SQ3.lo,2*SCALE) if is_tail else I(1)
    b,rho=profile_box(p,k)
    a,c=edge['a'],edge['c']
    # sqrt(h²+rho)<=h+sqrt(rho), c_E>0, and a_E+c_E<0.
    top=b[-1]+c[-1]*rho.sqrt()
    slope=-(a[-1]+c[-1]);assert slope.lo>0
    bound=top/slope
    H=max(1,(max(0,bound.hi)+SCALE-1)//SCALE+1)
    assert (top+(a[-1]+c[-1])*H).hi<0
    stack=[(k,I.raw(0,H*SCALE),0)]
    counts=Counter(); digest=sha256(); deepest=0
    # A compact exhaustive closed-box tree. Adjacent children overlap at midpoints.
    while stack:
        kk,hh,depth=stack.pop();counts['boxes']+=1;deepest=max(deepest,depth)
        assert counts['boxes']<=node_limit, (edge['edge_id'],p.get('contacts',p.get('wall')),counts)
        bb,rr=profile_box(p,kk)
        vE=evaluate(a[-1],bb[-1],c[-1],rr,hh)
        reason=None
        if vE.hi<0: reason='escape'
        else:
            for j in range(len(edge['outside'])):
                v=evaluate(a[j],bb[j],c[j],rr,hh)
                if forbidden(v): reason='wall:'+str(edge['outside'][j]);break
        if reason is not None:
            counts[reason]+=1
            digest.update(('%s:%d:%d:%d:%d;'%(reason,kk.lo,kk.hi,hh.lo,hh.hi)).encode())
            continue
        assert depth<180, (edge['edge_id'],p.get('contacts',p.get('wall')),kk.enc(),hh.enc())
        # Balance widths after scaling the two dimensions to the original box.
        split_k=is_tail and (kk.hi-kk.lo)*H*SCALE > (hh.hi-hh.lo)*(2*SCALE-SQ3.lo)
        if split_k:
            mid=(kk.lo+kk.hi)//2;assert kk.lo<mid<kk.hi
            stack.append((I.raw(mid,kk.hi),hh,depth+1));stack.append((I.raw(kk.lo,mid),hh,depth+1))
        else:
            mid=(hh.lo+hh.hi)//2;assert hh.lo<mid<hh.hi, (edge['edge_id'],p.get('contacts',p.get('wall')),kk.enc(),hh.enc(),vE.enc(),[(j,evaluate(a[j],bb[j],c[j],rr,hh).enc()) for j in range(len(edge['outside']))])
            stack.append((kk,I.raw(mid,hh.hi),depth+1));stack.append((kk,I.raw(hh.lo,mid),depth+1))
    out=dict(kind=p['kind'],h_bound=H,counts=dict(counts),maximum_depth=deepest,tree_sha256=digest.hexdigest())
    if is_tail:out['wall']=p['wall']
    else:out['contacts']=p['contacts']
    return out


def run(limit=None):
    edges,escape=coefficients();records=[];total=Counter()
    for edge in edges:
        rows=[]
        for p in edge['profiles']:
            r=certify(edge,p);rows.append(r);total.update(r['counts']);total[p['kind']+'_profiles']+=1
        records.append(dict(edge_id=edge['edge_id'],profiles=rows))
        print(edge['edge_id'],len(rows),'boxes',sum(r['counts']['boxes'] for r in rows),flush=True)
        if limit is not None and len(records)>=limit:break
    if limit is None:
        assert len(records)==52 and total['finite_profiles']==1825 and total['tail_profiles']==72
    sources=[BASE,PROFILES,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_ENDPOINT_PROFILES_AND_CONTINUOUS_TAILS_MISS_GAP',
        sources={s:sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
        dyadic_precision_bits=BITS,escape=escape,counts=dict(total),edges=records,
        parameterization='r=p+hU+sqrt(h^2+rho)V, h>=0; for each tail p=k*p0, rho=2-k^2/2, sqrt3<=k<=2',
        allowed_pairing_relaxation='{ -2cos(pi/m): m=2..15 } union (-infinity,-2cos(pi/16)]',
        cosine_enclosures=dict(pi=PI.enc(),isolated=[x.enc() for x in ALLOWED_POINTS],tail_threshold=HIGH_THRESHOLD.enc()),
        cosine_error_proof='Machin identity with alternating arctan next-term error<2^-200;degree128 cosine Taylor with absolute remainder<=2^130/130!<2^-140;all operations rounded outward to dyadic128',
        max_depth=max(r['maximum_depth'] for e in records for r in e['profiles']),
        terminal_boxes=sum(v for k,v in total.items() if k=='escape' or k.startswith('wall:')),
        scope='Every new finite endpoint on the saved old52 outgoing edges, retaining the12 old walls with Coxeter pairings, strictly misses the saved rational spacelike point. This is not an exclusion of other caps, later edges, or changed cores.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--limit',type=int)
    args=ap.parse_args();result=run(args.limit)
    path=ROOT/f'results/{TAG}.json'
    if args.limit is not None:
        print('PARTIAL',result['counts'])
    elif args.verify:
        assert result==json.loads(path.read_text());print('SAVED_REPLAY_PASS',result['counts'])
    else:
        path.write_text(json.dumps(result,indent=2)+'\n');print(result['status'],result['counts'])
