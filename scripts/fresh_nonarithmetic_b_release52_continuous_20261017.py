#!/usr/bin/env python3
"""Two-contact conics after release52; coupled finite-angle patch discovery."""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_b_release52_patch_20261017 as P
ROOT=P.ROOT;B=P.B;J=P.J;M=B.METRIC
TAG='fresh_nonarithmetic_b_release52_continuous_20261017'
ALPHABET=[2,3,4,5,6,8,10,12]
PAIRS=[[8,12],[12,47],[8,65],[64,65],[30,64],[7,30],[7,47]]

def roots(a,b,c):
    scale=max(abs(a),abs(b),abs(c),1.)
    if abs(a)<1e-12*scale:return []if abs(b)<1e-12*scale else[-c/b]
    disc=b*b-4*a*c
    if disc<-1e-10*scale*scale:return[]
    if disc<=0:return[-b/(2*a)]
    q=-(b+math.copysign(math.sqrt(disc),b))/2
    return sorted([q/a,c/q])

def intervals(polys):
    rr=sorted(set(round(t,11)for a,b,c in polys for t in roots(a,b,c)if t>1e-8 and t<1e8))
    cuts=[0.]+rr+[None];out=[]
    for lo,hi in zip(cuts,cuts[1:]):
        t=(lo+1 if hi is None else hi/2 if lo==0 else math.sqrt(lo*hi))
        if all(a*t+b+c/t < -1e-8 for a,b,c in polys):out.append([lo,hi])
    return out

def samples(iv):
    lo,hi=iv
    if hi is None:vv=[max(lo,1.)*f+1e-5 for f in[1.01,1.1,1.5,2,3,5,8,13,21,34]]
    elif lo==0:vv=[hi*f for f in[.01,.03,.1,.2,.35,.5,.65,.8,.9,.97,.99]]
    else:vv=[lo*(hi/lo)**f for f in[.01,.05,.15,.3,.5,.7,.85,.95,.99]]
    vv +=[float(k)for k in[1,2,3,4,5,7,10,16,25,40,64,100]if k>lo+1e-8 and(hi is None or k<hi-1e-8)]
    return sorted(set(vv))

def inside(t,ivs):return t>0 and any(t>lo+1e-8 and(hi is None or t<hi-1e-8)for lo,hi in ivs)
def evaluate(f,t):return np.array(f['p'])+t*np.array(f['a'])+np.array(f['b'])/t

def families(data):
    z,labels,ns,h,profiles=data;retained=[i for i in labels if i!=52];out=[];counts=Counter()
    for group,ids in enumerate(PAIRS):
        frame=np.array([ns[i]for i in ids]);g=frame@M@frame.T
        assert np.linalg.eigvalsh(g)[0]>0
        gi=np.linalg.inv(g);h0=h-(h@M@frame.T)@gi@frame;h0/=math.sqrt(-h0@M@h0)
        kk=np.linalg.svd(np.vstack([frame@M,h0@M]))[2][-1];kk/=math.sqrt(kk@M@kk)
        for pi,pr in enumerate(profiles):
            off=np.array(pr['offsets']);q=pr['square']
            for ms in product(ALPHABET,repeat=2):
                cc=-np.array([B.theta(m)for m in ms])-off[[labels.index(i)for i in ids]];p=cc@gi@frame;rho=q-p@M@p
                if rho<1e-8:continue
                for sign in[-1,1]:
                    a=math.sqrt(rho)*(h0+sign*kk)/2;b=math.sqrt(rho)*(-h0+sign*kk)/2
                    polys=[];bad=False;fixed=[]
                    for i in retained:
                        if i in ids:continue
                        n=ns[i];aa=a@M@n;bb=b@M@n;cc0=p@M@n+off[labels.index(i)]
                        if max(abs(aa),abs(bb))<1e-8:
                            rel=B.relation(cc0)
                            if rel is None or rel=='parallel_or_unresolved':bad=True;break
                            fixed.append([i,rel])
                        else:polys.append([float(aa),float(cc0+2),float(bb)])
                    if bad:counts['constant_gap']+=1;continue
                    polys.append([float(a@M@h),float(p@M@h+pr['center_offset']),float(b@M@h)])
                    ivs=intervals(polys);counts['conics']+=1
                    if not ivs:continue
                    out.append(dict(index=len(out),group=group,anchors=ids,labels=list(ms),profile=pi,branch=sign,
                        p=p.tolist(),a=a.tolist(),b=b.tolist(),intervals=ivs,fixed_contacts=fixed,
                        square=q,old_anchor_products=(-np.array([B.theta(m)for m in ms])).tolist()))
        print('FAMILY_GROUP',group,'count',sum(r['group']==group for r in out),flush=True)
    return out,counts

def gate(y,pi,data):
    z,labels,ns,h,profiles=data;pr=profiles[pi];off=pr['offsets']
    if y@M@h+pr['center_offset']>=-1e-7:return False
    return all((r:=B.relation(y@M@ns[i]+off[labels.index(i)]))is not None and r!='parallel_or_unresolved'for i in labels if i!=52)

def compute():
    data=P.setup();z,labels,ns,h,profiles=data;ff,counts=families(data)
    actual=[i for i in z['actual_sides']if i!=52];groups=[[f for f in ff if f['group']==g]for g in range(len(PAIRS))]
    pairs=[];seen=set()
    # Rebuild the entire12-face from its two exposed8/12 and12/47 ends.
    for fa in groups[0]:
        pa=profiles[fa['profile']];ca=np.array(pa['contacts'])
        for fb in groups[1]:
            pb=profiles[fb['profile']];cb=np.array(pb['contacts']);shift=ca@J@cb
            a12=-B.theta(fa['labels'][1]);b12=-B.theta(fb['labels'][0])
            for m in ALPHABET:
                ab=-B.theta(m)-shift
                corner=np.array([[2,a12,b12],[a12,pa['square'],ab],[b12,ab,pb['square']]])
                if np.linalg.eigvalsh(corner)[0]<1e-8:continue
                for iv in fa['intervals']:
                    for t in samples(iv):
                        a=evaluate(fa,t)
                        for u in roots(a@M@np.array(fb['a']),a@M@np.array(fb['p'])+shift+B.theta(m),a@M@np.array(fb['b'])):
                            if not inside(u,fb['intervals']):continue
                            b=evaluate(fb,u)
                            if not gate(a,fa['profile'],data)or not gate(b,fb['profile'],data):continue
                            key=tuple(np.round(np.r_[a,b],7))
                            if key in seen:continue
                            seen.add(key)
                            led=B.ledger(np.array([ns[i]for i in actual]+[a,b]),actual+[66,67])
                            row=dict(index=len(pairs),families=[fa['index'],fb['index']],parameters=[t,float(u)],mutual_label=m,
                                projections=[a.tolist(),b.tolist()],profiles=[fa['profile'],fb['profile']],boundary=led)
                            pairs.append(row)
        print('LEFT_DONE',fa['index'],'pairs',len(pairs),flush=True)
    pairs.sort(key=lambda r:(r['boundary']['open_ends'],r['boundary']['exterior'],-r['boundary']['finite']))
    return dict(status='BOUNDED_NUMERICAL_CONTINUOUS_RELEASE52_PATCH',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(P.__file__),P.C.PARENT,
            ROOT/'results/fresh_nonarithmetic_replace47_close12_20261017.json']},
        deleted_wall=52,profiles=profiles,pairs_of_old_anchors=PAIRS,alphabet=ALPHABET,families=ff,counts=dict(counts),pairs=pairs,
        scope='All seven A2 profiles at sampled anchor labels. Open all-other-dotted parameter intervals are obtained by every quadratic boundary event; isolated additional finite old contacts are not part of these open intervals. One conic parameter is sampled explicitly and the other solved from the prescribed mutual finite angle. Full numeric S6 boundary is scored. This is not an all-parameter exclusion or exact new chamber.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(compute()))
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['families']),len(out['pairs']),flush=True)
    if out['pairs']:print('BEST',out['pairs'][0],flush=True)
if __name__=='__main__':main()
