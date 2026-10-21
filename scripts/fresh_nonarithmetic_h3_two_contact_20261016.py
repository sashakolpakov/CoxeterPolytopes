#!/usr/bin/env python3
"""Continuous norm-conic pilot for two-contact caps on the actual H3 face.

Numerical discovery only: complete quadratic interval decomposition within
each listed bounded two-contact recipe, not an all-label/exact certificate.
"""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
from scipy.linalg import null_space
import fresh_nonarithmetic_h3_closure_20261016_scout as H3

P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_two_contact_20261016'
CAPS='results/fresh_nonarithmetic_nonright_s6_20261016_caps.json'
CONE='results/fresh_nonarithmetic_nonright_s6_20261016_cone.json'
S=H3.S;ETA=np.diag([-1.,1.,1.,1.]);TOL=2e-8

def clean(x):
    if isinstance(x,np.ndarray):return [clean(z)for z in x]
    if isinstance(x,(float,np.floating)):return float(format(float(x),'.13g'))
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,dict):return{k:clean(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return[clean(z)for z in x]
    return x

def compatible(g,tol=2e-7):
    if g < -2-tol:return True
    if abs(g)<tol:return True
    if not -2+tol<g<-1+tol:return False
    m=math.pi/math.acos(np.clip(-g/2,-1,1));return round(m)>=3 and abs(m-round(m))<tol

def setup():
    Q,N,H,E=P.geometry();old=P.P.B.read(H3.CAP);new=P.P.B.read(CAPS);cone=P.P.B.read(CONE)
    N.append(mp.matrix([P.mid(x)for x in old['normal_interval']]))
    N.extend(mp.matrix([P.mid(x)for x in r['normal_interval']])for r in new['caps']);assert len(N)==55
    assert all(abs((n.T*Q*n)[0]-2)<mp.mpf('1e-40')for n in N)
    W=mp.matrix([[N[j][i]for j in S]for i in range(10)]);J=(W.T*Q*W)**-1
    T=H-W*J*W.T*Q*H;T/=mp.sqrt(-(T.T*Q*T)[0]);basis=[T]
    for unused in range(3):
        choices=[]
        for n in N:
            v=n-W*J*W.T*Q*n;v+=((v.T*Q*T)[0])*T
            for z in basis[1:]:v-=((v.T*Q*z)[0])*z
            choices.append(((v.T*Q*v)[0],v))
        q,v=max(choices,key=lambda z:z[0]);assert q>mp.mpf('1e-30');basis.append(v/mp.sqrt(q))
    def coord(v):return np.array([float(-(v.T*Q*basis[0])[0])]+[float((v.T*Q*z)[0])for z in basis[1:]])
    NC=np.array([coord(n)for n in N]);rays=[]
    for kind,sets in cone['complete_census'].items():
        for ids in sets:
            if not set(S)<=set(ids):continue
            v=mp.matrix(P.P.null([list(Q*N[i])for i in ids]));v/=-(v.T*Q*T)[0]
            c=coord(v);assert abs(c[0]-1)<1e-10
            rays.append(dict(active=ids,kind=kind,point=c))
    assert Counter(r['kind']for r in rays)=={'finite_vertices':13,'hyperideal_rays':7}
    return NC,rays,[[mp.nstr(x,65)for x in v]for v in basis]

def roots(q):
    a,b,c=q;s=max(abs(a),abs(b),abs(c),1.)
    if abs(a)<1e-12*s:return[]if abs(b)<1e-12*s else[-c/b]
    d=b*b-4*a*c
    if d<-1e-11*s*s:return[]
    if d<0:d=0
    sd=math.sqrt(d);v=-.5*(b+math.copysign(sd,b));return[-b/(2*a)]if abs(v)<1e-16 else[v/a,c/v]

def intervals(polys,positive):
    endpoints=sorted(z for p in polys for z in roots(p)if math.isfinite(z)and(not positive or z>1e-11))
    unique=[]
    for z in endpoints:
        if not unique or abs(z-unique[-1])>1e-9*max(1,abs(z)):unique.append(z)
    cuts=([0.]if positive else[None])+unique+[None];out=[]
    for a,b in zip(cuts,cuts[1:]):
        if a is None:t=b-max(1,abs(b))
        elif b is None:t=a+max(1,abs(a))
        elif positive and a>0:t=math.sqrt(a*b)
        else:t=(a+b)/2
        if all(np.polyval(p,t)<1e-8*max(1,np.sum(np.abs(p)*np.array([t*t,abs(t),1])))for p in polys):out.append([a,b])
    return out

def evaluate(f,t):
    p,a,b=[np.array(f[k])for k in('p','a','b')]
    if f['kind']=='hyperbola':return p+a*t+b/t
    if f['kind']=='circle':return p+a*(1-t*t)/(1+t*t)+b*(2*t)/(1+t*t)
    return p+a*t

def poly(f,cov,offset=0.):
    p,a,b=[np.array(f[k])for k in('p','a','b')];c=p@cov+offset;aa=a@cov;bb=b@cov
    if f['kind']=='hyperbola':return np.array([aa,c,bb])
    if f['kind']=='circle':return np.array([c-aa,2*bb,c+aa])
    return np.array([0.,aa,c])

def parameter_samples(interval,positive):
    a,b=interval
    if a is None:return[b-d*max(1,abs(b))for d in(.01,.1,1,10)]
    if b is None:return[a+d*max(1,abs(a))for d in(.01,.1,1,10)]
    if positive and a>0:return[math.exp((1-s)*math.log(a)+s*math.log(b))for s in(.001,.1,.25,.5,.75,.9,.999)]
    return[(1-s)*a+s*b for s in(.001,.1,.25,.5,.75,.9,.999)]

def run(maxlabel):
    NC,rays,basis=setup();outside=[j for j in range(55)if j not in S];hyper=[r for r in rays if r['kind']=='hyperideal_rays']
    anchors=sorted({tuple(p)for r in hyper for p in itertools.combinations([j for j in r['active']if j not in S],2)})
    ideal_points=[np.r_[1.,r['point'][1:]/np.linalg.norm(r['point'][1:])]for r in hyper]
    families=[];counts=Counter();degenerate=[];recipes=[]
    for pair in anchors:
        V=NC[list(pair)];G=V@ETA@V.T;det=np.linalg.det(G)
        if abs(det)<1e-9:degenerate.append(dict(anchors=pair,reason='singular_contact_frame'));continue
        Ji=np.linalg.inv(G);K=null_space(V@ETA);ev,U=np.linalg.eigh(K.T@ETA@K);K=K@U@np.diag(1/np.sqrt(abs(ev)))
        if min(abs(ev))<1e-10:degenerate.append(dict(anchors=pair,reason='null_kernel'));continue
        for labels in itertools.product(range(2,maxlabel+1),repeat=2):
            counts['recipes']+=1;c=np.array([0. if m==2 else-2*math.cos(math.pi/m)for m in labels]);p=c@Ji@V;rho=2-p@ETA@p
            fs=[]
            if min(ev)>0:
                if rho>1e-9:fs=[dict(kind='circle',p=p,a=math.sqrt(rho)*K[:,0],b=math.sqrt(rho)*K[:,1],branch=0)]
                elif abs(rho)<1e-9:degenerate.append(dict(anchors=pair,labels=labels,reason='singleton_norm_branch'))
            else:
                if abs(rho)<1e-9:
                    for s in(-1,1):fs.append(dict(kind='line',p=p,a=K[:,0]+s*K[:,1],b=np.zeros(4),branch=s))
                else:
                    a=K[:,1]if rho>0 else K[:,0];b=K[:,0]if rho>0 else K[:,1]
                    for s in(-1,1):fs.append(dict(kind='hyperbola',p=p,a=s*math.sqrt(abs(rho))*(a+b)/2,b=s*math.sqrt(abs(rho))*(a-b)/2,branch=s))
            for f in fs:
                counts['real_conics']+=1;f.update(anchors=pair,labels=labels,projected_gram=G,particular_square=float(p@ETA@p))
                polys=[poly(f,np.array([-1.,0,0,0]))];fixed=[];reject=False
                for j in outside:
                    if j in pair:continue
                    cov=ETA@NC[j];vary=max(abs(np.array(f['a'])@cov),abs(np.array(f['b'])@cov));val=p@cov
                    if vary<1e-8:
                        if not compatible(val):reject=True;counts['constant_contact_gap']+=1;break
                        fixed.append(dict(wall=j,product=float(val)))
                    else:polys.append(poly(f,cov,2.))
                if reject:continue
                ints=intervals(polys,f['kind']=='hyperbola')
                good=[]
                for interval in ints:
                    samples=[]
                    left,right=interval;critical=[]
                    for point in ideal_points:
                        critical.extend(t for t in roots(poly(f,ETA@point))if(left is None or t>left+1e-9)and(right is None or t<right-1e-9))
                    critical=sorted(set(critical));cuts=[left]+critical+[right]
                    params=parameter_samples(interval,f['kind']=='hyperbola')
                    for a,b in zip(cuts,cuts[1:]):params.extend(parameter_samples([a,b],f['kind']=='hyperbola'))
                    for t in sorted(set(params)):
                        v=evaluate(f,t);g=NC@ETA@v
                        if v[0]<1e-7 or abs(v@ETA@v-2)>1e-6:continue
                        if not all(compatible(g[j],2e-5)for j in outside):continue
                        cover=[i for i,r in enumerate(hyper)if v@ETA@r['point']>1e-7]
                        samples.append(dict(parameter=t,coverage=cover,ideal_coverage=[i for i,point in enumerate(ideal_points)if v@ETA@point>1e-7],vector=v,minimum_dotted_margin=min([-2-g[j]for j in outside if j not in pair and all(x['wall']!=j for x in fixed)]+[99.])))
                    if samples:good.append(dict(interval=interval,ideal_pairing_critical_parameters=critical,samples=samples))
                if good:
                    f.update(fixed_contacts=fixed,intervals=good,index=len(families));families.append(f)
        print('pair',pair,'recipes',counts['recipes'],'families',len(families),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{H3.TAG}.py',f'scripts/{P.TAG}.py',H3.CAP,CAPS,CONE,
           'results/fresh_nonarithmetic_h3_joint_cover_20261016.json','scripts/fresh_nonarithmetic_h3_closure_20261016_repair_scout.py']
    return clean(dict(status='NUMERICAL_CONTINUOUS_TWO_CONTACT_H3_DOMAINS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                      support=S,retained_walls=55,max_anchor_label=maxlabel,anchors=anchors,counts=dict(counts),degeneracies=degenerate,
                      basis=basis,old_projected_normals=NC,actual_face_rays=rays,families=families,
                      scope='Numerical rational-conic interval pilot. All pairs of side anchors incident to the seven current exterior H3 rays and finite labels2..maxlabel are included. Remaining nonconstant old contacts are required strictly dotted; constant right/finite products are allowed. Additional isolated finite-contact parameter values, arbitrary labels, other anchor pairs and no/two/one-contact alternatives remain outside. Boundary equality/parallel cases are not accepted. No exact Coxeter cap or compact-face claim.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--maxlabel',type=int,default=10);ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run(args.maxlabel);p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text());print('PASS saved continuous-family numerical replay')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['families']),flush=True)
