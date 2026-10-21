#!/usr/bin/env python3
"""A changed four-exterior-ray chain based on alternative66 and a left cap."""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_alternative66_endcaps_20261017 as E

H,B,M,ROOT=E.H,E.B,E.M,E.ROOT
TAG='fresh_nonarithmetic_alternative66_chain_patch_20261017'
CHOICE=[4,2,3]
FRAMES=[[30,59,67],[66,58,67]]

def detailed_cone(ns,labels):
    cov=np.array([ns[i]for i in labels])@M;cov/=np.linalg.norm(cov,axis=1)[:,None];chart=cov.sum(axis=0)
    rays=[];seen=set()
    for tt in combinations(range(len(labels)),3):
        mat=np.vstack([cov[list(tt)],chart])
        if abs(np.linalg.det(mat))<1e-10:continue
        v=np.linalg.solve(mat,[0,0,0,-1]);v/=np.linalg.norm(v);g=cov@v
        if max(g)>2e-10:continue
        active=tuple(labels[j]for j in range(len(labels))if abs(g[j])<2e-10)
        if active in seen:continue
        seen.add(active);sq=float(v@M@v)
        rays.append(dict(index=len(rays),active=list(active),coordinates=v.tolist(),square=sq,kind='finite'if sq<-1e-12 else'exterior'if sq>1e-12 else'unresolved'))
    pp=defaultdict(list)
    for r in rays:
        for p in combinations(r['active'],2):pp[p].append(r['index'])
    edges=[];warnings=[]
    for p,ids in pp.items():
        if len(ids)!=2:warnings.append(dict(sides=p,rays=ids));continue
        a,b=[rays[i]for i in ids];ab=float(np.array(a['coordinates'])@M@np.array(b['coordinates']));f=sum(r['kind']=='finite'for r in[a,b])
        kind='compact'if f==2 else'one_open'if f==1 else'two_open'if ab<-1e-12 and a['square']*b['square']-ab*ab<-1e-15 else'outside'
        edges.append(dict(sides=list(p),rays=ids,kind=kind))
    return dict(rays=rays,edges=edges,warnings=warnings,counts=dict(Counter(r['kind']for r in rays)),open_ends=sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in edges))

def corner(frame,pi,ms,pair):
    a,b=pair;g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for i,m in enumerate(E.G.profile_labels(pi)):g[i][4]=g[4][i]=m
    g[2][3]=g[3][2]=5 if set(frame[j]for j in pair)=={30,59}else 4 if set(frame[j]for j in pair)=={66,58}else 2
    g[2][4]=g[4][2]=ms[a];g[3][4]=g[4][3]=ms[b]
    return E.components(tuple(q for r in g for q in r))

def domain(frame):
    low=[];tails=[]
    for pi in range(7):
        for ms in product(range(2,7),repeat=3):
            cc=[corner(frame,pi,ms,p)for p in[(0,1),(0,2)]]
            if any(x is None for x in cc):continue
            hi=[i for i,m in enumerate(ms)if m==6];assert len(hi)<=1
            r=dict(profile=pi,labels=list(ms),corner_types=cc)
            if hi:r['unbounded_label_position']=hi[0];tails.append(r)
            else:low.append(r)
    return dict(low=low,tails=tails)

def setup():
    z,labels,ns,h,profiles,primary=H.setup();old=H.A.setup()[0]
    gram=np.array([[B.mid(x)for x in row]for row in old['residual_gram']]);va,ve=np.linalg.eigh((gram+gram.T)/2);tr=ve@np.diag(1/np.sqrt(abs(va)))
    gate=json.loads((ROOT/'results'/f'{E.G.TAG}.json').read_text())
    pr=next(r for r in gate['records'][1]['exact_recipes']if r['profile']==0 and r['labels']==CHOICE)
    rr=next(r for r in pr['branches']if r['sign']==-1);ns[66]=np.linalg.solve(tr,np.array([B.mid(q)for q in rr['coordinates']]))
    parent=json.loads((ROOT/'results'/f'{E.TAG}.json').read_text());stage=next(s for s in parent['stages']if s['choice']==CHOICE)
    left=next(r for r in stage['accepted_low_caps']if r['profile']==0 and r['assigned_labels']==[4,2,2]and r['branch']==-1)
    assert left['index']==20 and not left['blockers']and left['center_product']<0
    ns[67]=np.array(left['projection']);labels=labels+[67]
    for p in profiles:p['offsets'].append(0.)
    actual=z['actual_sides']+[67];squares={r['wall']:B.mid(r['square'])for r in z['projections']};squares[67]=2.
    return z,labels,ns,h,profiles,actual,squares,left

def high_join_events(domains,caps,labels,ns,h,profiles):
    """All algebraic right-event candidates in each numerical high family.

    A high bond to67 is an isolated I2(m) in the required common corner.
    Thus the opposite cap has profile0, is right67, and must be right to
    this cap. Squaring the latter equality gives a quadratic in c.
    Coefficients and root tests are numerical; this is not a sign proof.
    """
    rows=[]
    for gi,frame in enumerate(FRAMES):
        f=np.array([ns[i]for i in frame]);inv=np.linalg.inv(f@M@f.T)
        k=np.linalg.svd(f@M)[2][-1];kk=float(k@M@k)
        for recipe in domains[gi]['tails']:
            assert recipe['profile']==0 and recipe['unbounded_label_position']==2
            rhs=-np.array([B.theta(m)for m in recipe['labels']]);rhs[2]=0
            p0=rhs@inv@f;p1=np.array([0.,0.,-1.])@inv@f
            rr=np.array([2-p0@M@p0,-2*p0@M@p1,-p1@M@p1])/kk
            for v in caps:
                if v['group']==gi or v['profile']!=0 or v['assigned_labels'][2]!=2:continue
                y=np.array(v['projection']);ab=np.array([p0@M@y,p1@M@y]);d=float(k@M@y)
                poly=np.polynomial.polynomial.polymul(ab,ab)-d*d*rr
                roots=np.roots(poly[::-1]);candidates=[]
                for root in roots:
                    if abs(root.imag)>1e-8 or not(np.sqrt(3.)-1e-8<=root.real<=2+1e-8):continue
                    c=float(root.real)
                    for branch in E.frame_solver(frame,labels,ns,h,profiles,recipe,c):
                        product_value=float(np.array(branch['projection'])@M@y)
                        if abs(product_value)<1e-5:candidates.append(dict(c=c,branch=branch['branch'],native_product=product_value,center_product=branch['center_product'],blockers=branch['blockers']))
                endpoints=[]
                for c in[np.sqrt(3.),2.]:
                    future=[r for r in E.frame_solver(frame,labels,ns,h,profiles,recipe,c)if r['center_product']<0]
                    assert len(future)==1
                    endpoints.append(float(np.array(future[0]['projection'])@M@y))
                rows.append(dict(group=gi,recipe=recipe,opposite_index=v['index'],quadratic_coefficients=poly.tolist(),
                    polynomial_roots=[[float(q.real),float(q.imag)]for q in roots],right_event_candidates=candidates,endpoint_products=endpoints))
    return rows

def run():
    z,labels,ns,h,profiles,actual,squares,left=setup();cone=detailed_cone(ns,actual)
    assert cone['counts']=={'finite':28,'exterior':6}and cone['open_ends']==10 and not cone['warnings']
    domains=[domain(f)for f in FRAMES];caps=[];rejected=[];groups=[];tails=[]
    for gi,frame in enumerate(FRAMES):
        group=[]
        for recipe in domains[gi]['low']:
            for r in E.frame_solver(frame,labels,ns,h,profiles,recipe):
                row=dict(**r,frame=frame,profile=recipe['profile'],assigned_labels=recipe['labels'],corner_types=recipe['corner_types'],group=gi)
                if r['center_product']>=-1e-8 or r['blockers']:rejected.append(row);continue
                y=np.array(row['projection']);pi=row['profile'];row['index']=len(caps)
                row['incidence']=H.incidence_score([ns[i]for i in actual]+[y],actual+[68],[squares[i]for i in actual]+[profiles[pi]['square']])
                row['boundary']=B.ledger(np.array([ns[i]for i in actual]+[y]),actual+[68]);group.append(len(caps));caps.append(row)
        groups.append(group)
        for recipe in domains[gi]['tails']:
            samples=[]
            for c in np.linspace(np.sqrt(3.),2.,33):
                for r in E.frame_solver(frame,labels,ns,h,profiles,recipe,float(c)):
                    samples.append(dict(c=float(c),branch=r['branch'],center_product=r['center_product'],blockers=r['blockers']))
            tails.append(dict(group=gi,recipe=recipe,samples=samples))
    pairs=[];pc=Counter();improved=[]
    for ia,ib in product(*groups):
        a,b=caps[ia],caps[ib];ya,yb=np.array(a['projection']),np.array(b['projection']);q=H.A.native(ya,a['profile'],yb,b['profile'],profiles);rel=B.relation(q);pc[str(rel)]+=1
        if not isinstance(rel,int):continue
        g=[[2]*5 for _ in range(5)]
        for i in range(5):g[i][i]=1
        g[0][1]=g[1][0]=3
        for k,r in[(3,a),(4,b)]:
            for j,m in enumerate(E.G.profile_labels(r['profile'])):g[j][k]=g[k][j]=m
            g[2][k]=g[k][2]=dict(r['finite_contacts'])[67]
        g[3][4]=g[4][3]=rel;ty=E.components(tuple(q for r in g for q in r))
        if ty is None:continue
        inc=H.incidence_score([ns[i]for i in actual]+[ya,yb],actual+[68,69],[squares[i]for i in actual]+[profiles[a['profile']]['square'],profiles[b['profile']]['square']])
        bd=B.ledger(np.array([ns[i]for i in actual]+[ya,yb]),actual+[68,69]);r=dict(indices=[ia,ib],mutual=q,label=rel,shared67_corner=ty,incidence=inc,boundary=bd);pairs.append(r)
        if inc.get('open_ends',99)<10:improved.append(r)
    events=high_join_events(domains,caps,labels,ns,h,profiles)
    assert len(events)==10
    assert all(len(r['finite_contacts'])==3 for r in caps)
    outer_products=[float(ns[f[1]]@M@ns[f[2]])for f in FRAMES]
    assert all(q<-2 for q in outer_products)
    paths=[Path(__file__),Path(E.__file__),ROOT/'results'/f'{E.TAG}.json',H.PARENT,H.BOUNDARY,ROOT/'results'/f'{E.G.TAG}.json']
    print('DOMAINS',[(len(d['low']),len(d['tails']))for d in domains],'CAPS',list(map(len,groups)),'PAIRS',dict(pc),'IMPROVEMENTS',len(improved),flush=True)
    return dict(status='BOUNDED_NUMERICAL_CHANGED66_CHAIN_TOPOLOGY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        choice66=CHOICE,left67_recipe=left,frames=FRAMES,control_cone=cone,domains=domains,groups=groups,accepted=caps,rejected=rejected,
        continuous_tail_samples=tails,high_label_right_event_quadratics=events,pair_relations=dict(pc),finite_pairs=pairs,improvements=improved,
        conditional_new_affine_exclusion=dict(outer_anchor_products=outer_products,
            argument='For any exact realization with these native relations, a new affine circuit uses only the three finite old neighbors and A2. The two outer neighbors are dotted to one another, so an affine circuit cannot contain both. It is therefore a principal subset of one of the two certified finite endpoint diagrams, a contradiction. This excludes new affine circuits, conditional on native realization; it does not replay the old affine inventory.'),
        scope='Changed native-coherent numerical64-wall branch; no further deletion. Complete finite endpoint shape enumeration with isolated high-label tails retained explicitly, low-label numerical roots,33-point continuous tail samples. Numerical full intrinsic boundary and separate PD incidence score. No exact promotion or all-parameter tail exclusion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
