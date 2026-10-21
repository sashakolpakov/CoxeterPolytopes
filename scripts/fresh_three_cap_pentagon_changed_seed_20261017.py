#!/usr/bin/env python3
"""Changed exterior-pentagon one/two/three-cap discovery; numerical geometry only."""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release12_47_52_pentagon_20261017 as S
import fresh_nonarithmetic_alternative66_chain_patch_20261017 as C

A,H,B,G,M,ROOT=S.A,S.H,S.B,S.G,S.M,S.ROOT
TAG='fresh_three_cap_pentagon_changed_seed_20261017'
DELETE=[12,47,52,64,65]
CYCLE=[7,8,44,63,30]
EDGE={frozenset((a,b)):m for a,b,m in[(7,8,5),(8,44,2),(44,63,2),(63,30,2),(30,7,5)]}

@lru_cache(None)
def finite_corner(a,b,pi,ma,mb):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for j,k in enumerate((a,b)):
        for i,m in enumerate(G.profile_labels(4 if k==7 else 0)):g[i][j+2]=g[j+2][i]=m
    for i,m in enumerate(G.profile_labels(pi)):g[i][4]=g[4][i]=m
    g[2][3]=g[3][2]=EDGE[frozenset((a,b))]
    g[2][4]=g[4][2]=ma;g[3][4]=g[4][3]=mb
    return A.components(tuple(x for r in g for x in r))

def domain(frame,closed=False):
    low=[];tails=[]
    pairs=list(zip(range(len(frame)-1),range(1,len(frame))))
    if closed:pairs.append((len(frame)-1,0))
    for pi in range(7):
        for ms in product(range(2,7),repeat=len(frame)):
            cc=[finite_corner(frame[i],frame[j],pi,ms[i],ms[j])for i,j in pairs]
            if any(x is None for x in cc):continue
            high=[i for i,m in enumerate(ms)if m==6]
            r=dict(profile=pi,labels=list(ms),corner_types=cc,high_positions=high)
            (tails if high else low).append(r)
    return dict(low=low,tails=tails)

def frame_roots(frame,labels,ns,profiles,recipe,cs=None):
    """Low roots, or norm roots in one continuous high parameter for >=4 contacts."""
    f=np.array([ns[i]for i in frame]);p=profiles[recipe['profile']]
    off=np.array(p['offsets'])[[labels.index(i)for i in frame]]
    native=-np.array([B.theta(m)if m!=6 else 0 for m in recipe['labels']]);high=recipe['high_positions']
    if cs is not None:
        for j,c in zip(high,cs):native[j]=-c
    rhs=native-off
    if len(frame)==3:
        gram=f@M@f.T;inv=np.linalg.inv(gram);k=np.linalg.svd(f@M)[2][-1]
        part=rhs@inv@f;ratio=(p['square']-part@M@part)/(k@M@k)
        if ratio<=1e-10:return []
        return [dict(projection=(part+s*np.sqrt(ratio)*k).tolist(),branch=s,tail_cosines=cs or[])for s in[-1,1]]
    # Full rank linear contacts determine the projected normal uniquely.
    mat=f@M;pinv=np.linalg.pinv(mat);base=pinv@rhs
    if high and cs is None:
        if len(high)!=1:return []
        j=high[0];d=np.zeros(len(frame));d[j]=-1;direction=pinv@d
        poly=[float(direction@M@direction),float(2*base@M@direction),float(base@M@base-p['square'])]
        roots=np.roots(poly if abs(poly[0])>1e-12 else poly[1:]);out=[]
        for r in roots:
            if abs(r.imag)>1e-8 or not np.sqrt(3)-1e-8<=r.real<=2+1e-8:continue
            c=float(r.real);y=base+c*direction
            if max(abs(mat@y-(rhs+c*d)))<1e-7:out.append(dict(projection=y.tolist(),branch=0,tail_cosines=[c],norm_polynomial=poly))
        return out
    if max(abs(mat@base-rhs))<1e-7 and abs(base@M@base-p['square'])<1e-7:
        return [dict(projection=base.tolist(),branch=0,tail_cosines=cs or[])]
    return []

def native_score(root,frame,recipe,labels,ns,h,profiles,actual,squares,with_boundary=True):
    y=np.array(root['projection']);pi=recipe['profile'];p=profiles[pi]
    values={i:float(y@M@ns[i]+p['offsets'][labels.index(i)])for i in labels}
    rel={i:B.relation(x)for i,x in values.items()}
    assigned=dict(zip(frame,recipe['labels']));tail_rels=[]
    for j,c in zip(recipe['high_positions'],root['tail_cosines']):
        q=B.relation(-c);assigned[frame[j]]=q if isinstance(q,int)else'continuous_unresolved';tail_rels.append(q)
    rel.update(assigned)
    bad=[dict(wall=i,product=values[i])for i in labels if i not in DELETE and (rel[i]is None or rel[i]=='parallel_or_unresolved')]
    center=float(y@M@h+p['center_offset'])
    out=dict(**root,frame=frame,recipe=recipe,center_product=center,blockers=bad,
             finite_contacts=[[i,m]for i,m in rel.items()if isinstance(m,int)and i not in DELETE],
             assigned_angles_coxeter=all(isinstance(m,int)for m in assigned.values()),tail_relations=tail_rels)
    if not bad and center<-1e-8 and with_boundary:
        nd=dict(ns);nd[67]=y;out['boundary']=C.detailed_cone(nd,actual+[67])
        out['incidence']=H.incidence_score([ns[i]for i in actual]+[y],actual+[67],[squares[i]for i in actual]+[p['square']])
    return out

def shared_corner(old,a,b,m):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for i,x in enumerate(G.profile_labels(4 if old==7 else 0)):g[i][2]=g[2][i]=x
    for j,r in[(3,a),(4,b)]:
        for i,x in enumerate(G.profile_labels(r['recipe']['profile'])):g[i][j]=g[j][i]=x
        g[2][j]=g[j][2]=dict(r['finite_contacts']).get(old,0)
        if not g[2][j]:return None
    g[3][4]=g[4][3]=m
    return A.components(tuple(x for r in g for x in r))

def run():
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in DELETE]
    squares={r['wall']:B.mid(r['square'])for r in z['projections']}
    cost=H.deletion(primary,z,DELETE);control=C.detailed_cone(ns,actual)
    assert cost['lost_finite_count']==11 and cost['compact_faces_by_dimension']=={2:6}
    assert control['counts']=={'finite':15,'exterior':7}and control['open_ends']==11 and not control['warnings']
    assert any(r['active']==[30,44,63]and r['kind']=='finite'for r in control['rays'])
    for (a,b),m in[(tuple(k),v)for k,v in EDGE.items()]:assert B.relation(float(ns[a]@M@ns[b]))==m
    stages=[];all_caps=[];cache={};unresolved=[]
    frames=[(tuple(CYCLE),True)]+[(tuple(CYCLE[(at+j)%5]for j in range(n)),False)for n in[3,4]for at in range(5)]
    for frame0,closed in frames:
        frame=list(frame0);dom=domain(frame,closed);accepted=[];low_trials=[];tail_trials=[];counts=Counter()
        for recipe in dom['low']+dom['tails']:
            high=recipe['high_positions'];samples=[]
            if high and len(frame)==3:
                if len(high)>1:
                    unresolved.append(dict(frame=frame,recipe=recipe,reason='multivariate_continuous_tail'));continue
                samples=[[float(c)]for c in np.linspace(np.sqrt(3.),2.,33)]
            elif len(high)>1:
                unresolved.append(dict(frame=frame,recipe=recipe,reason='multivariate_norm_locus'));continue
            else:samples=[None]
            for cs in samples:
                roots=frame_roots(frame,labels,ns,profiles,recipe,cs);counts['real_tail_branches'if high else'real_low_branches']+=len(roots)
                for root in roots:
                    row=native_score(root,frame,recipe,labels,ns,h,profiles,actual,squares,with_boundary=not high)
                    if not row['blockers']and row['center_product']<-1e-8:
                        if not high:row['index']=len(all_caps);all_caps.append(row);accepted.append(row['index'])
                        counts['coherent_tail_samples'if high else'coherent_low']+=1
                    (tail_trials if high else low_trials).append(row)
        stage=dict(frame=frame,closed=closed,domain=dom,low_trials=low_trials,tail_trials=tail_trials,accepted_low_indices=accepted,counts=dict(counts))
        stages.append(stage);cache[frame0]=stage
        print('CHANGED_PENTAGON_FRAME',frame,'domain',len(dom['low']),len(dom['tails']),'accepted',len(accepted),'counts',dict(counts),flush=True)
    two=[];three=[]
    for at in range(5):
        lf=tuple(CYCLE[(at+j)%5]for j in[-1,0,1]);rf=tuple(CYCLE[(at+j)%5]for j in[1,2,3,4])
        pairs=[];relcounts=Counter()
        for ia,ib in product(cache[lf]['accepted_low_indices'],cache[rf]['accepted_low_indices']):
            a,b=all_caps[ia],all_caps[ib];ya,yb=map(lambda r:np.array(r['projection']),(a,b))
            v=H.A.native(ya,a['recipe']['profile'],yb,b['recipe']['profile'],profiles);rel=B.relation(v);relcounts[str(rel)]+=1
            if not isinstance(rel,int):continue
            endpoints=[shared_corner(x,a,b,rel)for x in[lf[0],lf[-1]]]
            if any(x is None for x in endpoints):continue
            nd=dict(ns);nd.update({67:ya,68:yb});bd=C.detailed_cone(nd,actual+[67,68]);pairs.append(dict(indices=[ia,ib],mutual=v,label=rel,corner_types=endpoints,boundary=bd))
        two.append(dict(at=at,left_frame=lf,right_frame=rf,pair_relations=dict(relcounts),accepted_pairs=pairs))
        # Two two-edge caps in the central-triangle topology must meet finitely
        # at their common old boundary root before any third cap can help.
        af=tuple(CYCLE[(at+j)%5]for j in[1,2,3]);bf=tuple(CYCLE[(at+j)%5]for j in[3,4,5]);old=af[-1]
        joined=[];pc=Counter()
        for ia,ib in product(cache[af]['accepted_low_indices'],cache[bf]['accepted_low_indices']):
            a,b=all_caps[ia],all_caps[ib];v=H.A.native(np.array(a['projection']),a['recipe']['profile'],np.array(b['projection']),b['recipe']['profile'],profiles);rel=B.relation(v);pc[str(rel)]+=1
            if isinstance(rel,int):
                ty=shared_corner(old,a,b,rel)
                if ty is not None:joined.append(dict(indices=[ia,ib],mutual=v,label=rel,corner_types=ty))
        three.append(dict(at=at,one_edge_frame=[CYCLE[at],CYCLE[(at+1)%5]],two_edge_frames=[af,bf],pair_relations=dict(pc),finite_join_pairs=joined))
        print('CHANGED_PENTAGON_JOINS',at,'two',dict(relcounts),'three',dict(pc),'finite',len(pairs),len(joined),flush=True)
    paths=[Path(__file__),Path(S.__file__),Path(C.__file__),Path(A.__file__),Path(H.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_CHANGED_PENTAGON_ONE_TWO_THREE_CAP_DOMAINS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=DELETE,boundary_cycle=CYCLE,exact_old_deletion_cost=cost,deletion_boundary=control,retained_outside_finite_triangle=[30,44,63],
        frame_stages=stages,accepted_low_caps=all_caps,two_center_splits=two,three_center_splits=three,unresolved_tail_loci=unresolved,
        finite_domain_scope='Complete finite-corner labels2..5 plus symbolic6 for each isolated high I2 component; all seven A2 profiles. Five-contact one-cap and all five two/three-cap central-triangle rotations are represented. One-variable >=4-contact norm roots are solved numerically over the full cosine interval; three-contact tails have33 continuous samples. Multivariate tails are explicit gaps.',
        scope='Numerical root/contact/center and full intrinsic boundary diagnostic, not a new exact chamber or all-parameter exclusion. All native retained contacts are checked. The old chord30-44 and finite triangle30-44-63 lie outside the exterior patch and are retained. No third cap is attempted unless the two mandatory two-edge caps make an actual finite common corner.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/(TAG+'.json')
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'accepted',len(out['accepted_low_caps']),flush=True)
