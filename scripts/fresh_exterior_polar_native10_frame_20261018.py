#!/usr/bin/env python3
"""Right-polar diagnostic and changed native10 reconstruction of wall30."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
import fresh_compact10_edge_subdivision_20261018 as D

K,H,T,ROOT=D.K,D.H,D.T,D.ROOT
TAG='fresh_exterior_polar_native10_frame_20261018'
SEED=ROOT/'results/fresh_joint10_30_65_compact_pair_20261018.json'


def intervals(s,pairs):
    p=s.primary;dec=lambda x:tuple(map(F,x))
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    h=list(map(dec,p['common_center']));g=[[dec(q) for q in r] for r in p['computational_gram']]
    ns={}
    for r in p['facet_witnesses']:
        hp=dec(p['center_products'][str(r['wall'])]);assert hp[1]<0
        ns[r['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,r['point']))]
    out=[]
    for i,j in sorted(set(tuple(sorted(e)) for e in pairs)):
        q=(F(),F())
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);q=q[0]+z[0],q[1]+z[1]
        v=s.dot(s.old[i],s.ps(i,{}),s.old[j],s.ps(j,{}))
        assert float(q[0])-2e-5<=v<=float(q[1])+2e-5
        out.append(dict(walls=[i,j],native_product=v,relation=T.G.relation(v),exact_native_enclosure=list(map(str,q))))
    return out


def polar(s,seed):
    anchors=[8,11,33,44];n=s.old[10]
    assert s.ps(10,{})==1 and all(s.ps(i,{})==0 for i in anchors)
    assert max(abs(float(n@T.M@s.old[i])) for i in anchors)<2e-7
    assert abs(float(n@T.M@n)-4/3)<2e-7 and abs(float(n@T.M@s.old[7])+1/3)<2e-7
    gram=np.array([[float(s.old[i]@T.M@s.old[j]) for j in [8,11,44]] for i in [8,11,44]])
    assert np.linalg.det(gram)<-1 and np.linalg.matrix_rank(gram)==3
    # Only the seven inherited white-contact profiles are considered.
    # Native norm2 fixes the magnitude on the one-dimensional positive polar.
    base={i:v for i,v in s.old.items() if i not in (10,30,65)};rows=[]
    for pi in range(7):
        scale=np.sqrt(3*s.profiles[pi]['square']/4)
        for sign in (-1,1):
            v=sign*scale*n;values,bad,cp=s.native(v,pi,base,{})
            expected=s.off(pi,s.ps(7,{}))-sign*scale/3
            assert abs(values[7]-expected)<2e-7
            rows.append(dict(profile=pi,sign=sign,scale=float(scale),normal=v.tolist(),
                native7=values[7],native7_formula=float(expected),center_product=cp,
                blockers=[[i,values[i]] for i in bad],native_products=sorted(values.items()),
                restores_old10=pi==1 and sign==1))
    viable=[r for r in rows if not r['blockers'] and r['center_product']<-1e-8]
    assert len(viable)==1 and viable[0]['restores_old10']
    assert all(r['center_product']>0 for r in rows if r['sign']==-1)
    conflicts=[]
    for record in seed['records']:
        r=record['native'][0];v=np.array(r['normal'])
        val=s.dot(n,1,v,r['profile'])
        if not any(x['R_prefix']==record['selection']['R_prefix'] for x in conflicts):
            assert -1<val<0 and T.G.relation(val) is None
            conflicts.append(dict(R_prefix=record['selection']['R_prefix'],native10_product=val))
    assert len(conflicts)==2
    return dict(anchors=anchors,spanning_triple=[8,11,44],spanning_gram=gram.tolist(),
        old10_intrinsic_square=float(n@T.M@n),old10_projected7=float(n@T.M@s.old[7]),
        rows=rows,changed_R_conflicts=conflicts,
        exact_argument='Conditional on the inherited exact control relations,8/11/44 span a Lorentz3 space:8 is right to11/44 and11--44 is dotted. Its positive one-dimensional orthogonal complement is spanned by old10, of projected square4/3. The seven profile norms leave only plus/minus sqrt(3q/4)*old10. On the center-facing sheet, profile0 has native7 product -1/sqrt6 in(-1,0); profiles2..6 have positive product; profile1 is old10. The opposite sheet has positive full center product, including white-profile offset, as separately checked here. This is a fixed-scaffold, seven-profile obstruction, not an arbitrary-profile theorem or a new exact normal certificate.',
        numerical_scope='Saved original intervals check the native interface; scaling, center signs and changed-R products are numerical. No interval reconstruction of a new wall is claimed.')


def bank(s,base,ids):
    frame=[7,55,10];duties=[(0,1)];bonds=[2]*9
    m=T.G.relation(s.dot(base[7],s.ps(7,{}),base[55],s.ps(55,{})))
    assert isinstance(m,int);bonds[1]=m
    dom=T.E.recipes(tuple(s.ps(i,{}) for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms)) for pi,ms,high in dom if not high and all(2<=z<=5 for z in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
        if T.P.corner3(s.ps(7,{}),s.ps(55,{}),pi,m,ms[0],ms[1]) is not None]
    assert sorted(recipes)==sorted(direct)
    records=[];caps=[];hist=Counter();digest=sha256();seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]] for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([base[i] for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(records),profile=pi,labels=ms,norm=info,branches=[]);records.append(row)
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage='native10_norm',recipe=row['index'],norm=info));continue
        for sign,v in roots:
            vals,bad,cp=s.native(v,pi,base,{})
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@base[i])-q) for i,q in zip(frame,rhs)))
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad])
            row['branches'].append(br);hist[str(tuple(bad))]+=1
            digest.update(json.dumps([row['index'],sign,pi,ms,br['native_products'],cp],separators=(',',':')).encode())
            if err>2e-5:
                s.unresolved.append(dict(stage='native10_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nd={**base,67:v};tt={67:pi}
            assert s.corner([7,55,67],nd,tt) is not None
            before=len(s.unresolved);sid=s.score(nd,tt,ids+[67],'native10_frame')
            raw=s.ledgers[sid] if sid is not None else s.unresolved[before]['boundary']
            faces=[H.face(raw,i) for i in (8,10,11,33,44,55,63,64,65,66,67)]
            r=dict(index=len(caps),recipe=row['index'],branch=sign,profile=pi,labels=ms,normal=v.tolist(),score=sid,
                cone_diagnostic=None if sid is not None else s.unresolved[before:],whole_face_records=faces,
                restores_old30=pi==s.ps(30,{}) and float(np.max(abs(v-s.old[30])))<2e-6,
                finite_rays=raw['counts'].get('finite',0),open_ends=raw['open_ends'],
                compact10=faces[1]['all_rays_finite'],compact8=faces[0]['all_rays_finite'],
                exterior_supports=[r['active'] for r in raw['rays'] if r['kind']=='exterior'],
                open_edges=[e for e in raw['edges'] if e['kind'] in ('one_open','two_open')])
            br['cap']=r['index'];caps.append(r)
    return dict(frame=frame,duties=duties,low_labels=[2,3,4,5],recipes=records,caps=caps,
        symbolic_high_families_not_tested=sum(bool(h) for _,_,h in dom),branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def run():
    seed=json.loads(SEED.read_text())
    for p,h in seed['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s)
    pairs=list(combinations([7,55,10],2))+list(combinations([8,11,33,44],2))+[(10,i) for i in [7,8,11,33,44]]
    nativechecks=intervals(s,pairs);pol=polar(s,seed)
    debit=T.H.deletion(s.primary,s.cur,[30]);assert debit['lost_finite_count']==13 and debit['compact_faces_by_dimension']=={2:10,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    base={i:v for i,v in s.old.items() if i!=30};ids=[i for i in s.cur['actual_sides'] if i!=30]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');assert control is not None
    deletion=s.score(base,{},ids,'delete30')
    b=bank(s,base,ids)
    assert sum(c['restores_old30'] for c in b['caps'])==1
    paths=[Path(__file__),SEED,ROOT/'results/fresh_joint10_30_65_compact_pair_20261018.md',Path(D.__file__),Path(K.__file__),Path(H.__file__),
        Path(T.__file__),Path(T.H.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.S.__file__),Path(T.S.C.__file__),
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_EXTERIOR_POLAR_AND_NATIVE10_FRAME',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        decision_record='A right cap on the persistent8/11/33/44 exterior is forced back to removedold10 within the seven profiles and conflicts with each changedR. Therefore prescribe the missing native10 product directly in frame7/55/10; old33 and44 are now derived gates. This is different from frozen7/55/33 and10/44/7 banks. No finite7/10/new30 corner is imposed.',
        old_frame_interface=interface,original_native_checks=nativechecks,polar=pol,deleted_walls=[30],retained_old_walls=62,
        exact_deletion_debit=debit,control_score=control,deletion_score=deletion,bank=b,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(old_normals={str(i):v.tolist() for i,v in s.old.items()},center=s.center.tolist(),
            profiles=[dict(contacts=p['contacts'],square=p['square'],center_offset=p['center_offset'],native_old_offsets=sorted(zip(s.labels,p['offsets']))) for p in s.profiles]),
        scope='One changed bounded low-label frame, seven inherited profiles, both norm branches, all62 retained old walls including10/33/44 and full profile offsets/common-center gates. Every admitted cap has complete numericalS6 cone versus PD-incidence checks and11 whole-face records; provisional cones stay unaccepted. Exact13/10/1 oldsupport debit. No high-label/continuous completeness, exact newnormal, ambientfinite census, finite-volume or newdimension certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:assert not out.exists();out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],'recipes',len(d['bank']['recipes']),'branches',sum(len(r['branches']) for r in d['bank']['recipes']),'caps',len(d['bank']['caps']),'unresolved',len(d['unresolved']),flush=True)
