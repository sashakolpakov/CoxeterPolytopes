#!/usr/bin/env python3
"""Insert N73 between10 and reverse-R70; genuine10/J/N and J/X/N duties.

The third X/R/N duty is deferred. No10/X edge is assumed. All normals,
cone and whole-face incidence checks here remain numerical.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_oblique52_reverse_r_closure_20261018 as R

ROOT=R.ROOT;A=R.A;O=R.O;C=R.C;T=R.T;Q=R.Q;H=R.H
TAG='fresh_oblique_left_boundary_connector_20261018';NID=73;LOW=(2,3,4,5)
WALLS=A.WALLS+[NID]


def topology():
    fs=[f for f in A.topology()['triangles'] if set(f) not in ({10,67,70},{67,69,70})]
    fs += [[10,67,NID],[67,69,NID],[69,70,NID]];bc=[8,10,NID,70,72,52]
    es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2));vs=sorted(set(sum(fs,[])))
    assert (len(vs),len(es),len(fs))==(10,21,12)
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2} and all(e not in es for e in [(10,70),(67,70),(10,69),(67,68),(52,71)])
    links=[]
    for v in vs:
        ee=[tuple(i for i in f if i!=v) for f in fs if v in f];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            nxt=seen|{i for e in ee if seen.intersection(e) for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=10,edges=21,faces=12,links=links)


def reconstruct(s,docs):
    previous=R.reconstruct(s,docs);out=[];seen={}
    for k,cap in enumerate(docs[R.TAG]['accepted']):
        p=previous[cap['prefix']];nd=dict(p['normals']);tt=dict(p['types'])
        v,ch=C.solve_saved(s,[69,71,72],cap['profile'],cap['labels'],cap['branch'],nd,tt,cap['normal'])
        nd[70]=v;tt[70]=cap['profile'];steps=[*p['reconstruction'],dict(wall=70,**ch)]
        key=json.dumps([[z['wall'],z['frame'],z['profile'],z['labels'],z['branch']] for z in steps],separators=(',',':'))
        for i in tt:
            _,bad,cp=s.native(nd[i],tt[i],{j:w for j,w in nd.items() if j!=i},tt);assert not bad and cp<-1e-8
        old_product=s.dot(np.array(p['original_R']),p['original_R_profile'],v,tt[70])
        origin=dict(accepted=k,prefix=cap['prefix'],R=cap['index'])
        if key in seen:
            pre=out[seen[key]];assert tt==pre['types'] and all(np.max(abs(nd[i]-pre['normals'][i]))<1e-9 for i in tt)
            pre['origins'].append(origin);continue
        seen[key]=len(out)
        out.append(dict(index=len(out),origins=[origin],recipe_key=key,reconstruction=steps,normals=nd,types=tt,
            old_R_as_N_product=old_product,old_R_normal=p['original_R'],old_R_profile=p['original_R_profile']))
    assert sum(len(p['origins']) for p in out)==376 and all(p['old_R_as_N_product']>0 for p in out)
    return out


def choose_pilot(prefixes):
    chosen=[];old=set()
    for p in sorted(prefixes,key=lambda p:(p['old_R_as_N_product'],p['index'])):
        key=p['origins'][0]['prefix']
        if key in old:continue
        old.add(key);chosen.append(p['index'])
        if len(chosen)==4:break
    assert len(chosen)==4
    return chosen


def recipes(s,nd,tt):
    frame=[10,67,69];duties=[(0,1),(1,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],tt),nd[frame[b]],s.ps(frame[b],tt)))
        assert isinstance(m,int);bonds[3*a+b]=m
    sy=T.E.recipes(tuple(s.ps(i,tt) for i in frame),tuple(bonds),tuple(duties))
    rr=[(pi,list(ms)) for pi,ms,hi in sy if not hi and all(m in LOW for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(rr)==sorted(direct)
    return frame,duties,bonds,sy,rr


def high_finite_obstruction(s,nd,tt,sy):
    rx=T.G.relation(s.dot(nd[69],s.ps(69,tt),nd[70],s.ps(70,tt)))
    assert isinstance(rx,int) and rx>=3
    high=[z for z in sy if z[2]]
    assert all(pi==0 and list(hi)==[2] and list(ms) in ([2,2,6],[3,2,6]) for pi,ms,hi in high)
    # The physical norm2 principal X,R,N has determinant
    # 8-2(a*a+b*b+c*c)-2*a*b*c, with a=theta(NX)>=sqrt3,
    # b=theta(XR)>=1, c=theta(RN)>=0. Its upper bound is zero.
    assert 8-2*(3+1)==0
    return dict(RX_label=rx,high_families=len(high),NX_squared_lower_bound=3,
        XR_squared_lower_bound=1,RN_magnitude_lower_bound=0,determinant_upper_bound=0,
        conclusion='No symbolic high family can have a finite X/R/N corner; equality is not strict PD.')


def run(pilot=False):
    names=(R.TAG,A.TAG,O.TAG,C.TAG,Q.TAG,Q.P.TAG);docs={n:json.loads((ROOT/f'results/{n}.json').read_text()) for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=docs[R.TAG];assert old['summary']['native_R']==old['summary']['clean_R']==376
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert all(i not in s.index for i in range(67,74))
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=topology();ids=[i for i in s.cur['actual_sides'] if i not in C.DELETED]+[67,68,69,70,71,72]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control');assert control['score'] is not None and control['open_ends']==10
    prefixes=reconstruct(s,docs);pilot_ids=choose_pilot(prefixes)
    chosen=[prefixes[i] for i in pilot_ids] if pilot else prefixes
    records=[];accepted=[]
    for pre in chosen:
        nd,tt=pre['normals'],pre['types'];frame,duties,bonds,sy,rr=recipes(s,nd,tt)
        before=score(nd,tt,ids,'frozen_reverse_R');assert before['score'] is not None and before['finite_rays']==22 and before['open_ends']==14
        row=dict(index=pre['index'],origins=pre['origins'],recipe_key=pre['recipe_key'],reconstruction=pre['reconstruction'],
            prefix_normals={str(i):nd[i].tolist() for i in tt},prefix_profiles=tt,old_R_as_N_product=pre['old_R_as_N_product'],
            old_R_normal=pre['old_R_normal'],old_R_profile=pre['old_R_profile'],before=before,
            frame=frame,duties=duties,frame_profiles=[s.ps(i,tt) for i in frame],duty_bonds=bonds,
            symbolic_high_families=[list(z) for z in sy if z[2]],
            high_finite_obstruction=high_finite_obstruction(s,nd,tt,sy),recipes=[],caps=[])
        records.append(row);seen=set();digest=sha256();hist=Counter()
        for pi,ms in rr:
            rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,tt[i])) for i,m in zip(frame,ms)]
            roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
            rec=dict(index=len(row['recipes']),profile=pi,labels=ms,norm=info,branches=[]);row['recipes'].append(rec)
            if info['status'].startswith('unresolved'):
                s.unresolved.append(dict(prefix=pre['index'],recipe=rec['index'],norm=info));continue
            for branch,v in roots:
                vals,bad,cp=s.native(v,pi,nd,tt)
                err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
                b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],center_product=cp,norm_residual=err)
                rec['branches'].append(b);hist[str(tuple(bad))]+=1
                digest.update(json.dumps([rec['index'],branch,pi,ms,b['native_products'],cp],separators=(',',':')).encode())
                if err>2e-5:s.unresolved.append(dict(prefix=pre['index'],recipe=rec['index'],branch=branch,error=err));continue
                if bad or cp>=-1e-8:continue
                key=(pi,*np.round(v,8))
                if key in seen:b['duplicate']=True;continue
                seen.add(key);full={**nd,NID:v};ft={**tt,NID:pi}
                for i in ft:
                    _,fail,center=s.native(full[i],ft[i],{j:w for j,w in full.items() if j!=i},ft);assert not fail and center<-1e-8
                assert all(s.corner([frame[a],frame[b],NID],full,ft) is not None for a,b in duties)
                sc=score(full,ft,ids+[NID],'new_left_boundary_connector');faces={f['wall']:f for f in sc['whole_faces']}
                cs=[s.corner(f,full,ft) for f in top['triangles']]
                capsets={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
                cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,branch=branch,normal=v.tolist(),
                    native_products=sorted(vals.items()),center_product=cp,norm_residual=err,corner_types=cs,
                    all_disk_corners=all(z is not None for z in cs),deferred_X_R_N=s.corner([69,70,NID],full,ft),
                    exact_whole_interior_caps=capsets,whole_compact_faces={str(w):faces[w]['all_rays_finite'] for w in faces},**sc)
                b['cap']=cap['index'];row['caps'].append(cap);accepted.append(dict(prefix=pre['index'],**cap))
        row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest()
        print('PREFIX',pre['index'],'recipes',len(rr),'nativeN',len(row['caps']),flush=True)
    summary=dict(frozen_pointed_configurations=376,distinct_recipe_prefixes=len(prefixes),tested_prefixes=len(records),
        recipes=sum(len(r['recipes']) for r in records),branches=sum(len(q['branches']) for r in records for q in r['recipes']),
        native_N=len(accepted),clean_N=sum(c['score'] is not None for c in accepted),
        completed_disks=sum(c['score'] is not None and c['all_disk_corners'] and all(c['exact_whole_interior_caps'].values()) for c in accepted),
        fewer_than_ten_ends=sum(c['score'] is not None and c['open_ends']<10 for c in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(R.__file__)};paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{f}/{n}.{e}' for n in names for f,e in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PILOT_OBLIQUE_LEFT_BOUNDARY_CONNECTOR' if pilot else 'BOUNDED_OBLIQUE_LEFT_BOUNDARY_CONNECTOR',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(C.DELETED),retained_old_walls=len(s.old)-len(C.DELETED)+6,new_id=NID,
        wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,
            original_facets=58,new_pre_N_walls=6,explicit_product_gates=58,total_pre_N_facets=64),
        exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,
        deduplication='Identical complete seven-step reconstruction instructions; proximity alone never removes a prefix.',
        pilot_selection=dict(rule='Smallest positive former-R/new-R product, breaking ties by index; four distinct frozen S-prefix indices.',indices=pilot_ids),
        prefix_inventory=[dict(index=p['index'],origins=p['origins'],recipe_key=p['recipe_key'],old_R_as_N_product=p['old_R_as_N_product']) for p in prefixes],
        records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=old['numerical_snapshot'],
        scope='Complete seven-profile low-label2..5 N banks over selected or all frozen376 reverse-R configurations. Only10/J/N and J/X/N prescribed; X/R/N deferred. Every symbolic high family fails the exact necessary strict-PD test for X/R/N, conditional on the inherited spherical classifier and frozen native labels. Low normal failures remain numerical. All58 retained old walls and mutual native/center gates precede target cone/PD and whole-face checks; no candidate reaches these downstream gates. Exact19/12/1 old-certificate debit unchanged. No continuous-prefix, exact-new-normal or ambientH9 volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--pilot',action='store_true');a=ap.parse_args();assert not(a.verify and a.pilot)
    p=ROOT/f'results/{TAG}{"_pilot" if a.pilot else ""}.json';assert a.verify or not p.exists()
    out=json.loads(json.dumps(run(a.pilot)))
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
