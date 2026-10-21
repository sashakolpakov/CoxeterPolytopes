#!/usr/bin/env python3
"""Add a new boundary connector S72 after the oblique R/reverse M stage.

Required corners52/L/S and L/M/S; M/R/S is checked afterward.
No52/M,52/R or52/X edge is imposed. All native/center/cone gates retained.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_moving64_oblique_r_closure_20261018 as O

ROOT=O.ROOT;C=O.C;T=O.T;Q=O.Q;H=O.H
TAG='fresh_oblique52_boundary_connector_20261018';SID=72;LOW=(2,3,4,5)
WALLS=C.WALLS+[SID]


def topology():
    old=C.topology();fs=[f for f in old['triangles'] if set(f) not in ({52,71,68},{52,70,71})]
    fs += [[52,68,SID],[68,71,SID],[71,70,SID]];bc=[8,10,70,SID,52]
    es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2));vs=sorted(set(sum(fs,[])))
    assert (len(vs),len(es),len(fs))==(9,19,11)
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2} and all(e not in es for e in [(52,71),(52,70),(52,69),(67,68)])
    links=[]
    for v in vs:
        ee=[tuple(i for i in f if i!=v) for f in fs if v in f];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            nxt=seen|{i for e in ee if seen.intersection(e) for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=9,edges=19,faces=11,links=links)


def reconstruct(s,old,cdata,flip,seed):
    base={i:v for i,v in s.old.items() if i not in C.DELETED};out=[];groups={}
    for pointed in old['reverse_accepted']:
        rev=old['reverse_records'][pointed['reverse']];sel=rev['selection'];prior=cdata['records'][sel['prefix']]
        fr=flip['records'][prior['selection']['prefix']]
        j=seed['J_caps'][fr['J']];l=seed['L_bank']['caps'][fr['L']];x=fr['bank']['caps'][prior['selection']['X']]
        nd=dict(base);tt={};checks=[]
        def step(i,frame,pi,labels,branch,normal,anchors=False):
            v,ch=C.solve_saved(s,frame,pi,labels,branch,s.old if anchors else nd,{} if anchors else tt,normal)
            nd[i]=v;tt[i]=pi;checks.append(dict(wall=i,**ch))
        step(67,[8,10,63],0,[2,3,2],j['branch'],j['normal'],True)
        step(68,[8,52,64],l['profile'],l['labels'],l['branch'],l['normal'],True)
        step(69,[8,67,68],x['profile'],x['labels'],x['branch'],x['normal'])
        step(70,[10,67,69],sel['profile'],sel['labels'],sel['branch'],sel['normal'])
        m=rev['bank']['caps'][pointed['M']]
        step(71,[68,69,70],m['profile'],m['labels'],m['branch'],m['normal'])
        for i in tt:
            _,bad,cp=s.native(nd[i],tt[i],{k:v for k,v in nd.items() if k!=i},tt);assert not bad and cp<-1e-8
        # Equality of complete reconstruction instructions, not rounded-normal
        # proximity, is the criterion for dropping a duplicate prefix.
        key=json.dumps([[z['wall'],z['frame'],z['profile'],z['labels'],z['branch']] for z in checks],separators=(',',':'))
        origin=dict(reverse=pointed['reverse'],M=pointed['M'])
        if key in groups:
            rec=out[groups[key]];assert tt==rec['types']
            assert all(np.max(abs(nd[i]-rec['normals'][i]))<1e-9 for i in tt)
            rec['origins'].append(origin);continue
        groups[key]=len(out)
        out.append(dict(index=len(out),origins=[origin],reconstruction=checks,recipe_key=key,
            normals=nd,types=tt,saved_before=pointed))
    assert sum(len(r['origins']) for r in out)==52
    return out


def recipes(s,nd,tt):
    frame=[52,68,71];duties=[(0,1),(1,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],tt),nd[frame[b]],s.ps(frame[b],tt)))
        assert isinstance(m,int);bonds[3*a+b]=m
    sy=T.E.recipes(tuple(s.ps(i,tt) for i in frame),tuple(bonds),tuple(duties))
    rows=[(pi,list(ms)) for pi,ms,hi in sy if not hi and all(m in LOW for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(rows)==sorted(direct)
    return frame,duties,bonds,sy,rows


def run(pilot=False):
    names=(O.TAG,C.TAG,Q.TAG,Q.P.TAG);docs={n:json.loads((ROOT/f'results/{n}.json').read_text()) for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old,cdata,flip,seed=[docs[n] for n in names];assert old['summary']['reverse_native_M']==52
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert all(i not in s.index for i in range(67,73))
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=topology();ids=[i for i in s.cur['actual_sides'] if i not in C.DELETED]+[67,68,69,70,71]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control');assert control['score'] is not None and control['open_ends']==10
    prefixes=reconstruct(s,old,cdata,flip,seed);chosen=prefixes[:3] if pilot else prefixes
    records=[];accepted=[]
    for pre in chosen:
        nd,tt=pre['normals'],pre['types'];frame,duties,bonds,sy,rr=recipes(s,nd,tt)
        before=score(nd,tt,ids,'reconstructed_oblique_R_M');assert before['score'] is not None and before['finite_rays']==21 and before['open_ends']==13
        row=dict(index=pre['index'],origins=pre['origins'],recipe_key=pre['recipe_key'],reconstruction=pre['reconstruction'],
            prefix_normals={str(i):nd[i].tolist() for i in tt},prefix_profiles=tt,before=before,frame=frame,duties=duties,
            frame_profiles=[s.ps(i,tt) for i in frame],duty_bonds=bonds,symbolic_high_families=[list(z) for z in sy if z[2]],
            recipes=[],caps=[]);records.append(row);seen=set();digest=sha256();hist=Counter()
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
                seen.add(key);full={**nd,SID:v};ft={**tt,SID:pi}
                for i in ft:
                    _,fail,center=s.native(full[i],ft[i],{j:w for j,w in full.items() if j!=i},ft);assert not fail and center<-1e-8
                assert all(s.corner([frame[a],frame[b],SID],full,ft) is not None for a,b in duties)
                sc=score(full,ft,ids+[SID],'new52_boundary_connector');faces={f['wall']:f for f in sc['whole_faces']}
                cs=[s.corner(f,full,ft) for f in top['triangles']]
                caps={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
                cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,branch=branch,
                    normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,
                    corner_types=cs,all_disk_corners=all(z is not None for z in cs),deferred_M_R_S=s.corner([71,70,SID],full,ft),
                    exact_whole_interior_caps=caps,whole_compact_faces={str(w):faces[w]['all_rays_finite'] for w in faces},**sc)
                b['cap']=cap['index'];row['caps'].append(cap);accepted.append(dict(prefix=pre['index'],**cap))
        row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest()
        print('PREFIX',pre['index'],'recipes',len(rr),'nativeS',len(row['caps']),flush=True)
    summary=dict(frozen_pointed_configurations=52,distinct_recipe_prefixes=len(prefixes),tested_prefixes=len(records),
        recipes=sum(len(r['recipes']) for r in records),branches=sum(len(q['branches']) for r in records for q in r['recipes']),
        native_S=len(accepted),clean_S=sum(c['score'] is not None for c in accepted),
        completed_disks=sum(c['score'] is not None and c['all_disk_corners'] and all(c['exact_whole_interior_caps'].values()) for c in accepted),
        fewer_than_ten_ends=sum(c['score'] is not None and c['open_ends']<10 for c in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(O.__file__),Path(C.__file__),Path(Q.__file__),Path(H.__file__),Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__)}
    paths.update(ROOT/p for p in old['sources']);paths.update(ROOT/f'{f}/{n}.{e}' for n in names for f,e in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PILOT_OBLIQUE52_BOUNDARY_CONNECTOR' if pilot else 'BOUNDED_OBLIQUE52_BOUNDARY_CONNECTOR',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(C.DELETED),retained_old_walls=len(s.old)-len(C.DELETED)+6,new_id=SID,
        exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,
        deduplication='Identical complete five-step reconstruction instructions; numerical proximity alone never removes a prefix.',
        prefix_origins=[dict(index=p['index'],origins=p['origins'],recipe_key=p['recipe_key']) for p in prefixes],
        records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='Bounded seven-profile low-label2..5 connector with both real norm branches at the frozen52 native obliqueR/reverseM configurations. Newboundary8/10/R/S/52; finite52/L/S,L/M/S duties thenM/R/S andall disk/whole-face gates. All58retainedoldwalls plusmutualnative/strictcenter before fullnumericalcone/PD. Exactoriginal19/12/1certificate debit unchanged. No high-label/all-prefix/continuous completeness, exactnewnormal or ambientH9 volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--pilot',action='store_true');a=ap.parse_args();assert not(a.verify and a.pilot)
    p=ROOT/f'results/{TAG}.json';assert a.verify or a.pilot or not p.exists()
    out=json.loads(json.dumps(run(a.pilot)))
    if a.verify:assert out==json.loads(p.read_text())
    elif a.pilot:print(json.dumps(out['summary']),flush=True)
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
