#!/usr/bin/env python3
"""Move64 after the changed44/63/65 disk, with a fresh native profile identity.

Nine selected profile classes are reconstructed before native-first M71
banks in52/L/X. Old63 and64 occur only as explicit seed anchors, never as
retained facets. Every reached R70 closure restores all actual native gates.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_joint44_63_65_link_flip_20261018 as Q
T,H,ROOT=Q.T,Q.H,Q.ROOT
TAG='fresh_joint64_native_compact_closure_20261018'
DELETED=(30,44,63,64,65)
LOW=(2,3,4,5)
WALLS=[7,8,10,11,33,52,55,66,67,68,69,70,71]


def topology():
    fs=[[8,10,67],[10,67,70],[67,70,69],[70,71,69],[71,68,69],[8,68,69],[8,67,69],[8,52,68],[52,71,68],[52,70,71]]
    bc=[8,10,70,52];es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2));links=[]
    assert len(set(sum(fs,[])))==8 and len(es)==17 and len(fs)==10
    assert (67,68) not in es and (52,69) not in es
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])} and set(es.values())=={1,2}
    for v in sorted(set(sum(fs,[]))):
        ee=[tuple(i for i in f if i!=v) for f in fs if v in f];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            ss=seen|{i for e in ee if seen.intersection(e) for i in e}
            if ss==seen:break
            seen=ss
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=8,edges=17,faces=10,links=links)


def solve_saved(s,frame,pi,labels,branch,nd,tt,saved):
    rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,tt[i])) for i,m in zip(frame,labels)]
    roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
    assert not info['status'].startswith('unresolved')
    vv=[v for b,v in roots if b==branch];assert len(vv)==1
    v=vv[0];assert np.max(abs(v-np.array(saved)))<2e-6
    return v,dict(frame=frame,profile=pi,labels=labels,branch=branch,norm=info,max_saved_difference=float(np.max(abs(v-np.array(saved)))))

def bank(s,frame,duties,nd,types,ids,newid,stage):
    bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],types),nd[frame[b]],s.ps(frame[b],types)))
        assert isinstance(m,int);bonds[3*a+b]=m
    symbolic=T.E.recipes(tuple(s.ps(i,types) for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms)) for pi,ms,hi in symbolic if not hi and all(m in LOW for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],types),s.ps(frame[b],types),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    rows=[];caps=[];seen=set();hist=Counter();digest=sha256()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,types[i])) for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(rows),profile=pi,labels=ms,norm=info,branches=[]);rows.append(row)
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage=stage,recipe=row['index'],norm=info));continue
        for sign,v in roots:
            vals,bad,cp=s.native(v,pi,nd,types)
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad])
            row['branches'].append(br);hist[str(tuple(bad))]+=1
            digest.update(json.dumps([row['index'],sign,pi,ms,br['native_products'],cp],separators=(',',':')).encode())
            if err>2e-5:
                s.unresolved.append(dict(stage=stage,recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nn={**nd,newid:v};tt={**types,newid:pi}
            assert all(s.corner([frame[a],frame[b],newid],nn,tt) is not None for a,b in duties)
            # Recheck all previously added normals against the new one too.
            for i in tt:
                _,fail,center=s.native(nn[i],tt[i],{j:w for j,w in nn.items() if j!=i},tt)
                assert not fail and center<-1e-8
            walls=[i for i in WALLS if i in nn]
            rec=dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,branch=sign,normal=v.tolist(),
                native_products=sorted(vals.items()),center_product=cp,norm_residual=err,
                **Q.P.B.B.score(s,nn,tt,ids+[newid],stage,walls))
            br['cap']=rec['index'];caps.append(rec)
    return dict(frame=frame,duties=duties,frame_profiles=[s.ps(i,types) for i in frame],
        duty_bonds=bonds,symbolic_high_families=[list(r) for r in symbolic if r[2]],
        complete_low_recipe_count=len(recipes),recipes=rows,caps=caps,
        branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def run():
    prior=json.loads((ROOT/f'results/{Q.TAG}.json').read_text());seeds=json.loads((ROOT/f'results/{Q.P.TAG}.json').read_text())
    for doc in (prior,seeds):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);assert all(i not in s.index for i in (67,68,69,70,71))
    top=topology();debit=T.H.deletion(s.primary,s.cur,DELETED);olddebit=prior['exact_deletion_debit']
    assert json.loads(json.dumps({k:v for k,v in debit.items() if k!='deleted_walls'}))=={k:v for k,v in olddebit.items() if k!='deleted_walls'}
    assert debit['lost_finite_count']==19 and debit['compact_faces_by_dimension']=={2:12,3:1}
    base={i:v for i,v in s.old.items() if i not in DELETED};ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'control63');assert control['score'] is not None and control['open_ends']==10
    # High label27 is measured in the original control; it is not discarded
    # merely because the subsequent replacement banks have a bounded low cube.
    pi=s.ps(64,{});rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]] for i,m in zip([52,65,63],[2,2,27])]
    roots,info=T.N.norm_roots([s.old[i] for i in (52,65,63)],rhs,s.profiles[pi]['square'])
    controls=[]
    for sign,v in roots:
        vals,bad,cp=s.native(v,pi,{i:w for i,w in s.old.items() if i!=64},{})
        controls.append(dict(branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],center_product=cp,restores_old64=bool(np.max(abs(v-s.old[64]))<2e-6)))
    assert sum(c['restores_old64'] and not c['blockers'] and c['center_product']<-1e-8 for c in controls)==1
    control64=dict(frame=[52,65,63],profile=pi,labels=[2,2,27],norm=info,branches=controls)
    deletion=score(base,{},ids,'delete30_44_63_64_65')
    classes={}
    for rec in prior['records']:
        if seeds['J_caps'][rec['J']]['label10']!=3:continue
        for x in rec['bank']['caps']:
            key=(rec['L_profile'],x['profile']);value=dict(x['native_products'])[64]
            cand=dict(prefix=rec['index'],X=x['index'],L_profile=rec['L_profile'],X_profile=x['profile'],old_X64=value)
            order=lambda z:(abs(z['old_X64']+2),z['prefix'],z['X'])
            if key not in classes or order(cand)<order(classes[key]):classes[key]=cand
    selection=[classes[k] for k in sorted(classes)];assert len(selection)==9
    records=[];r_banks=[];accepted=[]
    for sel in selection:
        rec=prior['records'][sel['prefix']];j=seeds['J_caps'][rec['J']];l=seeds['L_bank']['caps'][rec['L']];x=rec['bank']['caps'][sel['X']]
        jv,jcheck=solve_saved(s,[8,10,63],0,[2,3,2],j['branch'],s.old,{},j['normal'])
        lv,lcheck=solve_saved(s,[8,52,64],l['profile'],l['labels'],l['branch'],s.old,{},l['normal'])
        initial={**base,67:jv,68:lv};ty={67:0,68:l['profile']}
        xv,xcheck=solve_saved(s,[8,67,68],x['profile'],x['labels'],x['branch'],initial,ty,x['normal'])
        nd={**initial,69:xv};tt={**ty,69:x['profile']}
        for i in (67,68,69):
            _,bad,cp=s.native(nd[i],tt[i],{a:v for a,v in nd.items() if a!=i},tt);assert not bad and cp<-1e-8
        prefix=score(nd,tt,ids+[67,68,69],'selected_JLX_without64')
        mb=bank(s,[52,68,69],[(0,1),(1,2)],nd,tt,ids+[67,68,69],71,'new64')
        row=dict(index=len(records),selection=sel,J=rec['J'],L=rec['L'],reconstruction=[jcheck,lcheck,xcheck],prefix=prefix,bank=mb);records.append(row)
        for m in mb['caps']:
            nn={**nd,71:np.array(m['normal'])};types={**tt,71:m['profile']}
            rb=bank(s,[10,33,67],[(0,1),(0,2)],nn,types,ids+[67,68,69,71],70,'new30')
            r_banks.append(dict(prefix=row['index'],M=m['index'],bank=rb))
            for r in rb['caps']:
                full={**nn,70:np.array(r['normal'])};ft={**types,70:r['profile']}
                cs=[s.corner(f,full,ft) for f in top['triangles']];extra=[[7,55,70],[55,59,70],[10,33,70]];ec=[s.corner(f,full,ft) for f in extra]
                faces={f['wall']:f for f in r['whole_faces']};whole={str(w):faces[w]['all_rays_finite'] for w in (8,10,55,67,68,69,71)}
                exactcaps={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69,71)}
                accepted.append(dict(prefix=row['index'],M=m['index'],R=r['index'],corner_types=cs,extra_duties=extra,extra_corner_types=ec,all_corners=all(c is not None for c in cs+ec),
                    whole_compact_faces=whole,exact_whole_new_caps=exactcaps,score=r['score'],clean=r['score'] is not None,open_ends=r['open_ends'],finite_rays=r['finite_rays'],whole_faces=r['whole_faces']))
        print('PREFIX',row['index'],'L/Xprofiles',sel['L_profile'],sel['X_profile'],'Mrecipes',len(mb['recipes']),'Mnative',len(mb['caps']),flush=True)
    summary=dict(selected_prefixes=len(selection),M_recipes=sum(len(r['bank']['recipes']) for r in records),M_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in records),
        M_native=sum(len(r['bank']['caps']) for r in records),R_banks=len(r_banks),R_recipes=sum(len(r['bank']['recipes']) for r in r_banks),R_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in r_banks),
        native_five_caps=len(accepted),clean_five_caps=sum(r['clean'] for r in accepted),completed_disks=sum(r['clean'] and r['all_corners'] and all(r['exact_whole_new_caps'].values()) for r in accepted),
        fewer_than_ten_ends=sum(r['clean'] and r['open_ends']<10 for r in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(Q.__file__),Path(Q.P.__file__),Path(H.__file__),Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY}
    paths.update(ROOT/f'{d}/{tag}.{e}' for tag in (Q.TAG,Q.P.TAG) for d,e in (('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/p for p in prior['sources'])
    return dict(status='BOUNDED_JOINT64_NATIVE_COMPACT_CLOSURE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,replacement64_id=71,replacement30_id=70,exact_deletion_debit=debit,old_frame_interface=interface,topology=top,
        control=control,original64_recipe=control64,deletion=deletion,selection_rule='FixJ10label3; oneclosestoldX64 representative per realizedL/Xwhite-profilepair, tiesbypriorprefixandXindices. Nineexplicitprefixes; notall270Xcaps.',
        selected=selection,records=records,R_banks=r_banks,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='Native-first replacement64 usesfreshID71 and genuine52/L/M,L/X/M duties withno52/Xedge. Nineexplicitfrozenrecipeprefixes, eachcomplete seven-profile2..5Mbank andallreachedRbanks. All58retainedoldwalls plusmutualnative/strictcenter precede fullnumericalcone/PD/whole-side acceptance. Exact19/12/1originalcertificate debit; provisionalrecords notchambers. No all-label/continuouscompleteness, exactnewnormal or ambientH9finitevolumecertificate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
