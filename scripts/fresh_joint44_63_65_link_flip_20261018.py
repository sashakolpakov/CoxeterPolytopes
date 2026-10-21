#!/usr/bin/env python3
"""Move63 and flip the dotted44/65 seam after the complete two-bank preflight.

The bounded domain is21 explicitly selected changed prefixes, one closest
pair in each L white profile for each of three changed J recipes. The
replacement63 has genuine8/J and8/L corner duties; J/L is a nonedge.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_joint44_65_compact_chain_20261018 as P
T,H,ROOT=P.T,P.H,P.ROOT
TAG='fresh_joint44_63_65_link_flip_20261018'
DELETED=(30,44,63,65)
LOW=(2,3,4,5)
WALLS=[7,8,10,11,33,52,55,64,66,67,68,69,70]


def topology():
    fs=[[8,10,67],[10,67,70],[67,70,69],[70,64,69],[64,68,69],[8,68,69],[8,67,69],[8,52,68],[52,64,68]]
    bc=[8,10,70,64,52];es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2));links=[]
    assert len(set(sum(fs,[])))==8 and len(es)==16 and len(fs)==9 and (67,68) not in es
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2}
    for v in sorted(set(sum(fs,[]))):
        ee=[tuple(i for i in f if i!=v) for f in fs if v in f];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            ss=seen|{a for e in ee if seen.intersection(e) for a in e}
            if ss==seen:break
            seen=ss
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=8,edges=16,faces=9,links=links)

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
                **P.B.B.score(s,nn,tt,ids+[newid],stage,walls))
            br['cap']=rec['index'];caps.append(rec)
    return dict(frame=frame,duties=duties,frame_profiles=[s.ps(i,types) for i in frame],
        duty_bonds=bonds,symbolic_high_families=[list(r) for r in symbolic if r[2]],
        complete_low_recipe_count=len(recipes),recipes=rows,caps=caps,
        branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def reconstruct(s,frame,profile,labels,branch,nd,types,saved):
    rhs=[-T.B.theta(m)-(s.profiles[profile]['offsets'][s.index[i]] if i in s.index else s.off(profile,types[i])) for i,m in zip(frame,labels)]
    roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[profile]['square'])
    assert not info['status'].startswith('unresolved')
    found=[v for b,v in roots if b==branch];assert len(found)==1
    v=found[0];assert np.max(abs(v-np.array(saved)))<2e-6
    vals,bad,cp=s.native(v,profile,nd,types);assert not bad and cp<-1e-8
    return v


def run():
    priorpath=ROOT/f'results/{P.TAG}.json';prior=json.loads(priorpath.read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    top=topology();s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s)
    debit=T.H.deletion(s.primary,s.cur,DELETED)
    assert debit==prior['additional63_deletion_debit'] or json.loads(json.dumps(debit))==prior['additional63_deletion_debit']
    assert debit['lost_finite_count']==19 and debit['compact_faces_by_dimension']=={2:12,3:1}
    base={i:v for i,v in s.old.items() if i not in DELETED};ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    score=lambda nd,tt,ii,stage:P.B.B.score(s,nd,tt,ii,stage,[i for i in WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'control63');assert control['score'] is not None and control['open_ends']==10
    deletion=score(base,{},ids,'delete30_44_63_65')
    chosen=[]
    for j in prior['J_caps']:
        if j['restores_old44']:continue
        for pi in range(7):
            choices=[x for x in prior['joins'] if x['J']==j['index'] and not x['blockers'] and prior['L_bank']['caps'][x['L']]['profile']==pi]
            assert choices
            x=min(choices,key=lambda x:(abs(x['native_product']+2),x['L']))
            assert x['relation']=='dotted';chosen.append(x)
    assert len(chosen)==21 and len({(x['J'],x['L']) for x in chosen})==21
    records=[];r_banks=[];accepted=[]
    for sel in chosen:
        j=prior['J_caps'][sel['J']];l=prior['L_bank']['caps'][sel['L']]
        # The first J recipe still refers to the deleted old63 anchor as an
        # explicit geometric seed; it does not retain that wall as a facet.
        rhs=[0.,-T.B.theta(j['label10']),0.]
        roots,info=T.N.norm_roots([s.old[i] for i in (8,10,63)],rhs,s.profiles[0]['square'])
        jv=next(v for b,v in roots if b==j['branch']);assert np.max(abs(jv-np.array(j['normal'])))<2e-6
        lv=reconstruct(s,[8,52,64],l['profile'],l['labels'],l['branch'],base,{},l['normal'])
        nd={**base,67:jv,68:lv};tt={67:0,68:l['profile']}
        for i in (67,68):
            _,bad,cp=s.native(nd[i],tt[i],{a:v for a,v in nd.items() if a!=i},tt);assert not bad and cp<-1e-8
        prefix=score(nd,tt,ids+[67,68],'selected_JL')
        xb=bank(s,[8,67,68],[(0,1),(0,2)],nd,tt,ids+[67,68],69,'new63_flip')
        rec=dict(index=len(records),selected_join=sel,J=j['index'],L=l['index'],J_profile=0,L_profile=l['profile'],prefix=prefix,bank=xb);records.append(rec)
        for x in xb['caps']:
            nn={**nd,69:np.array(x['normal'])};types={**tt,69:x['profile']}
            x['corner_L64']=s.corner([64,68,69],nn,types)
            if x['corner_L64'] is None:continue
            rb=bank(s,[10,33,67],[(0,1),(0,2)],nn,types,ids+[67,68,69],70,'new30_closure')
            r_banks.append(dict(prefix=rec['index'],X=x['index'],bank=rb))
            for r in rb['caps']:
                full={**nn,70:np.array(r['normal'])};ty={**types,70:r['profile']}
                cs=[s.corner(f,full,ty) for f in top['triangles']]
                extra=[[7,55,70],[55,59,70],[52,64,70],[10,33,70]];ec=[s.corner(f,full,ty) for f in extra]
                faces={f['wall']:f for f in r['whole_faces']};whole={str(w):faces[w]['all_rays_finite'] for w in (8,10,55,64,67,68,69)}
                exactcaps={str(w):faces[w]['all_rays_finite'] and sorted(sorted(q['active']) for q in faces[w]['rays'])==sorted(sorted(f) for f in top['triangles'] if w in f) for w in (67,68,69)}
                accepted.append(dict(prefix=rec['index'],X=x['index'],R=r['index'],corner_types=cs,extra_duties=extra,extra_corner_types=ec,all_corners=all(c is not None for c in cs+ec),
                    whole_compact_faces=whole,exact_whole_new_caps=exactcaps,score=r['score'],clean=r['score'] is not None,open_ends=r['open_ends'],finite_rays=r['finite_rays'],whole_faces=r['whole_faces']))
        print('PREFIX',rec['index'],'J',j['index'],'L',l['index'],'profile',l['profile'],'Xrecipes',len(xb['recipes']),'Xnative',len(xb['caps']),flush=True)
    summary=dict(selected_prefixes=len(chosen),X_recipes=sum(len(r['bank']['recipes']) for r in records),X_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in records),
        X_native=sum(len(r['bank']['caps']) for r in records),X_L64_finite=sum(x['corner_L64'] is not None for r in records for x in r['bank']['caps']),
        R_banks=len(r_banks),R_recipes=sum(len(r['bank']['recipes']) for r in r_banks),R_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in r_banks),
        native_quadruples=len(accepted),clean_quadruples=sum(r['clean'] for r in accepted),completed_disks=sum(r['clean'] and r['all_corners'] and all(r['exact_whole_new_caps'].values()) for r in accepted),
        fewer_than_ten_ends=sum(r['clean'] and r['open_ends']<10 for r in accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(P.__file__),priorpath,Path(H.__file__),Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY}
    paths.update(ROOT/p for p in prior['sources'])
    return dict(status='BOUNDED_JOINT44_63_65_LINK_FLIP',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,deletion=deletion,
        selection_rule='For each changed J recipe and each of seven L profiles choose the prior native pair with product closest to minus2, ties by L index. This is21 selected prefixes, not all prior345 pairs or all possible J/L normals.',
        selected=chosen,records=records,R_banks=r_banks,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='Changed disk removes J/L adjacency and replaces63 by X with link8-J-R-64-L-8. Twenty-one explicit native prefixes are reconstructed; each X and reached R bank covers its complete declared seven-profile2..5 cube and both norm branches. All retained old/mutual native and strict common-center gates precede full numerical cone/PD and whole-side tests. Provisional cones remain diagnostic. Exact19/12/1 original certificate loss; no all-label/continuous completeness, exact new-normal or ambient H9 volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
