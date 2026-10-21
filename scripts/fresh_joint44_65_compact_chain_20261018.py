#!/usr/bin/env python3
"""Replace44/65 together, retain10/33, and reconstruct30 if the new link closes.

Four native-first44 recipes change the obstructing old44/65 product. The
replacement65 bank is complete only for seven profiles and labels2..5 in
its stated two-corner frame. Exact original support debit is kept separate
from numerical new normals and target-face cone ledgers.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_moving33_compact10_bridge_20261018 as B

T,H,ROOT=B.T,B.H,B.ROOT
TAG='fresh_joint44_65_compact_chain_20261018'
LOW=(2,3,4,5)
DELETED=(30,44,65)
WALLS=[7,8,10,11,33,52,55,63,64,66,67,68,69]
PREVIOUS='fresh_compact10_three_chain_gate_20261018'

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
                **B.B.score(s,nn,tt,ids+[newid],stage,walls))
            br['cap']=rec['index'];caps.append(rec)
    return dict(frame=frame,duties=duties,frame_profiles=[s.ps(i,types) for i in frame],
        duty_bonds=bonds,symbolic_high_families=[list(r) for r in symbolic if r[2]],
        complete_low_recipe_count=len(recipes),recipes=rows,caps=caps,
        branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def run():
    prior=json.loads((ROOT/f'results/{PREVIOUS}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s)
    originals=B.B.F.intervals(s,list(combinations([8,10,33,44,52,63,64,65],2)))
    debit=T.H.deletion(s.primary,s.cur,DELETED)
    assert debit['lost_finite_count']==19 and debit['compact_faces_by_dimension']=={2:12,3:1}
    assert sum(8 in v for v in debit['lost_finite_supports'])==3
    score=lambda nd,tt,ii,stage:B.B.score(s,nd,tt,ii,stage,[i for i in WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'control63')
    assert control['score'] is not None and control['open_ends']==10
    base={i:v for i,v in s.old.items() if i not in DELETED};ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    deletion=score(base,{},ids,'delete30_44_65')
    first=[];caps=[]
    for m in LOW:
        frame=[8,10,63];rhs=[0.,-T.B.theta(m),0.]
        roots,info=T.N.norm_roots([base[i] for i in frame],rhs,s.profiles[0]['square'])
        row=dict(index=len(first),profile=0,labels=[2,m,2],norm=info,branches=[]);first.append(row)
        assert not info['status'].startswith('unresolved')
        for sign,v in roots:
            vals,bad,cp=s.native(v,0,base,{})
            err=max(abs(float(v@T.M@v)-2),max(abs(float(v@T.M@base[i])-q) for i,q in zip(frame,rhs)))
            assert err<2e-5
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],
                released_products={str(i):s.dot(v,0,s.old[i],s.ps(i,{})) for i in DELETED})
            row['branches'].append(br)
            if bad or cp>=-1e-8:continue
            nd={**base,67:v};tt={67:0}
            assert s.corner([8,10,67],nd,tt) is not None
            br['cap']=len(caps)
            caps.append(dict(index=len(caps),recipe=row['index'],label10=m,profile=0,branch=sign,normal=v.tolist(),restores_old44=bool(np.max(abs(v-s.old[44]))<2e-6),**score(nd,tt,ids+[67],'new44')))
    assert len(caps)==4 and sum(x['restores_old44'] for x in caps)==1
    print('NEW44',len(caps),'changed',sum(not x['restores_old44'] for x in caps),flush=True)
    lb=bank(s,[8,52,64],[(0,1),(1,2)],base,{},ids,68,'new65')
    for a in lb['caps']:a['restores_old65']=bool(a['profile']==s.ps(65,{}) and np.max(abs(np.array(a['normal'])-s.old[65]))<2e-6)
    assert any(a['restores_old65'] for a in lb['caps'])
    print('NEW65_BANK',len(lb['recipes']),sum(len(r['branches']) for r in lb['recipes']),len(lb['caps']),flush=True)
    duties=[[8,67,68],[63,67,68],[63,64,68],[52,64,68],[8,52,68]]
    joins=[];pairs=[];r_banks=[];final=[]
    for j in caps:
        for l in lb['caps']:
            nd={**base,67:np.array(j['normal']),68:np.array(l['normal'])};tt={67:0,68:l['profile']}
            vals,bad,cp=s.native(nd[68],tt[68],{i:v for i,v in nd.items() if i!=68},tt)
            corners=[s.corner(f,nd,tt) for f in duties]
            row=dict(index=len(joins),J=j['index'],L=l['index'],native_product=vals[67],relation=T.G.relation(vals[67]),blockers=[[i,vals[i]] for i in bad],center_product=cp,corner_types=corners)
            joins.append(row)
            if bad or cp>=-1e-8:continue
            for i in (67,68):
                _,fail,center=s.native(nd[i],tt[i],{a:v for a,v in nd.items() if a!=i},tt)
                assert not fail and center<-1e-8
            row['pair']=len(pairs)
            pair=dict(index=len(pairs),J=j['index'],L=l['index'],all_five_corners=all(c is not None for c in corners),both_controls=j['restores_old44'] and l['restores_old65'],**score(nd,tt,ids+[67,68],'new44_65'))
            pairs.append(pair)
            if not pair['all_five_corners']:continue
            rb=bank(s,[10,33,67],[(0,1),(0,2)],nd,tt,ids+[67,68],69,'new30')
            r_banks.append(dict(pair=pair['index'],bank=rb))
            for r in rb['caps']:
                nn={**nd,69:np.array(r['normal'])};types={**tt,69:r['profile']}
                faces={f['wall']:f for f in r['whole_faces']}
                extra=[[7,55,69],[55,59,69],[52,64,69],[63,64,69],[63,67,69]]
                ec=[s.corner(f,nn,types) for f in extra]
                compact={str(i):faces[i]['all_rays_finite'] for i in (8,10,55,63,64,67,68)}
                final.append(dict(index=len(final),pair=pair['index'],R=r['index'],restores_old30=bool(r['profile']==s.ps(30,{}) and np.max(abs(nn[69]-s.old[30]))<2e-6),both_controls=pair['both_controls'],
                    extra_corner_duties=extra,extra_corner_types=ec,all_extra_corners=all(c is not None for c in ec),whole_compact_faces=compact,
                    clean=r['score'] is not None,score=r['score'],open_ends=r['open_ends'],finite_rays=r['finite_rays'],whole_faces=r['whole_faces']))
            print('R_BANK','pair',pair['index'],'recipes',len(rb['recipes']),'native',len(rb['caps']),flush=True)
    necessary=bank(s,[63,64,8],[(0,1)],base,{},ids,68,'necessary63_64_first')
    for a in necessary['caps']:
        nd={**base,68:np.array(a['normal'])};tt={68:a['profile']}
        a['remaining8_52_corner']=s.corner([8,52,68],nd,tt)
        a['restores_old65']=bool(a['profile']==s.ps(65,{}) and np.max(abs(nd[68]-s.old[65]))<2e-6)
    assert [a['restores_old65'] for a in necessary['caps'] if a['remaining8_52_corner'] is not None]==[True]
    debit63=T.H.deletion(s.primary,s.cur,(30,44,63,65));assert {k:v for k,v in debit63.items() if k!="deleted_walls"}=={k:v for k,v in debit.items() if k!="deleted_walls"}
    summary=dict(J_recipes=len(first),J_branches=sum(len(r['branches']) for r in first),J_native=len(caps),J_changed=3,
        L_recipes=len(lb['recipes']),L_branches=sum(len(r['branches']) for r in lb['recipes']),L_native=len(lb['caps']),
        joins=len(joins),native_pairs=len(pairs),five_corner_pairs=sum(p['all_five_corners'] for p in pairs),
        changed_five_corner_pairs=sum(p['all_five_corners'] and not p['both_controls'] for p in pairs),
        R_banks=len(r_banks),R_recipes=sum(len(r['bank']['recipes']) for r in r_banks),R_branches=sum(sum(len(q['branches']) for q in r['bank']['recipes']) for r in r_banks),
        native_triples=len(final),changed_triples=sum(not(p['both_controls'] and p['restores_old30']) for p in final),clean_triples=sum(p['clean'] for p in final),
        restored_compact_faces=sum(p['clean'] and all(p['whole_compact_faces'].values()) for p in final),fewer_than_ten_ends=sum(p['clean'] and p['open_ends']<10 for p in final),necessary_recipes=len(necessary['recipes']),necessary_branches=sum(len(r['branches']) for r in necessary['recipes']),necessary_native=len(necessary['caps']),necessary_completed=sum(a['remaining8_52_corner'] is not None for a in necessary['caps']),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(B.__file__),Path(H.__file__),Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY}
    paths.update(ROOT/f'{d}/{tag}.{e}' for tag in (PREVIOUS,B.TAG) for d,e in (('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/p for p in prior['sources'])
    return dict(status='BOUNDED_JOINT44_65_COMPACT_CHAIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,exact_deletion_debit=debit,old_frame_interface=interface,original_native_intervals=originals,
        control=control,deletion=deletion,J_frame=[8,10,63],J_duty=[8,10,67],J_recipes=first,J_caps=caps,L_bank=lb,L_five_corner_duties=duties,joins=joins,pairs=pairs,R_banks=r_banks,final=final,summary=summary,
        necessary63_64_bank=necessary,additional63_deletion_debit=debit63,boundary_ledgers=s.ledgers,unresolved=s.unresolved,numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='Four explicit native-first profile0 replacements44; complete declared seven-profile low2..5 replacement65 cube, all pairs, and complete low replacement30 cubes only for native five-corner pairs. All retained old and mutual native products and strict center precede cone scoring; numerical cone/PD agreement and affected whole-face ledgers distinguish clean from provisional results. Original44/65/30 controls included. Exact19/12/1 deleted certificate debit, not exact new-normal or ambient H9 volume certification. High-label/other-profile44 families and continuous parameter completeness not claimed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
