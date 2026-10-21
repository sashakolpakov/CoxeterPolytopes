#!/usr/bin/env python3
"""Bounded joint10/30 replacement seeded by two sole10 native failures."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
import fresh_compact10_edge_subdivision_20261018 as V

K,H,T,ROOT=V.K,V.H,V.T,V.ROOT
TAG='fresh_joint10_30_native_replacement_20261018'
SEED=ROOT/f'results/{V.TAG}.json'


def changed_frame_checks(s):
    p=s.primary;dec=lambda q:tuple(map(F,q))
    def mul(a,b):
        q=[x*y for x in a for y in b];return min(q),max(q)
    h=list(map(dec,p['common_center']));g=[[dec(q) for q in r] for r in p['computational_gram']]
    ns={}
    for r in p['facet_witnesses']:
        hp=dec(p['center_products'][str(r['wall'])]);assert hp[1]<0
        ns[r['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,r['point']))]
    out=[]
    for i,j in combinations([8,44,63],2):
        q=(F(),F())
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);q=(q[0]+z[0],q[1]+z[1])
        val=s.dot(s.old[i],s.ps(i,{}),s.old[j],s.ps(j,{}))
        assert float(q[0])-2e-5<=val<=float(q[1])+2e-5
        out.append(dict(walls=[i,j],profiles=[s.ps(i,{}),s.ps(j,{})],native_product=val,
            relation=T.G.relation(val),exact_native_enclosure=list(map(str,q))))
    return out


def lowbank(s,nd,types,ids,iscontrol):
    frame=[44,33,67] if iscontrol else [8,44,63]
    duties=[(0,2),(1,2)] if iscontrol else [(0,1),(1,2)]
    bonds=[2]*9;checks=[]
    for a,b in combinations(range(3),2):
        val=s.dot(nd[frame[a]],s.ps(frame[a],types),nd[frame[b]],s.ps(frame[b],types))
        rel=T.G.relation(val);checks.append(dict(walls=[frame[a],frame[b]],native_product=val,relation=rel))
        if (a,b) in duties:
            assert isinstance(rel,int),(frame,a,b,rel)
            bonds[3*a+b]=rel
    dom=T.E.recipes(tuple(s.ps(i,types) for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms)) for pi,ms,high in dom if not high and all(2<=m<=5 for m in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
        if all(T.P.corner3(s.ps(frame[a],types),s.ps(frame[b],types),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    tests=[];caps=[];seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]] if i in s.index else s.off(pi,types[i])) for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i] for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(tests),profile=pi,labels=ms,norm=info,branches=[]);tests.append(row)
        s.counts['new10_norm_equations']+=1
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage='new10_norm',recipe=row['index'],info=info));continue
        for sign,v in roots:
            s.counts['new10_norm_branches']+=1
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q) for i,q in zip(frame,rhs)))
            vals,bad,cp=s.native(v,pi,nd,types)
            nn={**nd,68:v};tt={**types,68:pi}
            fail=[tri for tri in [[7,33,68],[33,67,68],[7,11,68],[8,11,68],[8,44,68]] if s.corner(tri,nn,tt) is None]
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,
                native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad],extra_corner_failures=fail)
            row['branches'].append(br)
            if err>2e-5:
                s.unresolved.append(dict(stage='new10_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            assert all(s.corner([frame[a],frame[b],68],nn,tt) is not None for a,b in duties)
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key)
            for new in (67,68):
                vv,bb,cc=s.native(nn[new],tt[new],{i:x for i,x in nn.items() if i!=new},tt)
                assert not bb and cc<-1e-8
            before=len(s.unresolved);sid=s.score(nn,tt,ids+[68],'joint10_30')
            br['cap']=len(caps)
            caps.append(dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,branch=sign,
                normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,
                extra_corner_failures=fail,score=sid,cone_validated=sid is not None,
                cone_diagnostic=None if sid is not None else s.unresolved[before:],
                restores_old10_numerically=pi==s.ps(10,{}) and float(np.max(abs(v-s.old[10])))<2e-6))
    return dict(frame=frame,duties=duties,frame_native_checks=checks,complete_low_recipe_count=len(recipes),
        symbolic_high_families_not_tested=sum(bool(h) for _,_,h in dom),recipes=tests,caps=caps)


def run():
    old=json.loads(SEED.read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);framechecks=V.frame_checks(s);newchecks=changed_frame_checks(s)
    debit=T.H.deletion(s.primary,s.cur,[10,30])
    assert debit['lost_finite_count']==75
    assert debit['compact_faces_by_dimension']=={2:313,3:325,4:244,5:129,6:46,7:10,8:1}
    assert debit['compact_S6_sides_touched']==[8,10,11,44,55,63,64]
    old8lost=[v for v in debit['lost_finite_supports'] if 8 in v]
    base={i:v for i,v in s.old.items() if i not in (10,30)}
    ids=[i for i in s.cur['actual_sides'] if i not in (10,30)]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');assert control is not None
    deleted=s.score(base,{},ids,'delete10_30')
    replacement={10:68,30:67}
    restore=[sorted(replacement.get(i,i) for i in r['active']) for r in s.ledgers[control]['rays']
        if r['kind']=='finite' and set(r['active'])&{10,30}]
    # The fixed two changed seeds are identified by the complete frozen bank's
    # sole-wall failure, not by a rescan with a larger label cutoff.
    selected=[]
    for row in old['banks'][0]['recipes']:
        for br in row['branches']:
            if len(br['blockers'])==1 and br['blockers'][0][0]==10 and br['center_product']<-1e-8:
                selected.append((row,br,False))
    controlcap=next(c for c in old['banks'][0]['caps'] if c['restores_old_numerically'])
    rr=old['banks'][0]['recipes'][controlcap['recipe']]
    bb=next(b for b in rr['branches'] if b['branch']==controlcap['branch'])
    selected.append((rr,bb,True))
    assert [(r['profile'],r['labels'],b['branch']) for r,b,c in selected]==[(0,[5,2,4],1),(0,[5,2,5],1),(0,[5,2,2],1)]
    prefixes=[];banks=[];coherent=[];digest=sha256()
    for row,br,iscontrol in selected:
        pi=row['profile'];ms=row['labels'];frame=[7,55,33]
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]] for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([s.old[i] for i in frame],rhs,s.profiles[pi]['square'])
        assert not info['status'].startswith('unresolved')
        v=next(v for sign,v in roots if sign==br['branch'])
        beforebase={i:x for i,x in s.old.items() if i!=30}
        before,bad,cp=s.native(v,pi,beforebase,{})
        assert bad==([] if iscontrol else [10])
        assert abs(cp-br['center_product'])<2e-8
        assert all(abs(before[i]-x)<2e-8 for i,x in br['blockers'])
        vals,bad,cp=s.native(v,pi,base,{})
        assert not bad and cp<-1e-8
        error=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@s.old[i])-q) for i,q in zip(frame,rhs)))
        assert error<2e-5
        nn={**base,67:v};tt={67:pi};actual=ids+[67]
        assert s.corner([7,55,67],nn,tt) is not None
        start=len(s.unresolved);sid=s.score(nn,tt,actual,'new30_prefix')
        pref=dict(index=len(prefixes),source_recipe=row['index'],profile=pi,labels=ms,branch=br['branch'],
            old30_control=iscontrol,normal=v.tolist(),norm_residual=error,center_product=cp,
            before_dropping10_native_products=sorted(before.items()),before_dropping10_blockers=br['blockers'],
            retained_native_products=sorted(vals.items()),score=sid,cone_diagnostic=None if sid is not None else s.unresolved[start:])
        prefixes.append(pref)
        r44=s.dot(v,pi,s.old[44],s.ps(44,{}))
        pref['native44_relation']=T.G.relation(r44)
        pref['native44_product']=r44
        assert pref['native44_relation']==(2 if iscontrol else 'dotted')
        bank=lowbank(s,nn,tt,actual,iscontrol);banks.append(bank)
        print('JOINT10_30',pref['index'],ms,'recipes',len(bank['recipes']),'native',len(bank['caps']),flush=True)
        for cap in bank['caps']:
            ci=cap['index'];digest.update(json.dumps([pref['index'],ci,cap['profile'],cap['labels'],cap['branch'],cap['native_products']],separators=(',',':')).encode())
            if cap['score'] is None:continue
            final={**nn,68:np.array(cap['normal'])};types={**tt,68:cap['profile']}
            ledger=s.ledgers[cap['score']]
            failed=[f for f in restore if s.corner(f,final,types) is None]
            faces=[H.face(ledger,i) for i in (8,11,33,44,55,63,64,67,68)]
            f44=next(f for f in faces if f['wall']==44)
            new44={tuple(sorted(t)) for t in [[8,44,68],[44,63,68],[8,44,65],[44,63,65]]}
            row=dict(prefix=pref['index'],cap=ci,score=cap['score'],old30_restored=iscontrol,
                old10_restored=cap['restores_old10_numerically'],finite_rays=ledger['counts'].get('finite',0),
                open_ends=ledger['open_ends'],failed_old_S6_corner_restorations=failed,whole_face_records=faces,
                all_old_compact_S6_faces_rebuilt=all(f['all_rays_finite'] for f in faces if f['wall'] in (8,11,44,55,63,64,68)),
                mandatory_7_33_new10=s.corner([7,33,68],final,types) is not None,
                changed44_quadrilateral_complete=f44['all_rays_finite'] and {tuple(sorted(r['active'])) for r in f44['rays']}==new44,
                exterior_supports=[r['active'] for r in ledger['rays'] if r['kind']=='exterior'],
                open_edges=[e for e in ledger['edges'] if e['kind'] in ('one_open','two_open')])
            coherent.append(row)
    assert any(r['old30_restored'] and r['old10_restored'] and r['open_ends']==10 and not r['failed_old_S6_corner_restorations'] for r in coherent)
    paths=[Path(__file__),SEED,ROOT/f'results/{V.TAG}.md',Path(V.__file__),Path(K.__file__),Path(H.__file__),
        Path(T.__file__),Path(T.H.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,T.H.PARENT,T.H.BOUNDARY,Path(T.S.__file__),Path(T.S.C.__file__)]
    return dict(status='BOUNDED_JOINT10_30_NATIVE_REPLACEMENT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        deleted_walls=[10,30],retained_old_wall_count=61,exact_deletion_debit=debit,lost_old8_finite_supports=old8lost,
        old_frame_interface=interface,old_frame_native_checks=framechecks,changed_frame_native_checks=newchecks,
        numerical_replay_snapshot=dict(old_normals=[[i,v.tolist()] for i,v in sorted(s.old.items())],center=s.center.tolist(),
            profiles=[dict(contacts=p['contacts'],square=p['square'],center_offset=p['center_offset'],
                native_old_offsets=[[i,p['offsets'][s.index[i]]] for i in sorted(s.old)]) for p in s.profiles]),
        control_score=control,deletion_score=deleted,
        selection='Exactly the two sole10 native failures in the frozen finite edge-subdivision R bank, plus its unique old30 restoration. Reconstruct before deleting10, then apply all61 retained old-wall gates. The changed R/44 products are dotted, so the proposed finiteR/44/new10 duty is rejected before search. Instead changed prefixes use frame8/44/63 to replace the old44 link segment8--10--30--63 by8--new10--63. The old30 control alone uses frame44/33/R.',
        low_labels=[2,3,4,5],prefixes=prefixes,new10_banks=banks,required_old_S6_corner_restorations=restore,
        complete_native_pair_count=sum(len(b['caps']) for b in banks),complete_native_pair_digest=digest.hexdigest(),
        coherent_pairs=coherent,boundary_ledgers=s.ledgers,counts=dict(s.counts),unresolved=s.unresolved,
        scope='Two selected changed R prefixes use all seven profiles and complete2..5 new10 labels in frame8/44/63 with finite8/44/new10 and44/63/new10; the old30 control uses finite44/R/new10 and33/R/new10 in frame44/33/R. No outer-anchor corner imposed. All retained native/center gates precede all admitted pair cone/PD checks; seven actual compact-face checks, changed whole44 quadrilateral and all old S6 finite-corner losses are recorded. Exact75 finite-support certificate debit with losses on retained8, not an ambient reconstruction claim. Numerical bounded scout; no exact newnormal, high/continuous exclusion, finite-volume or newdimension claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));target=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(target.read_text())
    else:assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'pairs',out['complete_native_pair_count'],'coherent',len(out['coherent_pairs']),'unresolved',len(out['unresolved']),flush=True)
