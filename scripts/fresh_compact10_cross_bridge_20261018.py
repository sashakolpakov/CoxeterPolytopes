#!/usr/bin/env python3
"""Replace30/44 by a bridge8--10--33 and a separately native new30.

The changed wall10 boundary omits its old44--30 chain. Complete low
bridge recipes are tested before joining frozen numerical new30 walls.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_exterior_polar_native10_frame_20261018 as F

T,H,ROOT=F.T,F.H,F.ROOT
TAG='fresh_compact10_cross_bridge_20261018'
SEED=ROOT/f'results/{F.TAG}.json'
LOW=(2,3,4,5)
DELETED=(30,44)


def whole(s,sid,diagnostic,walls):
    bd=s.ledgers[sid] if sid is not None else diagnostic[0]['boundary']
    return dict(finite_rays=bd['counts'].get('finite',0),open_ends=bd['open_ends'],
        exterior_supports=[r['active'] for r in bd['rays'] if r['kind']=='exterior'],
        open_edges=[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')],
        whole_faces=[H.face(bd,w) for w in walls])


def score(s,nd,types,ids,stage,walls):
    before=len(s.unresolved);sid=s.score(nd,types,ids,stage)
    ds=None if sid is not None else s.unresolved[before:]
    assert sid is not None or len(ds)==1
    return dict(score=sid,cone_diagnostic=ds,**whole(s,sid,ds,walls))


def bridge_bank(s,base,ids):
    frame=[8,10,33];duties=[(0,1),(1,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(base[frame[a]],s.ps(frame[a],{}),base[frame[b]],s.ps(frame[b],{})))
        assert m==2;bonds[3*a+b]=m
    symbolic=T.E.recipes(tuple(s.ps(i,{}) for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms)) for pi,ms,high in symbolic if not high and all(z in LOW for z in ms)]
    direct=[(pi,list(ms)) for pi in range(7) for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],{}),s.ps(frame[b],{}),pi,bonds[3*a+b],ms[a],ms[b]) is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    records=[];caps=[];seen=set();digest=sha256();hist=Counter()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]] for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([base[i] for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(records),profile=pi,labels=ms,norm=info,branches=[]);records.append(row)
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage='bridge_norm',recipe=row['index'],norm=info));continue
        for sign,v in roots:
            vals,bad,cp=s.native(v,pi,base,{})
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@base[i])-q) for i,q in zip(frame,rhs)))
            br=dict(branch=sign,normal=v.tolist(),norm_residual=err,center_product=cp,
                native_products=sorted(vals.items()),blockers=[[i,vals[i]] for i in bad])
            row['branches'].append(br);hist[str(tuple(bad))]+=1
            digest.update(json.dumps([row['index'],sign,pi,ms,br['native_products'],cp],separators=(',',':')).encode())
            if err>2e-5:
                s.unresolved.append(dict(stage='bridge_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nd={**base,67:v};tt={67:pi}
            assert all(s.corner([frame[a],frame[b],67],nd,tt) is not None for a,b in duties)
            rec=dict(index=len(caps),recipe=row['index'],branch=sign,profile=pi,labels=ms,normal=v.tolist(),
                **score(s,nd,tt,ids+[67],'bridge',[8,10,11,33,55,63,64,65,66,67]))
            br['cap']=rec['index'];caps.append(rec)
    return dict(frame=frame,duties=duties,nonedge=[8,33],recipes=records,caps=caps,
        symbolic_high_families_not_tested=sum(bool(h) for _,_,h in symbolic),
        branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def run():
    seed=json.loads(SEED.read_text())
    for p,h in seed['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    # The exact original frame-interface calculation is frozen. Reconstruct
    # its profile choices from the actual stored A2 offsets before using it.
    interface=[]
    for r in seed['old_frame_interface']['actual_side_profiles']:
        i=r['wall'];k=s.index[i];a=s.profiles[1]['offsets'][k];b=s.profiles[4]['offsets'][k]
        contacts=[-(2*a-b),-(2*b-a)];matches=[j for j,p in enumerate(s.profiles) if max(abs(x-y) for x,y in zip(contacts,p['contacts']))<2e-5]
        assert matches==[r['profile']];s.types[i]=matches[0]
        assert max(abs(x-y) for x,y in zip(contacts,r['stored_A2_contacts']))<2e-5
        assert np.max(abs(s.old[i]-np.array(seed['numerical_snapshot']['old_normals'][str(i)])))<2e-6
        interface.append(dict(wall=i,profile=matches[0],contacts=contacts))
    originals=F.intervals(s,list(combinations([8,10,33],2)))
    assert T.G.relation(s.dot(s.old[8],s.ps(8,{}),s.old[33],s.ps(33,{})))=='dotted'
    debit=T.H.deletion(s.primary,s.cur,DELETED)
    assert debit['lost_finite_count']==16 and debit['compact_faces_by_dimension']=={2:12,3:1}
    assert sum(8 in v for v in debit['lost_finite_supports'])==2
    control=score(s,s.old,{},s.cur['actual_sides'],'control63',[8,10,11,30,33,44,55,63,64,65,66])
    assert control['score'] is not None and control['open_ends']==10
    base={i:v for i,v in s.old.items() if i not in DELETED};ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    deletion=score(s,base,{},ids,'delete30_44',[8,10,11,33,55,63,64,65,66])
    bank=bridge_bank(s,base,ids)
    branches=[b for r in bank['recipes'] for b in r['branches']]
    assert len(bank['recipes'])==42 and len(branches)==84 and not bank['caps']
    shared=sorted(set.intersection(*(set(i for i,_ in b['blockers']) for b in branches)))
    assert shared==[7,55]
    obstruction=dict(shared_native_blockers=shared,
        native_product_ranges={str(i):[min(dict(b['native_products'])[i] for b in branches),
            max(dict(b['native_products'])[i] for b in branches)] for i in shared},
        positive_products={str(i):sum(dict(b['native_products'])[i]>2e-7 for b in branches) for i in shared},
        scope='Observed complete bounded low-cube obstruction only; high-label and continuous exclusions are not established.')
    print('BRIDGE_BANK','recipes',len(bank['recipes']),'branches',sum(len(r['branches']) for r in bank['recipes']),'native',len(bank['caps']),flush=True)
    joins=[];accepted=[];target10=sorted(map(sorted,([7,10,11],[8,10,11],[8,10,67],[10,33,67],[7,10,33])))
    for j in bank['caps']:
        nd={**base,67:np.array(j['normal'])};tt={67:j['profile']}
        for r in seed['bank']['caps']:
            v=np.array(r['normal']);vals,bad,cp=s.native(v,r['profile'],nd,tt)
            row=dict(bridge=j['index'],frozen_R=r['index'],native_product_JR=vals[67],relation=T.G.relation(vals[67]),
                blockers=[[i,vals[i]] for i in bad],center_product=cp,native_products=sorted(vals.items()))
            joins.append(row)
            if bad or cp>=-1e-8:continue
            nn={**nd,68:v};types={**tt,68:r['profile']}
            # Restore all native gates in both directions before cone scoring.
            for a in (67,68):
                aa,bb,cc=s.native(nn[a],types[a],{i:w for i,w in nn.items() if i!=a},types)
                assert not bb and cc<-1e-8
            report=score(s,nn,types,ids+[67,68],'bridge_and_R',[8,10,11,33,55,63,64,65,66,67,68])
            f10=next(f for f in report['whole_faces'] if f['wall']==10)
            f8=next(f for f in report['whole_faces'] if f['wall']==8)
            exact10=f10['all_rays_finite'] and sorted(sorted(q['active']) for q in f10['rays'])==target10
            row['accepted']=len(accepted)
            accepted.append(dict(bridge=j['index'],frozen_R=r['index'],restores_old30=r['restores_old30'],
                expected_compact10_rays=target10,complete_compact10=exact10,compact8=f8['all_rays_finite'],
                R7_55_corner=s.corner([7,55,68],nn,types),**report))
    paths={Path(__file__),SEED,ROOT/f'results/{F.TAG}.md',Path(F.__file__),Path(T.__file__),Path(H.__file__),
        Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.S.__file__),
        Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY}
    paths.update(ROOT/p for p in seed['sources'])
    return dict(status='BOUNDED_COMPACT10_CHANGED_CROSS_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        changed_incidence='Replace compact10 boundary chain8--44--30--33 by8--J--33. The new30/33 product may remain dotted; finite10/new30/33 is not required. Only two J corner duties are prescribed; no whole-J disk or ambient closure is assumed.',
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,retained_projected_old_walls=len(base),
        fixed_white_walls=[0,1,2,3,5,6],exact_deletion_debit=debit,actual_side_profiles=interface,
        original_frame_native_intervals=originals,control=control,deletion=deletion,bank=bank,joins=joins,accepted=accepted,
        frozen_R_domain_count=len(seed['bank']['caps']),bounded_obstruction=obstruction,
        boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        snapshot=seed['numerical_snapshot'],
        scope='Complete low2..5 bridge cube over seven inherited white profiles, both norm branches and every retained old native/center gate. If reached, every join with all35 frozen numerical R walls restores all mutual/old gates and full S6 cone versus PD support checks. Whole8/10 and all affected side ledgers are explicit. Exact16/12/1 old-certificate debit including2 wall8 supports; no interval newnormal, high-label/continuous completeness or H9 finite-volume certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()))
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'bridges',len(out['bank']['caps']),'joins',len(out['joins']),'accepted',len(out['accepted']),'unresolved',len(out['unresolved']),flush=True)
