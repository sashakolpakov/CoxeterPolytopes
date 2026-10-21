#!/usr/bin/env python3
"""Move33 jointly with the compact10 bridge; complete bounded two-stage banks.

Removing30/33/44 changes the previously exhausted fixed33 geometry. Every
native intermediate is retained, including provisional cones, until the
second cap is added. Whole-face and full cone/PD ledgers remain mandatory.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_compact10_cross_bridge_20261018 as B
import fresh_joint30_33_hinge_native_pilot_20261018 as H

T,ROOT=B.T,B.ROOT
TAG='fresh_moving33_compact10_bridge_20261018'
LOW=(2,3,4,5)
DELETED=(30,33,44)
WALLS=[7,8,10,11,55,63,64,65,66,67,68]
FROZEN='fresh_compact10_cross_bridge_high_interval_20261018'


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
                restores_old33=bool(pi==s.ps(33,{}) and np.max(abs(v-s.old[33]))<2e-6),
                **B.score(s,nn,tt,ids+[newid],stage,walls))
            br['cap']=rec['index'];caps.append(rec)
    return dict(frame=frame,duties=duties,frame_profiles=[s.ps(i,types) for i in frame],
        duty_bonds=bonds,symbolic_high_families=[list(r) for r in symbolic if r[2]],
        complete_low_recipe_count=len(recipes),recipes=rows,caps=caps,
        branch_histogram=dict(hist),complete_branch_digest=digest.hexdigest())


def run():
    prior=json.loads((ROOT/f'results/{FROZEN}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s)
    originals=B.F.intervals(s,list(combinations([7,8,10,33,66],2)))
    debit=T.H.deletion(s.primary,s.cur,DELETED)
    assert debit['lost_finite_count']==18 and debit['compact_faces_by_dimension']=={2:12,3:1}
    assert sum(8 in v for v in debit['lost_finite_supports'])==2
    control=B.score(s,s.old,{},s.cur['actual_sides'],'control63',[7,8,10,11,30,33,44,55,63,64,65,66])
    assert control['score'] is not None and control['open_ends']==10
    base={i:v for i,v in s.old.items() if i not in DELETED};ids=[i for i in s.cur['actual_sides'] if i not in DELETED]
    deletion=B.score(s,base,{},ids,'delete30_33_44',[i for i in WALLS if i in base])
    first=bank(s,[7,10,66],[(0,1),(0,2)],base,{},ids,67,'moving33')
    print('MOVING33_BANK',len(first['recipes']),sum(len(r['branches']) for r in first['recipes']),len(first['caps']),flush=True)
    assert any(c['restores_old33'] for c in first['caps'])
    records=[];accepted=[];target10=sorted(map(sorted,([7,10,11],[8,10,11],[8,10,68],[10,67,68],[7,10,67])))
    for a in first['caps']:
        nn={**base,67:np.array(a['normal'])};tt={67:a['profile']}
        second=bank(s,[8,10,67],[(0,1),(1,2)],nn,tt,ids+[67],68,'moving33_bridge')
        rec=dict(S=a['index'],S_restores_old33=a['restores_old33'],bank=second);records.append(rec)
        for b in second['caps']:
            nd={**nn,68:np.array(b['normal'])};types={**tt,68:b['profile']}
            faces={f['wall']:f for f in b['whole_faces']}
            finite10=faces[10]['all_rays_finite'] and sorted(sorted(r['active']) for r in faces[10]['rays'])==target10
            whole={str(i):faces[i]['all_rays_finite'] for i in (8,10,55,63,64,65,67,68)}
            accepted.append(dict(index=len(accepted),S=a['index'],J=b['index'],score=b['score'],
                cone_validated=b['score'] is not None,open_ends=b['open_ends'],finite_rays=b['finite_rays'],
                complete_compact10=finite10,expected_compact10_rays=target10,whole_faces_compact=whole,
                extra_corner_duties=dict(wall8_J65=s.corner([8,65,68],nd,types),J63_65=s.corner([63,65,68],nd,types)),
                new_wall_count=2,old_wall_count=len(base)+6))
        if a['restores_old33']:assert not second['caps'] # frozen exact obstruction, included as control
        print('BRIDGE',a['index'],'S',a['profile'],a['labels'],'recipes',len(second['recipes']),'native',len(second['caps']),flush=True)
    branches=lambda b:sum(len(r['branches']) for r in b['recipes'])
    summary=dict(first_recipes=len(first['recipes']),first_branches=branches(first),first_native=len(first['caps']),
        moving_S=sum(not a['restores_old33'] for a in first['caps']),second_banks=len(records),
        second_recipes=sum(len(r['bank']['recipes']) for r in records),second_branches=sum(branches(r['bank']) for r in records),
        native_pairs=len(accepted),clean_pairs=sum(a['cone_validated'] for a in accepted),
        whole_compact10=sum(a['cone_validated'] and a['complete_compact10'] for a in accepted),
        whole_new_faces=sum(a['cone_validated'] and a['whole_faces_compact']['67'] and a['whole_faces_compact']['68'] for a in accepted),
        fewer_than_ten_ends=sum(a['cone_validated'] and a['open_ends']<10 for a in accepted),
        unresolved=len(s.unresolved),complete_low_domains=True)
    paths={Path(__file__),Path(B.__file__),Path(H.__file__),Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),
        Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY}
    paths.update(ROOT/f'{d}/{n}.{e}' for n in(FROZEN,B.TAG,H.TAG) for d,e in(('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/p for p in prior['sources'])
    return dict(status='BOUNDED_MOVING33_COMPACT10_BRIDGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        changed_incidence='Replace30/33/44 by S,J: finite7/10/S and7/66/S first, then8/10/J and10/S/J. Prospective wall10 boundary7--11--8--J--S. This changes the old33 anchor in the exhausted bridge frame; no8/S or10/66 adjacency is imposed.',
        deleted_walls=list(DELETED),retained_old_walls=len(base)+6,exact_deletion_debit=debit,
        old_frame_interface=interface,original_native_intervals=originals,control=control,deletion=deletion,
        first_bank=first,records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        numerical_snapshot=dict(metric=T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist() for i,v in s.old.items()},
            labels=s.labels,profiles=s.profiles,actual_sides=s.cur['actual_sides']),
        scope='Complete declared seven-profile native2..5 two-stage banks, with every retained old/new native product, strict common-center and full numerical cone/PD agreement before clean acceptance. Provisional first caps are retained for second-stage closure. Exact18/12/1 old-certificate debit including2 wall8 supports. Whole target-side faces are explicit; high-label families, continuous parameter completeness, exact new normals and ambient H9 finite volume remain unproved.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
