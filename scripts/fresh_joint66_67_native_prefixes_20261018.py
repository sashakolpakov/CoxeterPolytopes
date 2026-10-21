#!/usr/bin/env python3
"""Move both old66 and the wall59 replacement; bounded native prefix search."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import fresh_wall66_seam_subdivision_20261018 as C

T=C.T; ROOT=C.ROOT
TAG='fresh_joint66_67_native_prefixes_20261018'
RELEASE=[59,60,61,66]


def run():
    s=T.Scout()
    s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    s.counts=Counter();s.blockers=Counter();s.near=[];s.improvements=[]
    debit=T.H.deletion(s.primary,s.cur,RELEASE)
    assert debit['lost_finite_count']==9 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    base={i:v for i,v in s.old.items() if i not in RELEASE}
    ids=[i for i in s.cur['actual_sides'] if i not in RELEASE]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    deletion=s.score(base,{},ids,'delete59_60_61_66')
    assert control is not None and deletion is not None
    specs=[dict(new_wall=67,frame=[55,30,58],duties=[(0,1),(0,2)]),
           dict(new_wall=68,frame=[7,33,58],duties=[(0,1),(0,2)])]
    banks=[]
    for spec in specs:
        b=C.lowbank(s,spec['frame'],spec['duties'],base,{},ids,spec['new_wall'],
                    'new%d'%spec['new_wall'])
        b['new_wall']=spec['new_wall']
        for cap in b['caps']:
            cap['restores_old66_numerically']=(cap['profile']==0 and
                float(np.max(abs(np.array(cap['normal'])-s.old[66])))<2e-6)
        banks.append(b)
        print('JOINT_BANK',spec['new_wall'],'recipes',b['complete_low_recipe_count'],
              'caps',len(b['caps']),'unresolved',len(s.unresolved),flush=True)
    pairs=[];hist=Counter();digest=sha256();accepted=[]
    for ix,x in enumerate(banks[0]['caps']):
        for iy,y in enumerate(banks[1]['caps']):
            value=s.dot(np.array(x['normal']),x['profile'],np.array(y['normal']),y['profile'])
            rel=T.G.relation(value);hist[str(rel)]+=1
            digest.update(json.dumps([ix,iy,value,rel],separators=(',',':')).encode())
            row=dict(index=len(pairs),x=ix,y=iy,native_product=value,relation=rel,score=None)
            pairs.append(row)
            if rel in (None,'parallel_or_unresolved'):continue
            nn={**base,67:np.array(x['normal']),68:np.array(y['normal'])}
            tt={67:x['profile'],68:y['profile']}
            # Both banks passed all59 old walls. Recheck the complete joint
            # support and center, including the actual mutual native product.
            for new in (67,68):
                rest={i:v for i,v in nn.items() if i!=new}
                vals,bad,cp=s.native(nn[new],tt[new],rest,tt)
                assert not bad and cp<-1e-8
            sid=s.score(nn,tt,ids+[67,68],'joint66_67')
            if sid is None:continue
            row['score']=sid;bd=s.ledgers[sid]
            faces=[]
            for wall in (55,67,68):
                rr=[r for r in bd['rays'] if wall in r['active']]
                faces.append(dict(wall=wall,rays=rr,all_rays_finite=all(r['kind']=='finite' for r in rr)))
            accepted.append(dict(pair=row['index'],score=sid,old66_restored=y['restores_old66_numerically'],
                open_ends=bd['open_ends'],finite_rays=bd['counts'].get('finite',0),
                exterior_supports=[r['active'] for r in bd['rays'] if r['kind']=='exterior'],
                open_edges=[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')],whole_face_records=faces))
        print('JOINT_ROW',ix,'pairs',len(pairs),'accepted',len(accepted),'unresolved',len(s.unresolved),flush=True)
    paths=[Path(__file__),Path(C.__file__),Path(C.P.__file__),Path(T.__file__),Path(T.H.__file__),
           Path(T.E.__file__),Path(T.N.__file__),Path(T.P.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_JOINT66_67_NATIVE_PREFIXES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        deleted_walls=RELEASE,retained_old_walls=59,exact_deletion_debit=debit,
        control_score=control,deletion_score=deletion,low_labels=[2,3,4,5],bank_specs=specs,banks=banks,
        complete_pair_count=len(pairs),pair_histogram=dict(hist),complete_pair_digest=digest.hexdigest(),
        pairs=pairs,accepted=accepted,boundary_ledgers=s.ledgers,counts=dict(s.counts),
        histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n) for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=list(b),count=n) for(a,b),n in s.blockers.most_common()],unresolved=s.unresolved,
        scope='Two independently rebuilt complete low2..5 seven-profile banks after joint59/60/61/66 deletion; both norm sheets and every mutual pair. All59 retained old native walls plus both new walls and commoncenter; every admitted cap/pair has complete numericalS6 cone and independentPD support agreement. Exact old9/6/1 support debit only. Old66 restorations flagged as controls. No companion search, high-label/continuous completeness, exact new normals or ambient finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:
        assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'banks',[len(b['caps']) for b in out['banks']],'pairs',out['complete_pair_count'],
          'accepted',len(out['accepted']),'unresolved',len(out['unresolved']),flush=True)
