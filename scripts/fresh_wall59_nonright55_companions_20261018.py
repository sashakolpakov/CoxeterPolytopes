#!/usr/bin/env python3
"""Test three changed nonright55 prefixes with complete bounded low cap banks.

All five wedge rotations; native/center gates and complete target-cone checks.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import fresh_wall59_low_bond_companions_20261018 as P
import fresh_wall59_nonright55_prefix_20261018 as N
import fresh_wall66_seam_subdivision_20261018 as C

T = P.T
ROOT = P.ROOT
TAG = 'fresh_wall59_nonright55_companions_20261018'
CYCLE = P.CYCLE


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def finite_edges(s, nodes, nd, types):
    return [list(e) for e in sorted({tuple(sorted(e)) for f in nodes for e in combinations(f,2)})
            if not isinstance(T.G.relation(s.dot(nd[e[0]],s.ps(e[0],types),
                                                nd[e[1]],s.ps(e[1],types))),int)]


def run():
    replacement = ROOT/f'results/{N.TAG}.json'
    old=json.loads(replacement.read_text())
    for name,h in old['sources'].items():assert digest(ROOT/name)==h,(name,h)
    wanted=[[3,3,4],[3,4,3],[3,5,5]]
    seeds=[next(r for r in old['accepted'] if r['labels']==labels) for labels in wanted]
    assert len(old['accepted'])==8 and len(seeds)==3
    s=T.Scout()
    s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    s.counts=Counter();s.blockers=Counter();s.near=[];s.improvements=[]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    debit=T.H.deletion(s.primary,s.cur,[59,60,61])
    assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    base={i:v for i,v in s.old.items() if i not in [59,60,61]}
    ids=[i for i in s.cur['actual_sides'] if i not in [59,60,61]]
    topologies=[P.topology(r,n) for r in range(5) for n in (2,3)]
    records=[];junctions=[];dynamic=[];completions=[];pairtotal=0
    for seed in seeds:
        nd={**base,67:np.array(seed['normal'])};types={67:0};actual=ids+[67]
        vals,bad,cp=s.native(nd[67],0,base,{})
        assert not bad and cp<-1e-8
        assert abs(float(nd[67]@T.M@nd[67])-s.profiles[0]['square'])<2e-5
        sid=s.score(nd,types,actual,'replacement_seed');assert sid is not None
        edges={tuple(sorted(e['sides'])):e['kind'] for e in s.ledgers[sid]['edges']}
        ce={tuple(sorted(e)) for e in zip(CYCLE,CYCLE[1:]+CYCLE[:1])}
        assert all(edges.get(e)=='one_open' for e in ce)
        rec=dict(source_index=seed['index'],labels=seed['labels'],normal=seed['normal'],
                 seed_score=sid,banks=[],pair_summaries=[])
        for rot in range(5):
            frame=[CYCLE[(rot+j)%5] for j in range(3)]
            # All seven profiles and complete declared low2..5 recipe cube.
            bank=C.lowbank(s,frame,[(0,1),(1,2)],nd,types,actual,68,
                           'seed%d_wedge%d'%(seed['index'],rot))
            for cap in bank['caps']:
                nn,tt=s.add(nd,types,[cap])
                assert s.corner([frame[0],frame[1],68],nn,tt) is not None
                assert s.corner([frame[1],frame[2],68],nn,tt) is not None
            rec['banks'].append(bank)
        for rot in range(5):
            bc=CYCLE[rot:]+CYCLE[:rot]
            bx=rec['banks'][rot]['caps'];by=rec['banks'][(rot+2)%5]['caps']
            hist=Counter();nearest=None;native_finite=0;pair_digest=sha256()
            for ix,x in enumerate(bx):
                for iy,y in enumerate(by):
                    value=s.dot(np.array(x['normal']),x['profile'],np.array(y['normal']),y['profile'])
                    rel=T.G.relation(value);hist[str(rel)]+=1;pairtotal+=1
                    pair_digest.update(json.dumps([ix,iy,value,rel],separators=(',',':')).encode())
                    if rel in [None,'parallel_or_unresolved']:continue
                    if nearest is None or value>nearest['native_product']:
                        nearest=dict(x=ix,y=iy,native_product=value,relation=rel)
                    if not isinstance(rel,int):continue
                    native_finite+=1
                    nn,tt=s.add(nd,types,[x,y])
                    corner=s.corner([bc[2],68,69],nn,tt)
                    if corner is None:continue
                    pairscore=s.score(nn,tt,actual+[68,69],'finite_junction')
                    if pairscore is None:continue
                    ji=len(junctions)
                    jr=dict(seed=len(records),rotation=rot,x=ix,y=iy,native_product=value,
                            label=rel,corner=corner,score=pairscore)
                    junctions.append(jr)
                    top=P.topology(rot,2)
                    missing=finite_edges(s,top['triangles'],nn,tt)
                    failed=[f for f in top['triangles'] if s.corner(f,nn,tt) is None]
                    whole=P.full_caps(s,pairscore,top)
                    jr['two_cap_test']=dict(missing_edges=missing,failed_corners=failed,whole_caps=whole)
                    if not missing and not failed and whole:
                        completions.append(dict(junction=ji,caps=2,score=pairscore))
                    zb=C.lowbank(s,[bc[4],bc[0],68],[(0,1),(1,2)],nn,tt,actual+[68,69],70,
                                 'junction%d_Z'%ji)
                    dr=dict(junction=ji,bank=zb,tests=[]);dynamic.append(dr)
                    top=P.topology(rot,3)
                    for zi,z in enumerate(zb['caps']):
                        n3={**nn,70:np.array(z['normal'])};t3={**tt,70:z['profile']}
                        missing=finite_edges(s,top['triangles'],n3,t3)
                        failed=[f for f in top['triangles'] if s.corner(f,n3,t3) is None]
                        test=dict(Z=zi,missing_edges=missing,failed_corners=failed,score=None,whole_caps=False)
                        if not missing and not failed:
                            zs=s.score(n3,t3,actual+[68,69,70],'all_finite_disk_corners')
                            if zs is not None:
                                test.update(score=zs,whole_caps=P.full_caps(s,zs,top))
                                if test['whole_caps']:
                                    completions.append(dict(junction=ji,caps=3,dynamic=len(dynamic)-1,Z=zi,score=zs))
                        dr['tests'].append(test)
            if nearest is not None:
                nn,tt=s.add(nd,types,[bx[nearest['x']],by[nearest['y']]])
                nearest['score']=s.score(nn,tt,actual+[68,69],'closest_native_pair')
            rec['pair_summaries'].append(dict(rotation=rot,domain_size=len(bx)*len(by),
                relation_histogram=dict(hist),native_finite=native_finite,closest_native_pair=nearest,
                complete_pair_value_digest=pair_digest.hexdigest()))
        records.append(rec)
        print('SEED',seed['index'],seed['labels'],'caps',[len(b['caps']) for b in rec['banks']],
              'finitejunctions',len(junctions),'completed',len(completions),'unresolved',len(s.unresolved),flush=True)
    paths=[Path(__file__),Path(P.__file__),Path(P.R.__file__),Path(N.__file__),Path(C.__file__),Path(T.__file__),
           Path(T.H.__file__),Path(T.E.__file__),Path(T.N.__file__),Path(T.P.__file__),Path(T.G.__file__),
           replacement,T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_NONRIGHT55_PREFIX_COMPANIONS',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        exact_deletion_debit=debit,control_score=control,
        selection=dict(all_new_prefixes=8,selected_labels=wanted,new_seed_count=3,
                       rule='Three distinct new55label3 prefixes: near each low asymmetric angle corner, then both remaining labels5.'),
        low_labels=[2,3,4,5],high_label_samples=[],topologies=topologies,
        records=records,complete_pair_tests=pairtotal,finite_junctions=junctions,
        dynamic_Z_banks=dynamic,completed_disks=completions,
        boundary_ledgers=s.ledgers,counts=dict(s.counts),
        histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)
                    for (a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=list(b),count=n) for (a,b),n in s.blockers.most_common()],
        unresolved=s.unresolved,
        scope='Three explicitly selected numerical nonright55 prefixes, all five wedge rotations, all seven profiles and complete low2..5 finite-corner recipes; no high slots sampled. Every accepted cap, finite junction, closest native pair and locally complete disk receives complete numericalS6 cone/PD comparison. All retained old/new native products and commoncenter required. Same7/6/1 deletion debit; no exact normal, continuous/high-label exclusion, old-prism restoration or ambient finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'seeds',len(out['records']),'pairs',out['complete_pair_tests'],
          'finitejoins',len(out['finite_junctions']),'completed',len(out['completed_disks']),flush=True)
