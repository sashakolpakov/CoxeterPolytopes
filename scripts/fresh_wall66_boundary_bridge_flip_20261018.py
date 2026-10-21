#!/usr/bin/env python3
"""Rebuild J against67, replacing X-Z by67-J on frozen X/Y prefixes."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
import fresh_wall66_seam_subdivision_20261018 as C

T=C.T;ROOT=C.ROOT;TAG='fresh_wall66_boundary_bridge_flip_20261018'


def topology():
    old=C.topology();removed=[{67,68,71},{68,70,71}]
    fs=[f for f in old['triangles']if set(f)not in removed]+[[67,68,70],[67,70,71]]
    edges=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]));bc=old['boundary']
    assert(len(vs),len(edges),len(fs))==(9,19,11)
    assert(68,71)not in edges and edges[67,70]==2
    assert{e for e,n in edges.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert set(edges.values())=={1,2};links=[]
    for v in sorted(vs):
        es=[tuple(q for q in f if q!=v)for f in fs if v in f];ds=Counter(q for e in es for q in e);seen={next(iter(ds))}
        while True:
            new=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(ds)and sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
        links.append(dict(vertex=v,edges=es))
    return dict(boundary=bc,caps=old['caps'],triangles=fs,vertices=9,edges=19,faces=11,links=links,
        removed_edge=[68,71],new_edge=[67,70])


def run():
    p=ROOT/f'results/{C.TAG}.json';data=json.loads(p.read_text())
    priorpath=ROOT/'results/fresh_wall59_low_bond_companions_20261018.json';prior=json.loads(priorpath.read_text())
    for name,h in data['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    debit=T.H.deletion(s.primary,s.cur,[59,60,61]);assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    base={i:v for i,v in s.old.items()if i not in[59,60,61]};ids=[i for i in s.cur['actual_sides']if i not in[59,60,61]]
    top=topology();records=[];joins=[];zrows=[];completion=[];tests=[];oldcontacts=[]
    for r in data['records']:
        for j in r['J_bank']['caps']:
            vals=dict(j['native_products']);oldcontacts.append([vals[30],vals[67]])
    assert len(oldcontacts)==95 and max(v for row in oldcontacts for v in row)<-2
    for ri,rec in enumerate(data['records']):
        pair=rec['selection']['pair'];seed=prior['seeds'][pair['seed']]
        x=seed['banks'][1]['caps'][pair['x']];y=seed['banks'][3]['caps'][pair['y']]
        nd={**base,67:np.array(seed['normal'])};tt={67:0};actual=ids+[67]
        for i,cap in[(68,x),(69,y)]:
            vals,bad,cp=s.native(np.array(cap['normal']),cap['profile'],nd,tt)
            assert not bad and cp<-1e-8
            nd[i]=np.array(cap['normal']);tt[i]=cap['profile'];actual.append(i)
        assert T.G.relation(s.dot(nd[68],tt[68],nd[69],tt[69]))=='dotted'
        score=s.score(nd,tt,actual,'prefix');assert score is not None
        # No fictitious66/67 adjacency is passed to the recipe generator.
        bank=C.lowbank(s,[66,67,68],[(0,2),(1,2)],nd,tt,actual,70,'new_J')
        records.append(dict(prefix=ri,selection=rec['selection'],prefix_score=score,J_bank=bank))
        for ji,j in enumerate(bank['caps']):
            nn={**nd,70:np.array(j['normal'])};types={**tt,70:j['profile']}
            value=s.dot(nn[69],types[69],nn[70],types[70]);rel=T.G.relation(value)
            corner=s.corner([66,69,70],nn,types)
            joins.append(dict(prefix=ri,J=ji,native_YJ=value,relation=rel,corner66_Y_J=corner))
            if corner is None:continue
            zbank=C.lowbank(s,[30,67,70],[(0,1),(1,2)],nn,types,actual+[70],71,'new_Z')
            zr=dict(prefix=ri,J=ji,bank=zbank);zrows.append(zr)
            for zi,z in enumerate(zbank['caps']):
                full={**nn,71:np.array(z['normal'])};types4={**types,71:z['profile']}
                bad=[f for f in top['triangles']if s.corner(f,full,types4)is None]
                whole=C.P.full_caps(s,z['score'],top)
                tests.append(dict(Z_bank=len(zrows)-1,Z=zi,nonfinite_corners=bad,whole_cap_ray_sets=whole,score=z['score']))
                if not bad and whole:completion.append(dict(Z_bank=len(zrows)-1,Z=zi,score=z['score']))
        print('BOUNDARY_BRIDGE',ri,'J',len(bank['caps']),'finiteYJ',sum(r['corner66_Y_J']is not None for r in joins),'Zbanks',len(zrows),'complete',len(completion),flush=True)
    paths=[Path(__file__),Path(C.__file__),Path(T.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),Path(T.H.__file__),p,priorpath,T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_CHANGED_BOUNDARY_BRIDGE_FLIP',sources={str(q.relative_to(ROOT)):sha256(q.read_bytes()).hexdigest()for q in paths},
        exact_old_deletion_debit=debit,topology=top,control_score=control,old_J_boundary_contacts=oldcontacts,
        prefix_count=16,complete_low_labels=C.LOW,records=records,YJ_tests=joins,Z_banks=zrows,completion_tests=tests,completions=completion,
        counts=dict(s.counts),boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in s.blockers.most_common()],unresolved=s.unresolved,
        scope='Actual changed diagonal XZ→67J; J is reconstructed, never reused from the95 dotted-boundary caps. Complete low2..5 recipes over seven inherited profiles at only16 frozen X/Y representatives. Every admitted stage has full retained old/new native and common-center gates, full numerical cone and independent separated PD incidence comparison. Full disk and whole-cap ray sets required for completion. Same exact7/6/1 debit; no high-label/continuous completeness, exact normal or ambient finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'J',len(out['YJ_tests']),'Zbanks',len(out['Z_banks']),'complete',len(out['completions']),'unresolved',len(out['unresolved']),flush=True)
