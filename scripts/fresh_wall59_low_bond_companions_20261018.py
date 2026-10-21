#!/usr/bin/env python3
"""Bounded changed pentagonal cavity after replacing the old59/60/61 prism."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_wall59_bond_replacement_20261018 as R
T=R.T;ROOT=R.ROOT;TAG='fresh_wall59_low_bond_companions_20261018'
CYCLE=[30,67,58,66,33]
SEEDS=[[2,2,3],[2,4,2]]

def topology(rotation,ncaps):
    bc=CYCLE[rotation:]+CYCLE[:rotation];caps=[68,69,70][:ncaps]
    arcs=[2,3]if ncaps==2 else[2,2,1];fs=[];pos=0;starts=[]
    for cap,n in zip(caps,arcs):
        starts.append(pos)
        for k in range(n):fs.append([bc[(pos+k)%5],bc[(pos+k+1)%5],cap])
        pos+=n
    for j in range(ncaps):fs.append([bc[starts[j]],caps[j-1],caps[j]])
    if ncaps==3:fs.append(caps)
    ee=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]))
    assert(len(vs),len(ee),len(fs))==((7,13,7)if ncaps==2 else(8,16,9))
    assert{e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert set(ee.values())=={1,2};links=[]
    for v in sorted(vs):
        es=[tuple(x for x in f if x!=v)for f in fs if v in f];ds=Counter(x for e in es for x in e);seen={next(iter(ds))}
        while True:
            ns=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
            if ns==seen:break
            seen=ns
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
        links.append(dict(vertex=v,edges=es))
    return dict(rotation=rotation,caps=caps,arc_lengths=arcs,boundary=bc,triangles=fs,vertices=len(vs),edges=len(ee),faces=len(fs),links=links)

def full_caps(s,score,top):
    bd=s.ledgers[score];fs={tuple(sorted(f))for f in top['triangles']};rays=bd['rays']
    return all({tuple(sorted(r['active']))for r in rays if c in r['active']}=={f for f in fs if c in f}
        and all(r['kind']=='finite'for r in rays if c in r['active'])for c in top['caps'])

def run():
    inp=ROOT/'results/fresh_wall59_bond_replacement_20261018.json';data=json.loads(inp.read_text())
    assert data['counts']['accepted_59,60,61']==89
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter();s.near=[];s.improvements=[]
    debit=T.H.deletion(s.primary,s.cur,[59,60,61]);assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    controls=s.score(s.old,{},s.cur['actual_sides'],'control63')
    base={i:v for i,v in s.old.items()if i not in[59,60,61]};ids=[i for i in s.cur['actual_sides']if i not in[59,60,61]]
    allseeds=[];tops=[topology(r,n)for r in range(5)for n in[2,3]];records=[];pairhist=Counter();pairs=[];completions=[];dynamic=[];unresolved=[]
    for seedlabels in SEEDS:
        matches=[r for r in data['accepted']if r['release']==[59,60,61]and r['labels']==seedlabels];assert len(matches)==1
        seed=matches[0];v=np.array(seed['normal']);vals,bad,cp=s.native(v,0,base,{})
        assert not bad and cp<-1e-8
        nd={**base,67:v};types={67:0};actual=ids+[67]
        sid=s.score(nd,types,actual,'seed');assert sid is not None
        bd=s.ledgers[sid];edges={tuple(sorted(e['sides'])):e['kind']for e in bd['edges']}
        ce={tuple(sorted(e))for e in zip(CYCLE,CYCLE[1:]+CYCLE[:1])}
        assert all(edges[e]=='one_open'for e in ce)
        assert{e for e,k in edges.items()if k in['one_open','two_open']}-ce=={(7,47),(7,30),(30,52),(47,52)}
        banks=[];si=len(allseeds)
        for r in range(5):
            frame=[CYCLE[(r+k)%5]for k in range(3)]
            bank=s.bank(frame,[(0,1),(1,2)],nd,types,actual,68,'s%d_wedge%d'%(si,r));banks.append(bank)
            print('WEDGE',si,r,len(bank['symbolic_recipes']),len(bank['caps']),'unresolved',len(s.unresolved),flush=True)
        allseeds.append(dict(labels=seedlabels,source_index=seed['index'],normal=seed['normal'],score=sid,banks=banks))
        for rot in range(5):
            bc=CYCLE[rot:]+CYCLE[:rot];bx=banks[rot]['caps'];by=banks[(rot+2)%5]['caps']
            for ix,x in enumerate(bx):
                for iy,y in enumerate(by):
                    value=s.dot(np.array(x['normal']),x['profile'],np.array(y['normal']),y['profile']);rel=T.G.relation(value)
                    pairhist[si,rot,str(rel)]+=1
                    if rel in[None,'parallel_or_unresolved']:continue
                    nn,tt=s.add(nd,types,[x,y])
                    sid=s.score(nn,tt,actual+[68,69],'compatible_xy');
                    if sid is None:continue
                    pairs.append(dict(seed=si,rotation=rot,x=ix,y=iy,native_product=value,relation=rel,score=sid))
                    if not isinstance(rel,int):continue
                    junction=s.corner([bc[2],68,69],nn,tt)
                    if junction is None:continue
                    rec=dict(seed=si,rotation=rot,x=ix,y=iy,native_product=value,label=rel,score=sid,two_complete=False)
                    records.append(rec)
                    top=topology(rot,2)
                    if all(s.corner(f,nn,tt)is not None for f in top['triangles']):
                        rec['two_complete']=full_caps(s,sid,top)
                        if rec['two_complete']:completions.append(dict(seed=si,rotation=rot,caps=2,pair_record=len(records)-1,score=sid))
                    # The third cap owns only b4-b0. Its second prescribed
                    # corner is b0-X-Z, so no unnecessary third oldwall duty.
                    db=s.bank([bc[4],bc[0],68],[(0,1),(1,2)],nn,tt,actual+[68,69],70,'s%d_r%d_z'%(si,rot))
                    dr=dict(seed=si,rotation=rot,x=ix,y=iy,bank=db);dynamic.append(dr)
                    top3=topology(rot,3)
                    for iz,z in enumerate(db['caps']):
                        n3={**nn,70:np.array(z['normal'])};t3={**tt,70:z['profile']}
                        if not all(s.corner(f,n3,t3)is not None for f in top3['triangles']):continue
                        complete=full_caps(s,z['score'],top3)
                        if complete:completions.append(dict(seed=si,rotation=rot,caps=3,dynamic=len(dynamic)-1,z=iz,score=z['score']))
            print('JOIN',si,rot,'finite_junctions',len(records),'dynamic',len(dynamic),'complete',len(completions),flush=True)
    paths=[Path(__file__),Path(R.__file__),inp,Path(T.__file__),Path(T.H.__file__),Path(T.P.__file__),Path(T.E.__file__),Path(T.N.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_LOW_BOND_PENTAGON_COMPANIONS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        exact_deletion_debit=debit,control_score=controls,seed_labels=SEEDS,cycle=CYCLE,topologies=tops,seeds=allseeds,
        pair_histograms=[dict(seed=a,rotation=b,relation=c,count=n)for(a,b,c),n in sorted(pairhist.items())],compatible_pairs=pairs,finite_junctions=records,dynamic_z_banks=dynamic,completions=completions,
        boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        counts=dict(s.counts),blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in s.blockers.most_common()],near=s.near,unresolved=unresolved+s.unresolved,
        scope='Only two explicit low-angle replacement seeds; finite-corner recipe low labels2..5 and inherited high samples6/8/12/24. All accepted caps and finite junction assemblies pass retained native/common-center gates and complete numericalS6 cone-versus-independentPD support tests. Exact disk/deletion bookkeeping is separate from numerical normals. Removed old prism certificates are not restored merely by a lower open-end count. No exact ambient chamber or finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'seeds',len(out['seeds']),'finitejoins',len(out['finite_junctions']),'complete',len(out['completions']),'unresolved',len(out['unresolved']),flush=True)
