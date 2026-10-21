#!/usr/bin/env python3
"""Bounded native three-cap replacement of the open wall66 hexagon.

Both two-arcs-per-cap rotations; fixed retained old walls, all seven
A2 profiles, low labels and four explicit high representatives.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_cu_native_kr_split7_20261017 as G

ROOT,H,E,N,S,P,B,M=G.ROOT,G.H,G.E,G.N,G.S,G.P,G.B,G.M
TAG='fresh_release66_three_cap_disk_20261017'
BOUNDARY=[7,33,30,59,60,58]
HIGH=[6,8,12,24]


def topology(rotation):
    bc=BOUNDARY[rotation:]+BOUNDARY[:rotation];frames=[bc[2*j:2*j+3]if j<2 else[bc[4],bc[5],bc[0]]for j in range(3)]
    caps=[67,68,69];fs=[]
    for cap,frame in zip(caps,frames):fs.extend([[frame[0],frame[1],cap],[frame[1],frame[2],cap]])
    fs += [[bc[2],67,68],[bc[4],68,69],[bc[0],69,67],caps]
    es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2));vs=set(sum(fs,[]))
    assert(len(vs),len(es),len(fs))==(9,18,10)and len(vs)-len(es)+len(fs)==1
    assert set(es.values())=={1,2}and{e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];dd=Counter(i for e in ee for i in e);seen={next(iter(dd))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(dd)and sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(rotation=rotation,boundary=bc,frames=frames,triangles=fs,vertices=9,edges=18,faces=10,links=links)


def run():
    cur,labels,old,h,profiles,primary=H.setup()
    deletion=H.deletion(primary,cur,[66]);assert deletion['lost_finite_count']==2 and not deletion['compact_faces_touched']
    # Rest frame of the inherited common center reduces cancellation near60.
    assert abs(float(h@M@h)+1)<1e-10 and h[0]>0
    boost=np.eye(4);boost[0,0]=h[0];boost[0,1:]=h[1:];boost[1:,0]=h[1:]
    boost[1:,1:] += np.outer(h[1:],h[1:])/(h[0]+1)
    assert np.max(abs(boost.T@M@boost-M))<1e-10
    inverse=M@boost.T@M;ns={i:inverse@v for i,v in old.items()};center=inverse@h
    assert np.max(abs(center-np.array([1.,0,0,0])))<1e-10
    oldindex={i:j for j,i in enumerate(labels)}
    off=lambda a,b:float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    ps=lambda i:4 if i==7 else 6 if i==60 else 0
    for i in BOUNDARY:assert abs(float(ns[i]@M@ns[i])-profiles[ps(i)]['square'])<1e-5
    bonds={tuple(sorted((i,j))):G.relation(float(ns[i]@M@ns[j])+off(ps(i),ps(j)))for i,j in zip(BOUNDARY,BOUNDARY[1:]+BOUNDARY[:1])}
    assert all(isinstance(m,int)for m in bonds.values())
    actual=[i for i in cur['actual_sides']if i!=66];base={i:v for i,v in ns.items()if i!=66}
    ledgers=[];ledgerkeys={}
    def score(nd,types,ids):
        bd=S.C.detailed_cone(nd,ids)
        sq=[profiles[types[i]]['square']if i in types else float(ns[i]@M@ns[i])for i in ids]
        ic=H.incidence_score([nd[i]for i in ids],ids,sq)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and{tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
        if not good:return None,dict(boundary=bd,PD_incidence=ic)
        result=G.F.D.compact_score(dict(boundary=bd,PD_incidence=ic));key=json.dumps(result,sort_keys=True)
        if key not in ledgerkeys:ledgerkeys[key]=len(ledgers);ledgers.append(result)
        return ledgerkeys[key],None
    control,problem=score(ns,{},cur['actual_sides']);assert problem is None
    assert ledgers[control]['counts']=={'finite':26,'exterior':6}and ledgers[control]['open_ends']==10
    deleted,problem=score(base,{},actual);assert problem is None
    def native(v,pi,nd,types):
        vals={i:float(v@M@w)+(profiles[pi]['offsets'][oldindex[i]]if i in oldindex else off(pi,types[i]))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if G.relation(x)in[None,'parallel_or_unresolved']]
        return vals,bad,float(v@M@center)+profiles[pi]['center_offset']
    counts=Counter();blockers=Counter();unresolved=[];rotations=[];improvements=[];hist=Counter();near=[]
    for rot in range(2):
        topo=topology(rot);pools=[];domains=[]
        for side,frame in enumerate(topo['frames']):
            bb=[2]*9;bb[1]=bonds[tuple(sorted(frame[:2]))];bb[5]=bonds[tuple(sorted(frame[1:]))]
            dom=E.recipes(tuple(ps(i)for i in frame),tuple(bb),((0,1),(1,2)))
            domains.append(dict(frame=frame,symbolic_recipes=[list(z)for z in dom]));pool=[];seen=set()
            for pi,low,hi in dom:
                for hh in product(HIGH,repeat=len(hi)):
                    ms=list(low)
                    for j,m in zip(hi,hh):ms[j]=m
                    rhs=[-B.theta(m)-profiles[pi]['offsets'][oldindex[i]]for i,m in zip(frame,ms)]
                    roots,info=N.norm_roots([base[i]for i in frame],rhs,profiles[pi]['square']);counts['norm_equations']+=1
                    if info['status'].startswith('unresolved'):
                        unresolved.append(dict(rotation=rot,side=side,profile=pi,labels=ms,stage='norm',info=info));continue
                    for sign,v in roots:
                        counts['norm_branches']+=1
                        err=max(abs(float(v@M@v)-profiles[pi]['square']),max(abs(float(v@M@base[i])-r)for i,r in zip(frame,rhs)))
                        if err>2e-5:unresolved.append(dict(rotation=rot,side=side,stage='residual',profile=pi,labels=ms,branch=sign,error=err));continue
                        vals,bad,cp=native(v,pi,base,{})
                        if bad or cp>=-1e-8:
                            blockers[rot,side,tuple(bad)]+=1
                            if len(bad)<=3:near.append(dict(rotation=rot,side=side,profile=pi,labels=ms,branch=sign,blockers=[[i,vals[i]]for i in bad],center_product=cp))
                            continue
                        key=(pi,)+tuple(np.round(v,8))
                        if key in seen:counts['duplicates']+=1;continue
                        seen.add(key);nn=dict(base);nn[67+side]=v;tt={67+side:pi}
                        sid,problem=score(nn,tt,actual+[67+side])
                        if problem:unresolved.append(dict(rotation=rot,side=side,stage='cone',profile=pi,labels=ms,**problem));continue
                        rec=dict(index=len(pool),profile=pi,labels=ms,branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,score=sid)
                        pool.append(rec);counts['native_single_caps']+=1
                        bd=ledgers[sid];hist['single',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                        if bd['open_ends']<10:improvements.append(dict(stage='single',rotation=rot,side=side,cap=rec))
            pools.append(pool)
            print('RELEASE66_POOL',rot,side,frame,'domain',len(dom),'pool',len(pool),'unresolved',len(unresolved),flush=True)
        pairmaps=[];pairrows=[];pairhist=Counter()
        for a,b,junction in[(0,1,topo['boundary'][2]),(1,2,topo['boundary'][4]),(0,2,topo['boundary'][0])]:
            finite={};rows=[]
            for ra,rb in product(pools[a],pools[b]):
                x=float(H.A.native(np.array(ra['normal']),ra['profile'],np.array(rb['normal']),rb['profile'],profiles));rel=G.relation(x)
                counts['mutual_pair_tests']+=1;pairhist[str(rel)]+=1
                if not isinstance(rel,int):continue
                ma=ra['labels'][topo['frames'][a].index(junction)];mb=rb['labels'][topo['frames'][b].index(junction)]
                corner=P.corner3(ps(junction),ra['profile'],rb['profile'],ma,mb,rel)
                if corner is None:counts['nonfinite_junction_pairs']+=1;continue
                finite[ra['index'],rb['index']]=rel;rows.append(dict(indices=[ra['index'],rb['index']],native_product=x,label=rel,corner=corner))
            pairmaps.append(finite);pairrows.append(dict(sides=[a,b],junction=junction,pairs=rows))
        triples=[]
        for(ia,ib),mab in pairmaps[0].items():
            for ic in range(len(pools[2])):
                if(ib,ic)not in pairmaps[1]or(ia,ic)not in pairmaps[2]:continue
                rs=[pools[0][ia],pools[1][ib],pools[2][ic]];mbc=pairmaps[1][ib,ic];mac=pairmaps[2][ia,ic]
                central=P.corner3(*(r['profile']for r in rs),mab,mac,mbc)
                counts['three_finite_junctions']+=1
                if central is None:counts['central_corner_nonfinite']+=1;continue
                nn=dict(base);tt={}
                for i,r in enumerate(rs):nn[67+i]=np.array(r['normal']);tt[67+i]=r['profile']
                for i in tt:
                    _,bad,cp=native(nn[i],tt[i],{j:v for j,v in nn.items()if j!=i},tt);assert not bad and cp<-1e-8
                sid,problem=score(nn,tt,actual+[67,68,69])
                if problem:unresolved.append(dict(rotation=rot,stage='triple_cone',indices=[ia,ib,ic],**problem));continue
                bd=ledgers[sid];rayset={tuple(sorted(z['active']))for z in bd['rays']if z['kind']=='finite'};whole={}
                for i in[67,68,69]:
                    expected={tuple(sorted(t))for t in topo['triangles']if i in t};rays=[z for z in bd['rays']if i in z['active']]
                    whole[i]=all(z['kind']=='finite'for z in rays)and{tuple(sorted(z['active']))for z in rays}==expected
                rec=dict(indices=[ia,ib,ic],mutual_labels=[mab,mbc,mac],central_corner=central,whole_faces=whole,completed_disk=all(tuple(sorted(t))in rayset for t in topo['triangles'])and all(whole.values()),score=sid)
                triples.append(rec);counts['coherent_three_caps']+=1
                hist['triple',bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']]+=1
                if bd['open_ends']<10:improvements.append(dict(stage='triple',rotation=rot,triple=rec,roots=rs))
        rotations.append(dict(topology=topo,domains=domains,pools=pools,pair_relations=dict(pairhist),finite_junction_pairs=pairrows,triples=triples))
        print('RELEASE66_ROTATION',rot,'pools',list(map(len,pools)),'pairs',dict(pairhist),'triples',len(triples),flush=True)
    paths=[Path(__file__),Path(G.__file__),Path(H.__file__),Path(E.__file__),Path(N.__file__),Path(S.__file__),Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_RELEASE66_THREE_CAP_DISK_SCOUT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=[66],retained_control_walls=62,exact_deletion_ledger=deletion,control_score=control,deletion_score=deleted,
        coordinate_change=dict(boost=boost.tolist(),inverse=inverse.tolist(),center=center.tolist()),old_boundary_bonds=[[list(k),v]for k,v in bonds.items()],
        high_label_samples=HIGH,counts=dict(counts),rotations=rotations,boundary_ledgers=ledgers,improvements=improvements,unresolved=unresolved,near_native_failures=near,
        histograms=[dict(stage=s,finite=f,exterior=e,open_ends=n,count=c)for(s,f,e,n),c in sorted(hist.items())],
        blockers=[dict(rotation=k[0],side=k[1],walls=list(k[2]),count=n)for k,n in sorted(blockers.items())],
        scope='Two distinct rotations of a three-cap two-arcs-per-cap triangulated hexagonal replacement for66. Seven A2profiles and complete low finite-corner recipe domain, but high labels sampled only6/8/12/24. All62 retained native gates and commoncenter; complete numericalS6cone/PD comparison for deletion, every coherent singlecap and every three-cap candidate with allfour inner finite corners. No exact cap certificate, all-label/parameter exclusion, full ambient census or finite-volume claim. Exact63wall199finite10end control remains frozen.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
