#!/usr/bin/env python3
"""Replace failed dotted X-Y seam by J, retaining the changed wall59 cavity."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_wall59_low_bond_companions_20261018 as P
T=P.T;ROOT=P.ROOT;TAG='fresh_wall66_seam_subdivision_20261018'
LOW=[2,3,4,5]

def topology():
    X,Y,J,Z=68,69,70,71;bc=[67,58,66,33,30]
    fs=[[67,58,X],[58,66,X],[66,33,Y],[33,30,Y],[30,67,Z],
        [67,X,Z],[30,Y,Z],[66,X,J],[66,Y,J],[X,J,Z],[Y,J,Z]]
    ee=Counter(tuple(sorted(e))for f in fs for e in combinations(f,2));vs=set(sum(fs,[]))
    assert(len(vs),len(ee),len(fs))==(9,19,11)and set(ee.values())=={1,2}
    assert{e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    assert(68,69)not in ee;links=[]
    for v in sorted(vs):
        es=[tuple(q for q in f if q!=v)for f in fs if v in f];dd=Counter(q for e in es for q in e);seen={next(iter(dd))}
        while True:
            ns=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
            if ns==seen:break
            seen=ns
        assert seen==set(dd)and sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=es))
    return dict(boundary=bc,caps=[X,Y,J,Z],triangles=fs,vertices=9,edges=19,faces=11,links=links)

def lowbank(s,frame,duties,nd,types,ids,newid,stage):
    bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],types),nd[frame[b]],s.ps(frame[b],types)))
        assert isinstance(m,int);bonds[3*a+b]=m
    symbolic=T.E.recipes(tuple(s.ps(i,types)for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms))for pi,ms,high in symbolic if not high and all(m in LOW for m in ms)]
    # Independent finite-corner enumeration over the complete low native cube.
    direct=[]
    for pi in range(len(s.profiles)):
        for ms in product(LOW,repeat=3):
            if all(T.P.corner3(s.ps(frame[a],types),s.ps(frame[b],types),pi,bonds[3*a+b],ms[a],ms[b])is not None for a,b in duties):direct.append((pi,list(ms)))
    assert sorted(recipes)==sorted(direct)
    caps=[];tests=[];seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,types[i]))for i,m in zip(frame,ms)]
        rr,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(tests),profile=pi,labels=ms,norm=info,branches=[]);tests.append(row)
        s.counts['norm_equations']+=1
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage=stage,recipe=row['index'],info=info));continue
        for sign,v in rr:
            s.counts['norm_branches']+=1
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-r)for i,r in zip(frame,rhs)))
            vals,bad,cp=s.native(v,pi,nd,types)
            br=dict(branch=sign,norm_residual=err,center_product=cp,blockers=[[i,vals[i]]for i in bad]);row['branches'].append(br)
            if err>2e-5:
                s.unresolved.append(dict(stage=stage,recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:
                s.blockers[stage,tuple(bad)]+=1;continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nn={**nd,newid:v};tt={**types,newid:pi}
            sid=s.score(nn,tt,ids+[newid],stage)
            if sid is None:continue
            assert all(s.corner([frame[a],frame[b],newid],nn,tt)is not None for a,b in duties)
            br['cap']=len(caps);caps.append(dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,score=sid))
    return dict(frame=frame,duties=duties,complete_low_recipe_count=len(recipes),symbolic_high_families_not_tested=sum(bool(h)for _,_,h in symbolic),recipes=tests,caps=caps)

def run():
    path=ROOT/'results/fresh_wall59_low_bond_companions_20261018.json';old=json.loads(path.read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    groups={};lowcount=0;prior66=[]
    for si,seed in enumerate(old['seeds']):
        for bank in(1,3):
            b=seed['banks'][bank];j=b['frame'].index(66)
            prior66.append(dict(seed=si,bank=bank,native66_label_counts=dict(Counter(r['labels'][j]for r in b['caps']))))
    for r in old['compatible_pairs']:
        seed=old['seeds'][r['seed']];x=seed['banks'][1]['caps'][r['x']];y=seed['banks'][3]['caps'][r['y']]
        assert r['rotation']==1 and r['relation']=='dotted'
        if max(x['labels']+y['labels'])>5:continue
        lowcount+=1;key=(r['seed'],x['labels'][2],y['labels'][0])
        if key not in groups or r['native_product']>groups[key]['native_product']:groups[key]=r
    chosen=[dict(group=list(k),pair=r)for k,r in sorted(groups.items())]
    assert lowcount==3808 and len(chosen)==16
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    debit=T.H.deletion(s.primary,s.cur,[59,60,61]);assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');base={i:v for i,v in s.old.items()if i not in[59,60,61]};ids=[i for i in s.cur['actual_sides']if i not in[59,60,61]]
    top=topology();records=[];completions=[];zs=[];completion_tests=[]
    for selection in chosen:
        r=selection['pair'];seed=old['seeds'][r['seed']];x=seed['banks'][1]['caps'][r['x']];y=seed['banks'][3]['caps'][r['y']]
        nd={**base,67:np.array(seed['normal'])};types={67:0};actual=ids+[67]
        for label,cap in[(68,x),(69,y)]:
            v=np.array(cap['normal']);vals,bad,cp=s.native(v,cap['profile'],nd,types);assert not bad and cp<-1e-8
            assert abs(float(v@T.M@v)-s.profiles[cap['profile']]['square'])<2e-5
            nd[label]=v;types[label]=cap['profile'];actual.append(label)
        nativeXY=s.dot(nd[68],types[68],nd[69],types[69]);assert abs(nativeXY-r['native_product'])<1e-8 and T.G.relation(nativeXY)=='dotted'
        sid=s.score(nd,types,actual,'prefix');assert sid is not None
        bank=lowbank(s,[66,68,69],[(0,1),(0,2)],nd,types,actual,70,'J')
        rec=dict(selection=selection,prefix_score=sid,J_bank=bank);records.append(rec)
        print('BRIDGE',len(records)-1,selection['group'],'recipes',bank['complete_low_recipe_count'],'caps',len(bank['caps']),'unresolved',len(s.unresolved),flush=True)
        for ji,j in enumerate(bank['caps']):
            nn={**nd,70:np.array(j['normal'])};tt={**types,70:j['profile']}
            zbank=lowbank(s,[30,67,68],[(0,1),(1,2)],nn,tt,actual+[70],71,'Z')
            zr=dict(prefix=len(records)-1,J=ji,bank=zbank);zs.append(zr)
            for zi,z in enumerate(zbank['caps']):
                n4={**nn,71:np.array(z['normal'])};t4={**tt,71:z['profile']}
                required=sorted({tuple(sorted(e))for f in top['triangles']for e in combinations(f,2)if 71 in e})
                missing=[list(e)for e in required if not isinstance(T.G.relation(s.dot(n4[e[0]],s.ps(e[0],t4),n4[e[1]],s.ps(e[1],t4))),int)]
                failures=[f for f in top['triangles']if s.corner(f,n4,t4)is None]
                whole=P.full_caps(s,z['score'],top)
                completion_tests.append(dict(Z_bank=len(zs)-1,Z=zi,missing_finite_edges=missing,nonfinite_required_corners=failures,whole_cap_ray_sets=whole,score=z['score']))
                if not failures and whole:
                    completions.append(dict(Z_bank=len(zs)-1,Z=zi,score=z['score']))
            print('CLOSURE',len(records)-1,ji,'Zcaps',len(zbank['caps']),'complete',len(completions),flush=True)
    paths=[Path(__file__),Path(P.__file__),path,Path(T.__file__),Path(T.E.__file__),Path(T.P.__file__),Path(T.N.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY]
    return dict(status='BOUNDED_CHANGED_WALL66_SEAM_SUBDIVISION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        exact_deletion_debit=debit,old66_contact_coverage=prior66,eligible_low_prefix_pairs=lowcount,selection_rule='Maximum native X/Y product in each(seed,label66X,label66Y) group among complete saved low2..5 prefix pairs; all16 groups retained.',
        complete_low_labels=LOW,selected_prefixes=chosen,topology=top,control_score=control,records=records,Z_banks=zs,completion_tests=completion_tests,completions=completions,
        counts=dict(s.counts),boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in s.blockers.most_common()],unresolved=s.unresolved,
        scope='Genuinely changed central disk: X-Y dotted seam is subdivided by J, removing direct X-Y finite requirement. Complete low2..5 J/Z native recipe cubes over seven inherited profiles are checked only for16 explicitly selected frozen numerical prefixes; no continuous/high-label completeness. Every admitted stage retains all old/new native walls and commoncenter and complete numericalS6 cone/PD support comparison. Same exact old7/6/1 debit; no deletion-count improvement, exact normal or ambient finite-volume claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'prefixes',len(out['records']),'J',sum(len(r['J_bank']['caps'])for r in out['records']),'Z',sum(len(r['bank']['caps'])for r in out['Z_banks']),'complete',len(out['completions']),'unresolved',len(out['unresolved']),flush=True)
