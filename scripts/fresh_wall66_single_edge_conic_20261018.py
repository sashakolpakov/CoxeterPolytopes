#!/usr/bin/env python3
"""Retain all63 walls; continue a cap on the single open59/60 edge.

Numerical conic inequalities and native-contact events. No exact new chamber.
"""
import argparse,json,math
from collections import Counter
from hashlib import sha256
from itertools import product,combinations
from pathlib import Path
import numpy as np
import fresh_release58_66_seven_boundary_20261017 as T
ROOT,M=T.ROOT,T.M
TAG='fresh_wall66_single_edge_conic_20261018'
ANCHORS=[59,60]
HIGH=[6,8,12,24]
EVENT_LABELS=[2,3,4,5,6,8,12,24]
COMPANION_SEEDS=[1,9,10]
SIX_CYCLE=[30,59,67,58,66,33]

def topology():
    bc=SIX_CYCLE;frames=[[bc[(2*j+k)%6]for k in range(3)]for j in range(3)]
    fs=[]
    for cap,frame in zip([68,69,70],frames):
        fs.extend([[frame[0],frame[1],cap],[frame[1],frame[2],cap]])
    fs += [[67,68,69],[66,69,70],[30,70,68],[68,69,70]]
    es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2));vs=set(sum(fs,[]))
    assert(len(vs),len(es),len(fs))==(9,18,10)
    assert set(es.values())=={1,2}and{e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];dd=Counter(i for e in ee for i in e);seen={next(iter(dd))}
        while True:
            z=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if z==seen:break
            seen=z
        assert seen==set(dd)
        assert sorted(dd.values())==([1,1]+[2]*(len(dd)-2)if v in bc else[2]*len(dd))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,frames=frames,triangles=fs,vertices=9,edges=18,faces=10,links=links)

def companions(s,caps,sheets):
    topo=topology();records=[];expected=[[[58,3],[59,2],[60,2]],[[58,2],[59,3],[60,2]],[[58,3],[59,3],[60,2]]]
    native9=s.dot(s.old[58],s.ps(58,{}),s.old[59],s.ps(59,{}))
    assert T.G.relation(native9)==9
    for seed,contacts in zip(COMPANION_SEEDS,expected):
        cap=caps[seed];assert cap['profile']==0 and cap['finite_old_contacts']==contacts
        row=sheets[cap['sheet_record']];t=cap['parameter']
        v=np.array(row['part'])+np.array(row['va'])*t+np.array(row['vb'])/t
        assert np.max(abs(v-np.array(cap['normal'])))<1e-10
        vals,bad,cp=s.native(v,0,s.old,{});assert not bad and cp<-1e-8
        assert abs(float(v@M@v)-2)<2e-6
        nd={**s.old,67:v};types={67:0};ids=s.cur['actual_sides']+[67]
        sid=s.score(nd,types,ids,'companion_seed');assert sid==cap['score']
        bd=s.ledgers[sid];edges={tuple(sorted(e['sides'])):e['kind']for e in bd['edges']}
        boundary={tuple(sorted(e))for e in zip(SIX_CYCLE,SIX_CYCLE[1:]+SIX_CYCLE[:1])}
        assert all(edges[e]=='one_open'for e in boundary)
        assert{e for e,k in edges.items()if k in['one_open','two_open']}-boundary=={(7,30),(7,47),(30,52),(47,52)}
        # All rays of old60 are finite and exactly the new triangle's vertices.
        rays60=[r for r in bd['rays']if 60 in r['active']]
        assert all(r['kind']=='finite'for r in rays60)
        assert{tuple(sorted(r['active']))for r in rays60}=={(58,59,60),(58,60,67),(59,60,67)}
        banks=[]
        for a,frame in enumerate(topo['frames']):
            bank=s.bank(frame,[(0,1),(1,2)],nd,types,ids,68+a,f'companion_{seed}_{a}')
            banks.append(bank)
            print('SIX_CYCLE_BANK',seed,a,frame,'recipes',len(bank['symbolic_recipes']),'caps',len(bank['caps']),flush=True)
        pairs=[];pairmaps=[];hist=Counter()
        for a,b,junction in[(0,1,67),(1,2,66),(0,2,30)]:
            finite={};prows=[]
            for ra,rb in product(banks[a]['caps'],banks[b]['caps']):
                value=s.dot(np.array(ra['normal']),ra['profile'],np.array(rb['normal']),rb['profile']);rel=T.G.relation(value)
                hist[str(rel)]+=1;s.counts['companion_pair_tests']+=1
                if rel in[None,'parallel_or_unresolved']:continue
                nn={**nd,68+a:np.array(ra['normal']),68+b:np.array(rb['normal'])};tt={**types,68+a:ra['profile'],68+b:rb['profile']}
                psid=s.score(nn,tt,ids+[68+a,68+b],'companion_pair')
                if psid is None:continue
                corner=s.corner([junction,68+a,68+b],nn,tt)
                rec=dict(indices=[ra['index'],rb['index']],native_product=value,relation=rel,junction_corner=corner,score=psid)
                prows.append(rec)
                if corner is not None:finite[ra['index'],rb['index']]=rel
                if s.ledgers[psid]['open_ends']<10:s.improvements.append(dict(stage='companion_pair',seed=seed,sides=[a,b],record=rec))
            pairs.append(dict(sides=[a,b],junction=junction,records=prows));pairmaps.append(finite)
        triples=[]
        for(ia,ib),mab in pairmaps[0].items():
            for ic in range(len(banks[2]['caps'])):
                if(ib,ic)not in pairmaps[1]or(ia,ic)not in pairmaps[2]:continue
                rs=[banks[0]['caps'][ia],banks[1]['caps'][ib],banks[2]['caps'][ic]]
                nn,tt=s.add(nd,types,rs)
                corners=[s.corner(t,nn,tt)for t in topo['triangles']]
                if any(c is None for c in corners):continue
                for i in[68,69,70]:
                    _,bad,cp=s.native(nn[i],tt[i],{j:w for j,w in nn.items()if j!=i},tt);assert not bad and cp<-1e-8
                tsid=s.score(nn,tt,ids+[68,69,70],'companion_triple')
                if tsid is None:continue
                rec=T.completed(topo,rs,tsid,s);rec['indices']=[ia,ib,ic];rec['corners']=corners;triples.append(rec)
                if s.ledgers[tsid]['open_ends']<10:s.improvements.append(dict(stage='companion_triple',seed=seed,record=rec))
        records.append(dict(seed=seed,contacts=contacts,seed_score=sid,old60_complete_triangle=rays60,banks=banks,pairs=pairs,pair_histogram=dict(hist),triples=triples))
        print('SIX_CYCLE_DONE',seed,'banks',[len(b['caps'])for b in banks],'pairtypes',dict(hist),'triples',len(triples),flush=True)
    return dict(selected_seed_indices=COMPANION_SEEDS,selection='Three distinct profile0 low-label contacts to58/59, each right to60; selection is bounded, not complete in all conic/native events.',topology=topo,native58_59_product=native9,native58_59_label=9,records=records,counts=dict(s.counts),
                histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
                blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],near_native_failures=s.near,improvements=s.improvements)

def roots(coeff):
    a,b,c=map(float,coeff);scale=max(1.,abs(a),abs(b),abs(c))
    if abs(a)<1e-13*scale:
        if abs(b)<1e-13*scale:return [],'constant_or_nearconstant'
        return [-c/b],'linear'
    delta=b*b-4*a*c
    if abs(delta)<1e-12*scale*scale:return [],'near_tangent'
    if delta<0:return [],'nonreal'
    q=-.5*(b+math.copysign(math.sqrt(delta),b));out=[q/a]
    if q:out.append(c/q)
    return sorted(set(out)),'quadratic'

def qset(coeff):
    rr,status=roots(coeff);xs=[0.]+sorted(x for x in rr if x>0)+[float('inf')];out=[]
    a,b,c=coeff
    for lo,hi in zip(xs,xs[1:]):
        t=(lo+hi)/2 if math.isfinite(hi)else max(1.,2*lo)
        if(a*t+b)*t+c<0:out.append([lo,hi])
    return out,status

def meet(aa,bb):return[[max(a,c),min(b,d)]for a,b in aa for c,d in bb if max(a,c)<min(b,d)-1e-10*max(1.,max(a,c))]

def encode(intervals):return[[a,None if not math.isfinite(b)else b]for a,b in intervals]

def run():
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    old=s.old;frame=np.array([old[i]for i in ANCHORS]);G=frame@M@frame.T
    assert np.linalg.eigvalsh(G)[0]>1e-8;J=np.linalg.inv(G)
    def project(v):return v-(v@M@frame.T)@J@frame
    time=project(s.center);time/=math.sqrt(-float(time@M@time))
    space=project(old[66]);space+=float(space@M@time)*time;space/=math.sqrt(float(space@M@space))
    assert max(abs(frame@M@time))<1e-8 and max(abs(frame@M@space))<1e-8
    assert abs(float(space@M@time))<1e-8 and abs(float(space@M@space)-1)<1e-8 and abs(float(time@M@time)+1)<1e-8
    bond=T.G.relation(s.dot(old[59],s.ps(59,{}),old[60],s.ps(60,{})));assert bond==2
    symbolic=[];recipes=[]
    for pi,a,b in product(range(7),range(2,7),range(2,7)):
        if T.P.corner3(s.ps(59,{}),s.ps(60,{}),pi,bond,a,b)is None:continue
        symbolic.append([pi,a,b])
        if 6 in(a,b):
            assert(pi,a,b)==(0,6,2)
            recipes.extend((pi,h,2)for h in HIGH)
        else:recipes.append((pi,a,b))
    assert symbolic==[[0,2,2],[0,3,2],[0,4,2],[0,5,2],[0,6,2],[1,2,2]]
    sheets=[];caps=[];counts=Counter();badcounts=Counter();unresolved=[];seen={};improvements=[]
    for pi,m59,m60 in recipes:
        rhs=np.array([-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]]for i,m in zip(ANCHORS,(m59,m60))])
        part=rhs@J@frame;remaining=s.profiles[pi]['square']-float(part@M@part);assert remaining>1e-9
        for sign in(-1,1):
            va=math.sqrt(remaining)/2*(sign*space+time);vb=math.sqrt(remaining)/2*(sign*space-time)
            coeff={i:[float(va@M@w),float(part@M@w)+s.profiles[pi]['offsets'][s.index[i]],float(vb@M@w)]for i,w in old.items()if i not in ANCHORS}
            center=[float(va@M@s.center),float(part@M@s.center)+s.profiles[pi]['center_offset'],float(vb@M@s.center)]
            row=dict(index=len(sheets),profile=pi,labels=[m59,m60],sheet=sign,part=part.tolist(),va=va.tolist(),vb=vb.tolist(),remaining_norm=remaining,
                     coefficients=sorted(coeff.items()),center_coefficients=center,events=[],continuous_samples=[])
            sheets.append(row);acute=[[0.,float('inf')]];dotted=[[0.,float('inf')]];steps=[]
            order=[i for i in[58,66,55,61,7]if i in coeff]+[i for i in sorted(coeff)if i not in[58,66,55,61,7]]+['center']
            for i in order:
                cc=center if i=='center'else coeff[i]
                intervals,status=qset(cc)
                if status=='near_tangent':unresolved.append(dict(stage='acute_inequality_near_tangent',sheet=row['index'],wall=i,coefficients=cc))
                # A true constant-zero old native product imposes no acuteness
                # restriction. This numerical identity convention is checked
                # again by all native gates for every accepted normal.
                if max(abs(x)for x in cc)<1e-10:intervals=[[0.,float('inf')]]
                acute=meet(acute,intervals);steps.append(dict(wall=i,after=encode(acute),status=status))
                if not acute:break
            row['acute_intersection_steps']=steps;row['acute_intervals']=encode(acute)
            # Strictly dotted to every unassigned old wall gives open native
            # families; exceptional constant native contacts are allowed.
            for i,cc in coeff.items():
                a,c,b=cc
                if max(abs(a),abs(b))<1e-10:
                    if T.G.relation(c)in(None,'parallel_or_unresolved'):dotted=[];break
                    continue
                intervals,status=qset([a,c+2,b]);dotted=meet(dotted,intervals)
                if status=='near_tangent':unresolved.append(dict(stage='dotted_inequality_near_tangent',sheet=row['index'],wall=i))
                if not dotted:break
            dotted=meet(dotted,qset(center)[0]);row['strict_dotted_intervals']=encode(dotted)
            tests=[]
            for lo,hi in dotted:
                t=math.sqrt(lo*hi)if lo>0 and math.isfinite(hi)else hi/2 if math.isfinite(hi)else max(1.,2*lo)
                tests.append(('dotted_interval',None,None,t))
            for wall,cc in coeff.items():
                for m in EVENT_LABELS:
                    a,c,b=cc;rr,status=roots([a,c+T.B.theta(m),b]);counts['event_equations']+=1
                    if status=='near_tangent':unresolved.append(dict(stage='native_event_near_tangent',sheet=row['index'],wall=wall,label=m));continue
                    for t in rr:
                        if t>0:tests.append(('native_event',wall,m,t))
            for kind,wall,label,t in tests:
                counts['positive_parameter_tests']+=1
                v=part+va*t+vb/t
                error=max(abs(float(v@M@v)-s.profiles[pi]['square']),max(abs(frame@M@v-rhs)))
                if not np.isfinite(error)or error>2e-6:
                    unresolved.append(dict(stage='normal_residual',sheet=row['index'],kind=kind,wall=wall,label=label,t=t,error=float(error)));continue
                vals,bad,cp=s.native(v,pi,old,{})
                record=dict(kind=kind,wall=wall,label=label,t=t,blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_error=error)
                if kind=='dotted_interval':row['continuous_samples'].append(record)
                if bad or cp>=-1e-8:
                    badcounts[tuple(bad)]+=1;counts['native_or_center_rejected']+=1
                    if len(bad)<=2:row['events'].append(record)
                    continue
                key=(pi,)+tuple(np.round(v,8));record['accepted']=True
                if key in seen:record['duplicate_cap']=seen[key];row['events'].append(record);counts['duplicate_caps']+=1;continue
                full={**old,67:v};sid=s.score(full,{67:pi},s.cur['actual_sides']+[67],'single59_60cap')
                if sid is None:continue
                cap=dict(index=len(caps),sheet_record=row['index'],profile=pi,parameter=t,normal=v.tolist(),native_products=sorted(vals.items()),
                    center_product=cp,norm_error=error,event=record,score=sid,finite_old_contacts=[[i,T.G.relation(x)]for i,x in vals.items()if isinstance(T.G.relation(x),int)])
                seen[key]=cap['index'];caps.append(cap);record['cap']=cap['index'];row['events'].append(record)
                if s.ledgers[sid]['open_ends']<10:improvements.append(cap['index'])
            print('SINGLE_EDGE_SHEET',row['index'],[pi,m59,m60],sign,'acute',row['acute_intervals'],'dotted',row['strict_dotted_intervals'],'caps',len(caps),flush=True)
    companion=companions(s,caps,sheets)
    paths=[Path(__file__),Path(T.__file__),Path(T.F.__file__),Path(T.G.__file__),Path(T.E.__file__),Path(T.N.__file__),Path(T.P.__file__),Path(T.H.__file__),T.H.PARENT,T.H.BOUNDARY,
           ROOT/'scripts/fresh_nonarithmetic_pi8_quad_completion_20261017.py',ROOT/'results/fresh_nonarithmetic_pi8_quad_completion_20261017.json']
    return dict(status='NUMERICAL_WALL66_SINGLE_EDGE_CONIC',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        retained_old_walls=63,deleted_walls=[],exact_deletion_debit=dict(finite_supports=0,compact_faces=0,wall8_supports=0),
        anchors=ANCHORS,anchor_native_label=bond,symbolic_finite_corner_domain=symbolic,high_anchor_label_samples=HIGH,event_labels=EVENT_LABELS,
        conic_basis=dict(frame=frame.tolist(),frame_gram=G.tolist(),time=time.tolist(),space=space.tolist()),
        sheets=sheets,caps=caps,counts=dict(counts),blockers=[dict(walls=list(k),count=v)for k,v in badcounts.most_common()],
        companion_experiment=companion,boundary_ledgers=s.ledgers,control_score=control,improvements=improvements,unresolved=unresolved+s.unresolved,
        scope='All63 oldwalls retained. Single mandatory finite59/60/cap corner; two conic sheets and numerical positive-parameter inequality domains. Low2..5 and four high samples for59 contact, finite third-wall labels2..5/6/8/12/24; no interval/all-label theorem. Every accepted cap passes all native offsets/strictcenter and complete numericalcone/independentPD support agreement. No exact newnormal or ambientvolume claim.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'caps',len(out['caps']),'improvements',out['improvements'],'unresolved',len(out['unresolved']),flush=True)
