#!/usr/bin/env python3
"""Changed seven-edge cavity after jointly releasing58/66 and adding67.

Numerical discovery with all retained native gates and complete S6 scoring.
Three-cap arc lengths2/2/3 and four-cap2/2/2/1, all seven rotations.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_release66_three_cap_disk_20261017 as F
H,E,N,S,P,B,M,G,ROOT=F.H,F.E,F.N,F.S,F.P,F.B,F.M,F.G,F.ROOT
TAG='fresh_release58_66_seven_boundary_20261017'
HIGH=[6,8,12,24]
CYCLE=[7,33,30,67,60,59,55]


def topology(rotation,ncaps,diagonal=0):
    bc=CYCLE[rotation:]+CYCLE[:rotation];cs=[68,69,70,71][:ncaps]
    lengths=[2,2,3]if ncaps==3 else[2,2,2,1]
    fs=[];pos=0;starts=[]
    for cap,count in zip(cs,lengths):
        starts.append(pos)
        for k in range(count):fs.append([bc[(pos+k)%7],bc[(pos+k+1)%7],cap])
        pos+=count
    for i,cap in enumerate(cs):fs.append([bc[starts[i]],cs[i-1],cap])
    if ncaps==3:fs.append(cs)
    elif diagonal==0:fs += [[68,69,70],[68,70,71]]
    else:fs += [[68,69,71],[69,70,71]]
    es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2));vs=set(sum(fs,[]))
    assert(len(vs),len(es),len(fs))==((10,20,11)if ncaps==3 else(11,23,13))
    assert len(vs)-len(es)+len(fs)==1 and set(es.values())=={1,2}
    assert{e for e,k in es.items()if k==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
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
    return dict(rotation=rotation,caps=ncaps,diagonal=diagonal,boundary=bc,arc_lengths=lengths,triangles=fs,vertices=len(vs),edges=len(es),faces=len(fs),links=links)


class Scout:
    def __init__(self):
        self.cur,self.labels,old,h,self.profiles,self.primary=H.setup()
        boost=np.eye(4);boost[0,0]=h[0];boost[0,1:]=h[1:];boost[1:,0]=h[1:]
        boost[1:,1:]+=np.outer(h[1:],h[1:])/(h[0]+1)
        assert np.max(abs(boost.T@M@boost-M))<1e-10
        inverse=M@boost.T@M
        self.old={i:inverse@v for i,v in old.items()};self.center=inverse@h
        self.change=dict(boost=boost.tolist(),inverse=inverse.tolist(),center=self.center.tolist())
        self.index={i:j for j,i in enumerate(self.labels)}
        self.types={7:4,60:6};self.counts=Counter();self.blockers=Counter();self.near=[]
        for i in set(CYCLE)-{67}:
            assert abs(float(self.old[i]@M@self.old[i])-self.profiles[self.types.get(i,0)]['square'])<1e-5
        self.unresolved=[];self.ledgers=[];self.ledgerkeys={};self.hist=Counter();self.improvements=[]
        self.control=self.score(self.old,{},self.cur['actual_sides'],'control')
        assert self.ledgers[self.control]['counts']==dict(finite=26,exterior=6)and self.ledgers[self.control]['open_ends']==10
        self.deletion=H.deletion(self.primary,self.cur,[58,66])
        assert self.deletion['lost_finite_count']==9 and self.deletion['compact_faces_by_dimension']=={2:6,3:1}
        self.base={i:v for i,v in self.old.items()if i not in[58,66]}
        self.actual=[i for i in self.cur['actual_sides']if i not in[58,66]]
        self.deleted=self.score(self.base,{},self.actual,'deletion')
        assert self.ledgers[self.deleted]['counts']==dict(finite=21,exterior=7)and self.ledgers[self.deleted]['open_ends']==11
    def ps(self,i,types):return types.get(i,self.types.get(i,0))
    def off(self,a,b):return float(H.A.native(np.zeros(4),a,np.zeros(4),b,self.profiles))
    def dot(self,va,pa,vb,pb):return float(va@M@vb)+self.off(pa,pb)
    def score(self,nd,types,ids,stage):
        bd=S.C.detailed_cone(nd,ids)
        squares=[self.profiles[types[i]]['square']if i in types else float(self.old[i]@M@self.old[i])for i in ids]
        ic=H.incidence_score([nd[i]for i in ids],ids,squares)
        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'
        good=good and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
        good=good and{tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(s))for s in ic['finite_supports']}
        if not good:
            self.unresolved.append(dict(stage=stage,boundary=bd,PD_incidence=ic));return None
        row=G.F.D.compact_score(dict(boundary=bd,PD_incidence=ic));key=json.dumps(row,sort_keys=True)
        if key not in self.ledgerkeys:self.ledgerkeys[key]=len(self.ledgers);self.ledgers.append(row)
        self.hist[stage,row['counts']['finite'],row['counts']['exterior'],row['open_ends']]+=1
        return self.ledgerkeys[key]
    def native(self,v,pi,nd,types):
        vals={i:float(v@M@w)+(self.profiles[pi]['offsets'][self.index[i]]if i in self.index else self.off(pi,types[i]))for i,w in nd.items()}
        bad=[i for i,x in vals.items()if G.relation(x)in[None,'parallel_or_unresolved']]
        return vals,bad,float(v@M@self.center)+self.profiles[pi]['center_offset']
    def bank(self,frame,duties,nd,types,ids,newid,stage,score=True):
        bonds=[2]*9
        for a,b in duties:
            m=G.relation(self.dot(nd[frame[a]],self.ps(frame[a],types),nd[frame[b]],self.ps(frame[b],types)))
            assert isinstance(m,int),(frame,duties,a,b,m)
            bonds[a*3+b]=m
        dom=E.recipes(tuple(self.ps(i,types)for i in frame),tuple(bonds),tuple(duties))
        rows=[];seen=set()
        for pi,low,high in dom:
            for hh in product(HIGH,repeat=len(high)):
                ms=list(low)
                for i,m in zip(high,hh):ms[i]=m
                rhs=[-B.theta(m)-(self.profiles[pi]['offsets'][self.index[i]]if i in self.index else self.off(pi,types[i]))for i,m in zip(frame,ms)]
                roots,info=N.norm_roots([nd[i]for i in frame],rhs,self.profiles[pi]['square']);self.counts['norm_equations']+=1
                if info['status'].startswith('unresolved'):
                    self.unresolved.append(dict(stage=stage,frame=frame,profile=pi,labels=ms,norm=info));continue
                for sign,v in roots:
                    self.counts['norm_branches']+=1
                    err=max(abs(float(v@M@v)-self.profiles[pi]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip(frame,rhs)))
                    if err>2e-5:
                        self.unresolved.append(dict(stage=stage,frame=frame,profile=pi,labels=ms,branch=sign,error=err));continue
                    vals,bad,cp=self.native(v,pi,nd,types)
                    if bad or cp>=-1e-8:
                        self.blockers[stage,tuple(bad)]+=1
                        if len(bad)<=2:self.near.append(dict(stage=stage,frame=frame,profile=pi,labels=ms,branch=sign,blockers=[[i,vals[i]]for i in bad],center_product=cp))
                        continue
                    key=(pi,)+tuple(np.round(v,8))
                    if key in seen:self.counts['duplicate_branches']+=1;continue
                    seen.add(key);rec=dict(index=len(rows),frame=frame,profile=pi,labels=ms,branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err)
                    if score:
                        nn=dict(nd);nn[newid]=v;tt=dict(types);tt[newid]=pi
                        sid=self.score(nn,tt,ids+[newid],stage)
                        if sid is None:continue
                        rec['score']=sid
                        if self.ledgers[sid]['open_ends']<10:self.improvements.append(dict(stage=stage,cap=rec))
                    rows.append(rec);self.counts['native_caps']+=1
        return dict(frame=frame,duties=duties,symbolic_recipes=[list(x)for x in dom],caps=rows)
    def pair(self,a,b):return G.relation(self.dot(np.array(a['normal']),a['profile'],np.array(b['normal']),b['profile']))
    def corner(self,nodes,nd,types):
        ms=[]
        for a,b in combinations(nodes,2):
            m=G.relation(self.dot(nd[a],self.ps(a,types),nd[b],self.ps(b,types)))
            if not isinstance(m,int):return None
            ms.append(m)
        return P.corner3(*(self.ps(i,types)for i in nodes),*ms)
    def add(self,nd,types,rs,start=68):
        nn=dict(nd);tt=dict(types)
        for i,r in enumerate(rs):nn[start+i]=np.array(r['normal']);tt[start+i]=r['profile']
        return nn,tt


def run():
    s=Scout();seedbank=s.bank([30,59,60],[(0,1),(1,2)],s.base,{},s.actual,67,'anchor_domain')
    anchors=[r for r in seedbank['caps']if r['profile']==0 and r['labels']==[3,2,2]]
    assert len(anchors)==1;anchor=anchors[0]
    nd=dict(s.base);nd[67]=np.array(anchor['normal']);types={67:0};actual=s.actual+[67]
    start=s.ledgers[anchor['score']];assert start['counts']==dict(finite=23,exterior=7)and start['open_ends']==11
    edges={tuple(sorted(e['sides'])):e['kind']for e in start['edges']}
    cycleedges={tuple(sorted(e))for e in zip(CYCLE,CYCLE[1:]+CYCLE[:1])}
    assert all(edges[e]=='one_open'for e in cycleedges)
    remaining={e for e,k in edges.items()if k in['one_open','two_open']}-cycleedges
    assert remaining=={(7,47),(7,30),(30,52),(47,52)}
    banks=[]
    for i in range(7):
        frame=[CYCLE[(i+j)%7]for j in range(3)]
        b=s.bank(frame,[(0,1),(1,2)],nd,types,actual,68,'wedge'+str(i));banks.append(b)
        print('WEDGE',i,frame,'domain',len(b['symbolic_recipes']),'caps',len(b['caps']),'unresolved',len(s.unresolved),flush=True)
    # The two nonempty banks can also coexist without fitting either disk.
    # Score every compatible pair, including dotted pairs, before interpreting
    # their finite-corner gains as any progress on the full boundary.
    free_pairs=[];free_hist=Counter()
    for i,j in combinations(range(7),2):
        for ra,rb in product(banks[i]['caps'],banks[j]['caps']):
            value=s.dot(np.array(ra['normal']),ra['profile'],np.array(rb['normal']),rb['profile'])
            rel=G.relation(value);free_hist[str(rel)]+=1;s.counts['unrestricted_bank_pair_tests']+=1
            if rel in[None,'parallel_or_unresolved']:continue
            nn,tt=s.add(nd,types,[ra,rb]);sid=s.score(nn,tt,actual+[68,69],'bank_pair')
            if sid is None:continue
            free_pairs.append(dict(banks=[i,j],indices=[ra['index'],rb['index']],native_product=value,relation=rel,score=sid))
            if s.ledgers[sid]['open_ends']<10:s.improvements.append(dict(stage='bank_pair',banks=[i,j],roots=[ra,rb],score=sid))
    print('ALL_BANK_PAIRS',dict(free_hist),'scored',len(free_pairs),'improvements',len(s.improvements),flush=True)
    joins=[];triple_rows=[];four_rows=[];pairhist=Counter();diagonalhist=Counter()
    for rot in range(7):
        bc=CYCLE[rot:]+CYCLE[:rot];bs=[banks[(rot+2*j)%7]['caps']for j in range(3)]
        pairmaps=[];jrows=[]
        for a,b,junction in[(0,1,bc[2]),(1,2,bc[4])]:
            accepted={};rows=[]
            for ia,ra in enumerate(bs[a]):
                for ib,rb in enumerate(bs[b]):
                    rel=s.pair(ra,rb);pairhist[str(rel)]+=1;s.counts['adjacent_pair_tests']+=1
                    if not isinstance(rel,int):continue
                    nn,tt=s.add(nd,types,[ra,rb])
                    co=s.corner([junction,68,69],nn,tt)
                    if co is None:continue
                    sid=s.score(nn,tt,actual+[68,69],'finite_pair')
                    if sid is None:continue
                    accepted[ia,ib]=rel;rows.append(dict(indices=[ia,ib],label=rel,corner=co,score=sid))
                    if s.ledgers[sid]['open_ends']<10:s.improvements.append(dict(stage='finite_pair',rotation=rot,sides=[a,b],roots=[ra,rb],score=sid))
            pairmaps.append(accepted);jrows.append(dict(sides=[a,b],junction=junction,pairs=rows))
        joins.append(dict(rotation=rot,joins=jrows));dcache={}
        for(ia,ib),mab in pairmaps[0].items():
            for ic in range(len(bs[2])):
                if(ib,ic)not in pairmaps[1]:continue
                rs=[bs[0][ia],bs[1][ib],bs[2][ic]];mac=s.pair(rs[0],rs[2]);diagonalhist[str(mac)]+=1
                if mac in[None,'parallel_or_unresolved']:continue
                nn,tt=s.add(nd,types,rs)
                sid=s.score(nn,tt,actual+[68,69,70],'three_chain')
                if sid is None:continue
                s.counts['coherent_three_chains']+=1
                if s.ledgers[sid]['open_ends']<10:s.improvements.append(dict(stage='three_chain',rotation=rot,roots=rs,score=sid))
                # The third cap additionally takes the seventh boundary edge.
                target=topology(rot,3)
                corners=[s.corner(t,nn,tt)for t in target['triangles']]
                if all(c is not None for c in corners):
                    row=completed(target,rs,sid,s);triple_rows.append(row)
                # Fourth cap takes that edge; both interior diagonals are tested.
                if ic not in dcache:
                    zn=dict(nd);zn[70]=np.array(rs[2]['normal']);zt=dict(types);zt[70]=rs[2]['profile']
                    dcache[ic]=s.bank([bc[6],bc[0],70],[(0,1),(0,2)],zn,zt,actual+[70],71,'fourth'+str(rot),score=True)
                for dr in dcache[ic]['caps']:
                    full,ft=s.add(nd,types,rs+[dr]);vals,bad,cp=s.native(full[71],ft[71],{i:v for i,v in full.items()if i!=71},ft)
                    if bad or cp>=-1e-8:continue
                    for diag in[0,1]:
                        target=topology(rot,4,diag)
                        corners=[s.corner(t,full,ft)for t in target['triangles']]
                        if not all(c is not None for c in corners):continue
                        sid4=s.score(full,ft,actual+[68,69,70,71],'four_caps')
                        if sid4 is None:continue
                        row=completed(target,rs+[dr],sid4,s);four_rows.append(row)
                        if s.ledgers[sid4]['open_ends']<10:s.improvements.append(dict(stage='four_caps',rotation=rot,roots=rs+[dr],score=sid4,completed=row['completed_disk']))
        joins[-1]['fourth_banks']=[dict(third_index=i,bank=b)for i,b in dcache.items()]
        print('ROTATION',rot,'banks',list(map(len,bs)),'joins',list(map(len,pairmaps)),'third',s.counts['coherent_three_chains'],'four',len(four_rows),'improvements',len(s.improvements),flush=True)
    paths=[Path(__file__),Path(F.__file__),ROOT/'results'/f'{F.TAG}.json',Path(H.__file__),Path(E.__file__),Path(N.__file__),Path(S.__file__),Path(P.__file__),Path(G.__file__),Path(G.F.D.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_RELEASE58_66_SEVEN_BOUNDARY_SCOUT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      exact_deletion_ledger=s.deletion,retained_control_walls=61,coordinate_change=s.change,control_score=s.control,deletion_score=s.deleted,
      anchor_domain=seedbank,chosen_anchor_index=anchor['index'],boundary_cycle=CYCLE,independent_B_open_edges=sorted(remaining),
      topologies=[topology(r,n,d)for r in range(7)for n in[3,4]for d in(range(2)if n==4 else[0])],high_label_samples=HIGH,
      wedge_banks=banks,all_bank_pair_histogram=dict(free_hist),all_coherent_bank_pairs=free_pairs,
      joins=joins,three_cap_records=triple_rows,four_cap_records=four_rows,counts=dict(s.counts),
      adjacent_pair_histogram=dict(pairhist),diagonal_histogram=dict(diagonalhist),boundary_ledgers=s.ledgers,improvements=s.improvements,unresolved=s.unresolved,
      histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
      blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],near_native_failures=s.near,
      scope='Changed actual geometry after releasing58/66 and adding the formerly58-only-blocked67 from30/59/60 with profile0 labels3/2/2. Exact inherited deletion debit; bounded numerical seven-boundary three/four-cap disk search, all seven rotations and both four-cap diagonals. SevenA2profiles, low2..5 and high6/8/12/24 only. All61 retained native products, all prior/new contacts and commoncenter; complete numericalS6cone/independentPD support checks for every accepted cap and coherent finite-pair/chain/completion. No exact new normal, all-label exclusion, ambient finite-census or finite-volume assertion; frozen63/199/10 infinite-volume control unchanged.')


def completed(topo,roots,sid,s):
    bd=s.ledgers[sid];finite={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'};whole={}
    for i in range(68,68+len(roots)):
        expected={tuple(sorted(t))for t in topo['triangles']if i in t};rays=[r for r in bd['rays']if i in r['active']]
        whole[i]=all(r['kind']=='finite'for r in rays)and{tuple(sorted(r['active']))for r in rays}==expected
    return dict(topology=topo,roots=roots,score=sid,whole_faces=whole,completed_disk=all(tuple(sorted(t))in finite for t in topo['triangles'])and all(whole.values()))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
