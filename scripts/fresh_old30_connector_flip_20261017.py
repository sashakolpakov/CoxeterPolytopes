#!/usr/bin/env python3
"""Subdivide BC by an old30 connector, then remove the independently blocked YB.

Exact disk combinatorics and bounded numerical geometry. No exact promotion.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_interior_connector_disk_20261017 as P
F,T,N,S,H,B,M,ROOT=P.F,P.T,P.N,P.S,P.H,P.B,P.M,P.ROOT
TAG='fresh_old30_connector_flip_20261017'
LABELS=range(2,7)

def topology(flip):
    old=P.disk();ts=[list(t)for t in old['triangles']]
    def replace(rem,add):
        nonlocal ts
        wanted=[frozenset(t)for t in rem]
        assert all(sum(frozenset(t)==w for t in ts)==1 for w in wanted)
        ts=[t for t in ts if frozenset(t)not in wanted]+add
    replace([['30','B','C'],['Z','B','C']],[['30','B','W'],['30','W','C'],['Z','B','W'],['Z','W','C']])
    if flip:
        assert not any({'R','Z'}<=set(t)for t in ts)
        replace([['Y','Z','B'],['Y','R','B']],[['Y','R','Z'],['R','B','Z']])
    boundary=old['boundary'];new=list('ARCXYBZW');edges=Counter(tuple(sorted(e))for t in ts for e in combinations(t,2));be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert {e for e,c in edges.items()if c==1}==be
    assert all(c==(1 if e in be else 2)for e,c in edges.items())
    assert len(ts)==len({frozenset(t)for t in ts})
    for v in boundary+new:
        lk=defaultdict(set)
        for t in ts:
            if v in t:
                a,b=[i for i in t if i!=v];lk[a].add(b);lk[b].add(a)
        todo=[next(iter(lk))];seen=set()
        while todo:
            a=todo.pop()
            if a not in seen:seen.add(a);todo.extend(lk[a]-seen)
        assert len(seen)==len(lk)
        assert sorted(map(len,lk.values()))==([1,1]+[2]*(len(lk)-2)if v in boundary else[2]*len(lk))
    xyz=Counter(sum(v in new for v in t)for t in ts)
    assert(len(boundary)+len(new),len(edges),len(ts))==(13,31,19)
    assert xyz=={1:5,2:7,3:7}
    return dict(boundary=boundary,triangles=ts,vertices=13,edges=31,triangles_count=19,x=5,y=7,z=7,conditional_delta_ends=-5,flip_YB_to_RZ=flip,
      scope='Exact abstract disk only; geometry must establish embedding, all stipulated corners, nonincidences and the retained outside30--44 face.')

def run():
    direct,flipped=topology(False),topology(True)
    current,labels,ns,h,profiles,primary=H.setup();actual=[i for i in current['actual_sides']if i not in S.DELETE]
    parent=ROOT/'results'/f'{P.TAG}.json';old=ROOT/'results/fresh_two_connector_coupled_disk_20261017.json';seeds=json.loads(parent.read_text());bank=json.loads(old.read_text())
    counts=Counter();blocks=Counter();unresolved=[];zbank=[];wbank=[];complete=[];direct_relations=Counter();prefix_scores=[];high_domains=[]
    def native(v,pi,assigned=()):
        if v@M@h+profiles[pi]['center_offset']>=-1e-8:return [('center',0.)]
        return T.gates(T.native_values(v,pi,labels,ns,profiles),assigned)
    def pair(v,pi,w,pj):return float(H.A.native(v,pi,w,pj,profiles))
    def label(v,pi,w,pj):return F.label(pair(v,pi,w,pj))
    def dotted(v,pi,w,pj):return v@M@w<-np.sqrt(profiles[pi]['square']*profiles[pj]['square'])-1e-7
    def nonincident(v,pi,allowed):
        return all(v@M@ns[i]<-np.sqrt(profiles[pi]['square']*(ns[i]@M@ns[i]))-1e-7 for i in actual if i not in allowed)
    def rhs(m,pi,pj):return -B.theta(m)-pair(np.zeros(4),pi,np.zeros(4),pj)
    for iseed,seed in enumerate(seeds['seeds']):
        a,c,x,y,r=map(np.array,seed['normals']);xp,yp=seed['X_profile'],seed['Y_profile'];m8,maxx,mcx=seed['X_labels'];my8,myr,myx=seed['Y_labels']
        matching=[]
        for ib,br in enumerate(bank['B_bank']):
            if br['R44']!=seed['R44']:continue
            bp=br['B_profile'];b=np.array(br['B_projection']);m30,m63,mrb=br['B_labels']
            assert not native(b,bp,[30,63])
            ryb=B.relation(pair(y,yp,b,bp));direct_relations[str(ryb)]+=1;counts['direct_YB_tests']+=1
            if isinstance(ryb,int):counts['direct_valid_finite_YB']+=1
            if not all(dotted(b,bp,v,pj)for v,pj in[(a,0),(c,0),(x,xp),(y,yp)]):continue
            if any(B.relation(pair(b,bp,v,pj))in[None,'parallel_or_unresolved']for v,pj in[(a,0),(c,0),(x,xp),(y,yp)]):continue
            if not nonincident(b,bp,[30,63]):continue
            counts['coherent_six_cap_prefixes']+=1;matching.append((ib,br,b,bp))
        if not matching:continue
        for zp,mzx,mzy,mzr in product(range(7),LABELS,LABELS,LABELS):
            if P.corner3(xp,yp,zp,myx,mzx,mzy)is None or P.corner3(yp,0,zp,myr,mzy,mzr)is None:continue
            if 6 in[mzx,mzy,mzr]:high_domains.append(dict(seed=iseed,profile=zp,labels=[mzx,mzy,mzr],scope='6 is sampled; corresponding isolated labels>=7 remain unsearched.'))
            zz,info=N.norm_roots([x,y,r],[rhs(mzx,zp,xp),rhs(mzy,zp,yp),rhs(mzr,zp,0)],profiles[zp]['square']);counts['Z_equations']+=1
            if info['status'].startswith('unresolved'):unresolved.append(['Z',iseed,zp,mzx,mzy,mzr,info]);continue
            for zsg,z in zz:
                counts['Z_branches']+=1;bad=native(z,zp)
                if bad:blocks[('Z',str(bad[0][0]))]+=1;continue
                if not nonincident(z,zp,[])or not dotted(z,zp,a,0):counts['Z_extra_intrinsic_contact']+=1;continue
                if B.relation(pair(z,zp,a,0))in[None,'parallel_or_unresolved']:counts['Z_A_native_gap']+=1;continue
                counts['native_Z']+=1;mzc=label(z,zp,c,0)
                if mzc is None:counts['Z_C_not_finite']+=1;continue
                if P.corner3(xp,0,zp,mcx,mzx,mzc)is None:counts['Z_C_corner_not_finite']+=1;continue
                iz=len(zbank);zbank.append(dict(seed=iseed,profile=zp,labels=[mzx,mzy,mzr,mzc],branch=zsg,projection=z.tolist()));counts['native_XYRCZ']+=1
                for ib,br,b,bp in matching:
                    m30,m63,mrb=br['B_labels'];mzb=label(z,zp,b,bp)
                    if mzb is None:counts['Z_B_not_finite']+=1;continue
                    if P.corner3(0,bp,zp,mrb,mzr,mzb)is None:continue
                    counts['coherent_seven_cap_prefixes']+=1
                    for wp,mw30,mwc,mwb in product(range(7),LABELS,LABELS,LABELS):
                        if P.corner3(0,0,wp,2,mw30,mwc)is None or P.corner3(0,bp,wp,m30,mw30,mwb)is None:continue
                        ww,wi=N.norm_roots([ns[30],c,b],[rhs(mw30,wp,0),rhs(mwc,wp,0),rhs(mwb,wp,bp)],profiles[wp]['square']);counts['W_equations']+=1
                        if wi['status'].startswith('unresolved'):unresolved.append(['W',iz,ib,wp,mw30,mwc,mwb,wi]);continue
                        for wsg,w in ww:
                            bad=native(w,wp,[30])
                            if bad:blocks[('W',str(bad[0][0]))]+=1;continue
                            if not nonincident(w,wp,[30])or not all(dotted(w,wp,v,pj)for v,pj in[(a,0),(r,0),(x,xp),(y,yp)]):continue
                            if any(B.relation(pair(w,wp,v,pj))in[None,'parallel_or_unresolved']for v,pj in[(a,0),(r,0),(x,xp),(y,yp)]):continue
                            counts['native_W']+=1;mwz=label(w,wp,z,zp)
                            if mwz is None:counts['W_Z_not_finite']+=1;continue
                            if P.corner3(bp,zp,wp,mzb,mwb,mwz)is None or P.corner3(0,zp,wp,mzc,mwc,mwz)is None:continue
                            new=[a,r,c,x,y,b,z,w];ps=[0,0,0,xp,yp,bp,zp,wp];newlabels=list(range(67,75));nd=dict(ns);nd.update(zip(newlabels,new));checks=[]
                            for i,j in combinations(range(8),2):
                                v=pair(new[i],ps[i],new[j],ps[j]);rel=B.relation(v);checks.append(dict(pair=[i,j],product=v,relation=rel))
                            if any(q['relation']is None or q['relation']=='parallel_or_unresolved'for q in checks):counts['full_pair_gate_failure']+=1;continue
                            bd=S.C.detailed_cone(nd,actual+newlabels);inc=H.incidence_score([nd[i]for i in actual+newlabels],actual+newlabels,[float(nd[i]@M@nd[i])for i in actual+newlabels])
                            if bd['warnings']or inc['status']!='numerically_separated_acute_incidence' or bd['open_ends']!=inc['open_ends']:unresolved.append(['complete_boundary',iz,ib,wp,mw30,mwc,mwb,bd,inc]);continue
                            # Every required triangle is an actual finite support in the PD census.
                            mp={str(i):i for i in[7,8,44,63,30]}|dict(zip(list('ARCXYBZW'),newlabels));supported={tuple(sorted(v))for v in inc['finite_supports']}
                            assert all(tuple(sorted(mp[v]for v in tri))in supported for tri in flipped['triangles'])
                            counts['coherent_complete_disks']+=1
                            complete.append(dict(Z_index=iz,B_index=ib,W_profile=wp,W_labels=[mw30,mwc,mwb,mwz],W=w.tolist(),native_pairs=checks,boundary=bd,PD_incidence=inc))
        print('OLD30_CONNECTOR',iseed+1,dict(counts),flush=True)
    assert counts['direct_YB_tests']==1985 and counts['direct_valid_finite_YB']==0
    paths=[Path(__file__),Path(P.__file__),Path(F.__file__),Path(T.__file__),Path(N.__file__),Path(H.__file__),parent,old,H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_OLD30_CONNECTOR_AND_YB_FLIP',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},direct_subdivision_topology=direct,flipped_topology=flipped,counts=dict(counts),direct_YB_relations=dict(direct_relations),blockers=[[list(k),v]for k,v in sorted(blocks.items())],Z_bank=zbank,complete=complete,unresolved=unresolved,symbolic_high_Z_domains=high_domains,
      deleted_walls=S.DELETE,exact_old_deletion_cost=seeds['exact_old_deletion_cost'],domain=dict(prefixes=20,B_bank=422,Z_W_A2_profiles=list(range(7)),assigned_labels=list(LABELS)),
      predecessor_prose_correction=dict(report=f'{P.TAG}.md',incorrect_sentence='Every accepted profile happens to have zero A2 component.',actual_profile_counts={'0,0':14,'0,1':2,'0,4':2,'1,0':2},scope='Saved source/data already used the actual profiles and remain unchanged.'),
      scope='Exact disk combinatorics, bounded numerical equations. Direct BC subdivision retains a YB contact absent in all saved matching pairs; the adjacent YB-to-RZ flip is tested separately. Higher isolated labels remain open; a symbolic6-domain record is not continuous coverage. Full native/common-center and intrinsic nonincidence gates precede any numerical cone score. No exact geometry or finite-volume promotion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],len(out['complete']),flush=True)
