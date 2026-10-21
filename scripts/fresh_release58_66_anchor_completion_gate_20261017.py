#!/usr/bin/env python3
"""Exact anchor67 provenance and complete nine-recipe continuation gate.

No saved numerical vector selects a branch or proves a native equality.
The frozen physical interval engine and finite-diagram classifier are reused.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_dotted_c10_domain_20261017 as D
import fresh_release58_66_seven_boundary_20261017 as N
import fresh_nonarithmetic_pi9_prism_completion_20261017 as C60
ROOT,P,U,I,Q,G=D.ROOT,D.P,D.U,D.I,D.Q,D.G
TAG='fresh_release58_66_anchor_completion_gate_20261017'
S=D.S

class Space:
    def __init__(self):
        self.d=D.data();self.z=self.d['z'];self.ids=self.d['ids']
        self.cp={i:[D.C.literal(self.z,i,j)for j in [5,6]]for i in self.ids if i not in S}
        # Old60 is never replaced downstream. Its original exact recipe fixes
        # this covector; later literal() routines only re-enclose it physically.
        assert {i:C60.PRESCRIBED[60][i]for i in S}=={0:2,1:2,2:2,3:2,5:2,6:5}
        assert all(G.overlap(a,b)for a,b in zip(self.cp[60],G.profile(6)))
        self.cp[60]=G.profile(6)
        self.g={(i,j):self.d['pair'](i,j)for i in self.cp for j in self.cp}
        ns=self.z['normals'];h=self.z['common_center'];pair=self.z['pair']
        wg=[[D.C.literal(self.z,i,j)for j in S]for i in S];wi,_=P.C53.inverse(wg)
        hc=[U.rnd(pair(h,ns[i]))for i in S];hcw=U.imv(wi,hc)
        self.hwhite=hcw[4:]
        self.projected={};self.hp={}
        for i in self.cp:
            cov=[D.C.literal(self.z,i,j)for j in S];wc=U.imv(wi,cov)
            self.projected[i]=[U.rnd(a-b)for a,b in zip(ns[i],P.combine(wc,[ns[j]for j in S]))]
            self.hp[i]=U.rnd(pair(h,ns[i])-U.idot(cov,hcw))
    def recipe(self,frame,pi,labels):
        cp=G.profile(pi);q=U.rnd(2-G.white(cp,cp))
        gram=[[self.g[i,j]for j in frame]for i in frame];inv,piv=P.C53.inverse(gram)
        rhs=[U.rnd(-Q.theta(m)-G.white(cp,self.cp[i]))for i,m in zip(frame,labels)]
        coeff=U.imv(inv,rhs);rem=U.rnd(q-U.idot(rhs,coeff))
        choices=[]
        for seed in self.cp:
            if seed in frame:continue
            cov=[self.g[seed,i]for i in frame];kc=U.imv(inv,cov)
            kq=U.rnd(self.g[seed,seed]-U.idot(cov,kc))
            if kq.l>0:choices.append((float((kq.h-kq.l)/kq.l),seed,kc,kq))
        assert choices,('positive complement unresolved',frame,pi,labels)
        _,seed,kc,kq=min(choices)
        proof=dict(frame=frame,profile=pi,labels=labels,gram=G.encmat(gram),inverse_pivots=piv,
                   projected_contacts=P.enc(rhs),kernel_seed=seed,kernel_square=kq.enc(),remaining_square=rem.enc())
        if rem.h<0:return proof,[]
        assert rem.l>0,('norm residual unresolved',proof)
        scale=Q.sqrt_closed(rem/kq);proof['positive_scale']=scale.enc()
        roots=[]
        for sign in[-1,1]:
            row={}
            for j in self.cp:
                col=[self.g[i,j]for i in frame]
                kb=kq if j==seed else U.rnd(self.g[seed,j]-U.idot(kc,col))
                row[j]=U.rnd(U.idot(coeff,col)+sign*scale*kb)
            for i,v in zip(frame,rhs):row[i]=v
            native={i:U.rnd(row[i]+G.white(cp,self.cp[i]))for i in row}
            for i,m in zip(frame,labels):native[i]=-Q.theta(m)
            hc=U.rnd(U.idot(coeff,[self.hp[i]for i in frame])+sign*scale*(self.hp[seed]-U.idot(kc,[self.hp[i]for i in frame])))
            center=U.rnd(hc+U.idot(cp,self.hwhite))
            root=P.combine(coeff,[self.projected[i]for i in frame])
            kv=[U.rnd(a-b)for a,b in zip(self.projected[seed],P.combine(kc,[self.projected[i]for i in frame]))]
            physical=[U.rnd(a+sign*scale*b)for a,b in zip(root,kv)]
            roots.append(dict(sign=sign,row=row,native=native,center=center,projected_center=hc,physical=physical,q=q,cp=cp))
        return proof,roots
    def add(self,name,root):
        for i,v in root['row'].items():self.g[name,i]=self.g[i,name]=v
        self.g[name,name]=root['q'];self.cp[name]=root['cp'];self.hp[name]=root['projected_center'];self.projected[name]=root['physical']


def rejection(v):
    if v.l>0:return dict(kind='strictly_obtuse',product=v.enc())
    if v.l>-1 and v.h<0:return dict(kind='gap_between_right_and_label3',product=v.enc())
    # Search for a witness gap only. Failure is unresolved, never accepted.
    for m in range(3,65):
        lo=-Q.theta(m+1);hi=-Q.theta(m)
        if lo.h<v.l and v.h<hi.l:return dict(kind='consecutive_native_angle_gap',labels=[m,m+1],lower=lo.enc(),upper=hi.enc(),product=v.enc())
    return None


def run():
    s=Space();retained=[i for i in s.ids if i not in[58,66]]
    proof,branches=s.recipe([30,59,60],0,[3,2,2]);assert len(branches)==2
    good=[];anchor_checks=[]
    prescribed={**{i:2 for i in S},30:3,59:2,60:2}
    for r in branches:
        gates=[];bad=[]
        for i in retained:
            if i in S:v=I()
            else:v=r['native'][i]
            if i in prescribed:
                assert v.l<=(-Q.theta(prescribed[i])).h and v.h>=(-Q.theta(prescribed[i])).l
                gates.append(dict(wall=i,label=prescribed[i],product=v.enc()))
            elif v.h<-2:gates.append(dict(wall=i,label='dotted',product=v.enc()))
            else:
                witness=rejection(v);assert witness is not None,('unresolved anchor native',r['sign'],i,v.enc())
                bad.append(dict(wall=i,**witness))
        anchor_checks.append(dict(sign=r['sign'],center=r['center'].enc(),gates=gates,rejections=bad))
        if not bad:
            assert r['center'].h<0;good.append(r)
    assert len(good)==1;anchor=good[0];s.add(67,anchor)
    removed={i:anchor['native'][i].enc()for i in[58,66]}
    assert rejection(anchor['native'][58]) and anchor['native'][66].l>0
    control=json.loads((ROOT/'results'/f'{D.C.TAG}.json').read_text())
    escape=G.decode(control['strict_escape']['point']);ev=U.rnd(s.z['pair'](escape,anchor['physical']))
    print('ANCHOR',anchor['sign'],'escape',float(ev.l),float(ev.h),flush=True)
    assert ev.h<0,'saved escape is cut; need a new complete boundary interpretation'
    # The old deletion ledger is exact incidence data, separately bound below.
    parent=json.loads((ROOT/'results'/f'{N.TAG}.json').read_text())
    debit=parent['exact_deletion_ledger'];assert debit['lost_finite_count']==9
    oldvertices=control['retained_old_finite_supports']+[r['support']for r in control['new_finite_vertices']]
    lost=[v for v in oldvertices if set(v)&{58,66}];assert len(lost)==9
    ns={i:s.z['normals'][i]for i in retained};ns[67]=anchor['physical']
    def literal(i,j):
        if i==j:return I(2)
        if 67 in(i,j):
            other=j if i==67 else i
            return I()if other in S else anchor['native'][other]
        return D.C.literal(s.z,i,j)
    newfinite=[];trials=[]
    for subset in combinations(sorted(prescribed),8):
        supp=sorted([67,*subset]);dots=[(i,j)for i,j in combinations(supp,2)if literal(i,j).h<-2]
        if dots:trials.append(dict(support=supp,dotted_pair=dots[0]));continue
        gram=[[literal(i,j)for j in supp]for i in supp];pivs=P.pd(gram);inv,_=P.C53.inverse(gram)
        hc=[U.rnd(s.z['pair'](s.z['common_center'],ns[i]))for i in supp];co=U.imv(inv,hc)
        pt=[U.rnd(a-b)for a,b in zip(s.z['common_center'],P.combine(co,[ns[i]for i in supp]))]
        sq=U.rnd(s.z['pair'](pt,pt));assert sq.h<0
        remaining={i:U.rnd(s.z['pair'](pt,n))for i,n in ns.items()if i not in supp};assert all(v.h<0 for v in remaining.values())
        newfinite.append(dict(support=supp,square=sq.enc(),remaining_products={str(i):v.enc()for i,v in remaining.items()}));trials.append(dict(support=supp,PD_pivots=pivs))
    assert len(newfinite)==2
    for wall,profile in[(30,0),(67,0),(60,6)]:
        assert [v.enc()for v in s.cp[wall]]==[v.enc()for v in G.profile(profile)],('literal profile',wall)
    bonds=[2]*9;bonds[1]=3;bonds[5]=2
    recipes=N.E.recipes((0,0,6),tuple(bonds),((0,1),(1,2)))
    expected=[(0,(2,2,2),()),(0,(2,3,2),()),(0,(2,4,2),()),(0,(2,5,2),()),(0,(3,2,2),()),(0,(4,2,2),()),(0,(5,2,2),()),(1,(2,2,2),()),(1,(3,2,2),())]
    assert recipes==expected
    records=[]
    for pi,ms,high in recipes:
        assert not high;pr,rr=s.recipe([30,67,60],pi,list(ms));out=[]
        for r in rr:
            witnesses=[]
            for i in retained:
                if i in S:continue
                w=rejection(r['native'][i])
                if w is not None:witnesses.append(dict(wall=i,**w))
            assert witnesses,('unresolved continuation',pi,ms,r['sign'])
            out.append(dict(sign=r['sign'],center=r['center'].enc(),rejections=witnesses))
        assert len(out)==2
        records.append(dict(recipe=pr,branches=out))
        print('FINITE_RECIPE',pi,ms,'branches',len(out),'blockers',[[x['wall']for x in b['rejections']]for b in out],flush=True)
    paths=[Path(__file__),Path(D.__file__),Path(D.C.__file__),Path(N.__file__),Path(N.E.__file__),Path(N.P.__file__),Path(C60.__file__),ROOT/'results'/f'{C60.TAG}.json',Path(Q.__file__),Path(G.__file__),ROOT/'results'/f'{N.TAG}.json',ROOT/'results'/f'{D.C.TAG}.json',N.H.BOUNDARY]
    return dict(status='PASS_EXACT_ANCHOR67_AND_NINE_RECIPE_OBSTRUCTION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      anchor_recipe=proof,anchor_branches=anchor_checks,chosen_sign=anchor['sign'],new_normal=P.enc(anchor['physical']),deleted_wall_products=removed,
      retained_walls=61,new_walls=1,total_walls=62,exact_deletion_ledger=debit,old_finite_count=199,lost_finite_supports=lost,retained_old_finite_count=190,new_finite_trials=trials,new_finite_vertices=newfinite,complete_finite_vertex_count=190+len(newfinite),
      strict_escape=dict(inherited_support=control['strict_escape']['support'],new_product=ev.enc(),conclusion='The inherited strict spacelike escape survives, proving infinite volume.'),
      finite_frame_profiles=[[30,0],[67,0],[60,6]],old60_profile_provenance=dict(source=C60.TAG,prescribed_S6_labels={str(i):C60.PRESCRIBED[60][i]for i in S}),symbolic_recipes=recipes,records=records,finite_recipe_count=len(records),norm_branches=sum(len(r['branches'])for r in records),
      scope='Exact norm2 anchor67 rightS6 with native contacts30:3,59:2,60:2; both algebraic signs checked against all61 retained walls and the strict center. The old58/66 deletion debit and complete new finite-vertex census are exact. The retained strict spacelike escape proves this62-wall seed has infinite volume. No saved numerical normal selects or validates a branch. Every finite-corner continuation in frame30/67/60 is rejected, conditional on the inherited complete seven-profile spherical classification. No exclusion of other incidence, other anchors or cap systems.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['complete_finite_vertex_count'],out['norm_branches'],flush=True)
