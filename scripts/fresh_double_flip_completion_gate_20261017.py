#!/usr/bin/env python3
"""Exact fixed-prefix gates for W completion and the WZ->XB disk flip.

All endpoint labels are covered on six stated algebraic prefix branches.
Physical outward rational interval and finite-diagram helpers are inherited.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_old8_triangle_double_flip_20261017 as SC
import fresh_old30_connector_native_gate_20261017 as E
Q,P,U,I,G,ROOT=E.Q,E.P,E.U,E.I,E.G,E.ROOT
TAG='fresh_double_flip_completion_gate_20261017'

def topology():
    old=SC.topology();tri=old['triangles'];remove=[frozenset(t)for t in[['X','Z','W'],['B','Z','W']]]
    assert all(sum(frozenset(t)==s for t in tri)==1 for s in remove)
    assert not any({'X','B'}<=set(t)for t in tri)
    tri=[t for t in tri if frozenset(t)not in remove]+[['X','B','Z'],['X','B','W']]
    boundary=old['boundary'];new=list('ARCXYBZW');vertices=boundary+new
    edges=Counter(tuple(sorted(e))for t in tri for e in combinations(t,2));be={tuple(sorted(e))for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert len(tri)==len({frozenset(t)for t in tri})==19 and len(edges)==31
    assert {e for e,n in edges.items()if n==1}==be and all(n==(1 if e in be else 2)for e,n in edges.items())
    for v in vertices:
        link=defaultdict(set)
        for t in tri:
            if v in t:
                a,b=[x for x in t if x!=v];link[a].add(b);link[b].add(a)
        todo=[next(iter(link))];seen=set()
        while todo:
            x=todo.pop()
            if x not in seen:seen.add(x);todo.extend(link[x]-seen)
        assert len(seen)==len(link)
        assert sorted(map(len,link.values()))==([1,1]+[2]*(len(link)-2)if v in boundary else[2]*len(link))
    assert Counter(sum(v in new for v in t)for t in tri)=={1:5,2:8,3:6}
    assert not any({'W','Z'}<=set(t)for t in tri)
    return dict(old,triangles=tri,third_flip='WZ->XB')

def run():
    top=topology();parentpath=ROOT/'results'/f'{SC.TAG}.json';doc=json.loads(parentpath.read_text());assert len(doc['prefixes'])==6
    base=Q.Engine();_,_,ns,_,_,_=SC.H.setup();records=[];counts=Counter();high=I(Q.theta(6).l,F(2))
    def rejected(g,frame,pi,ms,box=None,allowed=None):
        contacts=[-box if m==6 else-Q.theta(m)for m in ms]
        rows,info=g.norm_frame(frame,pi,contacts)
        if not rows:return dict(reason='negative_remaining_square',frame_certificate=info)
        cp=G.profile(pi);branches=[]
        for sg,row in rows:
            bad=[]
            for i in E.TARGETS:
                value=U.rnd(row[i]+G.white(cp,g.cp[i]));ex=Q.exclusion(value)
                if ex:bad.append(dict(wall=i,**ex))
            if bad:branches.append(dict(sign=sg,reason='old_native_contact',rejection=bad[0]));continue
            if allowed is not None:
                value=U.rnd(row['Z']+G.white(cp,g.cp['Z']));separations=[]
                for m in allowed:
                    target=-Q.theta(m)
                    if value.h<target.l:side='strictly_below'
                    elif value.l>target.h:side='strictly_above'
                    else:return None
                    separations.append(dict(label=m,threshold=target.enc(),side=side))
                assert separations
                branches.append(dict(sign=sg,reason='required_finite_WZ_angle_absent',product=value.enc(),allowed_labels=allowed,separations=separations));continue
            return None
        return dict(reason='both_norm_branches_excluded',frame_certificate=info,branches=branches)
    def cover(g,frame,pi,ms,box=None,allowed=None,depth=0):
        counts['tested_boxes']+=1;out=rejected(g,frame,pi,ms,box,allowed)
        if out is not None:
            counts['certified_leaves']+=1;counts['maximum_depth']=max(counts['maximum_depth'],depth)
            return dict(cosine_box=box.enc()if box is not None else None,certificate=out)
        assert box is not None and depth<24,('unresolved_exact_completion_gate',frame,pi,ms,allowed,depth,box.enc()if box else None)
        mid=(box.l+box.h)/2
        return dict(split=str(mid),children=[cover(g,frame,pi,ms,I(box.l,mid),allowed,depth+1),cover(g,frame,pi,ms,I(mid,box.h),allowed,depth+1)])
    for j,rec in enumerate(doc['prefixes']):
        a,r,c,x,y,z=map(np.array,rec['normals']);assert rec['profiles']==[0,0,0,0,0,0]
        m8z,mrz,mxz=rec['Z_labels'];assert m8z==3 and mrz in[8,12]and mxz in[3,4,5]and rec['R44']==15
        g=E.Gram(base,ns)
        g.add('C',[7,30,10],0,[2,2,9],c)
        g.add('A',[7,8,'C'],0,[2,2,3],a)
        g.add('R',[8,44,63],0,[2,15,2],r)
        g.add('X',[8,'A','C'],0,[2,3,2],x)
        g.add('Z',[8,'R','X'],0,[m8z,mrz,mxz],z)
        assert all(F(s['remaining_square'][0])>0 for s in g.recipes)
        wcases=[];bcases=[]
        # Original disk: mandatory30/C/W,C/X/W,X/Z/W corners. XZ>2
        # excludes any WX or WZ label>=6; the other two labels retain
        # their full isolated-I2 tails through the symbolic label6.
        for pi,m30,mc,mx in product(range(7),range(2,7),range(2,7),range(2,6)):
            if SC.P.corner3(0,0,pi,2,m30,mc)is None or SC.P.corner3(0,0,pi,2,mc,mx)is None:continue
            allowed=[mz for mz in range(2,6)if SC.P.corner3(0,0,pi,mxz,mx,mz)is not None]
            if not allowed:continue
            ms=[m30,mc,mx];assert ms.count(6)<=1
            counts['W_high_families'if 6 in ms else'W_low_recipes']+=1
            proof=cover(g,[30,'C','X'],pi,ms,high if 6 in ms else None,allowed)
            wcases.append(dict(profile=pi,labels=ms,high_positions=[i for i,m in enumerate(ms)if m==6],allowed_WZ=allowed,proof=proof))
        # Flipped disk: RZ>=6 forces BR=BZ=2 in its finite corner.
        # XZ>=3 forces BX<=5 in the other finite corner.
        for pi,mx in product(range(7),range(2,6)):
            if SC.P.corner3(0,0,pi,mrz,2,2)is None or SC.P.corner3(0,0,pi,mxz,mx,2)is None:continue
            counts['B_recipes']+=1
            bcases.append(dict(profile=pi,labels=[2,2,mx],proof=cover(g,['R','Z','X'],pi,[2,2,mx])))
        records.append(dict(prefix=j,exact_prefix_recipes=g.recipes,W_cases=wcases,B_cases=bcases))
        print('DOUBLE_FLIP_GATE',j+1,dict(counts),flush=True)
    assert (counts['W_low_recipes'],counts['W_high_families'],counts['B_recipes'])==(368,22,60)
    assert counts['tested_boxes']==counts['certified_leaves']==450 and counts['maximum_depth']==0
    paths=[Path(__file__),Path(SC.__file__),parentpath,Path(E.__file__),Path(Q.__file__),Path(G.__file__),Path(G.E.X.__file__),G.CONTROL,G.PARENT,SC.H.PARENT,SC.H.BOUNDARY]
    return dict(status='PASS_EXACT_ALL_LABEL_W_AND_FLIPPED_B_GATES_ON_SIX_PREFIXES',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},topology=top,counts=dict(counts),records=records,targets=E.TARGETS,
      completeness='An actual S6 facet in a compact target face is rightH4 and has one of seven A2 profiles. For W, all finite30/C/W andC/X/W diagrams are enumerated with symbolic6 standing for any isolated bond>=6. The nonrightXZ bond excludes WX,WZ>=6, so WZ has a complete finite label list2..5 filtered by its corner. There is at most one remaining high label and its magnitude is relaxed to the entire closed interval[sqrt3,2]. For the WZ->XB flip, finiteR/Z/B with RZ8 or12 forcesBR=BZ=2, while finiteX/Z/B withXZ3/4/5 forcesBX<=5. This yields60 literal B recipes and both norm branches. Every branch is excluded by an old native obtuse/strict Coxeter-gap witness, negative remaining norm after a positive-complement signature guard, or separation from every permitted finite WZ label.',
      scope='Exact necessary-condition obstruction for the original W duty and changed B duty on six specified algebraic prefixes. The numerical parent is not promoted to an exact chamber, nor are its full contacts or its bounded prefix search independently certified. No continuously deformed prefix or other disk incidence is excluded. The all-label scope applies to completion walls at these fixed prefixes only. Shared physical interval and finite-diagram implementations are inherited. No new finite-volume result or improvement over exact63/199/10 control is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
