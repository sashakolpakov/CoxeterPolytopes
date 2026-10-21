#!/usr/bin/env python3
"""Changed-boundary continuations of the deletion-free C3 companions.

Frames are derived from pairs of actual open intrinsic edges, not from a
fixed old patch. Complete finite diagrams restrict all contact recipes.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_c3_companion_caps_20261017 as C
J,H,B,M,N,S,P,ROOT=C.J,C.H,C.B,C.M,C.N,C.S,C.P,C.ROOT
TAG='fresh_c3_companion_continuation_20261017'
HIGH=[6,8,12]

@lru_cache(None)
def recipes(ps,bonds,duties):
    rows=[]
    for pi,ms in product(range(7),product(range(2,7),repeat=3)):
        if all(P.corner3(ps[a],ps[b],pi,bonds[a*3+b],ms[a],ms[b])is not None for a,b in duties):
            rows.append((pi,ms,tuple(i for i,m in enumerate(ms)if m==6)))
    return rows

def run():
    src=ROOT/'results'/f'{C.TAG}.json';parent=json.loads(src.read_text())
    cur,labels,base,h,profiles,primary=C.setup();actual=cur['actual_sides']+[67,68]
    counts=Counter();fail=Counter();rows=[];unresolved=[];hist=Counter();improvements=[]
    def off(a,b):return float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    for pr in parent['companions']:
        nd=dict(base);nd[68]=np.array(pr['normal']);prids={7:4,68:pr['profile']}
        ps=lambda i:prids.get(i,0)
        opens=[tuple(e['sides'])for e in pr['boundary']['edges']if e['kind']in['one_open','two_open']and set(e['sides'])<=set([7,30,47,52,67,68])]
        frames={}
        for e,f in combinations(opens,2):
            ids=tuple(sorted(set(e)|set(f)))
            if len(ids)!=3 or 68 not in ids:continue
            duties=tuple(sorted([tuple(sorted(ids.index(i)for i in edge))for edge in[e,f]]))
            frames[ids]=tuple(sorted(set(frames.get(ids,()))|set(duties)))
        local=[];seen=set();drecords=[]
        for ids,duties in sorted(frames.items()):
            bonds=[2]*9
            for a,b in duties:
                m=B.relation(float(nd[ids[a]]@M@nd[ids[b]])+off(ps(ids[a]),ps(ids[b])))
                assert isinstance(m,int);bonds[a*3+b]=m
            dom=recipes(tuple(ps(i)for i in ids),tuple(bonds),duties)
            drecords.append(dict(frame=list(ids),corner_pairs=[list(t)for t in duties],recipes=len(dom),high_recipes=sum(bool(r[2])for r in dom)))
            for pi,low,high in dom:
                for hh in product(HIGH,repeat=len(high)):
                    ms=list(low)
                    for i,m in zip(high,hh):ms[i]=m
                    rhs=[-B.theta(m)-off(pi,ps(i))for i,m in zip(ids,ms)]
                    roots,info=N.norm_roots([nd[i]for i in ids],rhs,profiles[pi]['square']);counts['norm_equations']+=1
                    if info['status'].startswith('unresolved'):
                        unresolved.append(dict(parent=pr['index'],frame=ids,profile=pi,labels=ms,norm=info));continue
                    for sg,v in roots:
                        counts['norm_branches']+=1
                        err=max(abs(float(v@M@v)-profiles[pi]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip(ids,rhs)))
                        if err>1e-5:
                            unresolved.append(dict(parent=pr['index'],frame=ids,profile=pi,labels=ms,branch=sg,error=err));continue
                        vals={i:float(v@M@base[i])+(profiles[pi]['offsets'][labels.index(i)]if i!=67 else 0.)for i in labels+[67]}
                        vals[68]=float(v@M@nd[68])+off(pi,ps(68))
                        bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                        if bad:fail[tuple(bad)]+=1;continue
                        cp=float(v@M@h)+profiles[pi]['center_offset']
                        if cp>=-1e-8:counts['center_rejected']+=1;continue
                        key=(pi,)+tuple(np.round(v,8))
                        if key in seen:counts['duplicates']+=1;continue
                        seen.add(key);nn=dict(nd);nn[69]=v;sel=actual+[69]
                        bd=S.C.detailed_cone(nn,sel);ic=H.incidence_score([nn[i]for i in sel],sel,[float(nn[i]@M@nn[i])for i in sel])
                        good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
                        good=good and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
                        if not good:unresolved.append(dict(parent=pr['index'],stage='boundary',boundary=bd,PD_incidence=ic));continue
                        rec=dict(index=len(local),frame=list(ids),profile=pi,labels=ms,branch=sg,normal=v.tolist(),center_product=cp,norm_residual=err,
                          native_products=sorted(vals.items()),boundary=bd,PD_incidence=ic)
                        local.append(rec);counts['coherent_continuations']+=1
                        score=(bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']);hist[score]+=1
                        if bd['open_ends']<10:improvements.append(dict(parent=pr['index'],continuation=rec))
        joined=[]
        for a,b in combinations(local,2):
            counts['continuation_pairs_tested']+=1
            va,vb=np.array(a['normal']),np.array(b['normal'])
            ab=float(H.A.native(va,a['profile'],vb,b['profile'],profiles));rel=B.relation(ab)
            if rel in[None,'parallel_or_unresolved']:continue
            nn=dict(nd);nn[69]=va;nn[70]=vb;sel=actual+[69,70]
            bd=S.C.detailed_cone(nn,sel);ic=H.incidence_score([nn[i]for i in sel],sel,[float(nn[i]@M@nn[i])for i in sel])
            good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
            good=good and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
            if not good:unresolved.append(dict(parent=pr['index'],stage='pair_boundary',indices=[a['index'],b['index']],boundary=bd,PD_incidence=ic));continue
            rec=dict(indices=[a['index'],b['index']],native_product=ab,relation=rel,boundary=bd,PD_incidence=ic)
            joined.append(rec);counts['coherent_continuation_pairs']+=1
            if bd['open_ends']<10:improvements.append(dict(parent=pr['index'],continuation_pair=rec))
        rows.append(dict(parent=pr['index'],parent_profile=pr['profile'],parent_labels=pr['labels'],open_edges=opens,domains=drecords,continuations=local,coherent_pairs=joined))
        print('C3_CHANGED_BOUNDARY',pr['index'],'frames',len(frames),'accepted',len(local),'pairs',len(joined),'improvements',len(improvements),flush=True)
    paths=[Path(__file__),Path(C.__file__),src,Path(P.__file__),H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_C3_CHANGED_BOUNDARY_CONTINUATION',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      deleted_walls=[],retained_control_walls=63,C3_retained=True,parent_companions=len(parent['companions']),
      high_label_samples=HIGH,counts=dict(counts),complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=i,count=n)for(f,e,i),n in sorted(hist.items())],
      blocker_histogram=[dict(walls=list(k),count=n)for k,n in sorted(fail.items())],records=rows,improvements=improvements,unresolved=unresolved,
      scope='Numerical changed-incidence discovery at all54 saved C3 companions. Frames come from pairs of actual open edges sharing a facet in the new B component; all finite A2 profiles and isolated high-label shapes are classified first. High6/8/12 are bounded representatives only. Every accepted continuation has all65 prior native products/common-center checks and complete target cone/independent PD incidence agreement. No finite-volume or exact new chamber claim, and no arbitrary high-label exclusion. All original63 walls are retained; oldfinite/compact faces persist in any exact acute realization with the stated strict center.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
