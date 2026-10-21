#!/usr/bin/env python3
"""Two-old-contact conics meeting one prior trial cap at a finite angle."""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import numpy as np
import mpmath as mp
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csr_matrix
import fresh_nonarithmetic_h3_two_contact_20261016 as C
import fresh_nonarithmetic_h3_two_contact_20261016_pairs as B
ROOT=C.ROOT;TAG=C.TAG+'_trial_extension';ETA=C.ETA
TRIAL='results/fresh_nonarithmetic_h3_closure_20261016_trial.json'

def run(iterations):
    d=json.loads((ROOT/f'results/{C.TAG}.json').read_text());tr=json.loads((ROOT/TRIAL).read_text())
    Q,N,H,E=C.P.geometry();basis=[mp.matrix([mp.mpf(x)for x in v])for v in d['basis']]
    def coord(v):return np.array([float(-(v.T*Q*basis[0])[0])]+[float((v.T*Q*z)[0])for z in basis[1:]])
    extra=np.array([coord(mp.matrix([mp.mpf(x)for x in v]))for v in tr['normal_intervals_or_proposals']])
    NC=np.array(d['old_projected_normals']);allN=np.vstack([NC,extra]);initial=[]
    assert len(extra)==6 and len(allN)==61
    for i,v in enumerate(extra):
        assert v[0]>0 and abs(v@ETA@v-2)<2e-5
        assert all(C.compatible(x,3e-5)for x in NC@ETA@v)
        assert all(C.compatible(w@ETA@v,3e-5)for w in extra[:i])
    for r in tr['rays']:
        v=coord(mp.matrix([mp.mpf(x)for x in r['point']]));v/=v[0];active=[i+2 if i>=53 else i for i in r['active']]
        assert max(allN@ETA@v)<3e-5
        initial.append(dict(mask=sum(1<<i for i in active),v=v,kind=r['kind']))
    pts=np.array([np.r_[1,r['v'][1:]/np.linalg.norm(r['v'][1:])]for r in initial if r['kind']=='hyperideal_rays'])
    stats=Counter();pool={}
    blocked=json.loads((ROOT/f'results/{C.TAG}_affine_gate.json').read_text())['excluded_family_indices']
    def admit(f,t,recipe):
        v=C.evaluate(f,t);g=allN@ETA@v
        if v[0]<1e-7 or abs(v@ETA@v-2)>2e-5 or not all(C.compatible(x,3e-5)for x in g):return
        cover=np.flatnonzero(pts@ETA@v>1e-7).tolist()
        if not cover:return
        key=tuple(np.round(v,7));pool.setdefault(key,dict(family=f['index'],parameter=t,recipe=recipe,vector=v,coverage=cover,all61_products=g))
    for fi,f in enumerate(d['families']):
        if f['index']in blocked:continue
        intervals=[r['interval']for r in f['intervals']]
        for j,n in enumerate(extra):
            for m in range(2,11):
                stats['finite_equations']+=1;goal=0. if m==2 else-2*math.cos(math.pi/m)
                qq=C.poly(f,ETA@n,-goal)
                if max(abs(qq))<1e-10:stats['identical_contact']+=1;continue
                for t in C.roots(qq):
                    if B.inside(t,intervals,f['kind']=='hyperbola'):admit(f,t,dict(trial_wall=55+j,label=m))
        # Generic all-trial-dotted intervals are also represented in every
        # initial replacement-ideal coverage cell; no trial finite contact.
        polys=[C.poly(f,ETA@n,2.)for n in extra]
        for interval in C.intervals(polys,f['kind']=='hyperbola'):
            for old in intervals:
                lo=max([x for x in(interval[0],old[0])if x is not None],default=None)
                hi=min([x for x in(interval[1],old[1])if x is not None],default=None)
                if lo is not None and hi is not None and lo>=hi-1e-9:continue
                cuts=sorted({t for p in pts for t in C.roots(C.poly(f,ETA@p))if(lo is None or t>lo+1e-8)and(hi is None or t<hi-1e-8)})
                cuts=[lo]+cuts+[hi]
                for a,b in zip(cuts,cuts[1:]):
                    samples=C.parameter_samples([a,b],f['kind']=='hyperbola')
                    for t in samples[::max(1,len(samples)//3)]:admit(f,t,dict(all_trial_contacts='strictly_dotted'))
        if(fi+1)%100==0:print('families',fi+1,'pool',len(pool),flush=True)
    stats['new_two_contact_candidates']=len(pool)
    oldbank=json.loads((ROOT/'results/fresh_nonarithmetic_h3_joint_cover_20261016.json').read_text())
    for index,r in enumerate(oldbank['pool']):
        anchors=[i+2 if i>=53 else i for i in r['anchors']];A=allN[anchors];J=np.linalg.inv(A@ETA@A.T)
        c=np.array([0. if m==2 else-2*math.cos(math.pi/m)for m in r['labels']]);p=c@J@A
        h=np.array([1.,0,0,0]);k=h-(h@ETA@A.T)@J@A;rho=(2-p@ETA@p)/(k@ETA@k)
        assert rho>0;v=p+r['branch']*math.sqrt(rho)*k;g=allN@ETA@v
        assert v[0]>0 and abs(v@ETA@v-2)<2e-5 and all(C.compatible(x,3e-5)for x in g)
        cover=np.flatnonzero(pts@ETA@v>1e-7).tolist()
        if not cover:continue
        pool.setdefault(tuple(np.round(v,7)),dict(old_pool_index=index,recipe=dict(anchors=anchors,labels=r['labels'],branch=r['branch']),vector=v,coverage=cover,all61_products=g))
    stats['union_candidates']=len(pool)
    rows=list(pool.values());V=np.array([r['vector']for r in rows]);samples=list(pts);conflicts=set();history=[];proposal=None
    def fits(g):
        lab=np.pi/np.arccos(np.clip(-g/2,-1+1e-15,1-1e-15))
        return(g<-2-3e-5)|(abs(g)<3e-5)|((g<-1+3e-5)&(g>-2+3e-5)&(abs(lab-np.rint(lab))<3e-5))
    # Partition each old-ideal coverage group into literal incompatibility
    # cliques, then compress every cross-clique bad biclique exactly.
    groups={}
    for i,r in enumerate(rows):groups.setdefault(tuple(r['coverage']),[]).append(i)
    cliques=[]
    for ids in groups.values():
        ok=fits(V[ids]@ETA@V[ids].T);local=[]
        for j,i in enumerate(ids):
            target=next((cl for cl in local if not ok[j,cl].any()),None)
            if target is None:local.append([j])
            else:target.append(j)
        cliques.extend([[ids[j]for j in cl]for cl in local])
    fixed={tuple(sorted(cl))for cl in cliques}
    for aa,bb in itertools.combinations(cliques,2):
        bad=~fits(V[aa]@ETA@V[bb].T);patterns={}
        for i,r in zip(aa,bad):
            if r.any():patterns.setdefault(np.packbits(r).tobytes(),[]).append(i)
        for packed,ii in patterns.items():
            js=np.flatnonzero(np.unpackbits(np.frombuffer(packed,dtype=np.uint8))[:len(bb)])
            fixed.add(tuple(sorted(ii+[bb[j]for j in js])))
    fixed=sorted(fixed);stats['complete_pair_tests']=len(rows)*(len(rows)-1)//2;stats['incompatibility_cliques']=len(cliques);stats['compressed_pair_rows']=len(fixed)
    print('COMPLETE_PAIR_CONSTRAINTS',stats['complete_pair_tests'],len(cliques),len(fixed),flush=True)
    for it in range(iterations):
        matrix=[{j:1 for j in row}for row in fixed];lo=[-np.inf]*len(matrix);hi=[1]*len(matrix)
        for p in samples:
            ids=np.flatnonzero(V@ETA@p>1e-7);matrix.append({int(i):1 for i in ids});lo.append(1);hi.append(np.inf)
        for i,j in sorted(conflicts):matrix.append({i:1,j:1});lo.append(-np.inf);hi.append(1)
        rr=[];cc=[];vv=[]
        for i,r in enumerate(matrix):
            for j,x in r.items():rr.append(i);cc.append(j);vv.append(x)
        A=csr_matrix((vv,(rr,cc)),shape=(len(matrix),len(rows)))
        ans=milp(c=np.ones(len(rows)),integrality=np.ones(len(rows)),bounds=Bounds(0,1),constraints=LinearConstraint(A,lo,hi),options=dict(time_limit=30,mip_rel_gap=0))
        rec=dict(iteration=it,samples=len(samples),conflicts=len(conflicts),status=int(ans.status),message=ans.message)
        if ans.x is None:history.append(rec);print(rec,flush=True);break
        chosen=np.flatnonzero(ans.x>.5).tolist();bad=[]
        for i,j in itertools.combinations(chosen,2):
            if not C.compatible(V[i]@ETA@V[j],3e-5):conflicts.add((i,j));bad.append([i,j])
        rec.update(selected=chosen,new_incompatible_pairs=bad)
        if bad:history.append(rec);print('selection',it,'newconflicts',len(bad),flush=True);continue
        current=initial
        try:
            for wall,i in enumerate(chosen,61):current=B.clip(current,V[i],wall)
        except ArithmeticError as e:rec['unresolved_clip']=str(e);history.append(rec);break
        counts=dict(Counter(r['kind']for r in current));rec['face_counts']=counts;history.append(rec);print('coherent',it,'selected',len(chosen),'counts',counts,flush=True)
        if not counts.get('hyperideal_rays',0)and not counts.get('unresolved_null',0):proposal=dict(selected=chosen,counts=counts);break
        added=0
        for r in current:
            if r['kind']!='hyperideal_rays':continue
            p=np.r_[1,r['v'][1:]/np.linalg.norm(r['v'][1:])]
            if all(np.linalg.norm(p-q)>1e-9 for q in samples):samples.append(p);added+=1
        if not added:break
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{B.TAG}.py',f'results/{C.TAG}.json',f'results/{C.TAG}_affine_gate.json',TRIAL,
           'results/fresh_nonarithmetic_h3_joint_cover_20261016.json']
    return C.clean(dict(status='BOUNDED_NUMERICAL_TWO_CONTACT_TRIAL_EXTENSION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                        counts=dict(stats),candidate_count=len(rows),pool=rows,history=history,compact_face_proposal=proposal,
                        scope='Two original55-wall finite contacts, one finite contact to a prior numerical trial side or continuous all-trial-dotted intervals, combined with the old4329 bank; all61 contact gates checked numerically. The exact25 null-line exclusions are removed from the new families. Complete finite-bank incompatibility constraints and full numerical face clipping do not constitute an exact or all-family nonexistence proof.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--iterations',type=int,default=50);ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run(a.iterations);p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved trial-extension replay')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidate_count'],flush=True)
