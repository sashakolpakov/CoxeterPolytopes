#!/usr/bin/env python3
"""Release the obsolete prism blockers; rebuild from the actual open boundary.

Exact finite-incidence deletion debit; bounded numerical replacement/closure.
The exact old63 and anchor62 chambers are controls, not promoted numerics.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import sympy as sp
import fresh_release58_66_seven_boundary_20261017 as N
ROOT,M,G=N.ROOT,N.M,N.G
TAG='fresh_release55_59_61_prism_replacement_20261017'
ANCHOR='fresh_release58_66_anchor_completion_gate_20261017'
DELETED=[55,58,59,61,66]
PAIR_POOL=24


def compact_incidence(vs):
    edges=defaultdict(list);faces=defaultdict(list)
    for i,v in enumerate(vs):
        for e in combinations(v,8):edges[e].append(i)
        for k in range(8):
            for w in combinations(v,k):faces[w].append(i)
    compact=[]
    for w,ids in faces.items():
        if len(ids)<10-len(w):continue
        if all(all(len(edges[tuple(j for j in vs[i]if j!=k)])==2 for k in vs[i]if k not in w)for i in ids):compact.append((w,ids))
    return compact


def debit(vs,compact,deleted):
    lost={i for i,v in enumerate(vs)if set(v)&set(deleted)}
    touched=[dict(support=w,dimension=9-len(w),finite_supports=[vs[i]for i in ids])for w,ids in compact if lost.intersection(ids)]
    return dict(deleted=deleted,lost_finite_supports=[vs[i]for i in sorted(lost)],lost_finite_count=len(lost),compact_faces_touched=touched,
        compact_faces_by_dimension=dict(Counter(r['dimension']for r in touched)),lost_supports_containing8=sum(8 in vs[i]for i in lost))


def annotate_unresolved(s,start,nd,types,b):
    """Keep tangent equations unresolved, with their literal B context.

    The projected vector below is only a floating-point diagnostic. The
    separate exact_tangencies() certificate, not this calculation, resolves
    the five saved literal equations.
    """
    for row in s.unresolved[start:]:
        row['B_context']=dict(index=b['index'],profile=b['profile'],labels=b['labels'],branch=b['branch'])
        if row.get('norm',{}).get('status')!='unresolved_tangent_norm':continue
        frame=row['frame'];pi=row['profile'];ms=row['labels'];F=np.array([nd[i]for i in frame])
        rhs=np.array([-N.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,types[i]))for i,m in zip(frame,ms)])
        v=np.linalg.solve(F@M@F.T,rhs)@F
        vals,bad,cp=s.native(v,pi,nd,types)
        row['projected_vector_diagnostic']=dict(normal=v.tolist(),norm_residual=float(v@M@v)-s.profiles[pi]['square'],
          blockers=[[i,vals[i]]for i in bad],center_product=cp,
          scope='Numerical projected-vector diagnostic only; the separate symbolic recipe certificate decides the exact tangent equation.')


def exact_tangencies(unresolved,anchor):
    """Symbolic frame certificates, without saved floating-point normals."""
    assert anchor['finite_frame_profiles']==[[30,0],[67,0],[60,6]]
    assert anchor['old60_profile_provenance']['prescribed_S6_labels']=={'0':2,'1':2,'2':2,'3':2,'5':2,'6':5}
    assert anchor['anchor_recipe']['profile']==0 and anchor['anchor_recipe']['labels']==[3,2,2]
    phi=(1+sp.sqrt(5))/2;q60=sp.simplify(2-2*phi**2/3)
    assert q60==(3-sp.sqrt(5))/3 and q60>0
    theta=lambda m:sp.simplify(2*sp.cos(sp.pi/m))
    out=[]
    for r in unresolved:
        assert r['norm']['status']=='unresolved_tangent_norm'
        assert r['profile']==r['B_context']['profile']==0
        bm=r['B_context']['labels'];frame=r['frame']
        if frame==[60,67,68]:
            gram=sp.Matrix([[q60,0,-theta(bm[2])],[0,2,-theta(bm[1])],[-theta(bm[2]),-theta(bm[1]),2]])
        else:
            assert frame==[30,67,68]
            gram=sp.Matrix([[2,-1,-theta(bm[0])],[-1,2,-theta(bm[1])],[-theta(bm[0]),-theta(bm[1]),2]])
        rhs=sp.Matrix([-theta(m)for m in r['labels']])
        minors=[sp.simplify(gram[:k,:k].det())for k in[1,2,3]]
        assert all(x>0 for x in minors)
        coeff=(gram.inv()*rhs).applyfunc(sp.simplify);square=sp.simplify((coeff.T*rhs)[0])
        assert square==2 and all(x<=0 for x in coeff)and any(x<0 for x in coeff)
        # Only norm2, zero-profile frame roots occur with nonzero coefficient.
        assert all(i!=60 or x==0 for i,x in zip(frame,coeff))
        out.append(dict(B_context=r['B_context'],frame=frame,labels=r['labels'],
          gram=[[str(x)for x in row]for row in gram.tolist()],rhs=list(map(str,rhs)),
          positive_principal_minors=list(map(str,minors)),projection_coefficients=list(map(str,coeff)),
          projection_square=str(square),complement_signature='one_negative',norm_residual='0',
          conclusion='The unique norm2 root is the displayed nonpositive combination of actual zero-profile frame walls, so its product with any strict chamber center is positive.'))
    assert len(out)==5
    return out


def run():
    parent=json.loads((ROOT/'results'/f'{N.TAG}.json').read_text());anchor_exact=json.loads((ROOT/'results'/f'{ANCHOR}.json').read_text())
    for d in[parent,anchor_exact]:
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=N.Scout();oldvs=[tuple(v)for v in s.primary['retained_old_finite_supports']]+[tuple(v['support'])for v in s.primary['new_finite_vertices']]
    avs=[v for v in oldvs if not set(v)&{58,66}]+[tuple(v['support'])for v in anchor_exact['new_finite_vertices']]
    assert len(oldvs)==199 and len(avs)==192
    oc=compact_incidence(oldvs);ac=compact_incidence(avs);comparisons=[]
    for extra in[[55],[61],[55,61],[55,61,7],[55,61,59]]:
        comparisons.append(dict(extra=extra,relative_control=debit(oldvs,oc,[58,66,*extra]),relative_anchor=debit(avs,ac,extra)))
    chosen=comparisons[-1];assert chosen['relative_control']['lost_finite_count']==15 and chosen['relative_anchor']['lost_finite_count']==8
    assert chosen['relative_control']['compact_faces_by_dimension']=={2:10,3:2}and chosen['relative_control']['lost_supports_containing8']==0
    ar=parent['anchor_domain']['caps'][parent['chosen_anchor_index']];assert ar['profile']==0 and ar['labels']==[3,2,2]
    base={i:v for i,v in s.old.items()if i not in DELETED};base[67]=np.array(ar['normal']);ids=[i for i in s.cur['actual_sides']if i not in DELETED]+[67]
    types={67:0}
    # Recover the complete old-side A2 profile from two independent white
    # pairings. Native evaluation itself always uses the frozen full offsets.
    for i in ids:
        if i==67:continue
        ix=s.index[i];target=[s.profiles[p]['offsets'][ix]for p in[1,4]]
        matches=[p for p in range(7)if max(abs(s.off(q,p)-v)for q,v in zip([1,4],target))<1e-7]
        assert len(matches)==1,(i,target,matches);types[i]=matches[0]
        assert abs(float(base[i]@M@base[i])-s.profiles[types[i]]['square'])<1e-5
    vals,bad,h=s.native(base[67],0,{i:v for i,v in base.items()if i!=67},types);assert not bad and h<-1e-8
    before=s.score(base,types,ids,'released_region');assert before is not None
    assert s.ledgers[before]['counts']==dict(finite=19,exterior=7)and s.ledgers[before]['open_ends']==11
    seeds=s.bank([30,67,60],[(0,1),(1,2)],base,types,ids,68,'replacement_B')
    assert json.loads(json.dumps(seeds['symbolic_recipes']))==anchor_exact['symbolic_recipes']
    print('REPLACEMENT_B','recipes',len(seeds['symbolic_recipes']),'accepted',len(seeds['caps']),'unresolved',len(s.unresolved),flush=True)
    rows=[];pairhist=Counter();escapes=[]
    for b in seeds['caps']:
        nb,tb=s.add(base,types,[b],68);bids=ids+[68];ledger=s.ledgers[b['score']]
        b['numerical_deleted_wall_matches']=[dict(wall=i,max_coordinate_difference=float(np.max(abs(nb[68]-s.old[i]))))
          for i in DELETED if np.max(abs(nb[68]-s.old[i]))<1e-7 and
          max(abs(s.off(q,b['profile'])-s.profiles[q]['offsets'][s.index[i]])for q in[1,4])<1e-7]
        opens=[tuple(e['sides'])for e in ledger['edges']if e['kind']in['one_open','two_open']]
        frames={}
        for e,f in combinations(opens,2):
            frame=tuple(sorted(set(e)|set(f)))
            if len(frame)!=3 or 68 not in frame:continue
            duties=tuple(sorted({tuple(sorted(frame.index(i)for i in x))for x in[e,f]}))
            frames[frame]=tuple(sorted(set(frames.get(frame,()))|set(duties)))
        banks=[];caps=[];seen={}
        for frame,duties in sorted(frames.items()):
            unresolved_start=len(s.unresolved)
            bank=s.bank(list(frame),list(duties),nb,tb,bids,69,'boundary_C')
            annotate_unresolved(s,unresolved_start,nb,tb,b)
            banks.append(bank)
            for c in bank['caps']:
                key=(c['profile'],*tuple(np.round(c['normal'],8)))
                if key in seen:
                    caps[seen[key]]['additional_recipes'].append(dict(frame=c['frame'],labels=c['labels'],branch=c['branch']));continue
                seen[key]=len(caps);caps.append(dict(c,index=len(caps),additional_recipes=[]))
        pool=sorted(caps,key=lambda c:(s.ledgers[c['score']]['open_ends'],-s.ledgers[c['score']]['counts']['finite'],max(c['labels']),c['index']))[:PAIR_POOL]
        pairs=[]
        for c,d in combinations(pool,2):
            nn,tt=s.add(nb,tb,[c,d],69);nv,bad,cp=s.native(nn[70],tt[70],{i:v for i,v in nn.items()if i!=70},tt)
            relation=G.relation(nv[69]);pairhist[str(relation)]+=1
            if bad or cp>=-1e-8:continue
            score=s.score(nn,tt,bids+[69,70],'boundary_C_D')
            if score is None:continue
            rr=dict(indices=[c['index'],d['index']],native_product=nv[69],relation=relation,score=score);pairs.append(rr)
            if s.ledgers[score]['open_ends']<10:s.improvements.append(dict(stage='boundary_C_D',B=b,roots=[c,d],record=rr))
        # A complete target cone with exterior rays is still unclosed. Save
        # their literal supports explicitly rather than infer volume from
        # isolated new compact faces.
        local_scores=sorted({b['score'],*[c['score']for c in caps],*[x['score']for x in pairs]})
        for score in local_scores:
            supports=[r['active']for r in s.ledgers[score]['rays']if r['kind']=='exterior']
            assert {(7,30,47),(30,47,52)}<={tuple(sorted(x))for x in supports}
            escapes.append(dict(B=b['index'],score=score,exterior_supports=supports))
        rows.append(dict(B=b,open_edges=opens,derived_frames=[dict(frame=k,duties=v)for k,v in sorted(frames.items())],banks=banks,continuations=caps,pair_pool=[c['index']for c in pool],coherent_pairs=pairs))
        print('BOUNDARY_B',b['index'],'frames',len(frames),'C',len(caps),'pairs',len(pairs),'scores',[(s.ledgers[j]['counts'],s.ledgers[j]['open_ends'])for j in local_scores],'unresolved',len(s.unresolved),flush=True)
    tangencies=exact_tangencies(s.unresolved,anchor_exact)
    paths=[Path(__file__),Path(N.__file__),ROOT/'results'/f'{N.TAG}.json',ROOT/'scripts'/f'{ANCHOR}.py',ROOT/'results'/f'{ANCHOR}.json',Path(N.H.__file__),Path(N.E.__file__),Path(N.N.__file__),Path(N.S.__file__),Path(N.P.__file__),Path(N.G.__file__),N.H.PARENT,N.H.BOUNDARY]
    return dict(status='BOUNDED_RELEASED_PRISM_BLOCKER_REPLACEMENT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      deleted_old_walls=DELETED,retained_old_wall_count=58,initial_anchor_wall_count=59,exact_deletion_comparisons=comparisons,chosen_debit=chosen,
      old_side_profiles=sorted(types.items()),control_score=s.control,released_score=before,anchor=ar,replacement_seeds=seeds,
      high_label_samples=N.HIGH,continuation_pair_pool_limit=PAIR_POOL,records=rows,pair_relation_histogram=dict(pairhist),exterior_ledgers=escapes,
      persistent_numerical_exterior_supports=[[7,30,47],[30,47,52]],
      counts=dict(s.counts),boundary_ledgers=s.ledgers,improvements=s.improvements,unresolved=s.unresolved,
      exact_tangent_gates=tangencies,remaining_unresolved=[],
      histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
      blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],near_native_failures=s.near,
      scope='Exact incidence debit from the complete199/192 finite censuses; no wall8 finite vertex is deleted. Actual blockers55/61 plus cheap59 are released along with prior58/66, yielding58 retained old walls plus anchor67. The frozen nine recipe roots are reevaluated in this genuinely changed chamber, then new continuation frames are derived from actual complete open-edge ledgers. All seven native A2 profiles, complete low2..5 recipes, isolated high samples6/8/12/24, both norm signs; at most24 ranked continuation roots per B participate in paired tests. All retained native products/commoncenter and complete numericalS6cone/PD support agreement are required. Five numerical tangent flags are excluded by separate exact symbolic recipe/common-center certificates. No exact new cap, ambient finite-vertex census, complete label-domain exclusion, or finite-volume claim; exact controls unchanged.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
