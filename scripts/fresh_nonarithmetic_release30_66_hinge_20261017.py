#!/usr/bin/env python3
"""Changed7/30 hinge using new66, with full native gates/deletion costs.

One bounded numerical event experiment; no numerical geometry promotion.
The old paired-B all-label gate is checked to survive its deleted49/62.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import sympy as sy
import fresh_nonarithmetic_release47_cycle_patch_20261017 as A

ROOT=A.ROOT;B=A.B;M=A.M
TAG='fresh_nonarithmetic_release30_66_hinge_20261017'
PARENT=ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json'
BOUNDARY=ROOT/'results/fresh_nonarithmetic_release49_62_boundary_20261017.json'
ALPHABET=[2,3,4,5]
FRAMES=[[7,66,k]for k in[33,47,52,55,59,63,64,44]]+[[7,47,52],[66,47,52]]

def setup():
    old,labels,ns,h,profiles=A.setup()
    primary=json.loads(PARENT.read_text());current=json.loads(BOUNDARY.read_text())['after63']
    gram=np.array([[B.mid(x)for x in row]for row in old['residual_gram']])
    vals,vecs=np.linalg.eigh((gram+gram.T)/2)
    transform=vecs@np.diag(1/np.sqrt(abs(vals)))
    root=np.array([B.mid(x)for x in primary['replacement_recipe']['residual_coordinates']])
    ns[66]=np.linalg.solve(transform,root)
    oldlabels=labels;labels=sorted(set(labels)-{49,62}|{66})
    ns={i:ns[i]for i in labels}
    for p in profiles:
        p['offsets']=[p['offsets'][oldlabels.index(i)]if i!=66 else 0. for i in labels]
    assert abs(ns[66]@M@ns[66]-2)<1e-7
    return current,labels,ns,h,profiles,primary

def surviving_pair_gate():
    path=ROOT/'results/fresh_nonarithmetic_replace47_b_pair_gate_20261017.json'
    x=json.loads(path.read_text());bad=[r for r in x['records']if'rejection'in r]
    assert len(bad)==16 and all(r['rejection']['wall']in[4,7]for r in bad)
    assert {tuple(p['indices'])for p in x['low_pairs']}==set(product(*x['accepted_groups']))
    assert all(p['rejection']['kind']in['obtuse','dotted_no_finite_pair','strict_consecutive_angle_gap']for p in x['low_pairs'])
    assert all(t['coefficient_square_positive']and t['linear_nonnegative']for t in x['high_label_tails'])
    return dict(source=str(path.relative_to(ROOT)),unchanged_old_walls=[4,7,30,47,52],
        low_rejected_branches=len(bad),mutual_rejected_pairs=len(x['low_pairs']),infinite_tails=len(x['high_label_tails']),
        scope='All old native branch blockers survive; the joins and tail norm equations use only unchanged walls. Adding66 can only further restrict this specified two-cap family, so the old gate persists.')

def deletion(primary,current,deleted=(30,)):
    vs=[tuple(v)for v in primary['retained_old_finite_supports']]+[tuple(v['support'])for v in primary['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    edges=defaultdict(list);faces=defaultdict(list)
    for i,v in enumerate(vs):
        for e in combinations(v,8):edges[e].append(i)
        for k in range(8):
            for w in combinations(v,k):faces[w].append(i)
    deleted=set(deleted)
    lost={i for i,v in enumerate(vs)if set(v)&deleted};touched=[]
    for w,ids in faces.items():
        if not lost.intersection(ids)or len(ids)<10-len(w):continue
        if all(all(len(edges[tuple(j for j in vs[i]if j!=k)])==2 for k in vs[i]if k not in w)for i in ids):
            touched.append(dict(support=list(w),dimension=9-len(w),finite_supports=[list(vs[i])for i in ids]))
    return dict(deleted_walls=sorted(deleted),lost_finite_supports=[list(vs[i])for i in sorted(lost)],lost_finite_count=len(lost),
        compact_faces_touched=touched,compact_faces_by_dimension=dict(Counter(r['dimension']for r in touched)),
        compact_S6_sides_touched=[r['wall']for r in current['facets']if r['compact']and any(set(v)&deleted for v in r['vertex_supports'])],
        scope='Exact loss of old certificates, reconstructed from complete199 finite incidence; replacement geometry may rebuild them but is not assumed to do so.')

def incidence_score(normals,labels,squares):
    """Numerical acute-face support check, flagging possible affine supports.

    Positive definite supports are actual by strict-center projection.
    With no parabolic2/3 support an acute H3 corner is finite or open;
    its finite vertices are simple and I=2E-3Vf. No approximate-zero
    determinant is accepted as a sign or as an ideal classification.
    """
    normals=np.array(normals);g=normals@M@normals.T
    error=float(max(abs(np.diag(g)-squares)))
    if error>1e-5:return dict(status='unresolved_normalization',max_norm_error=error)
    np.fill_diagonal(g,squares);eps=1e-8;edges=[];vertices=[];ambiguous=[]
    for a,b in combinations(range(len(labels)),2):
        d=g[a,a]*g[b,b]-g[a,b]**2
        if abs(d)<=eps*g[a,a]*g[b,b]:ambiguous.append([labels[a],labels[b]])
        elif d>0:edges.append([labels[a],labels[b]])
    for t in combinations(range(len(labels)),3):
        vals=np.linalg.eigvalsh(g[np.ix_(t,t)])
        if vals[0]>eps:vertices.append([labels[i]for i in t])
        elif abs(vals[0])<=eps:ambiguous.append([labels[i]for i in t])
    if ambiguous:return dict(status='unresolved_possible_affine_supports',max_norm_error=error,ambiguous_supports=ambiguous,
        definite_finite_supports=vertices,definite_geodesic_edges=edges)
    return dict(status='numerically_separated_acute_incidence',max_norm_error=error,finite=len(vertices),
        geodesic_edges=len(edges),open_ends=2*len(edges)-3*len(vertices),finite_supports=vertices,edge_supports=edges)

def literal_affine_obstruction(row,profiles):
    """An exact assigned-contact obstruction, independent of numeric roots."""
    assigned=dict(zip(row['frame'],row['labels']));pr=profiles[row['profile']]
    theta={2:sy.Integer(0),3:sy.Integer(1),4:sy.sqrt(2),5:(1+sy.sqrt(5))/2}
    trials=[]
    for support in row['acute_incidence_crosscheck']['ambiguous_supports']:
        assert set(support)<={7,55,66,67} and 7 in support and 67 in support
        assert set(support)-{67}<=set(assigned)
        ids=[5,6]+support;g=2*sy.eye(len(ids));edges=[]
        for a,b in combinations(range(len(ids)),2):
            i,j=ids[a],ids[b]
            if 67 in(i,j):
                k=j if i==67 else i
                m=(pr['label']if k==pr['node']else 2)if k in[5,6]else assigned[k]
            else:m=3 if{i,j}in[{5,6},{6,7}]else 2
            g[a,b]=g[b,a]=-theta[m]
            if m!=2:edges.append([i,j,m])
        for size in range(3,len(ids)+1):
            for ix in combinations(range(len(ids)),size):
                if ids.index(67)not in ix:continue
                h=g.extract(ix,ix)
                if sy.simplify(h.det())!=0:continue
                null=h.nullspace()
                if len(null)!=1:continue
                k=null[0]
                if k[0].is_negative:k=-k
                if not all(q.is_positive for q in k):continue
                pd=True
                for j in range(size):
                    sub=h.minor_submatrix(j,j)
                    if not all(sy.simplify(sub[:r,:r].det()).is_positive for r in range(1,size)):pd=False;break
                if not pd:continue
                vertices=[ids[i]for i in ix]
                return dict(affine_roots=vertices,gram=[[str(q)for q in r]for r in h.tolist()],
                    positive_kernel=[str(q)for q in k],finite_proper_deletions=True,
                    assigned_nonright_edges=[e for e in edges if set(e[:2])<=set(vertices)],
                    orthogonal_retained_H4=[0,1,2,3],
                    consequence='Every coherent exact realization of this assigned recipe has a persistent affine null direction incident to the retained noncrystallographic H4; it admits no finite-volume acute completion retaining these walls. No numerical root equality is used.')
    raise AssertionError(('unproved near-affine recipe',row['index'],row['frame'],row['labels'],row['profile']))

def run():
    z,labels,ns,h,profiles,primary=setup();actual=[i for i in z['actual_sides']if i!=30]
    retained=[i for i in labels if i!=30];old=np.array([ns[i]for i in retained]);poly=np.array([ns[i]for i in actual])
    baseline=B.ledger(np.array([ns[i]for i in z['actual_sides']]),z['actual_sides'])
    assert (baseline['finite'],baseline['exterior'],baseline['open_ends'])==(26,6,10),baseline
    removed=B.ledger(poly,actual);cost=deletion(primary,z);gate=surviving_pair_gate()
    print('CURRENT_CONTROL',baseline['finite'],baseline['exterior'],baseline['open_ends'],'DELETE30',removed['finite'],removed['exterior'],removed['open_ends'],'ambient_lost',cost['lost_finite_count'],flush=True)
    ms=np.array(list(product(ALPHABET,repeat=3)));desired=-np.array([[B.theta(m)for m in row]for row in ms])
    counts=Counter();scores=Counter();blockers=Counter();records=[];improvements=[];seen=set()
    joint_records=[];joint_seen=set();joint_scores=Counter()
    old_squares={r['wall']:B.mid(r['square'])for r in z['projections']}
    for ids in FRAMES:
        f=np.array([ns[i]for i in ids]);gram=f@M@f.T
        if abs(np.linalg.det(gram))<1e-10:counts['singular_frame']+=1;continue
        inv=np.linalg.inv(gram);k=np.linalg.svd(f@M)[2][-1];kq=k@M@k
        if abs(kq)<1e-10:counts['null_kernel']+=1;continue
        for pi,p in enumerate(profiles):
            off=np.array(p['offsets']);cc=desired-off[[labels.index(i)for i in ids]];pp=cc@inv@f
            rem=(p['square']-np.einsum('ij,jk,ik->i',pp,M,pp))/kq;good=np.where(rem>1e-9)[0]
            counts['real_norm_branches']+=2*len(good)
            for sign in[-1,1]:
                yy=pp[good]+sign*np.sqrt(rem[good])[:,None]*k
                products=yy@M@old.T+off[[labels.index(i)for i in retained]]
                for j,i in enumerate(ids):products[:,retained.index(i)]=desired[good,j]
                centers=yy@M@h+p['center_offset']
                for local in np.where(centers<-1e-8)[0]:
                    rel=[B.relation(x)for x in products[local]]
                    bad=tuple(i for i,r in zip(retained,rel)if r is None or r=='parallel_or_unresolved')
                    if bad:
                        blockers[bad]+=1
                        if bad in [(51,),(58,)]:
                            y=yy[local];key=bad+(pi,)+tuple(np.round(y,7))
                            if key not in joint_seen:
                                joint_seen.add(key)
                                ja=[i for i in actual if i not in bad]
                                jb=B.ledger(np.array([ns[i]for i in ja]+[y]),ja+[67])
                                ic=incidence_score([ns[i]for i in ja]+[y],ja+[67],[old_squares[i]for i in ja]+[p['square']])
                                joint_scores[bad[0],jb['finite'],jb['exterior'],jb['open_ends']]+=1
                                joint_records.append(dict(index=len(joint_records),additional_deleted=list(bad),
                                    frame=ids,labels=ms[good[local]].tolist(),profile=pi,branch=sign,projection=y.tolist(),
                                    center_product=float(centers[local]),
                                    removed_blocker_product=float(products[local,retained.index(bad[0])]),
                                    finite_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)and i not in bad],
                                    boundary=jb,acute_incidence_crosscheck=ic))
                        continue
                    y=yy[local];key=(pi,)+tuple(np.round(y,7))
                    if key in seen:counts['duplicate']+=1;continue
                    seen.add(key);counts['coherent_candidates']+=1
                    boundary=B.ledger(np.vstack([poly,y]),actual+[67])
                    ic=incidence_score([ns[i]for i in actual]+[y],actual+[67],[old_squares[i]for i in actual]+[p['square']])
                    score=(boundary['finite'],boundary['exterior'],boundary['open_ends']);scores[score]+=1
                    rec=dict(index=len(records),frame=ids,labels=ms[good[local]].tolist(),profile=pi,branch=sign,
                        projection=y.tolist(),center_product=float(centers[local]),finite_contacts=[[i,r]for i,r in zip(retained,rel)if isinstance(r,int)],
                        deleted30_product=float(y@M@ns[30]+off[labels.index(30)]),boundary=boundary,acute_incidence_crosscheck=ic)
                    records.append(rec)
                    if boundary['open_ends']<10:improvements.append(rec);print('IMPROVEMENT',rec['index'],ids,rec['labels'],score,flush=True)
        print('HINGE_FRAME',ids,'accepted',len(records),'scores',dict(scores),flush=True)
    joint_costs=[deletion(primary,z,[30,i])for i in[51,58]]
    assert not improvements and all(r['boundary']['open_ends']>=10 for r in joint_records)
    for r in records+joint_records:
        if r['acute_incidence_crosscheck']['status']=='unresolved_possible_affine_supports':
            r['literal_affine_obstruction']=literal_affine_obstruction(r,profiles)
    def crosscheck_summary(rows):
        hist=Counter();status=Counter();warnings=0;disagree=[];improved=[]
        for r in rows:
            c=r['acute_incidence_crosscheck'];status[c['status']]+=1
            warnings+=bool(r['boundary']['edges'].get('degenerate_face',0)or r['boundary']['rays']!=r['boundary']['finite']+r['boundary']['exterior'])
            if c['status']=='numerically_separated_acute_incidence':
                hist[c['finite'],c['open_ends']]+=1
                if c['open_ends']!=r['boundary']['open_ends']:disagree.append(r['index'])
                if c['open_ends']<10:improved.append(r['index'])
        return dict(status_counts=dict(status),numerical_projective_warning_cases=warnings,disagreement_indices=disagree,
            incidence_histogram=[dict(finite=f,open_ends=e,count=n)for(f,e),n in sorted(hist.items())],improvement_indices=improved,
            exact_affine_recipe_exclusions=sum('literal_affine_obstruction'in r for r in rows))
    paths=[Path(__file__),Path(A.__file__),Path(A.C.__file__),Path(B.__file__),PARENT,BOUNDARY,A.INPUT,A.C.PARENT,
        ROOT/'results/fresh_nonarithmetic_replace47_b_pair_gate_20261017.json']
    paths.append(ROOT/'results/fresh_h9_fixed_core_cusp_obstruction_20261013.md')
    return dict(status='BOUNDED_NUMERICAL_RELEASE30_NEW66_HINGE',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        current_control=baseline,after_deletion30=removed,deletion_cost=cost,surviving_pair_obstruction=gate,
        finite_anchor_labels=ALPHABET,frames=FRAMES,white_profiles=profiles,retained_labels=retained,
        counts=dict(counts),complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=n,count=c)for(f,e,n),c in sorted(scores.items())],
        blocker_histogram=[dict(walls=list(k),count=c)for k,c in sorted(blockers.items(),key=lambda z:(len(z[0]),z[0]))],
        records=records,improvements=improvements,joint_deletion_costs=joint_costs,joint_records=joint_records,
        joint_complete_boundary_histogram=[dict(additional_deleted=i,finite=f,exterior=e,open_ends=n,count=c)for(i,f,e,n),c in sorted(joint_scores.items())],
        numerical_boundary_crosschecks=dict(single=crosscheck_summary(records),joint=crosscheck_summary(joint_records)),
        scope='Exact old incidence/deletion audit, inherited all-label paired-B obstruction, and literal assigned-contact affine obstructions for56 near-affine recipes. New norm construction/native recognition/projective boundary scoring are bounded floating-point diagnostics:10 specified frames, labels2..5, all7 A2 profiles/both branches. Projective warnings and disagreements are retained; separated acute pair/triple incidence gives a distinct numerical score. Additional deletion51 or58 is tested only for sole native blockers. No high-label/continuous exclusion or exact promoted chamber; no finite-volume example.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['complete_boundary_histogram'],flush=True)
