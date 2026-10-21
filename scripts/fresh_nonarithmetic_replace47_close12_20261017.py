#!/usr/bin/env python3
"""Exact replacement47: right S6+12+52, pi/3 to7; compact new12 face."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_nonarithmetic_release54_replace12_20261017 as C

ROOT=C.ROOT
TAG='fresh_nonarithmetic_replace47_close12_20261017'
P,U,I,D=C.P,C.U,C.I,C.D
S=C.S
PARENT=ROOT/'results'/(C.TAG+'.json')
PRESCRIBED={**{i:2 for i in S+[12,52]},7:3}


def geometry():
    z=C.geometry();old=z['normals'];pair=z['pair'];h=z['common_center']
    projected7=[U.rnd(a+(b+2*c)/3)for a,b,c in zip(old[7],old[5],old[6])]
    frame=[projected7,old[12],old[52]]
    c=U.rnd(pair(old[7],old[52]))
    gram=[[I(F(4,3)),I(),c],[I(),I(2),-D.theta(6)],[c,-D.theta(6),I(2)]]
    inv,pivots=P.C53.inverse(gram)
    # Old8 supplies a well-conditioned kernel with entirely literal frame
    # products. Using old47 here would repeatedly expand dotted products
    # involving new12 and needlessly widen the exact interval enclosure.
    old_products=[-D.theta(5),I(),I()]
    old_coeff=U.imv(inv,old_products)
    projection=P.combine(old_coeff,frame)
    kernel=[U.rnd(a-b)for a,b in zip(old[8],projection)]
    kernel_square=U.rnd(2-U.idot(old_products,old_coeff))
    contacts=[I(-1),I(),I()];coeff=U.imv(inv,contacts)
    particular=P.combine(coeff,frame);residual=U.rnd(2-U.idot(contacts,coeff))
    assert kernel_square.l>0 and residual.l>0
    scale=U.sqrt(residual/kernel_square)
    candidates=[];branch_tests=[]
    for sign in [-1,1]:
        n=[U.rnd(a+sign*scale*b)for a,b in zip(particular,kernel)]
        tests={i:U.rnd(pair(n,v))for i,v in old.items()if i!=47 and i not in PRESCRIBED}
        bad=[i for i,v in tests.items()if v.l>0]
        if bad:
            branch_tests.append(dict(sign=sign,rejected_obtuse_wall=bad[0],product=tests[bad[0]].enc()))
            continue
        assert all(q.h<-2 for q in tests.values())
        assert pair(n,h).h<0
        candidates.append((sign,n));branch_tests.append(dict(sign=sign,accepted=True))
    assert len(candidates)==1
    sign,n=candidates[0];normals={**old,47:n}
    z.update(normals=normals,old_normals=old,replacement=n,
        replacement_recipe=dict(projected7=P.enc(projected7),frame_gram=[P.enc(r)for r in gram],
            inverse_pivots=pivots,kernel_seed=8,old_products=P.enc(old_products),contacts=P.enc(contacts),
            particular=P.enc(particular),kernel=P.enc(kernel),kernel_square=kernel_square.enc(),
            residual=residual.enc(),positive_scale=scale.enc(),chosen_sign=sign,branch_tests=branch_tests))
    return z


def literal(z,i,j):
    if i==j:return I(2)
    ns=z['normals']
    prescribed={47:PRESCRIBED,12:{**{k:2 for k in S+[7,8]},52:6},
        52:{**{k:2 for k in S+[8]},30:4},
        63:{**{k:2 for k in S+[30,44,65]},64:27},
        64:{**{k:2 for k in S+[30,52,65]},63:27},
        65:{**{k:2 for k in S+[8,63,64]},44:15,52:4}}
    m=prescribed.get(i,{}).get(j,prescribed.get(j,{}).get(i))
    if m is not None:return -D.theta(m)
    if i in[12,47]or j in[12,47]:return U.rnd(z['pair'](ns[i],ns[j]))
    labels=D.E.labels()
    if max(i,j)<55 and labels[i][j]:return -D.theta(labels[i][j])
    return U.rnd(z['pair'](ns[i],ns[j]))


def run():
    z=geometry();ns=z['normals'];pair=z['pair'];h=z['common_center'];r=z['replacement'];labels=sorted(ns)
    source=json.loads(PARENT.read_text());norm=U.rnd(pair(r,r));assert norm.l<=2<=norm.h
    gates=[]
    for i in labels:
        if i==47:continue
        q=U.rnd(pair(r,ns[i]));m=PRESCRIBED.get(i)
        if m is None:assert q.h<-2,(i,q.enc())
        else:
            delta=U.rnd(q+D.theta(m));assert delta.l<=0<=delta.h
        gates.append(dict(wall=i,label=m or'dotted',product=q.enc()))
    assert len(labels)==64 and len(gates)==63
    hh=U.rnd(pair(h,h));hp={i:U.rnd(pair(h,ns[i]))for i in labels}
    assert hh.h<0 and all(q.h<0 for q in hp.values())
    facets=[]
    for i in labels:
        point=[U.rnd(a-hp[i]*b/2)for a,b in zip(h,ns[i])]
        square=U.rnd(hh-hp[i]*hp[i]/2);assert square.h<0
        facets.append(dict(wall=i,point=P.enc(point),square=square.enc()))
    old_supports=source['retained_old_finite_supports']+[v['support']for v in source['new_finite_vertices']]
    assert len(old_supports)==197
    lost=[s for s in old_supports if 47 in s];retained=[s for s in old_supports if 47 not in s]
    assert lost==[S+[7,30,47]] and len(retained)==196
    old_compact_edges=[]
    for s in retained:
        common=sorted(set(s)&set(lost[0]))
        if len(common)==8:old_compact_edges.append(dict(support=common,other_vertex=s))
    assert old_compact_edges==[dict(support=S+[7,30],other_vertex=S+[7,30,55])]
    trials=[];vertices=[]
    for subset in combinations(sorted(PRESCRIBED),8):
        support=sorted([47]+list(subset))
        if 7 in support and 52 in support:
            q=literal(z,7,52);assert q.h<-2
            trials.append(dict(support=support,rejected_dotted_pair=[7,52],product=q.enc()));continue
        gram=[[literal(z,i,j)for j in support]for i in support]
        pivots=P.pd(gram);inv,ip=P.C53.inverse(gram);rhs=[hp[i]for i in support]
        coeff=U.imv(inv,rhs)
        point=[U.rnd(a-b)for a,b in zip(h,P.combine(coeff,[ns[i]for i in support]))]
        square=U.rnd(hh-U.idot(rhs,coeff));assert square.h<0
        strict=[]
        for i in labels:
            if i in support:continue
            q=U.rnd(pair(point,ns[i]));assert q.h<0,(support,i,q.enc());strict.append([i,q.enc()])
        vertices.append(dict(support=support,point=P.enc(point),square=square.enc(),strict_other_products=strict))
        trials.append(dict(support=support,PD_pivots=pivots,inverse_pivots=ip,vertex=len(vertices)-1))
    assert len(trials)==9 and len(vertices)==2
    assert {tuple(v['support'])for v in vertices}=={tuple(S+[7,12,47]),tuple(S+[12,47,52])}
    escape=C.decode(source['strict_escape']['point']);en=U.rnd(pair(escape,escape));ef=U.rnd(pair(escape,h))
    ep={i:U.rnd(pair(escape,ns[i]))for i in labels}
    assert en.l>0 and ef.h<0 and all(q.h<0 for q in ep.values())
    return dict(status='PASS_EXACT_64_WALL_REPLACEMENT47_CLOSE12',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(C.__file__),PARENT]},
        deleted_old_walls=[47],reused_label_for_new_wall=47,total_wall_count=64,
        computational_root_basis=z['basis_ids'],computational_gram=[P.enc(row)for row in z['exact_gram']],
        replacement_recipe=z['replacement_recipe'],new_normal=P.enc(r),norm_consistency=norm.enc(),new_contact_checks=gates,
        old47_product=pair(r,z['old_normals'][47]).enc(),common_center=P.enc(h),common_center_norm=hh.enc(),
        common_center_products={str(i):v.enc()for i,v in hp.items()},facet_witnesses=facets,
        old_finite_count=197,lost_old_finite_supports=lost,retained_old_finite_supports=retained,
        new_wall_finite_support_trials=trials,new_finite_vertices=vertices,complete_finite_vertex_count=198,
        old_compact_edges_at_lost_vertex=old_compact_edges,compact_faces_lost_dimension_at_least2=[],
        compact_face_preservation_proof='Old47 has exactly one finite vertex, with exactly one compact incident edge. Therefore no compact face of dimension>=2 meets old47. All other finite vertices persist strictly against new47 by acute projection. The old compact edge S6+7+30 becomes open; this is explicitly not monotone containment.',
        finite_census_proof='All unchanged old-only PD9 supports were already actual vertices of the parent acute chamber. Every finite vertex involving new47 uses8 of its9 finite neighbors; all9 possibilities are exhausted,7 rejected by dotted7/52 and2 actual.',
        strict_escape=dict(point=P.enc(escape),norm=en.enc(),future=ef.enc(),products={str(i):v.enc()for i,v in ep.items()}),
        scope='All64 actual facets, complete198 finitevertices and strict escape. Compact S6+12 quadrilateral and complete intrinsic boundary are certified separately; no complete ambient exterior cone or finite volume claim.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    path=ROOT/'results'/(TAG+'.json')
    if a.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_finite_vertex_count'],flush=True)

if __name__=='__main__':main()
