#!/usr/bin/env python3
"""Exact replacement12 at angle pi/6 to52, after releasing vertex-free54."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_nonarithmetic_s6_release34_exact_20261017 as X

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_release54_replace12_20261017'
P,U,I,D=X.P,X.U,X.I,X.D
S=X.S
PARENT=ROOT/'results'/(X.TAG+'.json')


def decode(row):return [I(F(a),F(b))for a,b in row]


def geometry():
    z=X.geometry();old=z['normals'];pair=z['pair']
    # Exact A2 inverse: old7 has its sole white contact -1 at white6.
    projected7=[U.rnd(a+(b+2*c)/3)for a,b,c in zip(old[7],old[5],old[6])]
    frame=[projected7,old[8],old[52]]
    c=U.rnd(pair(old[7],old[52]))
    gram=[[I(F(4,3)),-D.theta(5),c],[-D.theta(5),I(2),I()],[c,I(),I(2)]]
    inv,pivots=P.C53.inverse(gram)
    old_products=[I(),I(),U.rnd(pair(old[12],old[52]))]
    old_coeff=U.imv(inv,old_products)
    projection=P.combine(old_coeff,frame)
    kernel=[U.rnd(a-b)for a,b in zip(old[12],projection)]
    kernel_square=U.rnd(2-U.idot(old_products,old_coeff))
    contacts=[I(),I(),-D.theta(6)]
    coeff=U.imv(inv,contacts)
    particular=P.combine(coeff,frame)
    residual=U.rnd(2-U.idot(contacts,coeff))
    assert kernel_square.l>0 and residual.l>0
    scale=U.sqrt(residual/kernel_square)
    normal=[U.rnd(a+scale*b)for a,b in zip(particular,kernel)]
    normals={i:n for i,n in old.items()if i not in[12,54]};normals[12]=normal
    z.update(normals=normals,old_normals=old,
        replacement_recipe=dict(projected7=P.enc(projected7),frame_gram=[P.enc(r)for r in gram],
            inverse_pivots=pivots,old_products=P.enc(old_products),contacts=P.enc(contacts),
            particular=P.enc(particular),kernel=P.enc(kernel),kernel_square=kernel_square.enc(),
            residual=residual.enc(),positive_scale=scale.enc()),replacement=normal)
    return z


def run():
    z=geometry();ns=z['normals'];pair=z['pair'];h=z['common_center'];r=z['replacement']
    source=json.loads(PARENT.read_text());labels=sorted(ns)
    norm=U.rnd(pair(r,r));assert norm.l<=2<=norm.h
    prescribed={**{i:2 for i in S+[7,8]},52:6}
    gates=[]
    for i in labels:
        if i==12:continue
        value=U.rnd(pair(r,ns[i]))
        m=prescribed.get(i)
        if m is None:assert value.h<-2,(i,value.enc())
        else:
            # These equalities follow from the projection/kernel recipe.
            delta=U.rnd(value+D.theta(m));assert delta.l<=0<=delta.h
        gates.append(dict(wall=i,label=m or 'dotted',product=value.enc()))
    assert len(labels)==64 and len(gates)==63
    hh=U.rnd(pair(h,h));hp={i:U.rnd(pair(h,ns[i]))for i in labels}
    assert hh.h<0 and all(q.h<0 for q in hp.values())
    facets=[]
    for i in labels:
        point=[U.rnd(a-hp[i]*b/2)for a,b in zip(h,ns[i])]
        square=U.rnd(hh-hp[i]*hp[i]/2);assert square.h<0
        facets.append(dict(wall=i,point=P.enc(point),square=square.enc()))
    old_supports=[row['support']for key in ['preserved_old_finite_vertices','finite_corners']for row in source[key]]
    assert len(old_supports)==196 and len({tuple(s)for s in old_supports})==196
    lost=[s for s in old_supports if 12 in s or 54 in s]
    retained=[s for s in old_supports if s not in lost]
    assert lost==[S+[7,8,12]] and all(54 not in s for s in old_supports)
    # A compact two-face incident to lost v would have two compact edges.
    # All edges except the one omitting12 retain12 in their active support;
    # the complete old census has no second vertex containing12.
    old_compact_edges_at_lost=[]
    for s in retained:
        common=sorted(set(s)&set(lost[0]))
        if len(common)==8:old_compact_edges_at_lost.append(dict(support=common,other_vertex=s))
    assert len(old_compact_edges_at_lost)==1 and old_compact_edges_at_lost[0]['support']==S+[7,8]

    old_labels=D.E.labels()
    old_new={52:{**{i:2 for i in S+[8]},30:4},
             63:{**{i:2 for i in S+[30,44,65]},64:27},
             64:{**{i:2 for i in S+[30,52,65]},63:27},
             65:{**{i:2 for i in S+[8,63,64]},44:15,52:4}}
    def literal(i,j):
        if i==j:return I(2)
        if 12 in[i,j]:
            other=j if i==12 else i
            if other in prescribed:return -D.theta(prescribed[other])
        m=old_new.get(i,{}).get(j,old_new.get(j,{}).get(i))
        if m is not None:return -D.theta(m)
        if max(i,j)<55 and old_labels[i][j]:return -D.theta(old_labels[i][j])
        return U.rnd(pair(ns[i],ns[j]))
    trials=[];vertices=[]
    for subset in combinations(sorted(prescribed),8):
        support=sorted([12]+list(subset))
        if 7 in support and 52 in support:
            witness=literal(7,52);assert witness.h<-2
            trials.append(dict(support=support,rejected_dotted_pair=[7,52],product=witness.enc()));continue
        gram=[[literal(i,j)for j in support]for i in support]
        pivots=P.pd(gram);inv,ip=P.C53.inverse(gram)
        rhs=[hp[i]for i in support]
        coeff=U.imv(inv,rhs)
        point=[U.rnd(a-b)for a,b in zip(h,P.combine(coeff,[ns[i]for i in support]))]
        square=U.rnd(hh-U.idot(rhs,coeff));assert square.h<0
        strict=[]
        for i in labels:
            if i in support:continue
            v=U.rnd(pair(point,ns[i]));assert v.h<0,(support,i,v.enc());strict.append([i,v.enc()])
        vertices.append(dict(support=support,point=P.enc(point),square=square.enc(),strict_other_products=strict))
        trials.append(dict(support=support,PD_pivots=pivots,inverse_pivots=ip,vertex=len(vertices)-1))
    assert len(trials)==9 and len(vertices)==2
    assert {tuple(v['support'])for v in vertices}=={tuple(S+[7,8,12]),tuple(S+[8,12,52])}
    # Every retained old vertex is strict for the new root by the acute
    # projection lemma, since its active n-root set is unchanged and PD.
    assert len(retained)==195
    escape=decode(source['strict_escape']['computational_point'])
    en=U.rnd(pair(escape,escape));ef=U.rnd(pair(escape,h));ep={i:U.rnd(pair(escape,ns[i]))for i in labels}
    assert en.l>0 and ef.h<0 and all(q.h<0 for q in ep.values())
    return dict(status='PASS_EXACT_64_WALL_RELEASE54_REPLACEMENT12',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in[Path(__file__),Path(X.__file__),PARENT]},
        deleted_old_walls=[12,54],reused_label_for_new_wall=12,total_wall_count=64,
        computational_root_basis=z['basis_ids'],computational_gram=[P.enc(row)for row in z['exact_gram']],
        replacement_recipe=z['replacement_recipe'],new_normal=P.enc(r),norm_consistency=norm.enc(),new_contact_checks=gates,
        old54_obtuse_product=pair(r,z['old_normals'][54]).enc(),old12_product=pair(r,z['old_normals'][12]).enc(),
        common_center=P.enc(h),common_center_norm=hh.enc(),common_center_products={str(i):v.enc()for i,v in hp.items()},
        facet_witnesses=facets,old_finite_count=196,lost_old_finite_supports=lost,
        retained_old_finite_supports=retained,new_wall_finite_support_trials=trials,new_finite_vertices=vertices,
        complete_finite_vertex_count=len(retained)+len(vertices),old_compact_edges_at_lost_vertex=old_compact_edges_at_lost,
        compact_faces_lost_dimension_at_least2=[],
        compact_face_preservation_proof='Old54 has no finite vertex. Old12 has only one finite vertex v, and v has only one compact incident edge. Hence no old compact face of dimension>=2 meets either deleted wall; all such faces remain, since their vertices are retained and strict against new12. The old compact edge S6+7+8 acquires a new endpoint.',
        finite_census_proof='Every retained old PD n-root support was already an actual old vertex by acute projection; deleting walls cannot reveal another old-only finite vertex. New12 has exactly9 finite neighbors, and all9 eight-subsets are exhausted, with7 dotted-pair rejections and2 actual vertices.',
        strict_escape=dict(point=P.enc(escape),norm=en.enc(),future=ef.enc(),products={str(i):v.enc()for i,v in ep.items()}),
        scope='Complete finite-vertex census and all64 facets; infinite volume certified. No complete ambient exterior cone is asserted.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_finite_vertex_count'],flush=True)


if __name__=='__main__':main()
