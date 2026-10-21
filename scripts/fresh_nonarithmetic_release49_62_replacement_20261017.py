#!/usr/bin/env python3
"""Exact joint49+62 replacement by one endpoint wall; scoped construction."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_nonarithmetic_release49_endpoint_gate_20261017 as G

E=G.E;OLD=E.X
ROOT,P,U,I,D=G.ROOT,G.P,G.U,G.I,G.D
TAG='fresh_nonarithmetic_release49_62_replacement_20261017'
S=E.S
PRESCRIBED={**{i:2 for i in S+[58,7]},33:4}
PARENT=ROOT/'results'/f'{OLD.TAG}.json'
GATE=ROOT/'results'/f'{G.TAG}.json'


def geometry():
    z=OLD.geometry();old=z['normals'];pair=z['pair']
    control=json.loads(G.CONTROL.read_text())['after64'];sp=G.Space(control)
    assert G.encmat([[OLD.literal(z,i,j)for j in E.FRAME]for i in E.FRAME])==control['residual_gram']
    proof,branches=sp.recipe([58,7,33],0,[2,2,4])
    chosen=next(b for b in branches if b['sign']==-1)
    frozen=json.loads(GATE.read_text())
    saved=next(a for a in frozen['records'][1]['exact_recipes']if a['profile']==0 and a['labels']==[2,2,4])
    assert proof['kernel_seed']==saved['kernel_seed']
    previous=next(b for b in saved['branches']if b['sign']==-1)
    assert all(G.overlap(a,b)for a,b in zip(G.decode(chosen['coordinates']),G.decode(previous['coordinates'])))
    root=P.combine(G.decode(chosen['coordinates']),[old[i]for i in E.FRAME])
    z.update(old_normals=old,normals={**{i:n for i,n in old.items()if i not in(49,62)},66:root},
             replacement=root,replacement_recipe=dict(**proof,chosen_sign=-1,residual_coordinates=chosen['coordinates']))
    return z


def literal(z,i,j):
    if i==j:return I(2)
    if 66 in (i,j):
        other=j if i==66 else i
        if other in PRESCRIBED:return -D.theta(PRESCRIBED[other])
        return U.rnd(z['pair'](z['normals'][i],z['normals'][j]))
    return OLD.literal(z,i,j)


def run():
    z=geometry();ns=z['normals'];pair=z['pair'];h=z['common_center'];n=z['replacement']
    old=json.loads(PARENT.read_text());hh=U.rnd(pair(h,h));assert hh.h<0
    norm=U.rnd(pair(n,n));assert norm.l<=2<=norm.h
    products={i:U.rnd(pair(h,r))for i,r in ns.items()};assert all(q.h<0 for q in products.values())
    gates=[]
    for i in sorted(ns):
        if i==66:continue
        q=U.rnd(pair(n,ns[i]));m=PRESCRIBED.get(i)
        if m is None:assert q.h<-2,(i,q.enc())
        else:assert G.overlap(q,-D.theta(m))
        gates.append(dict(wall=i,label=m or'dotted',product=q.enc()))
    assert len(ns)==63 and len(gates)==62
    facets=[]
    for i in sorted(ns):
        point=[U.rnd(a-products[i]*b/2)for a,b in zip(h,ns[i])]
        square=U.rnd(hh-products[i]*products[i]/2);assert square.h<0
        facets.append(dict(wall=i,point=P.enc(point),square=square.enc()))
    oldvertices=old['retained_old_finite_supports']+[r['support']for r in old['new_finite_vertices']]
    assert len(oldvertices)==198
    lost=[v for v in oldvertices if set(v)&{49,62}];retained=[v for v in oldvertices if not set(v)&{49,62}]
    assert lost==[S+[7,49,58]] and len(retained)==197
    checks=[]
    for support in retained:
        gram=[[OLD.literal(dict(normals=z['old_normals'],pair=pair),i,j)for j in support]for i in support]
        inv,_=P.C53.inverse(gram);rhs=[products[i]for i in support]
        coefficients=U.imv(inv,rhs)
        point=[U.rnd(a-b)for a,b in zip(h,P.combine(coefficients,[ns[i]for i in support]))]
        # Acute projection gives strict persistence algebraically as well.
        val=U.rnd(pair(point,n));assert val.h<0,(support,val.enc())
        checks.append(dict(support=support,new_product=val.enc()))
    trials=[];vertices=[]
    for subset in combinations(sorted(PRESCRIBED),8):
        support=sorted([66,*subset])
        if 58 in support and 33 in support:
            q=literal(z,58,33);assert q.h<-2
            trials.append(dict(support=support,dotted_pair=[33,58],product=q.enc()));continue
        gram=[[literal(z,i,j)for j in support]for i in support];pivots=P.pd(gram)
        inv,ip=P.C53.inverse(gram);rhs=[products[i]for i in support];coeff=U.imv(inv,rhs)
        point=[U.rnd(a-b)for a,b in zip(h,P.combine(coeff,[ns[i]for i in support]))]
        square=U.rnd(hh-U.idot(rhs,coeff));assert square.h<0
        other={i:U.rnd(pair(point,r))for i,r in ns.items()if i not in support}
        assert all(q.h<0 for q in other.values()),support
        vertices.append(dict(support=support,point=P.enc(point),square=square.enc(),
            other_products=[[i,q.enc()]for i,q in sorted(other.items())]))
        trials.append(dict(support=support,PD_pivots=pivots,inverse_pivots=ip,vertex=len(vertices)-1))
    assert len(trials)==9 and len(vertices)==2
    assert {tuple(v['support'])for v in vertices}=={tuple(S+[7,33,66]),tuple(S+[7,58,66])}
    inherited=G.decode(old['strict_escape']['point']);inherited_new=U.rnd(pair(inherited,n))
    assert inherited_new.l>0
    control=json.loads(G.CONTROL.read_text())['after64'];escape=None
    for row in control['rays']:
        if row['kind']!='exterior_ray' or 49 in row['active']:continue
        ray=P.combine(G.decode(row['coordinates']),[ns[i]for i in E.FRAME])
        if pair(ray,n).h>=0:continue
        for power in range(1,50):
            epsilon=G.Q(1,2**power);point=[U.rnd(a+epsilon*b)for a,b in zip(ray,h)]
            en=U.rnd(pair(point,point));ef=U.rnd(pair(point,h));ep={i:U.rnd(pair(point,r))for i,r in ns.items()}
            if en.l>0 and ef.h<0 and all(q.h<0 for q in ep.values()):
                escape=point;escape_seed=dict(parent_ray=row['index'],support=row['active'],epsilon=str(epsilon));break
        if escape is not None:break
    assert escape is not None,'need a strict renewed spacelike escape'
    cutter=json.loads((ROOT/'results/fresh_nonarithmetic_escape62_20261017_cap.json').read_text())
    oldescape=G.decode(cutter['old60_escape']['computational'])
    oldep={i:U.rnd(pair(oldescape,r))for i,r in ns.items()}
    oldnorm=U.rnd(pair(oldescape,oldescape));oldfuture=U.rnd(pair(oldescape,h))
    deleted={i:U.rnd(pair(n,z['old_normals'][i]))for i in(49,62)}
    assert all(q.l>0 for q in deleted.values())
    paths=[Path(__file__),Path(OLD.__file__),Path(G.__file__),Path(E.__file__),PARENT,GATE,G.CONTROL,
           ROOT/'results/fresh_nonarithmetic_release49_hexagon_20261017.json',
           ROOT/'results/fresh_nonarithmetic_escape62_20261017_cap.json']
    return dict(status='PASS_EXACT_JOINT49_62_REPLACEMENT_63_WALLS_199_FINITE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_old_walls=[49,62],new_wall=66,total_wall_count=63,
        computational_root_basis=z['basis_ids'],computational_gram=G.encmat(z['exact_gram']),
        residual_basis=E.FRAME,replacement_recipe=z['replacement_recipe'],new_normal=P.enc(n),new_norm=norm.enc(),
        native_gates=gates,deleted_wall_products={str(i):q.enc()for i,q in deleted.items()},
        common_center=P.enc(h),center_norm=hh.enc(),center_products={str(i):q.enc()for i,q in products.items()},
        facet_witnesses=facets,old_finite_count=198,lost_old_finite_supports=lost,
        retained_old_finite_supports=retained,retention_checks=checks,
        new_finite_support_trials=trials,new_finite_vertices=vertices,complete_finite_vertex_count=199,
        compact_faces_lost_dimension_at_least2=[],
        deletion_scope='The frozen release49 incidence audit proves49 touches no compact face of dimension>=2. No old finite support contains62, so62 touches no compact face. Exactly one finite vertex S6+7+49+58 is lost; its compact incident edges are removed, with the complete intrinsic ledger in the successor boundary package.',
        finite_census_proof='Every remaining old-only PD9 support is an old actual finite vertex. Exactly197 avoid deleted49/62 and are checked strictly against66. Every new finite vertex uses8 of66s9 finite neighbors; all9 supports are exhausted,7 contain dotted33/58 and2 are finite actual corners. Thus199 is complete.',
        inherited64_escape_new_product=inherited_new.enc(),
        strict_escape=dict(**escape_seed,point=P.enc(escape),norm=en.enc(),future=ef.enc(),products={str(i):q.enc()for i,q in ep.items()}),
        old60_escape_after_replacement=dict(point=P.enc(oldescape),norm=oldnorm.enc(),future=oldfuture.enc(),
            products={str(i):q.enc()for i,q in oldep.items()},strictly_restored=all(q.h<0 for q in oldep.values())),
        scope='Exact changed63-wall acute Coxeter chamber, complete199 finitevertices, all actual facets and a strict spacelike escape. The complete S6 boundary is a separate successor. Both earlier strict escape witnesses are tested explicitly; no blanket containment assertion follows from deleting49/62. No full ambient exterior cone or finite-volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_finite_vertex_count'],flush=True)
