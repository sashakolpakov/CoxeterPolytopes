#!/usr/bin/env python3
"""Exact payment separation on complete selected-endpoint double-H4 triangles."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_double_h4_triangle_separation_20261017'
TRI='fresh_sqrt5_rank13_double_h4_endpoints_triangles_20261017'
LP='fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017'
K,Z=A.K,A.Z

def literal(target,occurrence):
    g={i:{} for i in range(14)}
    for i,j,m in target['edges']:g[i][j]=g[j][i]=m
    port=occurrence['port'];other=occurrence['other']
    attachment=target['ports'][port]['attachments'][occurrence['attachment']]
    for i,m in attachment['contacts']:g[i][13]=g[13][i]=m
    if occurrence['label']!=2:g[port][13]=g[13][port]=occurrence['label']
    return g,[port,other,13]

def geometry(g,black,corners):
    assert all(m in (3,5) for row in g.values() for m in row.values())
    phi=K(F(1,2),F(1,2));mag={3:K(1),5:phi}
    gram=[[K(2) if i==j else -mag[g[i][j]] if j in g[i] else Z for j in range(14)]for i in range(14)]
    assert A.inertia(gram)==[13,1,0]
    assert A.inertia([[v.conj()for v in row]for row in gram])==[14,0,0]
    assert all(x.integral()for row in gram for x in row)
    inverse=A.inverse(gram);white=sorted(set(range(14))-set(black));rays=[]
    for removed in black:
        active=[i for i in range(14)if i!=removed]
        assert A.inertia(A.sub(gram,active))==[13,0,0]
        assert A.inertia([[x.conj()for x in row]for row in A.sub(gram,active)])==[13,0,0]
        v=[-inverse[i][removed]for i in range(14)];products=A.mv(gram,v)
        assert products==[K(-1)if i==removed else Z for i in range(14)]
        assert A.pair(gram,v,v).sign()<0
        rays.append(v)
    assert all(A.pair(gram,v,w).sign()<0 for v in rays for w in rays)
    h=[sum((v[i]for v in rays),Z)for i in range(14)]
    push=A.mv(inverse,[K(-1)if i in white else Z for i in range(14)])
    epsilon=F(1)
    while True:
        center=[h[i]+epsilon*push[i]for i in range(14)]
        if A.pair(gram,center,center).sign()<0 and all(A.pair(gram,center,v).sign()<0 for v in rays):break
        epsilon/=2
    assert all(x.sign()<0 for x in A.mv(gram,center))
    face_tests=[]
    for removed,ray in zip(black,rays):
        active=[i for i in range(14)if i!=removed]
        kernel=A.kernel([gram[i]for i in active]);assert len(kernel)==1
        products=A.mv(gram,kernel[0]);assert products[removed]
        normalized=[-v/products[removed]for v in kernel[0]]
        assert normalized==ray
        face_tests.append(dict(removed=removed,active=active,ray=[v.enc()for v in ray],norm=A.pair(gram,ray,ray).enc()))
    return dict(gram=A.enc(gram),physical_inertia=[13,1,0],conjugate_inertia=[14,0,0],
        integral_cartan_entries=196,black=black,white=white,
        edges=[[i,j,m]for i,row in g.items()for j,m in row.items()if i<j],
        center_epsilon=str(epsilon),strict_center=[v.enc()for v in center],
        center_norm=A.pair(gram,center,center).enc(),complete_finite_face_rays=face_tests,
        corners=corners,total_payment=str(sum((F(c['payment'])for c in corners),F(0))),
        cone_proof='The full Gram is invertible. On the common white-wall subspace, every vector satisfying the three black inequalities is uniquely a nonnegative combination of the three displayed dual rays. All three rays are future timelike and mutually future-oriented, so their projectivized cone is exactly a compact triangle. Every normal occurs at a finite corner with13 independent active normals and a strict remaining inequality, hence every wall is an actual local facet.')

def run():
    triangles=B.read(TRI);lp=B.read(LP)
    primal=next(c for c in lp['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    weights={i:F(v)for i,v in primal['weights']}
    lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']for v in lp['variables']}
    finite={B.C.code(g):name for name,g in B.C.models()if len(g)<=13}
    records=[];counts=Counter();witness=None
    for ti,target in enumerate(triangles['targets']):
        for oi,occurrence in enumerate(target['triangle_occurrences']):
            g,black=literal(target,occurrence);corners=[]
            for removed in [13,occurrence['port'],occurrence['other']]:
                pair=[i for i in black if i!=removed]
                end={i:{j:m for j,m in row.items()if j!=removed}for i,row in g.items()if i!=removed}
                pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)]for p in pieces))
                code=sorted(B.C.code(p,pair)for p in pieces);index=lookup.get((forest,B.dump(code)))
                corners.append(dict(removed=removed,black=pair,forest=list(forest),marked_forest_code=code,
                    variable=index,payment=None if index is None else str(weights.get(index,F(0)))))
            assert corners[0]['forest']==target['forest']
            assert corners[1]['forest']==occurrence['destination']
            assert corners[2]['forest']==occurrence['third_forest']
            assert B.dump(corners[2]['marked_forest_code'])==B.dump(occurrence['third_marked_forest_code'])
            known=all(c['variable']is not None for c in corners)
            total=sum((F(c['payment'])for c in corners),F(0))if known else None
            record=dict(index=len(records),target=ti,occurrence=oi,framing_label=occurrence['label'],
                corners=corners,known=known,payment=None if total is None else str(total),
                unpaid=known and total<2)
            records.append(record);counts['occurrences']+=1
            counts['known'if known else'undefined_at_absent_forest']+=1
            counts['unpaid']+=record['unpaid']
            if record['unpaid'] and witness is None and occurrence['label']in(2,3,5):
                witness=dict(record=record['index'],target=ti,occurrence=oi,geometry=geometry(g,black,corners))
    assert witness is not None
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'scripts/{B.TAG}.py',
        f'scripts/{B.C.TAG}.py',f'results/{TRI}.json',f'results/{LP}.json']
    return dict(status='PASS_EXACT_DOUBLE_H4_TRIANGLE_SEPARATION',dimension=13,
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),records=records,witness=witness,
        minimum_known_payment=str(min(F(r['payment'])for r in records if r['known'])),
        scope='All1395 occurrences in the complete three-target selected-endpoint triangle domain are evaluated exactly. The displayed witness is an actual compact local triangle in an admissible integral-reflection chamber, not a compact finite-volume ambient polytope. An unpaid witness disproves coverage of the saved partial payment vector, without invalidating its certificate on its stated restricted rows or proving that no repaired rule exists.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['minimum_known_payment'],d['witness']['record'],flush=True)
