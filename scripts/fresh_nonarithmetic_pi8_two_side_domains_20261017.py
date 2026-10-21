#!/usr/bin/env python3
"""Complete continuous endpoint domains for the two new r8/S6 open edges."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_pi8_tetrahedron_completion_20261017 as C

P, U, I, D = C.P, C.U, C.I, C.D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_pi8_two_side_domains_20261017'
GATE='fresh_nonarithmetic_pi8_triangle_facet_gate_20261017'


def run():
    z=C.geometry();pair=z['pair'];normals=z['normals'][:]
    cap=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    gate=json.loads((ROOT/f'results/{GATE}.json').read_text())
    for data in(cap,gate):
        for path,digest in data['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    normals += [[I(*x)for x in row['computational']]for row in cap['new_normals'][1:]]
    assert len(normals)==58
    # The already certified finite point is the common old endpoint.
    physical_endpoint=[U.rnd(z['t8']*a-b/z['t8'])for a,b in zip(z['a'],z['b'])]
    endpoint=z['coordinates'](physical_endpoint)
    endpoint_incident=[0,1,2,3,5,6,7,30,55]
    endpoint_products=[U.rnd(pair(endpoint,n))for n in normals]
    assert pair(endpoint,endpoint).h<0 and pair(endpoint,z['h']).h<0
    for i,x in enumerate(endpoint_products):
        if i in endpoint_incident:assert x.l<=0<=x.h
        else:assert x.h<0
    records=[]
    for domain in gate['endpoint_domains']:
        side=domain['side'];ids=[0,1,2,3]+domain['nodes'];frame=[normals[i]for i in ids]
        gram=[[pair(a,b)for b in frame]for a in frame];pd=P.pd(gram)
        inverse,pivots=P.C53.inverse(gram)
        def project(v):
            return P.sub(v,P.combine(U.imv(inverse,[pair(v,n)for n in frame]),frame))
        k0=project(normals[51]);q0=pair(k0,k0);assert q0.h<0
        if pair(k0,z['h']).l>0:k0=[-x for x in k0]
        assert pair(k0,z['h']).h<0
        k0=[U.rnd(x/U.sqrt(-q0))for x in k0]
        for coordinate in range(10):
            seed=[I(int(i==coordinate))for i in range(10)];v=project(seed);a=pair(v,k0)
            k1=[U.rnd(x+a*y)for x,y in zip(v,k0)];q1=pair(k1,k1)
            if q1.l>0:break
        else:raise AssertionError('positive kernel direction')
        k1=[U.rnd(x/U.sqrt(q1))for x in k1]
        duals=[P.combine([inverse[j][4+i]for j in range(8)],frame)for i in range(4)]
        response=[[U.rnd(pair(d,n))for d in duals]for n in normals]
        orientations=[]
        for sign in(-1,1):
            aa=[U.rnd((x+sign*y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
            bb=[U.rnd((-x+sign*y)/U.sqrt(I(2)))for x,y in zip(k0,k1)]
            a_products=[U.rnd(pair(aa,n))for n in normals]
            b_products=[U.rnd(pair(bb,n))for n in normals]
            orientations.append(dict(kernel_side=sign,A=P.enc(aa),B=P.enc(bb),
                A_products=P.enc(a_products),B_products=P.enc(b_products),
                vertex_domain_quadratics=[dict(wall=i,quadratic=P.enc([a_products[i],I(0),-b_products[i]]))for i in range(58)if i not in ids],
                common_vertex_A=pair(aa,z['vertex']).enc(),common_vertex_B=pair(bb,z['vertex']).enc()))
        profiles=[]
        for row in domain['profiles']['finite_label_profiles']:
            contacts=[D.theta(m)for m in row['labels']]
            p=P.combine([-c for c in contacts],duals);rho=U.rnd(2-pair(p,p));assert rho.l>0
            profiles.append(dict(labels=row['labels'],magnitudes=P.enc(contacts),
                projection=P.enc(p),positive_kernel_square=rho.enc(),
                projection_products=P.enc([pair(p,n)for n in normals]),
                common_vertex_projection=pair(p,z['vertex']).enc()))
        records.append(dict(edge_sides=[side,55],full_support=ids,residual_nodes=domain['nodes'],
            finite_extension_profiles=domain['profiles'],finite_profiles=profiles,
            positive_gram_pivots=pd,inverse_pivots=pivots,dual_vectors=[P.enc(d)for d in duals],
            kernel_future_unit=P.enc(k0),kernel_positive_unit=P.enc(k1),
            kernel_future_seed=51,kernel_positive_seed=coordinate,
            native_response_matrix=[P.enc(row)for row in response],orientations=orientations,
            common_vertex_dual_products=P.enc([pair(d,z['vertex'])for d in duals])))
    assert [len(r['finite_profiles'])for r in records]==[16,36]
    assert [len(r['finite_extension_profiles']['tails'])for r in records]==[1,2]
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',
        f'scripts/{GATE}.py',f'results/{GATE}.json',f'scripts/{D.TAG}.py',f'results/{D.TAG}.json',
        f'scripts/{D.E.TAG}.py',f'results/{D.E.TAG}.json']
    return dict(status='PASS_COMPLETE_CONTINUOUS_R8_TWO_EDGE_ENDPOINT_DOMAINS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},retained_walls=58,
        computational_root_basis=z['basis_ids'],computational_gram=[P.enc(row)for row in z['exact_gram']],
        old_finite_endpoint=P.enc(endpoint),old_finite_endpoint_incident=endpoint_incident,
        old_finite_endpoint_products=P.enc(endpoint_products),
        common_persistent_vertex=P.enc(z['vertex']),domains=records,
        finite_profile_count=52,high_label_tail_count=3,orientation_branches_per_profile=2,
        normal_formula='For contact magnitudes c, p=-sum c_i*d_i, rho=2-p^2>0. For either sign s let A=(k0+s*k1)/sqrt2,B=(-k0+s*k1)/sqrt2. Every normal is r(t)=p+sqrt(rho/2)*(t*A+B/t), t>0. A.B=1, A²=B²=0.',
        vertex_formula='The corresponding future finite endpoint is x(t)=t*A-B/t; x²=-2 and x.r=0 exactly. The saved quadratic a_i*t²-b_i<=0 is precisely its old-wall inequality.',
        gate_formula='Each old product is p.n_i+sqrt(rho/2)*(a_i*t+b_i/t). Native Coxeter coherence requires it to be0,-2cos(pi/m),or<=-2. Acuteness and the actual-endpoint inequalities are saved necessary outer constraints. The complete58 offsets are retained, including56/57 and49/50/51.',
        common_interior='The saved V is an actual old finite vertex. Requiring r.V<0 is sufficient for a single or finite family of acute additions to share a strict future interior, by a sufficiently small common perturbation of V toward the old center. Mutual acuteness and all Coxeter contacts still need separate checking.',
        completeness='H4 cannot be enlarged inside a finite endpoint diagram. The other four-root supports are A3+A1 and A2+A1+A1. Their complete inverse-cost contact domains are inherited from the all-label triangular-facet certificate:16 and36 finite profiles, plus1 and2 isolated-A1 high-label tails. Positive residual gives a Lorentz2 hyperbola with exactly the two displayed orientations. No root field or high-label cutoff is imposed.',
        scope='A complete continuous family domain, not a solved mutual-contact system or a new wall addition. Finite-angle events to the remaining58 walls and coordinated polygon closure remain to be solved. The three high-label tails remain explicit continuous-label domains.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data==json.loads(out.read_text());print('PASS saved two-edge continuous-domain replay')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
