#!/usr/bin/env python3
"""Cross-coface cap capacities with core orbits and arithmetic two-cycles."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_orbit_integral_cap_capacity_20261017'
TAGS=dict(slices='fresh_sqrt5_labelled_coface_base_packing_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    actual='fresh_sqrt5_actual_provenance_dual_20261017',
    quad='fresh_sqrt5_rank13_packing_quad_separation_20261017',
    raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017',
    old='fresh_sqrt5_cross_coface_cap_capacity_20261017',
    model='fresh_sqrt5_cross_coface_cap_full_model_20261017')


def sign(z):
    a,b=z
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))


def square(z):
    a,b=z;return a*a+5*b*b,2*a*b


def integral(z):
    a,b=z
    by_basis=(a-b).denominator==1 and (2*b).denominator==1
    by_polynomial=(2*a).denominator==1 and (a*a-5*b*b).denominator==1
    assert by_basis==by_polynomial
    return by_polynomial


def edge(m,a='0',b='0'):
    return (m,str(F(a)),str(F(b))) if m==0 else (m,'0','0')


def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for obj in (d,d.get('prepared',{})):
        for p,h in obj.get('sources',{}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d


def cyclotomic_orders():
    # phi(N)<=4 forces prime divisors in{2,3,5} and the shown exponents.
    out=[]
    for a,b,c in product(range(4),range(2),range(2)):
        n=2**a*3**b*5**c
        phi=(2**(a-1) if a else 1)*(2 if b else 1)*(4 if c else 1)
        if phi<=4:out.append(n)
    assert sorted(out)==[1,2,3,4,5,6,8,10,12]
    assert sorted(n//2 for n in out if n>=6 and n%2==0)==[3,4,5,6]
    return sorted(out)


def run():
    D={k:read(t) for k,t in TAGS.items()}
    S,I,A,Q,R,old,model=(D[k] for k in ('slices','incidence','actual','quad','raw','old','model'))
    assert model['solver_invocations']==1
    dual=next(c for c in model['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(a) for i,a in dual['face_multipliers']}
    profiles={p['base_type']:p for p in S['coface_slice_demands']}
    cores={c['core']:c for c in S['core_symmetry_and_subset_orbits']}
    templates={t['canonical_core']:t for t in I['prism_templates']}
    counts=Counter();groups=defaultdict(list);occurrences=[];core_checks={}
    for B in S['actual_bases']:
        if B['coface_location']!=['finite','A1',0]:continue
        counts['A1_base_occurrences']+=1
        t=B['base_type'];p=profiles[t];c=cores[p['core']]
        ci=A['provenance_masses'][B['provenance']]['canonical_prism_core']
        si=Q['canonical_cores'][ci]['source_core'];r=R['records'][si]
        T=next(z for z in templates[ci]['bases'] if z['deleted_normal']==B['deleted'])
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism'])
        if set(range(r['nodes']))!={b,q}|set(mp):continue
        counts['sole_A1_extra']+=1
        assert b!=q and b not in mp and q not in mp
        assert sorted(mp.values())==list(range(c['nodes']))
        assert T['base_core']==p['core'] and T['opposite_base_facet']==q and T['finite_extra']==['A1']
        raw_edges={ij:edge(m,*r['parameters'][f'{ij[0]},{ij[1]}']) if m==0 else edge(m)
                   for ij,m in zip(combinations(range(r['nodes']),2),r['row']['labels'])}
        def E(i,j):return raw_edges[tuple(sorted((i,j)))]
        assert all(E(q,i)==edge(2) for i in mp)
        if E(b,q)[0]!=0:continue
        tau=tuple(map(F,E(b,q)[1:]));assert sign((tau[0]-4,tau[1]))>0 and integral(tau)
        h=(2-tau[0]/2,-tau[1]/2);w=square(h)
        gap=sign((h[0]+1,h[1]))>0
        finite=sign((h[0]+2,h[1]))>0 and h not in ((F(-1),F()),(F(-1,2),F(-1,2)))
        noninteger=not integral(w)
        if not(finite or noninteger):continue
        assert sign(h)<0
        counts['eligible_occurrences']+=1;counts['old_open_gap']+=gap
        counts['forbidden_finite_angle']+=finite;counts['nonintegral_two_cycle']+=noninteger
        cg=I['cores'][p['core']]
        ce={tuple(sorted((i,j))):edge(*e) for i,j,e in cg['edges']}
        def C(i,j):return ce.get(tuple(sorted((i,j))),edge(2))
        assert all(E(i,j)==C(mp[i],mp[j]) for i,j in combinations(mp,2))
        contacts=tuple(E(b,i) for i in sorted(mp,key=mp.get));vertices=tuple(B['canonical_vertices'])
        assert len(vertices)==len(set(vertices)) and set(vertices)<=set(cg['spherical_deletions'])
        orbit=set()
        for perm in c['automorphisms']:
            assert sorted(perm)==list(range(c['nodes']))
            assert {perm[i] for i in c['black']}==set(c['black'])
            assert all(C(i,j)==C(perm[i],perm[j]) for i,j in combinations(range(c['nodes']),2))
            moved=[None]*c['nodes']
            for i in range(c['nodes']):moved[perm[i]]=contacts[i]
            orbit.add((tuple(moved),tuple(sorted(perm[i] for i in vertices))))
        assert (contacts,vertices) in orbit
        assert len(c['automorphisms'])%len(orbit)==0
        key=(t,min(orbit),tuple(map(str,tau)))
        o=dict(base_index=B['index'],base_type=t,provenance=B['provenance'],face=B['face'],
            core=p['core'],canonical_prism_core=ci,source_core=si,deleted_normal=b,coface_root=q,
            core_isomorphism=sorted(map(list,mp.items())),contacts=contacts,base_vertices=vertices,
            automorphism_count=len(c['automorphisms']),orientation_orbit=sorted(orbit),
            capacity=len(orbit),dotted_squared_product=list(map(str,tau)),
            repeated_wall_product=list(map(str,h)),two_cycle=list(map(str,w)),
            two_cycle_trace=str(2*w[0]),two_cycle_norm=str(w[0]*w[0]-5*w[1]*w[1]),
            old_open_gap=gap,forbidden_finite_angle=finite,nonintegral_two_cycle=noninteger,
            arithmetic_needed=not finite)
        occurrences.append(o);groups[key].append(o)
        core_checks[p['core']]=dict(core=p['core'],nodes=c['nodes'],automorphism_count=len(c['automorphisms']))
    columns=[];required=defaultdict(F)
    for key,rows in sorted(groups.items()):
        p=profiles[key[0]];k=rows[0]['capacity']
        assert all(z['capacity']==k and z['orientation_orbit']==rows[0]['orientation_orbit'] for z in rows)
        co=Counter(z['face'] for z in rows);pr=Counter(z['provenance'] for z in rows)
        load=sum((m*alpha.get(fi,F()) for fi,m in co.items()),F());mass=alpha.get(p['simplex_face'],F())
        required[p['simplex_face']]=max(required[p['simplex_face']],load/k)
        columns.append(dict(index=len(columns),base_type=key[0],simplex_face=p['simplex_face'],
            simplex_dimension=p['simplex_dimension'],core=p['core'],capacity=k,
            source_bases=sorted(z['base_index'] for z in rows),
            positive_face_coefficients=sorted(map(list,co.items())),provenance_coefficients=sorted(map(list,pr.items())),
            orientation_orbit=rows[0]['orientation_orbit'],dotted_squared_product=rows[0]['dotted_squared_product'],
            repeated_wall_product=rows[0]['repeated_wall_product'],two_cycle=rows[0]['two_cycle'],
            forbidden_finite_angle=rows[0]['forbidden_finite_angle'],nonintegral_two_cycle=rows[0]['nonintegral_two_cycle'],
            load=str(load),initial_simplex_mass=str(mass),slack=str(k*mass-load)))
    indexed={tuple(z['source_bases']):z for z in columns}
    for z in old['columns']:
        n=indexed[tuple(z['source_bases'])]
        assert n['capacity']==z['capacity']==1 and n['simplex_face']==z['simplex_face']
        assert n['positive_face_coefficients']==z['positive_face_coefficients']
    violations=[z['index'] for z in columns if F(z['slack'])<0]
    assert len(occurrences)==len(columns)==815 and len(violations)==19
    assert sum(F(z['load'])>0 for z in columns)==167
    paths=[f'scripts/{TAG}.py']+[f'results/{t}.json' for t in TAGS.values()]
    paths+=['results/fresh_sqrt5_cross_coface_cap_capacity_20261017.md',
            'results/fresh_sqrt5_cross_coface_cap_capacity_20261017_review.md']
    return dict(status='PASS_EXACT_ORBIT_AND_INTEGRALITY_CAP_CAPACITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(counts,columns=len(columns),positive_columns=167,violations=19,
                    inherited_columns=len(old['columns']),new_columns=len(columns)-len(old['columns'])),
        cyclotomic_orders_degree_at_most_two=cyclotomic_orders(),core_checks=list(core_checks.values()),
        occurrences=occurrences,columns=columns,violating_columns=violations,
        capacity_histogram=sorted(Counter(c['capacity'] for c in columns).items()),
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)] for i,a in sorted(required.items())],
        input_bound=dual['bound'],new_repaired_bound=None,optimizer_invocations=0,
        scope='Selected arithmetic cap family with sole isolated internal A1 after base deletion. Colored-core orbit size bounds orientations; a repeated orientation has a forbidden Coxeter product or nonintegral Vinberg two-cycle. All old241columns preserved. This invalidates the current exact dual in the enlarged815-cap model, but computes neither a repair nor an optimum nor a dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    r=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert r==json.loads(out.read_text())
    else:out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['counts'],flush=True)
