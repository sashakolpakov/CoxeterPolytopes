#!/usr/bin/env python3
"""Simultaneous coface/subset base packing: exact scan, no optimization."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_labelled_coface_base_packing_20261017'
FLAGS = 'fresh_sqrt5_flag_refined_base_capacity_20261017'
POINTED = 'fresh_sqrt5_pointed_simplex_base_capacity_20261017'
INCIDENCE = 'fresh_sqrt5_prism_simplex_base_incidence_20261017'
SCOUT = 'fresh_sqrt5_coface_integer_base_scout_20261017'
BASE = 'fresh_sqrt5_joint_five_port_diagnostic_20261017'
MODEL = 'fresh_sqrt5_flag_local_packing_full_model_20261017'


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def read(stem):
    return json.loads((ROOT / 'results' / (stem + '.json')).read_text())


def freeze(x):
    return json.loads(json.dumps(x))


def automorphisms(core):
    n = core['nodes']; black = set(core['spherical_deletions'])
    g = {i: {} for i in range(n)}
    for i, j, label in core['edges']:
        assert i < j
        g[i][j] = g[j][i] = tuple(label)
    colors = {i: int(i in black) for i in g}
    for _ in range(n):
        codes = {i: (colors[i], tuple(sorted((m, colors[j]) for j, m in g[i].items()))) for i in g}
        ids = {x: j for j, x in enumerate(sorted(set(codes.values())))}
        colors = {i: ids[codes[i]] for i in g}
    classes = [[i for i in g if colors[i] == c] for c in sorted(set(colors.values()))]
    autos = []; tested = 0
    for targets in product(*(permutations(c) for c in classes)):
        mp = dict(pair for left, right in zip(classes, targets) for pair in zip(left, right))
        tested += 1
        if not all(g[i].get(j, 2) == g[mp[i]].get(mp[j], 2) for i, j in combinations(g, 2)):
            continue
        assert {mp[i] for i in black} == black
        autos.append([mp[i] for i in range(n)])
    assert autos and len(autos) == len(set(map(tuple, autos)))
    assert list(range(n)) in autos
    # Exhaustive color-class permutations are complete; closure is an extra check.
    aset = set(map(tuple, autos))
    assert all(tuple(a[b[i]] for i in range(n)) in aset for a in autos for b in autos)
    return dict(core=core['index'], nodes=n, black=sorted(black),
                white=sorted(set(g)-black), color_classes=classes,
                tested_permutations=tested, automorphisms=autos)


def orbits(elements, action):
    todo = set(elements); result = []
    while todo:
        x = min(todo); orbit = sorted(set(action(x)))
        assert x in orbit and set(orbit) <= todo
        result.append(orbit); todo -= set(orbit)
    return result


def run():
    flags, pointed, incidence, scout, base, model = map(read, [FLAGS, POINTED, INCIDENCE, SCOUT, BASE, MODEL])
    for data in [flags, pointed, incidence, scout, base]:
        for p, h in data['sources'].items():
            assert sha(p) == h, p
    alpha = {i: F(a) for i, a in base['face_multipliers']}
    assert len(alpha) == 2219 and all(a > 0 for a in alpha.values())
    pf = dict(model['prepared']['provenance_to_face'])
    assert pf == dict(pointed['provenance_to_face'])
    assert len(pf) == len(set(pf.values())) == 1858
    assert len(set(pf.values()) & set(alpha)) == 355
    profiles = {p['base_type']: p for p in pointed['simplex_profiles']}
    cores = {ci: automorphisms(incidence['cores'][ci]) for ci in sorted({p['core'] for p in profiles.values()})}
    assert len(cores) == 287
    core_orbits, subset_orbits, q_orbits = {}, {}, {}
    for ci, c in cores.items():
        A, black, white = c['automorphisms'], c['black'], c['white']
        subsets = [tuple(s) for k in range(2, len(black)+1) for s in combinations(black, k)]
        subset_orbits[ci] = orbits(subsets, lambda T: [tuple(sorted(a[i] for i in T)) for a in A])
        q_orbits[ci] = orbits(white, lambda q: [a[q] for a in A])
        core_orbits[ci] = orbits([(q, T) for q in white for T in subsets],
            lambda x: [(a[x[0]], tuple(sorted(a[i] for i in x[1]))) for a in A])
        assert sum(map(len, core_orbits[ci])) == len(white)*len(subsets)
        c['subset_orbits'] = subset_orbits[ci]
        c['coface_root_orbits'] = q_orbits[ci]
        c['simultaneous_coface_subset_orbits'] = core_orbits[ci]
    grouped = defaultdict(list)
    for o in flags['literal_flag_occurrences']:
        grouped[o['provenance'], o['deleted']].append(o)
    saved_bases = {(b['provenance'], b['deleted']): b for b in scout['actual_bases']}
    bases, by_profile = [], defaultdict(list)
    for key, rows in sorted(grouped.items()):
        j, deleted = key; first = rows[0]; t = first['base_type']; p = profiles[t]
        vertices = tuple(sorted(o['canonical_vertex'] for o in rows))
        loc = tuple(first['coface_location']); q = first['coface_root']
        assert all(o['base_type'] == t and tuple(o['coface_location']) == loc and o['coface_root'] == q for o in rows)
        old = saved_bases[key]
        assert list(vertices) == old['canonical_vertices'] and len(vertices) == old['prism_dimension']
        assert (old['base_type'], old['simplex_face'], old['simplex_dimension'], old['coface_root']) == (t,p['simplex_face'],p['simplex_dimension'],q)
        assert len(set(vertices)) == len(vertices) and set(vertices) <= set(cores[p['core']]['black'])
        if loc[0] == 'core':
            assert loc[1] in cores[p['core']]['white']
        else:
            assert loc[0] == 'finite' and Counter(p['spectator_forest'])[loc[1]] > 0
            assert loc[2] < len(flags['finite_root_orbits'][loc[1]])
        b = dict(index=len(bases), provenance=j, face=pf[j], deleted=deleted,
                 base_type=t, coface_root=q, coface_location=list(loc),
                 canonical_vertices=list(vertices), prism_dimension=len(vertices))
        bases.append(b); by_profile[t].append(b)
    assert set(grouped) == set(saved_bases) and len(bases) == 6016

    columns, violations, slices, slice_checks = [], [], [], []
    positive_loads = 0; source_count = 0
    identity_checks = []; required = {}; histogram = Counter()
    for serial, (t, p) in enumerate(sorted(profiles.items())):
        ci, si, d = p['core'], p['simplex_face'], p['simplex_dimension']
        c = cores[ci]; m = len(c['black']); assert m == d+1
        copies = Counter(p['spectator_forest'])
        finite_count = sum(num*sum(len(o['roots']) for o in flags['finite_root_orbits'][name]) for name,num in copies.items())
        assert len(c['white']) + finite_count == 13-d
        entries = []
        for oi, orbit in enumerate(core_orbits[ci]):
            entries.append((('core', oi), len(orbit), len(orbit[0][1]), orbit))
        for name, num in sorted(copies.items()):
            for roi, rr in enumerate(flags['finite_root_orbits'][name]):
                for oi, orbit in enumerate(subset_orbits[ci]):
                    entries.append((('finite', name, roi, oi), num*len(rr['roots'])*len(orbit), len(orbit[0]), orbit))
        totals = defaultdict(Counter); cap_totals = Counter()
        for label, orbit_size, u, orbit in entries:
            for divisor in range(2,u+1):
                cap = orbit_size*(u//divisor); terms = Counter()
                for b in by_profile[t]:
                    B = set(b['canonical_vertices']); loc = b['coface_location']
                    if label[0] == 'core':
                        if loc[0] != 'core':continue
                        co = sum(len(B & set(T))//divisor for q,T in orbit if q == loc[1])
                    else:
                        if tuple(loc) != label[:3]:continue
                        co = sum(len(B & set(T))//divisor for T in orbit)
                    if co:terms[b['provenance']] += co
                totals[u,divisor].update(terms);cap_totals[u,divisor] += cap
                load = sum((n*alpha.get(pf[j],F()) for j,n in terms.items()),F())
                old = alpha.get(si,F());slack = cap*old-load
                row = dict(index=len(columns),base_type=t,simplex_face=si,core=ci,
                    simplex_dimension=d,orbit_label=list(label),subset_size=u,
                    divisor=divisor,orbit_size=orbit_size,capacity=cap,
                    load=str(load),simplex_mass=str(old),slack=str(slack))
                columns.append(row);histogram[d,u,divisor] += 1
                positive_loads += load>0
                if slack < 0:
                    row = dict(row,provenance_coefficients=[[j,n] for j,n in sorted(terms.items())],
                        positive_face_coefficients=[[pf[j],n] for j,n in sorted(terms.items())])
                    violations.append(row);required[si] = max(required.get(si,old),load/cap)
        for u in range(2,m+1):
            for divisor in range(2,u+1):
                expected = Counter()
                for b in by_profile[t]:
                    k=b['prism_dimension']
                    co=sum((h//divisor)*comb(k,h)*comb(m-k,u-h)
                           for h in range(max(0,u-(m-k)),min(k,u)+1))
                    if co:expected[b['provenance']] += co
                assert totals[u,divisor] == expected
                assert cap_totals[u,divisor] == (13-d)*comb(m,u)*(u//divisor)
                identity_checks.append([t,u,divisor,cap_totals[u,divisor],sum(expected.values())])

        # Export a complete orbit-averaged subset demand for each coface class.
        # Its mass is summed over that coface orbit; cap = |orbit(q)| alpha_S.
        slice_defs = []
        for oi, oq in enumerate(q_orbits[ci]):
            slice_defs.append((('core',oi),len(oq),oq))
        for name,num in sorted(copies.items()):
            for roi,rr in enumerate(flags['finite_root_orbits'][name]):
                slice_defs.append((('finite',name,roi),num*len(rr['roots']),rr['roots']))
        root_bits={v:1<<i for i,v in enumerate(c['black'])}
        for label, mult, oq in slice_defs:
            block_coeff=defaultdict(lambda:defaultdict(F));source_bases=[]
            rep=oq[0]
            for b in by_profile[t]:
                loc=b['coface_location'];B=b['canonical_vertices']
                if label[0]=='core':
                    if loc[0]!='core' or loc[1] not in oq:continue
                    maps=[a for a in c['automorphisms'] if a[loc[1]]==rep]
                    assert maps and len(maps)*mult==len(c['automorphisms'])
                else:
                    if tuple(loc)!=label:continue
                    maps=c['automorphisms']
                counts=Counter(sum(root_bits[a[v]] for v in B) for a in maps)
                assert sum(counts.values())==len(maps)
                source_bases.append(b['index'])
                for mask,n in counts.items():block_coeff[mask][b['provenance']] += F(n,len(maps))
            marginals={mask:sum((co*alpha.get(pf[j],F()) for j,co in terms.items()),F())
                       for mask,terms in block_coeff.items()}
            mass=sum(marginals.values(),F())
            assert mass==sum((alpha.get(bases[z]['face'],F()) for z in source_bases),F())
            # Slices include zero-mass types to specify the complete saved domain.
            slices.append(dict(index=len(slices),base_type=t,simplex_face=si,core=ci,
                simplex_dimension=d,coface_label=list(label),coface_multiplicity=mult,
                representative_coface_root=rep,vertex_roots=c['black'],
                capacity=str(mult*alpha.get(si,F())),simplex_mass=str(alpha.get(si,F())),
                source_bases=source_bases,
                block_coefficients=[[mask,[[j,str(co)] for j,co in sorted(terms.items())]] for mask,terms in sorted(block_coeff.items())],
                block_marginals=[[mask,str(a)] for mask,a in sorted(marginals.items()) if a],
                total_block_mass=str(mass)))
            source_count += len(source_bases)
        if serial%500==0:print('LABELLED_COFACE',serial+1,'/',len(profiles),'columns',len(columns),'violations',len(violations),flush=True)
    assert source_count==len(bases)
    assert sum((F(s['total_block_mass']) for s in slices),F())==sum((alpha.get(b['face'],F()) for b in bases),F())
    paths=[]
    for stem in [FLAGS,POINTED,INCIDENCE,SCOUT,BASE,MODEL]:
        paths.extend([f'scripts/{stem}.py',f'results/{stem}.json'])
    paths.extend([f'scripts/{TAG}.py','results/fresh_prism_simplex_shared_capacity_20261017.md',
                  'results/fresh_vertex_weighted_prism_base_capacity_20261017.md'])
    return dict(status='PASS_EXACT_LABELLED_COFACE_SUBSET_PACKING_DIAGNOSTIC',
        sources={p:sha(p) for p in paths},dimension=13,optimizer_invocations=0,
        input_bound=base['bound'],alpha_changed=False,mu_changed=False,
        counts=dict(simplex_profiles=len(profiles),core_types=len(cores),actual_provenances=len(pf),
            actual_bases=len(bases),columns=len(columns),positive_load_columns=positive_loads,
            violations=len(violations),simplex_rows_requiring_increase=len(required),
            coface_slices=len(slices),positive_mass_slices=sum(bool(s['block_marginals']) for s in slices),
            total_automorphisms=sum(len(c['automorphisms']) for c in cores.values()),
            tested_permutations=sum(c['tested_permutations'] for c in cores.values())),
        core_symmetry_and_subset_orbits=list(cores.values()),actual_bases=bases,
        column_checks=columns,violating_columns=violations,
        necessary_simplex_mass_lower_bounds=[[si,str(alpha.get(si,F())),str(a)] for si,a in sorted(required.items())],
        histogram=[[d,u,r,n] for (d,u,r),n in sorted(histogram.items())],
        sum_over_subset_identities=identity_checks,coface_slice_demands=slices,
        face_multipliers=base['face_multipliers'],
        scope='Exact simultaneous(q,T) orbit floor(|B intersection T|/r) cuts on every represented saved base type, with every identical finite spectator copy. All 1858 actual provenances use a common fixed alpha;355 have positive mass. No alpha/mu repair or global LP. Slice demands are orbit-total averages for optional future full disjoint-base packing tests; this scan alone does not claim that full cone is exhausted or give an ambient dimension bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=freeze(run());path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert json.loads(path.read_text())==result
    else:path.write_text(json.dumps(result,separators=(',',':'),sort_keys=True)+'\n')
    print(result['status'],result['counts'],flush=True)
