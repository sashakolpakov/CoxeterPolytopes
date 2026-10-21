#!/usr/bin/env python3
"""Exact finite-rank capacities for caps on three isolated A1 copies."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import fresh_sqrt5_multispectator_cap_capacity_20261017 as ARITH
V=ARITH.V
K=ARITH.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_triple_a1_cap_capacity_20261017'
TAGS=dict(ARITH.TAGS,model='fresh_sqrt5_multispectator_cap_joint_ports_20261017')
PAIR_TYPES=[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2),(1,2)]
PATTERN_NAMES=['j0_c0','j1_c0','j2_c0','j0_c1','j1_c1','j0_c2','j1_c2']
read=ARITH.read
integral=ARITH.integral
ksqrt=ARITH.ksqrt
CLIQUE_CACHE={}


def pair_gate(h,c,tau):
    """Necessary acuteness, integrality AND the fixed ground-field test.

    In norm2 normalization the wall two-cycle is the squared product.
    If sqrt(tau) is not in K and its coefficient in that square is nonzero,
    the square is outside K even when it is an algebraic integer.
    """
    row=ARITH.pair_gate(h,c,tau)
    row['passes_acuteness_and_integrality']=row['allowed_by_necessary_tests']
    row['two_cycle_in_Qsqrt5']=row['branch']!='relative_quadratic_norm'
    row['allowed_by_necessary_tests']=row['passes_acuteness_and_integrality'] and row['two_cycle_in_Qsqrt5']
    return row


def cap(pattern,k):
    """Exact maximum on the complete finite labelled spectator domain.

    Sym(k) is transitive on positions (q,{p,r}), so a maximum clique may
    be assumed to contain position0. The saved include/exclude recurrence
    is complete, without a heuristic cutoff or floating optimization.
    """
    key=(tuple(pattern),k)
    if key in CLIQUE_CACHE:return CLIQUE_CACHE[key]
    assert k>=3
    positions=[(q,tuple(P))for q in range(k)for P in combinations([i for i in range(k)if i!=q],2)]
    assert positions[0]==(0,(1,2))
    adj=[];pair_counts=Counter()
    for i,(q,P0) in enumerate(positions):
        P=set(P0);mask=0
        for j,(r,R0) in enumerate(positions):
            if i==j:continue
            if q==r:
                if i<j:pair_counts['same_q']+=1
                continue
            R=set(R0);pair=(len(P&R),int(r in P)+int(q in R))
            assert pair in PAIR_TYPES
            z=PAIR_TYPES.index(pair)
            if i<j:pair_counts[PATTERN_NAMES[z]]+=1
            if pattern[z]:mask|=1<<j
        adj.append(mask)
    assert all(bool(adj[i]>>j&1)==bool(adj[j]>>i&1)for i in range(len(adj))for j in range(len(adj)))
    assert sum(pair_counts.values())==len(positions)*(len(positions)-1)//2
    proof={0:(0,())}
    def search(mask):
        if mask in proof:return proof[mask]
        bit=mask&-mask;v=bit.bit_length()-1;rest=mask^bit
        a,aw=search(rest);b,bw=search(rest&adj[v]);b+=1
        answer=(b,(v,)+bw)if b>a else(a,aw)
        proof[mask]=answer;return answer
    value,witness=search(adj[0]);value+=1;witness=(0,)+witness
    assert len(witness)==value and all(adj[i]>>j&1 for i,j in combinations(witness,2))
    # Every actual position is mapped to0 by one literal spectator permutation.
    maps=[]
    for q,P in positions:
        order=[q,*P]+[i for i in range(k)if i not in {q,*P}]
        perm=[None]*k
        for j,i in enumerate(order):perm[i]=j
        assert perm[q]==0 and sorted(perm[p]for p in P)==[1,2]
        assert sorted(perm)==list(range(k));maps.append(perm)
    # Recheck each saved recurrence from its children and fixed literal graph.
    for mask,(best,_) in sorted(proof.items()):
        if not mask:assert best==0;continue
        bit=mask&-mask;v=bit.bit_length()-1;rest=mask^bit
        assert best==max(proof[rest][0],1+proof[rest&adj[v]][0])
    row=dict(pattern=list(pattern),A1_count=k,position_count=len(positions),
        positions=positions,adjacency_masks_hex=[hex(m)for m in adj],
        pair_type_counts=dict(pair_counts),fixed_position=0,
        transitivity_maps_to_fixed=maps,maximum=value,
        witness_indices=witness,witness_positions=[positions[i]for i in witness],
        recurrence=[dict(mask_hex=hex(m),maximum=z[0],witness=list(z[1]))for m,z in sorted(proof.items())],
        scope='Complete exact finite spectator count only; graph uses necessary pair tests, so its clique maximum is a valid geometric upper bound, not a realization claim.')
    CLIQUE_CACHE[key]=row;return row


def run():
    D={k:read(t) for k,t in TAGS.items()}
    S,I,A,Q,R,M=(D[k] for k in ('slices','incidence','actual','quad','raw','model'))
    assert M['status'].startswith('PASS')
    alpha={i:F(a)for i,a in M['face_multipliers']};assert all(a>=0 for a in alpha.values())
    profiles={p['base_type']:p for p in S['coface_slice_demands']}
    savedcores={c['core']:c for c in S['core_symmetry_and_subset_orbits']}
    templates={t['canonical_core']:t for t in I['prism_templates']}
    counts=Counter();groups=defaultdict(list);occurrences=[];cores={};tau_checks={};seen=set()
    for B in S['actual_bases']:
        if B['coface_location']!=['finite','A1',0]:continue
        counts['A1_coface_base_occurrences']+=1
        ci=A['provenance_masses'][B['provenance']]['canonical_prism_core']
        si=Q['canonical_cores'][ci]['source_core'];r=R['records'][si]
        T=next(z for z in templates[ci]['bases'] if z['deleted_normal']==B['deleted'])
        if T['finite_extra']!=['A1','A1','A1']:continue
        counts['three_A1_internal_occurrences']+=1
        key0=(B['provenance'],B['deleted']);assert key0 not in seen;seen.add(key0)
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism'])
        ptype=profiles[B['base_type']];typ=I['base_types'][B['base_type']];core=ptype['core']
        assert core==typ['core']==T['base_core'] and T['opposite_base_facet']==q
        assert ptype['simplex_dimension']==typ['simplex_dimension']
        extra=set(range(r['nodes']))-set(mp)-{b};assert q in extra and len(extra)==3
        ordinary=sorted(extra-{q});assert len(ordinary)==2 and b not in extra
        graph=V.raw_graph(r);cg=V.core_graph(I['cores'][core])
        assert sorted(mp.values())==list(range(len(cg)))
        assert all(graph[i][j]==cg[mp[i]][mp[j]] for i,j in combinations(mp,2))
        assert all(graph[x][i]==V.RIGHT for x in extra for i in set(range(r['nodes']))-{b,x})
        assert all(graph[b][p]==V.edge(3) for p in ordinary) and graph[b][q][0] in (0,5,10)
        tau=V.squared(graph[b][q]);assert tau.sign()>0 and integral(tau)
        assert tau!=1
        contacts=tuple(graph[b][i]for i in sorted(mp,key=mp.get))
        assert all(e[0] in (2,3)for e in contacts)
        ports=[i for i,e in enumerate(contacts)if e[0]==3];assert ports
        vertices=set();base_pairs=set()
        for left,right in Q['canonical_cores'][ci]['quad_partitions']:
            if b in left:side,other=left,right
            elif b in right:side,other=right,left
            else:continue
            assert set(side)=={b,q}
            vertices.update(other);base_pairs.update(tuple(sorted((b,x)))for x in other)
        assert vertices and vertices<=set(mp)
        assert base_pairs<={tuple(sorted(z))for z in Q['canonical_cores'][ci]['pd_deletions']}
        canonical=tuple(sorted(mp[x]for x in vertices))
        assert list(canonical)==B['canonical_vertices']
        assert len(canonical)==B['prism_dimension']==templates[ci]['prism_dimension']
        assert set(canonical)<=set(I['cores'][core]['spherical_deletions'])
        assert sorted(mp[x]for x in T['full_spherical_deletions'])==I['cores'][core]['spherical_deletions']
        k=Counter(typ['spectator_forest'])['A1'];assert k>=3
        rank=sum(map(V.rank,typ['spectator_forest']))
        assert len(cg)+rank==14 and r['nodes']+rank-3==15
        if core not in cores:
            aut=V.automorphisms(cg)
            saved=savedcores[core]
            assert sorted(aut)==sorted(saved['automorphisms'])
            assert all({perm[i]for i in saved['black']}==set(saved['black'])for perm in aut)
            H,scales=V.physical_core(cg);J,det=V.inverse(H);assert det.sign()<0
            assert all(sum((H[i][k0]*J[k0][j]for k0 in range(len(H))),K())==int(i==j)
                       for i in range(len(H))for j in range(len(H)))
            pd_roots=[]
            for root in range(len(cg)):
                ids=[i for i in range(len(cg))if i!=root]
                if V.pd([[H[i][j]for j in ids]for i in ids]):pd_roots.append(root)
            assert pd_roots==I['cores'][core]['spherical_deletions'] and pd_roots
            cores[core]=dict(aut=aut,H=H,scales=scales,J=J,det=det,pd_roots=pd_roots)
        C=cores[core];J,scales=C['J'],C['scales']
        if scales is None:
            norm=sum((J[i][j]for i in ports for j in ports),K());projection_gauge=None
        else:
            ref=scales[ports[0]];ratios={i:ksqrt(scales[i]/ref)for i in ports}
            assert all(x is not None and x.sign()>0 for x in ratios.values())
            norm=ref*sum((ratios[i]*ratios[j]*J[i][j]for i in ports for j in ports),K())
            projection_gauge=dict(reference_scale=ref.out(),root_scale_ratios=[[i,ratios[i].out()]for i in ports])
        h=1-tau/2;assert norm==h and h+tau/2==1
        tau_key=tuple(tau.out())
        if tau_key not in tau_checks:
            pair_tests={name:pair_gate(h+F(j,2),F(c,2),tau) for name,(j,c) in zip(PATTERN_NAMES,PAIR_TYPES)}
            pattern=tuple(pair_tests[x]['allowed_by_necessary_tests']for x in PATTERN_NAMES)
            tau_checks[tau_key]=dict(tau=tau.out(),h=h.out(),pairs=pair_tests,pattern=pattern)
        tests=tau_checks[tau_key];clique=cap(tests['pattern'],k);single=clique['maximum'];proofkind='exact_actual_rank13_position_graph_clique'
        orbit=set()
        for perm in C['aut']:
            moved=[None]*len(cg)
            for i in range(len(cg)):moved[perm[i]]=contacts[i]
            orbit.add((tuple(moved),tuple(sorted(perm[i]for i in canonical))))
        assert (contacts,canonical)in orbit and len(C['aut'])%len(orbit)==0
        o=dict(base_index=B['index'],base_type=B['base_type'],provenance=B['provenance'],face=B['face'],
            simplex_face=ptype['simplex_face'],core=core,canonical_prism_core=ci,source_core=si,
            deleted_normal=b,coface_root=q,ordinary_A1_roots=ordinary,core_isomorphism=sorted(map(list,mp.items())),
            core_contacts=contacts,canonical_vertices=canonical,raw_vertex_omitted_pairs=sorted(map(list,base_pairs)),
            tau=tau.out(),projection_norm=norm.out(),projection_gauge=projection_gauge,
            A1_count=k,automorphism_count=len(C['aut']),orientation_orbit=sorted(orbit),
            orientation_count=len(orbit),pair_pattern=list(tests['pattern']),single_orientation_capacity=single,
            proof_kind=proofkind,capacity=single*len(orbit))
        occurrences.append(o)
        groups[(B['base_type'],min(orbit),tau_key)].append(o)
    columns=[];required=defaultdict(F)
    for key,rows in sorted(groups.items()):
        r=rows[0];p=profiles[r['base_type']];capacity=r['capacity']
        assert all(z['capacity']==capacity and z['orientation_orbit']==r['orientation_orbit']for z in rows)
        co=Counter(z['face']for z in rows);pr=Counter(z['provenance']for z in rows)
        assert len({z['base_index']for z in rows})==len(rows)
        load=sum((m*alpha.get(fi,F())for fi,m in co.items()),F());mass=alpha.get(p['simplex_face'],F())
        required[p['simplex_face']]=max(required[p['simplex_face']],load/capacity)
        columns.append(dict(index=len(columns),kind='three_isolated_A1_cap_profile',
            base_type=r['base_type'],simplex_face=p['simplex_face'],simplex_dimension=p['simplex_dimension'],core=r['core'],
            A1_count=r['A1_count'],tau=r['tau'],orientation_orbit=r['orientation_orbit'],
            orientation_count=r['orientation_count'],single_orientation_capacity=r['single_orientation_capacity'],
            pair_pattern=r['pair_pattern'],proof_kind=r['proof_kind'],capacity=capacity,
            source_bases=sorted(z['base_index']for z in rows),
            positive_face_coefficients=sorted(map(list,co.items())),provenance_coefficients=sorted(map(list,pr.items())),
            load=str(load),initial_simplex_mass=str(mass),slack=str(capacity*mass-load)))
    violations=[c['index']for c in columns if F(c['slack'])<0]
    assert counts['three_A1_internal_occurrences']==len(occurrences)==93
    assert len(tau_checks)==5
    counts.update(columns=len(columns),positive_columns=sum(F(c['load'])>0 for c in columns),
                  violations=len(violations),cores=len(cores),tau_values=len(tau_checks))
    paths=[Path(__file__),Path(V.__file__),Path(ARITH.__file__)]+[ROOT/'results'/f'{t}.json'for t in TAGS.values()]
    paths+=[ROOT/'results/fresh_prism_simplex_shared_capacity_20261017.md',
            ROOT/'results/fresh_sqrt5_cross_coface_cap_capacity_20261017_review.md']
    return dict(status='PASS_EXACT_THREE_A1_MULTISPECTATOR_CAP_CAPACITIES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),pair_pattern_order=PATTERN_NAMES,
        tau_pair_tests=[tau_checks[k]for k in sorted(tau_checks)],
        core_checks=[dict(core=c,nodes=len(z['H']),gram=[[x.out()for x in row]for row in z['H']],
                         gauge_scales=None if z['scales']is None else[x.out()for x in z['scales']],
                         inverse=[[x.out()for x in row]for row in z['J']],determinant=z['det'].out(),
                         spherical_deletions=z['pd_roots'],automorphisms=z['aut'])for c,z in sorted(cores.items())],
        complete_finite_position_cliques=[CLIQUE_CACHE[k] for k in sorted(CLIQUE_CACHE)],
        occurrences=occurrences,columns=columns,violating_columns=violations,
        orientation_histogram=sorted(Counter(c['orientation_count']for c in columns).items()),
        capacity_histogram=sorted(Counter(c['capacity']for c in columns).items()),
        occurrence_pattern_histogram=[dict(pattern=list(k),count=v)for k,v in sorted(Counter(tuple(o['pair_pattern'])for o in occurrences).items())],
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)]for i,a in sorted(required.items())],
        input_bound=M['bound'],new_repaired_bound=None,optimizer_invocations=0,
        arithmetic_ground_field='Q(sqrt5)',
        scope='Selected actual dotted-prism bases with three internal isolated A1 components after base deletion, two simple3 contacts and one distinct finite5/10 or dotted contact. All93 raw maps, projections, colored core automorphisms and profile multiplicities are reconstructed. Pair compatibility uses necessary acuteness, arithmetic two-cycle integrality and membership in the fixed ground field Q(sqrt5). Complete exact clique recurrences cover every actual rank13 spectator count, with no general all-rank claim. All capacities include colored-core orientation factors. No optimizer, witness repair, global optimum or dimension theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'orientations',out['orientation_histogram'],flush=True)
