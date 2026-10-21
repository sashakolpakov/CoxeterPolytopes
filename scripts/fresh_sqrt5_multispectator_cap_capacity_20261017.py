#!/usr/bin/env python3
"""Exact arithmetic capacities for caps supported on two isolated A1 copies.

Reconstructs the selected actual base maps, normalized projections and full
colored core orientations.  Uses only necessary acuteness/integrality pair
gates, followed by uniform combinatorial directed-pair capacities. No solve.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
import json
from pathlib import Path

import fresh_sqrt5_cross_coface_cap_capacity_20261017_review as V

K=V.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_multispectator_cap_capacity_20261017'
TAGS=dict(slices='fresh_sqrt5_labelled_coface_base_packing_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    actual='fresh_sqrt5_actual_provenance_dual_20261017',
    quad='fresh_sqrt5_rank13_packing_quad_separation_20261017',
    raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017',
    model='fresh_sqrt5_profile_clique_joint_ports_20261017')
PATTERN_NAMES=['disjoint','common_head','chain','swapped']


def read(tag):
    d=json.loads((ROOT/'results'/f'{tag}.json').read_text())
    for obj in (d,d.get('prepared',{})):
        for p,h in obj.get('sources',{}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d


def rsqrt(x):
    x=F(x)
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b) if a*a==x.numerator and b*b==x.denominator else None


def ksqrt(x):
    """Complete positive square-root test in Q(sqrt5), not a numerical test."""
    x=K(x);assert x.sign()>0
    if not x.b:
        a=rsqrt(x.a)
        if a is not None:return K(a)
        b=rsqrt(x.a/5)
        return None if b is None else K(0,b)
    d=rsqrt(x.a*x.a-5*x.b*x.b)
    if d is None:return None
    for sg in (1,-1):
        a=rsqrt((x.a+sg*d)/2)
        if a:
            z=K(a,x.b/(2*a));assert z*z==x
            return z if z.sign()>0 else -z
    return None


def integral(x):
    by_basis=(x.a-x.b).denominator==1 and (2*x.b).denominator==1
    by_polynomial=(2*x.a).denominator==1 and (x.a*x.a-5*x.b*x.b).denominator==1
    assert by_basis==by_polynomial
    return by_polynomial


def signrad(a,b,t):
    assert t.sign()>0
    if not b:return a.sign()
    if not a:return b.sign()
    if a.sign()==b.sign():return a.sign()
    return a.sign()*(a*a-b*b*t).sign()


def evalK(poly,x):
    out=K()
    for z in reversed(poly):out=out*x+z
    return out


def mulrad(x,y,t):return (x[0]*y[0]+t*x[1]*y[1],x[0]*y[1]+x[1]*y[0])


def evalrad(poly,x,t):
    out=(K(),K())
    for z in reversed(poly):
        out=mulrad(out,x,t);out=(out[0]+z,out[1])
    return out


def pair_gate(h,c,t):
    """Exact product sign and algebraic integrality of its native two-cycle."""
    c=F(c);sg=signrad(h,K(c),t);root=ksqrt(t)
    if root is not None:
        assert root.sign()>0 and root*root==t
        v=h+c*root;assert v.sign()==sg
        w=v*v;poly=[w.a*w.a-5*w.b*w.b,-2*w.a,F(1)]
        assert not evalK(poly,w)
        ok=integral(w);assert ok==all(x.denominator==1 for x in poly)
        rec=dict(branch='square_root_in_Qsqrt5',square_root=root.out(),
                 product=v.out(),two_cycle=w.out(),monic_Q_polynomial=list(map(str,poly)))
    else:
        A=h*h+c*c*t;B=2*c*h
        # w=A+B sqrt(t), so Norm_{K(sqrt(t))/K}(X-w)=X²-2AX+A²-B²t.
        if not B:
            poly=[A.a*A.a-5*A.b*A.b,-2*A.a,F(1)]
            assert not evalK(poly,A)
            ok=integral(A);assert ok==all(x.denominator==1 for x in poly)
            rec=dict(branch='square_root_not_in_K_but_two_cycle_in_K',
                     two_cycle=A.out(),monic_Q_polynomial=list(map(str,poly)))
        else:
            rel=[A*A-B*B*t,-2*A,K(1)]
            assert evalrad(rel,(A,B),t)==(K(),K())
            conj=[K(x.a,-x.b) for x in rel];norm=[K() for _ in range(5)]
            for i,a in enumerate(rel):
                for j,b in enumerate(conj):norm[i+j]=norm[i+j]+a*b
            assert all(not z.b for z in norm) and norm[-1]==1
            poly=[z.a for z in norm]
            assert evalrad(poly,(A,B),t)==(K(),K())
            ok=all(z.denominator==1 for z in poly)
            rec=dict(branch='relative_quadratic_norm',two_cycle_radical=[A.out(),B.out()],
                     relative_K_polynomial=[z.out() for z in rel],
                     monic_Q_polynomial=list(map(str,poly)))
    rec.update(product_sign=sg,two_cycle_integral=ok,allowed_by_necessary_tests=(sg<=0 and ok))
    return rec


def cap(pattern,k):
    assert k>=2
    if pattern==(False,False,False,False):return 1,'pairwise_incompatible'
    if pattern==(True,False,False,True):return 2*(k//2),'disjoint_arcs_and_two_cycles'
    if pattern==(True,False,False,False):return k//2,'oriented_matching'
    if pattern==(False,False,True,False):return (3 if k>=3 else 1),'directed_triangle_or_shorter_chain'
    if pattern==(False,True,False,False):return k-1,'common_head_star'
    raise AssertionError(('unclassified compatibility pattern',pattern,k))


def pair_kind(x,y):
    q,p=x;r,s=y
    assert q!=p and r!=s and x!=y
    if q==r:return 'same_tail'
    if p==s:return 'common_head'
    if q==s and p==r:return 'swapped'
    if q==s or p==r:return 'chain'
    assert len({q,p,r,s})==4
    return 'disjoint'


def small_combinatorial_checks():
    # Independent exhaustive subset recursion only for k<=5; all-k proofs
    # are written in the report.  This never controls a geometric exclusion.
    patterns=[(False,False,False,False),(True,False,False,True),
              (True,False,False,False),(False,False,True,False),(False,True,False,False)]
    out=[]
    for pattern in patterns:
        allowed=dict(zip(PATTERN_NAMES,pattern));allowed['same_tail']=False
        for k in range(2,6):
            arcs=[(q,p)for q in range(k)for p in range(k)if q!=p]
            masks=[]
            for i,x in enumerate(arcs):
                masks.append(sum(1<<j for j,y in enumerate(arcs) if i!=j and allowed[pair_kind(x,y)]))
            best=0
            def visit(candidates,size):
                nonlocal best
                if size+candidates.bit_count()<=best:return
                if not candidates:best=max(best,size);return
                while candidates:
                    bit=candidates & -candidates;i=bit.bit_length()-1;candidates^=bit
                    visit(candidates & masks[i],size+1)
                best=max(best,size)
            visit((1<<len(arcs))-1,0)
            expected,kind=cap(pattern,k);assert best==expected
            out.append(dict(pattern=list(pattern),A1_count=k,maximum=best,proof_kind=kind))
    return out


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
        if T['finite_extra']!=['A1','A1']:continue
        counts['two_A1_internal_occurrences']+=1
        key0=(B['provenance'],B['deleted']);assert key0 not in seen;seen.add(key0)
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism'])
        ptype=profiles[B['base_type']];typ=I['base_types'][B['base_type']];core=ptype['core']
        assert core==typ['core']==T['base_core'] and T['opposite_base_facet']==q
        assert ptype['simplex_dimension']==typ['simplex_dimension']
        extra=set(range(r['nodes']))-set(mp)-{b};assert q in extra and len(extra)==2
        p=next(iter(extra-{q}));assert p!=q and b not in extra
        graph=V.raw_graph(r);cg=V.core_graph(I['cores'][core])
        assert sorted(mp.values())==list(range(len(cg)))
        assert all(graph[i][j]==cg[mp[i]][mp[j]] for i,j in combinations(mp,2))
        assert all(graph[x][i]==V.RIGHT for x in (q,p)for i in set(range(r['nodes']))-{b,x})
        assert graph[b][p]==V.edge(3) and graph[b][q][0] in (0,5,10)
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
        k=Counter(typ['spectator_forest'])['A1'];assert k>=2
        rank=sum(map(V.rank,typ['spectator_forest']))
        assert len(cg)+rank==14 and r['nodes']+rank-2==15
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
        h=(3-tau)/2;assert norm==h and h+tau/2==K(F(3,2))
        tau_key=tuple(tau.out())
        if tau_key not in tau_checks:
            pair_tests=dict(disjoint=pair_gate(h,0,tau),common_head=pair_gate(h+F(1,2),0,tau),
                            chain=pair_gate(h,F(1,2),tau),swapped=pair_gate(h,1,tau))
            pattern=tuple(pair_tests[x]['allowed_by_necessary_tests']for x in PATTERN_NAMES)
            tau_checks[tau_key]=dict(tau=tau.out(),h=h.out(),pairs=pair_tests,pattern=pattern)
        tests=tau_checks[tau_key];single,proofkind=cap(tests['pattern'],k)
        orbit=set()
        for perm in C['aut']:
            moved=[None]*len(cg)
            for i in range(len(cg)):moved[perm[i]]=contacts[i]
            orbit.add((tuple(moved),tuple(sorted(perm[i]for i in canonical))))
        assert (contacts,canonical)in orbit and len(C['aut'])%len(orbit)==0
        o=dict(base_index=B['index'],base_type=B['base_type'],provenance=B['provenance'],face=B['face'],
            simplex_face=ptype['simplex_face'],core=core,canonical_prism_core=ci,source_core=si,
            deleted_normal=b,coface_root=q,other_A1_root=p,core_isomorphism=sorted(map(list,mp.items())),
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
        columns.append(dict(index=len(columns),kind='two_isolated_A1_cap_profile',
            base_type=r['base_type'],simplex_face=p['simplex_face'],simplex_dimension=p['simplex_dimension'],core=r['core'],
            A1_count=r['A1_count'],tau=r['tau'],orientation_orbit=r['orientation_orbit'],
            orientation_count=r['orientation_count'],single_orientation_capacity=r['single_orientation_capacity'],
            pair_pattern=r['pair_pattern'],proof_kind=r['proof_kind'],capacity=capacity,
            source_bases=sorted(z['base_index']for z in rows),
            positive_face_coefficients=sorted(map(list,co.items())),provenance_coefficients=sorted(map(list,pr.items())),
            load=str(load),initial_simplex_mass=str(mass),slack=str(capacity*mass-load)))
    violations=[c['index']for c in columns if F(c['slack'])<0]
    assert counts['two_A1_internal_occurrences']==len(occurrences)==555
    assert len(cores)==25 and len(tau_checks)==11
    counts.update(columns=len(columns),positive_columns=sum(F(c['load'])>0 for c in columns),
                  violations=len(violations),cores=len(cores),tau_values=len(tau_checks))
    paths=[Path(__file__),Path(V.__file__)]+[ROOT/'results'/f'{t}.json'for t in TAGS.values()]
    paths+=[ROOT/'results/fresh_prism_simplex_shared_capacity_20261017.md',
            ROOT/'results/fresh_sqrt5_cross_coface_cap_capacity_20261017_review.md']
    return dict(status='PASS_EXACT_TWO_A1_MULTISPECTATOR_CAP_CAPACITIES',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(counts),pair_pattern_order=PATTERN_NAMES,
        tau_pair_tests=[tau_checks[k]for k in sorted(tau_checks)],
        core_checks=[dict(core=c,nodes=len(z['H']),gram=[[x.out()for x in row]for row in z['H']],
                         gauge_scales=None if z['scales']is None else[x.out()for x in z['scales']],
                         inverse=[[x.out()for x in row]for row in z['J']],determinant=z['det'].out(),
                         spherical_deletions=z['pd_roots'],automorphisms=z['aut'])for c,z in sorted(cores.items())],
        combinatorial_small_checks=small_combinatorial_checks(),
        occurrences=occurrences,columns=columns,violating_columns=violations,
        orientation_histogram=sorted(Counter(c['orientation_count']for c in columns).items()),
        capacity_histogram=sorted(Counter(c['capacity']for c in columns).items()),
        occurrence_pattern_histogram=[dict(pattern=list(k),count=v)for k,v in sorted(Counter(tuple(o['pair_pattern'])for o in occurrences).items())],
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)]for i,a in sorted(required.items())],
        input_bound=M['bound'],new_repaired_bound=None,optimizer_invocations=0,
        scope='Selected actual dotted-prism bases with two internal isolated A1 components after base deletion, one simple3 contact and one distinct finite5/10 or dotted contact. All555 raw maps, projections, core automorphisms and profile multiplicities are reconstructed. Compatibility tests use only necessary acuteness and arithmetic two-cycle integrality, and all capacities include full colored-core orientation factors. Other spectator domains are omitted safely. No new optimizer, repair, global optimum or dimension bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],'orientations',out['orientation_histogram'],flush=True)
