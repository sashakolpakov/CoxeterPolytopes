#!/usr/bin/env python3
"""Exact prism/base-simplex incidence separation of one saved combined dual."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_simplex_quad_combined_20261017 as P
import fresh_sqrt5_rank14_core_spectator_domain_20261017 as I
import fresh_sqrt5_rank13_packing_quad_separation_20261017 as Q

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_prism_simplex_base_incidence_20261017'
LP='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
DISJOINT='fresh_simplex_disjoint_ports_20261017'


def run():
    lp=P.read(LP);faces=lp['faces']
    dual=next(c for c in lp['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    alpha={i:F(q) for i,q in dual['face_multipliers']}
    star_to_face={f['provenance'][0]['star']:f['index'] for f in faces[:10185]}
    assert len(star_to_face)==10185
    del lp
    old=P.read(I.TAG);sep=P.read(Q.TAG);source=P.read(Q.CORES)
    aggregation=P.read(Q.Q.TAG)['dotted_prism_census']['records']
    spectator=P.read(Q.SPECTATORS)['contexts']
    models=dict(B.C.models());ranks={name:len(g) for name,g in models.items()}
    finite={B.C.code(g):name for name,g in models.items()}
    def component_name(g):
        if sum(map(len,g.values()))!=2*len(g)-2:return None
        if any(m[0]==0 for row in g.values() for m in row.values()):return None
        h={i:{j:m[0] for j,m in row.items()} for i,row in g.items()}
        return finite.get(B.C.code(h))
    def restrict(g,ids):
        ids=set(ids)
        return {i:{j:m for j,m in g[i].items() if j in ids} for i in sorted(ids)}
    def spherical_deletions(g):
        answer=[]
        for p in g:
            pieces=B.C.pieces(restrict(g,set(g)-{p}))
            if all(component_name(h) is not None for h in pieces):answer.append(p)
        return answer
    graphs=[];colors=[];buckets=defaultdict(list);cores=[]
    def register(g,expected=None):
        ids=sorted(g);local={i:j for j,i in enumerate(ids)}
        h={local[i]:{local[j]:m for j,m in row.items()} for i,row in g.items()}
        cg,key=I.refinement(h)
        for ci in buckets[key]:
            mp=I.isomorphism(h,graphs[ci],cg,colors[ci])
            if mp is not None:
                if expected is not None:assert sorted(mp[local[p]] for p in expected)==cores[ci]['spherical_deletions']
                return ci,{p:mp[local[p]] for p in ids}
        ci=len(cores);S=spherical_deletions(h)
        if expected is not None:assert sorted(local[p] for p in expected)==S
        assert len(S)>=2
        graphs.append(h);colors.append(cg);buckets[key].append(ci)
        cores.append(dict(index=ci,nodes=len(h),edges=[[i,j,list(m)] for i,row in h.items() for j,m in sorted(row.items()) if i<j],
                          spherical_deletions=S,simplex_dimension=len(S)-1,old_catalogue_core=ci if ci<310 else None))
        return ci,local
    for c in old['cores']:
        g={i:{} for i in range(c['nodes'])}
        for i,j,m in c['edges']:g[i][j]=g[j][i]=(m,'0','0')
        ci,mp=register(g,c['spherical_deletions']);assert ci==c['index']
    assert len(cores)==310
    old_stars={}
    for r in old['inherited_star_maps']:
        key=r['core'],tuple(r['old_spectator_forest'])
        assert key not in old_stars;old_stars[key]=r['star']
    del old
    templates=[];base_dimension_counts=Counter();iso_count=0
    for canonical in sep['canonical_cores']:
        si=canonical['source_core'];record=source['records'][si];agg=aggregation[si]
        g,_=Q.graph(record);k=agg['prism_dimension'];factors=agg['maximal_prism_factors']
        assert canonical['required']==comb(k,2)
        if k==2:
            assert sorted(map(len,factors))==[2,2]
            deleted=sorted(set(factors[0]+factors[1]))
        else:
            deleted=agg['dotted_base'];assert len(deleted)==2
        bases=[]
        for p in deleted:
            # Removing a facet normal exposes the star of the OPPOSITE base.
            same=next(side for side in factors if p in side)
            assert len(same)==2
            base_facet=next(i for i in same if i!=p)
            expected=sorted(next(i for i in e['omitted'] if i!=p)
                            for e in agg['positive_definite_two_deletions'] if p in e['omitted'])
            assert len(expected)>=k and base_facet not in expected
            h=restrict(g,set(g)-{p});pieces=B.C.pieces(h)
            nonfinite=[piece for piece in pieces if component_name(piece) is None]
            assert len(nonfinite)==1
            core=nonfinite[0];assert set(expected)<=set(core)
            finite_extra=sorted(component_name(piece) for piece in pieces if piece is not core)
            assert None not in finite_extra
            ci,mp=register(core,expected);d=cores[ci]['simplex_dimension']
            assert d>=k-1
            # Checking every one-root deletion in the entire reduced core is
            # cheap and independently validates the inherited PD incidence map.
            actual=spherical_deletions(h);assert actual==expected
            bases.append(dict(deleted_normal=p,opposite_base_facet=base_facet,
                base_core=ci,finite_extra=finite_extra,simplex_dimension=d,
                full_spherical_deletions=expected,core_isomorphism=[[i,mp[i]] for i in sorted(mp)]))
            base_dimension_counts[k,d]+=1
        templates.append(dict(canonical_core=canonical['index'],source_core=si,prism_dimension=k,bases=bases))
        if len(templates)%300==0:print('PRISM_BASE_TEMPLATES',len(templates),len(cores),flush=True)
    assert len(templates)==1335
    types=[];type_lookup={};contexts=[];context_map={};counts=Counter()
    for fi,face in enumerate(faces):
        if face['kind']!='canonical_dotted_quad_aggregate':continue
        for ctx_index,ci in face['source_contexts']:
            assert ctx_index not in context_map
            ctx=spectator[ctx_index];template=templates[ci];sp=ctx['spectator_forest']
            assert ci==sep['source_core_maps'][ctx['core']]['canonical_core']
            assert sum(ranks[n] for n in sp)+sep['canonical_cores'][ci]['nodes']==15
            cc=Counter()
            for base in template['bases']:
                core=cores[base['base_core']];forest=tuple(sorted(sp+base['finite_extra']))
                assert core['nodes']+sum(ranks[n] for n in forest)==14
                fixed=core['nodes']%2==0 and all(ranks[n]%2==0 for n in forest)
                assert not fixed # parity inherited from the flexible parent star
                key=base['base_core'],forest
                if key not in type_lookup:
                    j=len(types);type_lookup[key]=j
                    si=old_stars.get(key);known_face=star_to_face.get(si)
                    types.append(dict(index=j,core=key[0],spectator_forest=list(forest),
                        simplex_dimension=core['simplex_dimension'],old_star=si,old_face=known_face,
                        catalogued=si is not None,represented_in_lp=known_face is not None))
                cc[type_lookup[key]]+=1
            assert sum(cc.values())==(4 if template['prism_dimension']==2 else 2)
            rec=[ctx_index,fi,ci,template['prism_dimension'],[list(x) for x in sorted(cc.items())]]
            context_map[ctx_index]=rec;contexts.append(rec)
    assert len(contexts)==81183
    loads=defaultdict(F);terms=defaultdict(list);upper_loads=defaultdict(F)
    row_diagnostics=[];missing_total_lower=F(0);all_missing_lower=F(0)
    for fi,a in sorted(alpha.items()):
        face=faces[fi]
        if face['kind']!='canonical_dotted_quad_aggregate':continue
        records=[context_map[ctx] for ctx,ci in face['source_contexts']]
        ks={r[3] for r in records};assert len(ks)==1;k=next(iter(ks))
        alternatives=[dict(r[4]) for r in records]
        keys=set().union(*(set(x) for x in alternatives));lower=[];upper=[]
        for t in sorted(keys):
            lo=min(x.get(t,0) for x in alternatives);hi=max(x.get(t,0) for x in alternatives)
            if lo:
                loads[t,k]+=lo*a;terms[t,k].append([fi,lo,str(a)]);lower.append([t,lo])
            upper_loads[t,k]+=hi*a;upper.append([t,hi])
        omitted=min(sum(m for t,m in x.items() if not types[t]['represented_in_lp']) for x in alternatives)
        missing_total_lower+=omitted*a
        all_missing_lower+=min(sum(m for t,m in x.items() if not types[t]['catalogued']) for x in alternatives)*a
        invariant=all(alternatives[0]==x for x in alternatives[1:])
        counts['positive_quad_rows']+=1;counts['positive_quad_provenances']+=len(records)
        counts['invariant_positive_rows' if invariant else 'varying_positive_rows']+=1
        row_diagnostics.append(dict(face=fi,prism_dimension=k,alpha=str(a),provenances=len(records),
            all_base_coefficients_identical=invariant,coefficientwise_minima=lower,
            coefficientwise_maxima=upper,minimum_unrepresented_bases=omitted))
    inequalities=[];missing=[];violations=[]
    for (t,k),load in sorted(loads.items()):
        typ=types[t];d=typ['simplex_dimension'];assert d>=k-1
        capacity=comb(d+1,k)*(13-d);assert capacity>0
        rec=dict(base_type=t,prism_dimension=k,simplex_dimension=d,capacity_per_simplex=capacity,
                 robust_prism_base_load=str(load),upper_prism_base_load=str(upper_loads[t,k]),
                 quad_terms=terms[t,k],old_star=typ['old_star'],old_face=typ['old_face'])
        if typ['represented_in_lp']:
            mass=alpha.get(typ['old_face'],F(0));rhs=capacity*mass
            rec.update(simplex_mass=str(mass),capacity_mass=str(rhs),excess=str(load-rhs),violated=load>rhs)
            inequalities.append(rec)
            if load>rhs:violations.append(rec)
        else:
            rec.update(status='MISSING_BASE_SIMPLEX_VARIABLE_NO_ZERO_MASS_ASSUMPTION')
            missing.append(rec)
    violations.sort(key=lambda x:(-F(x['excess']),x['base_type'],x['prism_dimension']))
    counts.update(canonical_prism_templates=len(templates),base_core_types=len(cores),
        inherited_base_core_types=310,new_base_core_types=len(cores)-310,flexible_quad_provenances=len(contexts),
        full_base_types=len(types),base_types_in_old_catalogue=sum(t['catalogued'] for t in types),
        base_types_represented_in_lp=sum(t['represented_in_lp'] for t in types),
        missing_base_types=sum(not t['catalogued'] for t in types),
        missing_edge_base_types=sum(t['simplex_dimension']==1 and not t['catalogued'] for t in types),
        checked_represented_inequalities=len(inequalities),robust_violations=len(violations),
        positive_missing_base_requirements=len(missing))
    witness=None
    if violations:
        v=min(violations,key=lambda r:(len(r['quad_terms']),r['base_type'],r['prism_dimension']))
        typ=types[v['base_type']]
        witness=dict(inequality=v,base_type=typ,base_core=cores[typ['core']],
                     quad_rows=[dict(face=fi,source_contexts=faces[fi]['source_contexts']) for fi,m,a in v['quad_terms']])
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{I.TAG}.py',f'scripts/{Q.TAG}.py',
        f'results/{LP}.json',f'results/{I.TAG}.json',f'results/{Q.TAG}.json',f'results/{Q.CORES}.json',
        f'results/{Q.SPECTATORS}.json',f'results/{Q.Q.TAG}.json',f'results/{DISJOINT}.json',
        f'results/{DISJOINT}.md']
    return dict(status='PASS_EXACT_PRISM_BASE_INCIDENCE_DUAL_SEPARATION',sources={p:P.sha(p) for p in paths},
        dimension=13,counts=dict(counts),base_dimension_counts=[[k,d,n] for (k,d),n in sorted(base_dimension_counts.items())],
        cores=cores,prism_templates=templates,base_types=types,
        context_schema=['source_context','combined_face','canonical_prism_core','prism_dimension','base_type_multiplicities'],
        contexts=contexts,positive_quad_rows=row_diagnostics,represented_inequalities=inequalities,
        violations=violations,missing_base_requirements=missing,minimal_violation=witness,
        forced_unrepresented_base_incidence_lower=str(missing_total_lower),
        forced_uncatalogued_base_incidence_lower=str(all_missing_lower),
        theorem='For a maximal d-simplex S and k>=2, each(k-1)-simplex subface has n-d k-cofaces outside S, since n-k+1 total minus d-k+1 inside. A prism is not a simplex face. Summing over C(d+1,k) subfaces yields sum_Q m(Q,S,k) f_Q <= C(d+1,k)(n-d) f_S. For k>=3 use the two simplex bases; for k=2 count all four edges invariantly. Deleting a base-facet normal gives the full star of the OPPOSITE base. Its complete PD one-deletion set gives the containing unique maximal simplex, allowing dimension1 for maximal edges.',
        provenance_rule='Each merged quad row retains every source context. For each base type separately, use the minimum base multiplicity across all its represented provenances. Therefore a reported excess persists under every nonnegative split of that row mass. Missing base types are recorded with unknown mass and never set to zero.',
        scope='No optimizer or new face enumeration. This exactly separates one saved dual using forced base stars of every81183 canonical flexible dotted prism provenance. It does not establish a feasible successor or a dimension bound, nor simultaneous compact ambient realization. Observed simplex dimensions are computed, not assumed bounded by4.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(B.dump(d)+'\n')
    print(d['status'],d['counts'],flush=True)
