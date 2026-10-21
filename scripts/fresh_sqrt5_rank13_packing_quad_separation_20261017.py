#!/usr/bin/env python3
"""Exact generic-rule separation on canonical dotted full-star quad aggregates."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank14_core_spectator_domain_20261017 as I
import fresh_quad_product_aggregation_20261017 as Q

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_packing_quad_separation_20261017'
RULE='fresh_sqrt5_rank13_simplex_packing_repair_20261017'
CORES='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017'
SPECTATORS='fresh_sqrt5_rank13_dotted_b1_spectators_20261017'


def graph(record):
    n=record['nodes'];g={i:{} for i in range(n)};edges={}
    for (i,j),m in zip(combinations(range(n),2),record['row']['labels']):
        edges[i,j]=m
        if m==2:continue
        if m==0:
            a,b=record['parameters'][f'{i},{j}'];weight=(0,str(F(a)),str(F(b)))
        else:weight=(m,'0','0')
        g[i][j]=g[j][i]=weight
    return g,edges


def partition(a,b):
    return tuple(sorted((tuple(sorted(a)),tuple(sorted(b)))))


def run():
    rule=B.read(RULE);source=B.read(CORES);spectators=B.read(SPECTATORS);aggregation=B.read(Q.TAG)
    rawagg=aggregation['dotted_prism_census']['records']
    assert len(source['records'])==len(rawagg)==1399
    cert=next(c for c in rule['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    B.dual_bound=D.sparse_dual_bound;B.verify_certificate(rule,cert)
    weights={i:F(a) for i,a in cert['weights']}
    models=dict(B.C.models());ranks={name:len(g) for name,g in models.items()}
    finite={B.C.code(g):name for name,g in models.items()}
    codes={name:json.loads(B.dump(B.C.code(g))) for name,g in models.items()}
    pair_lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']
        for v in rule['variables'] if v.get('kind')!='maximal_simplex_clique_charge'}
    known_forests={tuple(v['forest']) for v in rule['vertices']}
    canonical=[];graphs=[];colors=[];buckets=defaultdict(list);maps=[]
    for si,record in enumerate(source['records']):
        assert record['index']==rawagg[si]['index']==si
        g,edges=graph(record);assert len(B.C.pieces(g))==1
        cg,key=I.refinement(g);ci=None;mp=None
        for candidate in buckets[key]:
            mp=I.isomorphism(g,graphs[candidate],cg,colors[candidate])
            if mp is not None:ci=candidate;break
        agg=rawagg[si]
        if ci is None:
            ci=len(canonical);occurrences=[];counter=Counter()
            for pos,neg in agg['quadrilateral_black_partitions']:
                assert len(pos)==len(neg)==2
                for i,j in product(pos,neg):
                    pair=[next(a for a in pos if a!=i),next(b for b in neg if b!=j)]
                    ids=[a for a in range(record['nodes']) if a not in(i,j)]
                    pieces=B.C.pieces(B.C.graph(ids,edges));forest=tuple(sorted(finite[B.C.code(p)] for p in pieces))
                    code=B.dump(sorted(B.C.code(p,pair) for p in pieces))
                    counter[forest,code]+=1
                    occurrences.append(dict(vertex_omitted=[i,j],marked_corner_roots=pair,forest=list(forest),marked_forest_code=json.loads(code)))
            assert len(occurrences)==4*agg['quadrilateral_count']==agg['aggregate_corner_count']
            canonical.append(dict(index=ci,source_core=si,nodes=record['nodes'],quotient_rank=record['nodes']-1,
                quotient_even=(record['nodes']-1)%2==0,required=agg['quadrilateral_count'],
                pd_deletions=[x['omitted'] for x in agg['positive_definite_two_deletions']],
                quad_partitions=agg['quadrilateral_black_partitions'],
                aggregate_corner_types=[dict(forest=list(f),marked_forest_code=json.loads(c),multiplicity=m)
                                        for (f,c),m in sorted(counter.items())],
                aggregate_corner_occurrences=occurrences))
            graphs.append(g);colors.append(cg);buckets[key].append(ci);mp={i:i for i in g}
        target=canonical[ci]
        assert {tuple(sorted(mp[i] for i in x['omitted'])) for x in agg['positive_definite_two_deletions']}=={
            tuple(sorted(x)) for x in target['pd_deletions']}
        assert {partition([mp[i] for i in a],[mp[j] for j in b]) for a,b in agg['quadrilateral_black_partitions']}=={
            partition(a,b) for a,b in target['quad_partitions']}
        assert target['required']==agg['quadrilateral_count']
        maps.append(dict(source_core=si,canonical_core=ci,isomorphism=[mp[i] for i in range(len(g))]))
    assert len(maps)==1399
    seen=set();representatives=[];unpaid=[];missing=set();counts=Counter();bycore=defaultdict(Counter)
    trace=hashlib.sha256();witnesses=[];forest_seen=set();flexible_forests=set()
    for ctx in spectators['contexts']:
        ci=maps[ctx['core']]['canonical_core'];core=canonical[ci];sp=tuple(ctx['spectator_forest'])
        assert sum(ranks[n] for n in sp)+core['nodes']==15
        fixed=core['quotient_even'] and all(ranks[n]%2==0 for n in sp)
        counts['raw_fixed' if fixed else 'raw_flexible']+=1
        key=ci,sp
        if key in seen:continue
        seen.add(key);representatives.append(ctx['index'])
        counts['canonical_fixed' if fixed else 'canonical_flexible']+=1
        bycore[ci]['fixed' if fixed else 'flexible']+=1
        if fixed:continue
        coeff=Counter();absent=set()
        for corner in core['aggregate_corner_types']:
            forest=tuple(sorted(corner['forest']+list(sp)))
            assert sum(ranks[n] for n in forest)==13
            forest_seen.add(forest);flexible_forests.add(forest)
            code=B.dump(sorted(corner['marked_forest_code']+[codes[n] for n in sp]))
            if forest in known_forests:
                assert (forest,code) in pair_lookup
                coeff[pair_lookup[forest,code]]+=corner['multiplicity']
            else:absent.add(forest)
        missing.update(absent)
        total=sum((m*weights.get(j,F(0)) for j,m in coeff.items()),F(0))
        required=core['required'];passed=total>=required
        counts['flexible_paid' if passed else 'flexible_unpaid']+=1
        bycore[ci]['paid' if passed else 'unpaid']+=1
        if absent:counts['flexible_using_zero_extended_forests']+=1
        trace.update((B.dump([ctx['index'],ci,sp,required,str(total),sorted(coeff.items()),sorted(absent)])+'\n').encode())
        if not passed:
            unpaid.append([ctx['index'],str(total),required])
            witnesses.append((core['nodes'],len(sp),ctx['index'],dict(source_context=ctx['index'],canonical_core=ci,
                source_core=core['source_core'],spectator_forest=list(sp),required=required,payment=str(total),
                pair_coefficients=[list(x) for x in sorted(coeff.items())],simplex_packing_payment='0',
                zero_extended_forests=[list(x) for x in sorted(absent)])))
    assert counts['raw_fixed']==30369 and counts['raw_flexible']==89329
    assert counts['flexible_paid']+counts['flexible_unpaid']==counts['canonical_flexible']
    def attach_witness(rows):
        if not rows:return None
        w=min(rows,key=lambda x:x[:3])[3]
        return dict(w,core_record=source['records'][w['source_core']],
                    aggregate_corner_occurrences=canonical[w['canonical_core']]['aggregate_corner_occurrences'])
    represented=[r for r in witnesses if not r[3]['zero_extended_forests']]
    def only_ade_h(r):
        w=r[3];names=list(w['spectator_forest'])
        for corner in canonical[w['canonical_core']]['aggregate_corner_types']:names+=corner['forest']
        return all(not(name.startswith('B') or name in ('F4','I2(6)','I2(10)')) for name in names)
    represented_ade_h=[r for r in represented if only_ade_h(r)]
    counts['unpaid_all_forests_represented']=len(represented)
    witness=attach_witness(witnesses)
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{I.TAG}.py',f'scripts/{Q.TAG}.py',
        f'results/{RULE}.json',f'results/{CORES}.json',f'results/{SPECTATORS}.json',f'results/{Q.TAG}.json']
    return dict(status='PASS_EXACT_FLEXIBLE_FULL_STAR_QUAD_SEPARATION',sources={p:I.sha(p) for p in paths},
        counts=dict(counts,source_cores=1399,canonical_connected_cores=len(canonical),canonical_full_stars=len(seen),
            flexible_corner_forests=len(flexible_forests),zero_extended_forests_seen=len(missing)),
        determinant_rule='The connected corank1 core has nonsingular quotient rank n-1. Under any K-valued gauge with the same normalized Gram it changes by one common totally-positive scale and K-diagonal congruence. Its quotient determinant squareclass is invariant when n-1 is even. Each finite spectator contributes its actual rank. Thus the full ambient rank14 squareclass is fixed exactly in the conservative parity test: n odd and every spectator component even; otherwise at least one odd-rank connected factor leaves the determinant class flexible.',
        aggregation_rule='All PD two-root deletions and all K2,2 partitions are inherited from the exact full-star census and checked under every weighted graph isomorphism. Both deleted roots belong to the unique Lorentz core: deleting at most one core root cannot leave a PD block of rank n-1. Thus finite spectators do not create new deletion edges or quadrilaterals. Every listed K2,2 contributes exactly its four literal corner pairs once. Simplex packing charges contribute zero to these quad aggregates.',
        canonical_cores=canonical,source_core_maps=maps,representative_context_indices=representatives,
        evaluated_flexible_trace_sha256=trace.hexdigest(),unpaid_contexts=unpaid,
        core_counts=[dict(canonical_core=i,**dict(c)) for i,c in sorted(bycore.items())],
        zero_extended_forests=[list(x) for x in sorted(missing)],minimal_witness=witness,
        minimal_represented_witness=attach_witness(represented),
        minimal_ade_h_represented_witness=attach_witness(represented_ade_h),
        scope='No optimizer or new face census. This evaluates only the flexible canonical full-star aggregates of the1399 certified dotted cores and their frozen finite spectator contexts. Fixed-determinant aggregates are counted but not evaluated. Source context indices reference the existing spectator package rather than duplicating it. An unpaid aggregate refutes coverage by this particular generic rule; no compact ambient realization or dimension conclusion is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
