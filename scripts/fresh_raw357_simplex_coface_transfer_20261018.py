#!/usr/bin/env python3
"""Complete raw357 typed coface incidence; bounded existing-row export.

The report proves injectivity into actual (maximal simplex, A2 facet) pairs.
No native pair census, missing-simplex matrix patch or optimizer is run.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw357_simplex_coface_transfer_20261018'
INC='fresh_raw357_extra_vertex_incidence_pilot_20261018'
BASE='fresh_sqrt5_prism_simplex_base_incidence_20261017'
FORCED='fresh_sqrt5_forced_base_simplex_domain_20261017'
MODEL='fresh_all_marker_families_dual_continuation_20261018'
BARRIER='fresh_raw357_inherited_charge_barrier_20261018'


def run():
    tags=(INC,BASE,FORCED,MODEL);docs={t:json.loads((ROOT/f'results/{t}.json').read_text())for t in tags}
    paths={str(Path(__file__).relative_to(ROOT))}
    for t,d in docs.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(t,p)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{t}.{ext}'
            if(ROOT/p).exists():paths.add(p)
    for p in('results/fresh_all_marker_families_full_model_20261018.md',
             'results/fresh_prism_simplex_shared_capacity_20261017.md'):
        paths.add(p)
    inc,base,forced,model=(docs[t]for t in tags);geom=inc['geometry']
    assert(geom['raw_family'],geom['canonical_core'],geom['required_quads'])==(357,318,1)
    g=geom['native_gram'];right=[2,'0','0'];assert len(g)==10
    assert geom['native_automorphisms']==[list(range(10))]
    S=geom['compact_simplex_faces'];assert len(S)==1
    S=S[0];assert(S['excluded_wall'],S['fixed'],S['sides'],S['dimension'])==(2,[0,7,8,9],[1,3,4,5,6],4)
    tpl=next(t for t in base['prism_templates']if t['canonical_core']==318)
    assert tpl['source_core']==357 and tpl['prism_dimension']==2
    chosen=[b for b in tpl['bases']if b['simplex_dimension']==4]
    assert len(chosen)==1;chosen=chosen[0]
    assert(chosen['deleted_normal'],chosen['opposite_base_facet'],chosen['base_core'])==(2,0,1)
    assert chosen['finite_extra']==['A2','A2']and chosen['full_spherical_deletions']==S['sides']
    mp=dict(chosen['core_isomorphism']);assert set(mp)==set(S['sides'])
    core=base['cores'][1];assert core['nodes']==5 and core['spherical_deletions']==list(range(5))
    cg=[[right[:]for j in range(5)]for i in range(5)]
    for i in range(5):cg[i][i]=[1,'0','0']
    for i,j,e in core['edges']:cg[i][j]=cg[j][i]=e
    assert all(g[i][j]==cg[mp[i]][mp[j]]for i,j in combinations(mp,2))
    autos=[list(p)for p in permutations(range(5))if all(cg[i][j]==cg[p[i]][p[j]]for i,j in combinations(range(5),2))]
    assert autos==[list(range(5))]
    # Removing cap2 leaves the intrinsic Lorentz core plus two whole A2s.
    a2=[[0,9],[7,8]];deleted=set(range(10))-{2}
    for part in a2:
        assert g[part[0]][part[1]]==[3,'0','0']
        assert all(g[i][j]==right for i in part for j in deleted-set(part))
    assert set(mp)|set(sum(a2,[]))==deleted
    quad=geom['quad_partitions'];assert quad==[[[1,3],[0,2]]]
    edge_omissions=[[1,2],[2,3]]
    assert all(x in geom['counted_vertex_omissions']for x in edge_omissions)
    assert all(set(range(10))-set(x)in[set(t)for t in S['vertex_supports']]for x in edge_omissions)
    edge_core=sorted(mp[i]for i in(1,3));assert edge_core==[0,4]
    # Choose one endpoint intrinsically; there is no factor for its mate.
    endpoint_omitted=min(edge_core);assert endpoint_omitted==0
    byctx={r[0]:r for r in base['contexts']}
    existing=dict(forced['base_to_new_face'])
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']}
    rows=[];columns=[];missing=[];refs=[];seen=set()
    for src in inc['complete_family_inventory']:
        ctx=byctx[src['source_context']]
        assert ctx[1:4]==[src['face'],318,2]
        targets=[(t,m)for t,m in ctx[4]if base['base_types'][t]['core']==1]
        assert len(targets)==1;ti,m=targets[0];assert m==1
        typ=base['base_types'][ti];sp=sorted(src['spectator_forest']+['A2','A2'])
        assert typ['spectator_forest']==sp and typ['simplex_dimension']==4
        assert ti not in seen;seen.add(ti)
        b=Counter(sp)['A2'];assert b in(2,3,4)
        cap=2*b
        fi=typ['old_face']if typ['represented_in_lp']else existing.get(ti)
        assert typ['represented_in_lp']==(typ['old_face']is not None)
        source_mass=alpha.get(src['face'],F());target_mass=alpha.get(fi,F())
        assert source_mass==F(src['current_mass'])
        row=dict(source_face=src['face'],source_context=src['source_context'],source_required_quads=1,
            base_type=ti,base_core=1,passive_forest=sp,whole_A2_count=b,
            base_occurrences_per_quad=1,capacity=cap,existing_simplex_face=fi,
            simplex_required_triangle_payment=20,source_mass=str(source_mass),
            current_simplex_mass=str(target_mass)if fi is not None else None,
            coefficient_interpretation='one source quad per assigned actual simplex; at most one per actual A2 facet q',
            model_installed=False)
        if fi is None:
            assert not source_mass
            missing.append(dict(base_type=ti,source_face=src['face'],passive_forest=sp,
                required=20,reason='Target simplex row absent from the current finite model; no signed column installed.'))
        else:
            load=source_mass-cap*target_mass
            c=dict(source_face=src['face'],simplex_face=fi,base_type=ti,
                integer_face_coefficients=[[src['face'],1],[fi,-cap]],
                integer_vertex_capacities=[],vertex_cost=0,capacity=cap,
                current_exact_dual_load=str(load),strict_current_violation=load>0)
            columns.append(c)
            # Compare against the older universal capacity45 constraint.
            inherited=[z for z in forced['new_shared_transfer_columns']if z['simplex_face']==fi]
            for z in inherited:
                assert z['negative_simplex_coefficient']==[fi,-45]
                old=dict(z['positive_quad_coefficients']);assert old[src['face']]==2
                refs.append(dict(source_face=src['face'],simplex_face=fi,source=FORCED,
                    old_source_coefficient=2,old_simplex_capacity=45,
                    new_source_coefficient=1,new_simplex_capacity=cap))
        rows.append(row)
    assert len(rows)==len(seen)==47 and len(columns)==16 and len(missing)==31
    assert len({c['simplex_face']for c in columns})==16
    bad=[c for c in columns if c['strict_current_violation']]
    assert[c['source_face']for c in bad]==[70404,70406,70408,70414,70420,70422]
    assert sum((F(c['current_exact_dual_load'])for c in bad),F())==F(37155891,10**10)
    return dict(status='PASS_RAW357_COMPLETE_TYPED_SIMPLEX_COFACE_TRANSFER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        geometry=dict(source_core=318,raw_family=357,native_gram=g,source_quad=quad,
            deleted_cap=2,coface_facet=0,coface_whole_A2=a2[0],other_whole_A2=a2[1],
            simplex=S,canonical_core=core,core_isomorphism=[list(x)for x in sorted(mp.items())],
            core_automorphisms=autos,distinguished_edge_omissions=edge_core,
            canonical_endpoint_omitted=endpoint_omitted,base_occurrence_multiplicity=1),
        family=rows,existing_target_transfer_columns=columns,missing_simplex_rows=missing,
        inherited_capacity_comparisons=refs,
        counts=dict(source_classes=47,actual_simplex_types=47,current_target_rows=16,
            missing_target_rows=31,current_strict_violations=len(bad)),
        total_positive_load=str(sum((F(c['current_exact_dual_load'])for c in bad),F())),
        optimizer_invocations=0,full_matrix_reconstructed=False,model_changed=False,
        scope='Complete47-type actual coface theorem and16 valid existing-row signed transfers, not installed.31 absent simplex rows remain unpatched. Exact current-dual violations only; no repaired certificate, optimum change, chamber or dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists(),path
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['total_positive_load'],flush=True)
