#!/usr/bin/env python3
"""Order15 B/D totals92: forced B completion and typed whole D coupons.

The complete geometric censuses and residual/no-bridge theorem are inherited.
This reconstructs every eligible B map and every affected forced completion,
and rechecks all20 sealed D variants with the narrower upgraded coupon roles.
No finite forest bound or ambient dimension improvement is asserted.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_order15_affine_d_half_coupon_transfer_20261017 as H
import fresh_order15_affine_b_residual_terminal_transfer_20261017 as B

D=H.D; X=H.X; P=H.P; Q=H.Q; L=H.L; ROOT=H.ROOT
TAG='fresh_order15_bd_terminal_coupon_transfer_20261018'
CURRENT='fresh_order11_order16_dimension60_20261018'
DEPTH='fresh_depth10_finite7_zero_transfer_20261018'
BUDGET='fresh_rank61_complete_finite420_budget_20261018'
BDROP={(9,13),(9,14)}
EXPECTED=[366,369,373,383,386,395,398,535,540,611,662,665,729,766,769,
          835,986,989,994,997,1001,1031,1037,1051,1118,1126,1159,1164,
          1177,1503,1614,1619,1642,1648,1668,1796,1798,1801,1834,1907,
          1912,1935,1966,2189,2191,2194,2227,2300,2305,2328,2359,2423,2431]
STEMS=tuple(dict.fromkeys((CURRENT,DEPTH,BUDGET,B.TAG,H.TAG,H.C15,D.TAG,
    X.TAG,X.SMALL,X.OVER,X.COMBINED,
    'fresh_distance9_supporting_extensions_review_20261017',
    'fresh_order11_bc_zero_transfer_20261018',
    'fresh_order16_d_distance9_half_transfer_20261018',
    'fresh_order17_e6_unit_coupon_zero_transfer_20261018',
    'fresh_finite7_ideal5_triangle_payments_20261017',
    'fresh_finite_rank15_existing_half_transfer_20261018')))


def b_proof(c,core,maps,target):
    g,J,du,dv,ends=X.constraints(c,core);base=c['marked_component'];n=len(g)
    rows=[];control=None
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)}
        missing=sorted(set(range(15))-set(image));entries=[]
        assert all((int(2-g[a][b]) if a!=b and g[a][b] else 0)==target[mp[i]][mp[j]]
                   for i,a in enumerate(base) for j,b in enumerate(base))
        for new in missing:
            neighbors=[(old,node,target[new][old]) for old,node in image.items() if target[new][old]]
            assert len(neighbors)<=1
            for old,i,label in neighbors:
                assert i in c['possible_contacts'] and label in (3,4)
                square=F(1 if label==3 else 2)
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=i,
                    native_label=label,physical_contact_square=str(square),
                    inverse_diagonal=str(J[i][i]),projection_lower_bound=str(square*J[i][i])))
        assert entries
        row=dict(map=mp,missing_target_nodes=missing,forced_entries=entries)
        bad=next((r for r in entries if F(r['projection_lower_bound'])>2),None)
        if bad:
            row.update(outcome='forced_entry_projection_exclusion',witness=bad);rows.append(row);continue
        if control is None:
            gg,control=D.core_constraints(c,core);assert gg==g
        assert len(entries)==1 and entries[0]['native_label']==3
        assert missing in ([11,12],[10,11,12]) and 12 not in image
        first=entries[0]['new_target_node'];assert first==min(missing)
        assert all(target[t][t+1]==(4 if t==11 else 3) for t in missing[:-1])
        exterior=None
        if c['method']=='relevant_component':
            assert control['joint_component_order_bound']==14<15
            exterior=dict(actual_factor_order=15,whole_order_with_exterior_root=14)
        else:assert control['kind']=='at_most_one_core_contact_root'
        # The unique relevant first entry exhausts possible core contacts;
        # the optional exterior root is excluded by the actual factor order.
        # Thus every later tail root is outside the old core and has zero
        # old-core contacts. These zeros are forced before this construction.
        image.update({t:n+k for k,t in enumerate(missing)});nn=n+len(missing)
        h=[r+[F(0)]*len(missing) for r in g]+[[F(0)]*nn for _ in missing]
        for t in missing:
            # Physical roots all have norm2. Only the NEW terminal12 is
            # rescaled by1/sqrt2: norm1, native4 bond -1, all other bonds -1.
            h[image[t]][image[t]]=F(1 if t==12 else 2)
            for old in range(15):
                if target[t][old]:
                    assert target[t][old]==(4 if {t,old}=={11,12} else 3)
                    h[image[t]][image[old]]=h[image[old]][image[t]]=F(-1)
        ids=[image[t] for t in range(15)]
        affine=L.sub(h,ids);assert Q.inertia(affine)==[14,0,1]
        kernel=Q.nullspace(affine);assert len(kernel)==1
        z=kernel[0];z=[v if z[0]>0 else -v for v in z];assert all(v>0 for v in z)
        sig=Q.inertia(h);assert sig in ([16,2,0],[17,2,0])
        row.update(outcome='forced_completed_index_two',residual_contact_control=control,
            exterior_root_order_exclusion=exterior,terminal_congruence=dict(target_node=12,
                physical_norm='2',scale='1/sqrt2',scaled_norm='1',scaled_native4_bond='-1'),
            target_to_literal=ids,gram=L.encode(h),signature=sig,
            positive_affine_kernel=list(map(str,z)))
        rows.append(row)
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
        end_projection_proof=ends,maps=rows)


def path_role(g,role):
    ns=role['small_nodes'];mark=role['small_mark']
    assert role['small_type'] in ('A4','A5') and len(ns)==int(role['small_type'][1:])
    adj={i:[j for j in ns if i!=j and g[i][j]] for i in ns}
    assert sum(map(len,adj.values()))==2*(len(ns)-1)
    assert all(g[i][j]==-1 for i in ns for j in adj[i])
    leaves=[i for i in ns if len(adj[i])==1];assert len(leaves)==2
    distances={mark:0};queue=deque([mark])
    while queue:
        i=queue.popleft()
        for j in adj[i]:
            if j not in distances:distances[j]=distances[i]+1;queue.append(j)
    assert set(distances)==set(ns)
    lengths=sorted(distances[i] for i in leaves)
    assert lengths==([0,3] if len(ns)==4 else [1,3])
    ds=role['large_nodes'];dm=role['large_mark'];assert role['large_type']=='D11'
    da={i:[j for j in ds if i!=j and g[i][j]] for i in ds}
    forks=[i for i in ds if len(da[i])==3];assert len(forks)==1
    eligible=sorted(i for i in da[forks[0]] if len(da[i])==1)
    assert len(eligible)==2 and dm in eligible and role['whole_components']
    assert set(ns).isdisjoint(ds) and all(not g[i][j] for i in ns for j in ds)
    return dict(small_type=role['small_type'],small_nodes=ns,small_mark=mark,
        small_role='endpoint' if len(ns)==4 else 'second_from_endpoint',
        distances_to_endpoints=lengths,large_type='D11',large_nodes=ds,large_mark=dm,
        large_role='short_fork_leaf',large_eligible_marks=eligible,whole_components=True)


def d_units(doc,over):
    rows=[];roles=Counter()
    for old in doc['cases']:
        if old['index'] not in H.FULL:continue
        c=over['isolation']['cases'][old['index']];core=over['isolation']['cores'][c['core_key']]
        proof=H.full_case(c,core)
        assert proof=={k:v for k,v in old['proof'].items() if k!='maps'}
        for v in proof['variants']:
            g=[[F(x) for x in row] for row in v['gram']];corners=[];selected=None;other=None;role=None
            for corner in v['corners']:
                a,b=F(corner['used_ideal_payment']),F(corner['used_finite_coupon'])
                assert a in (0,F(1,2)) and b in (0,F(1,2))
                if a:selected=corner['marks'];assert not b
                if b:
                    assert len(corner['coupon_roles'])==1
                    role=path_role(g,corner['coupon_roles'][0]);other=corner['marks']
                    roles[role['small_type']]+=1
                corners.append(dict(marks=corner['marks'],ideal_payment='0',
                    finite_coupon=str(int(bool(b))),total=str(int(bool(b)))))
            assert role is not None and selected!=other and sum(F(r['total']) for r in corners)==1
            assert v['used_same_component_finite_units'] is False
            rows.append(dict(case=old['index'],optional_core_contact=v['optional_core_contact'],
                original_ideal_marks=selected,distinct_unit_marks=other,upgraded_role=role,
                corners=corners,total='1'))
    assert len(rows)==20 and roles==dict(A4=12,A5=8)
    return rows,dict(roles)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[CURRENT]['dimension_bound']==60 and docs[DEPTH]['feature_depth']==10
    oldcosts=docs[H.TAG]['simultaneous_order15_costs'];assert oldcosts==dict(A=75,B=94,C=91,D=94)
    small=docs[X.SMALL];over=docs[X.OVER];rec=next(r for r in small['records'] if r['kind']=='B' and r['order']==15)
    target=rec['target'];assert target==docs[B.TAG]['target']
    paid=set(map(tuple,docs[B.TAG]['retained_unit_pairs']));assert len(paid)==94 and BDROP<=paid
    bweights={p:F(int(p in paid-BDROP)) for p in combinations(range(15),2)}
    dweights={tuple(r['pair']):F(r['weight']) for r in docs[H.TAG]['weighted_pairs']}
    for p in H.HALVES:assert dweights[p]==F(1,2);dweights[p]=F()
    assert sum(bweights.values())==sum(dweights.values())==92
    rules=[]
    for kind,t,changes,weights,old in [('B',target,BDROP,bweights,docs[B.TAG]),
                                     ('D',docs[H.TAG]['target'],H.HALVES,dweights,docs[H.TAG])]:
        _,dist,_=P.distance_data(t);autos=P.embeddings(t,t)
        assert all({tuple(sorted(mp[i] for i in pair)) for pair in changes}==changes for mp in autos)
        assert all(dist[i][j]==10 for i,j in changes)
        assert all(weights[tuple(p)]==1 for p in old['basic_pairs_retained']+old['triangle_pairs_retained'])
        assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=8)
        rules.append(dict(kind=kind,target=t,automorphisms=autos,old_cost='94',new_cost='92',
            zero_pairs=sorted(changes),weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
            retained_basic=old['basic_pairs_retained'],retained_triangle=old['triangle_pairs_retained'],
            retained_short8=[p for p in weights if dist[p[0]][p[1]]<=8]))
    patterns=[p for p in small['patterns'] if p['order_bound']>=15 and len(p['nodes'])<=15]
    assert len(patterns)==536
    saved={r['case']:r for r in rec['case_embeddings']};support=[];affected=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings'] and pairs==list(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        chosen=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in BDROP]
        if chosen:affected.append((p['index'],chosen))
    assert len(support)==len(saved) and [i for i,m in affected]==EXPECTED
    cases=[];outcomes=Counter();signatures=Counter();bounds=Counter()
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        proof=b_proof(c,core,maps,target)
        for r in proof['maps']:
            outcomes[r['outcome']]+=1
            if r['outcome']=='forced_completed_index_two':signatures[str(r['signature'])]+=1
            else:bounds[r['witness']['projection_lower_bound']]+=1
        cases.append(dict(index=index,method=c['method'],core_key=c['core_key'],proof=proof))
    assert outcomes==dict(forced_entry_projection_exclusion=68,forced_completed_index_two=40)
    assert signatures=={'[16, 2, 0]':24,'[17, 2, 0]':16}
    units,rolecounts=d_units(docs[H.TAG],over)
    cs=over['isolation']['cases'];assert len(cs)==2614
    groups=dict(low=[c for c in cs if 0<(c.get('component_order_bound') or 0)<=14],
        protected=[c for c in cs if c.get('component_order_bound') is None],
        too_large=[c for c in cs if (c.get('component_order_bound') or 0)>=15 and len(c['marked_component'])>15],
        eligible=[c for c in cs if (c.get('component_order_bound') or 0)>=15 and len(c['marked_component'])<=15])
    routing={k:len(v) for k,v in groups.items()};assert routing==dict(low=1691,protected=83,too_large=304,eligible=536)
    assert Counter(c['method'] for c in groups['protected'])==dict(finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    coupons=[dict(small_type='whole '+s,small_role=role,small_eligible_marks=marks,
        large_type='whole D11',large_role='short_fork_leaf',large_eligible_mark_count=2,
        old_marked_pair_weight='1/2',new_marked_pair_weight='1',changed_marked_pairs=4,
        component_pair_cost_increment='2')
        for s,role,marks in [('A4','endpoint',[0,3]),('A5','second_from_endpoint',[1,3])]]
    paths=[Path(__file__),Path(H.__file__),Path(B.__file__),Path(D.__file__),Path(X.__file__),
           Path(P.__file__),Path(Q.__file__),Path(L.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
              (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER15_BD92_TYPED_UNIT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        baseline_dimension_bound=60,old_family_costs=oldcosts,new_family_costs=dict(A=75,B=92,C=91,D=92),
        new_order15_maximum='92',new_order15_envelope='92/15',rules=rules,
        complete_high_overlap_routing=routing,complete_eligible_patterns=536,
        B_reconstructed_support=support,B_affected_contexts=53,B_affected_maps=108,B_cases=cases,
        B_outcomes=dict(outcomes),B_signature_counts=dict(signatures),B_projection_bounds=dict(bounds),
        D_inherited_affected_maps=408,D_inherited_proper_exclusions=312,D_inherited_whole_maps=96,
        D_whole_variants=units,D_role_counts=rolecounts,new_cross_component_coupons=coupons,
        uses_finite_same_component_units=False,all_other_finite_weights_unchanged=True,
        finite_forest_bound_proved=False,finite_budget_requires_new_coupon_increments=True,
        compatible_depth10_finite7_zero=True,triangle_minimum=2,quadrilateral_minimum=1,
        actual_corner_priority=['current independent short8 route',
            'changed order15 D whole-factor coupon or inherited proper exclusion',
            'changed order15 B selected completion exclusion',
            'unchanged actual-face inheritance with complete nonselected interfaces'],
        new_dimension_bound=None,
        scope='Complete local weighted refinement, with inherited all-label census and residual/no-bridge inputs. Native norm2 contact bounds precede terminal congruence and completed Gram zeros. Twenty sealed D variants use distinct full typed finite coupons. New +2/+2 finite costs are explicit and not yet globally budgeted; no ambient improvement claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists(),'Refuse to overwrite an existing certificate.'
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['B_outcomes'],out['D_role_counts'],out['new_family_costs'],flush=True)
