#!/usr/bin/env python3
"""Four affine-D14 weights become halves using existing whole A4/A5-D11 coupons."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_d16_partial_weight_transfer_20261017 as D
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=D.ROOT;L=D.L;Q=D.Q;P=D.P
TAG='fresh_order15_affine_d_half_coupon_transfer_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
C15='fresh_order15_affine_c_terminal_projection_transfer_20261017'
HALVES={(1,13),(1,14),(9,11),(9,12)}
FULL=[691,795,825,827,1874,1976,2267,2369,2481,2485,2506,2513]
STEMS=(BASE,C15,D.TAG,X.TAG,D.SMALL,D.OVER,X.COMBINED,
       'fresh_distance9_supporting_extensions_review_20261017',
       'fresh_finite8_a17_d11_coupon_budget_20261017',
       'fresh_a2_d11_full_coupon_forest_budget_20261017')


def proper_case(c,core,maps,target):
    g,J,du,dv,ends=X.constraints(c,core);base=c['marked_component'];n=len(g)
    contact_sums=[]
    for end in c['end_sets']:
        for i,j in combinations(end,2):
            z=J[i][i]+J[j][j]+2*J[i][j];assert z>2
            contact_sums.append(dict(nodes=[i,j],two_unit_projection_lower=str(z)))
    rows=[]
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)}
        missing=sorted(set(range(15))-set(mp));entries=[]
        assert all((int(2-g[a][b]) if a!=b and g[a][b] else 0)==target[mp[i]][mp[j]]
                   for i,a in enumerate(base) for j,b in enumerate(base))
        for new in missing:
            ns=[image[t] for t in image if target[new][t]]
            assert len(ns)<=1
            if ns:
                i=ns[0];assert i in c['possible_contacts']
                entries.append(dict(new_target_node=new,core_node=i,
                    native_label=3,projection_square=str(J[i][i]),residual_square=str(2-J[i][i])))
        assert entries
        row=dict(map=mp,missing_target_nodes=missing,forced_entries=entries)
        bad=next((z for z in entries if F(z['projection_square'])>2),None)
        if bad:
            row.update(outcome='projection_exclusion',witness=bad);rows.append(row);continue
        assert all(F(z['projection_square'])>1 for z in entries)
        if len(entries)>1:
            witnesses=[]
            for a,b in combinations(entries,2):
                i,j=a['core_node'],b['core_node'];v=J[i][j]
                residual=(2-J[i][i])*(2-J[j][j])
                shortcut=min(du[i]+dv[j],du[j]+dv[i])+3
                if v*v>residual and (v>0 or shortcut<=8):
                    witnesses.append(dict(new_target_nodes=[a['new_target_node'],b['new_target_node']],
                        core_nodes=[i,j],projection_product=str(v),residual_product_bound=str(residual),
                        conclusion='positive_product' if v>0 else 'retained_short_path',
                        shortcut=shortcut if v<0 else None))
            assert witnesses
            row.update(outcome='two_entry_exclusion_or_shortcut',witness=witnesses[0]);rows.append(row);continue
        # The exact old core contacts of every remaining missing target root
        # must be fixed before forming this completed Gram.
        gg,control=D.core_constraints(c,core);assert gg==g
        if c['method']=='relevant_component':
            assert control['joint_component_order_bound']<15
            exterior=dict(actual_factor_order=15,
                maximum_whole_order_with_exterior_root=control['joint_component_order_bound'])
        else:exterior=None
        image.update({t:n+i for i,t in enumerate(missing)});nn=n+len(missing)
        h=[r+[F(0)]*len(missing) for r in g]+[[F(0)]*nn for _ in missing]
        for a in missing:
            h[image[a]][image[a]]=F(2)
            for b in range(15):
                if target[a][b]:h[image[a]][image[b]]=h[image[b]][image[a]]=F(-1)
        ids=[image[t] for t in range(15)]
        assert Q.inertia(L.sub(h,ids))==[14,0,1]
        sig=Q.inertia(h);assert sig==[16,2,0]
        row.update(outcome='forced_lorentz_index_exclusion',residual_contact_control=control,
            exterior_root_order_exclusion=exterior,target_to_literal=ids,gram=L.encode(h),signature=sig)
        rows.append(row)
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
        preferred_corner=c['preferred_corner'],end_projection_proof=ends,
        two_contact_projection_tests=contact_sums,maps=rows)


def full_case(c,core):
    g=[[F(x) for x in row] for row in core['gram']];n=len(g)
    assert core['gauge']==[1]*n and all(g[i][i]==2 for i in range(n))
    assert set(x for row in g for x in row)<={F(-1),F(0),F(2)}
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(n))-black);assert len(black)==4
    assert Q.inertia(L.sub(g,white))==[len(white),0,0]
    for b in black:
        ids=sorted(set(white)|{b});assert Q.inertia(L.sub(g,ids))==[len(ids),0,0]
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    affine=c['marked_component'];assert len(affine)==15
    assert Q.inertia(L.sub(g,affine))==[14,0,1]
    ab=next(b for b in c['blocks'] if b['nodes']==affine)
    kernel=list(map(F,ab['positive_kernel']));assert all(x>0 for x in kernel)
    assert all(sum(g[i][j]*kernel[k] for k,j in enumerate(affine))==0 for i in affine)
    contacts=[None]
    if c['method']=='radical':
        assert Q.inertia(g)==[n-2,1,1]
        z=list(map(F,core['radical'][0]));assert any(z)
        assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
        allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
        pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
        assert pos and neg
        du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);shortcuts=[]
        for i,j in product(pos,neg):
            length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
            shortcuts.append(dict(opposite_sign_nodes=[i,j],length=length))
        control=dict(kind='radical_contact_sealing',radical=list(map(str,z)),
            allowed_contacts=allowed,positive_nodes=pos,negative_nodes=neg,forced_shortcuts=shortcuts,
            conclusion='Every additional supporting white is orthogonal to the entire core unless short8 pays.')
    else:
        assert c['method']=='lorentz_residual'
        gg,control=D.core_constraints(c,core);assert gg==g
        assert c['possible_contacts']==c['omitted_black']
        assert len(c['single_contact_candidates'])==2
        assert all(z['projection_square']=='2' and z['residual_square']=='0'
                   and z['maximum_attached_orthogonal_rank']==0 for z in c['single_contact_candidates'])
        contacts += [z['contact_node'] for z in c['single_contact_candidates']]
        control.update(conclusion='At most one simple core-contact root, of residual zero; all later whites are orthogonal to it and the core.')
    variants=[]
    for entry in contacts:
        h=[r[:] for r in g]
        if entry is not None:
            h=[r+[F(0)] for r in h]+[[F(0)]*(n+1)]
            h[n][n]=F(2);h[n][entry]=h[entry][n]=F(-1)
        nn=len(h);ww=sorted(set(range(nn))-black)
        sig=Q.inertia(h);assert sig[1]==1
        assert Q.inertia(L.sub(h,ww))==[len(ww),0,0]
        for b in black:
            ids=sorted(set(ww)|{b});assert Q.inertia(L.sub(h,ids))==[len(ids),0,0]
        corners=[];coupon_total=F()
        for old in c['all_corners']:
            marks=old['preferred_corner'];ids=sorted(set(ww)|set(marks));sigc=Q.inertia(L.sub(h,ids))
            if marks==c['preferred_corner']:
                assert sigc==[len(ids)-1,0,1] and affine in D.components(h,ids)
                corners.append(dict(marks=marks,signature=sigc,kind='affineD14_order15',
                    used_ideal_payment='1/2',used_finite_coupon='0'))
                continue
            assert sigc==[len(ids),0,0]
            blocks=D.finite_blocks(h,ids)
            selected=[next(b for b in blocks if m in b['nodes']) for m in marks]
            coupons=[]
            for a,b,am,bm in ((selected[0],selected[1],marks[0],marks[1]),
                              (selected[1],selected[0],marks[1],marks[0])):
                if a['type'] in ('A4','A5') and b['type']=='D11' and b['roles'].get(str(bm))=='short_fork_leaf':
                    coupons.append(dict(small_type=a['type'],small_nodes=a['nodes'],small_mark=am,
                        large_type=b['type'],large_nodes=b['nodes'],large_mark=bm,
                        large_mark_role='short_fork_leaf',weight='1/2',whole_components=True))
            assert len(coupons)<=1
            value=F(1,2) if coupons else F();coupon_total+=value
            unused=[]
            if selected[0] is selected[1]:
                unused.append(dict(type=selected[0]['type'],nodes=selected[0]['nodes'],marks=marks,
                    marked_distance=Q.graph_dist(h,[marks[0]])[marks[1]],used_in_this_proof=False))
            corners.append(dict(marks=marks,signature=sigc,kind='finite',components=blocks,
                marked_component_indices=[blocks.index(b) for b in selected],used_ideal_payment='0',
                used_finite_coupon=str(value),coupon_roles=coupons,
                unused_same_component_payments=unused))
        assert coupon_total==F(1,2)
        assert sum(F(q['used_ideal_payment']) for q in corners)==F(1,2)
        variants.append(dict(optional_core_contact=entry,gram=L.encode(h),signature=sig,white=ww,
            corners=corners,used_finite_coupon_total='1/2',used_ideal_total='1/2',used_total='1',
            used_same_component_finite_units=False))
    return dict(core_gram=core['gram'],affine_nodes=affine,positive_affine_kernel=list(map(str,kernel)),
        all_extra_root_control=control,variants=variants)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    small=data[D.SMALL];over=data[D.OVER]
    rec=next(r for r in small['records'] if r['order']==15 and r['kind']=='D')
    target=rec['target'];prior=set(map(tuple,rec['retained_paid_pairs']));assert len(prior)==96
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    basic=X.L.S.basic_pairs(target);assert len(basic)==92 and prior-basic==HALVES
    assert all(dist[i][j]==10 for i,j in HALVES)
    assert all({tuple(sorted(a[i] for i in p)) for p in HALVES}==HALVES for a in autos)
    weights={p:F(1,2) if p in HALVES else F(1) if p in prior else F(0)
             for p in combinations(range(15),2)}
    assert sum(weights.values(),F())==94
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if r['kind']=='D' and r['order']==15)
    assert all(weights[tuple(p)]==1 for p in triangle['required_pairs'])
    assert all(weights[p]==1 for p in basic)
    short8={p for p in weights if dist[p[0]][p[1]]<=8};assert all(weights[p]==1 for p in short8)
    patterns=[p for p in small['patterns'] if p['order_bound']>=15 and len(p['nodes'])<=15]
    assert len(patterns)==536
    saved={r['case']:r for r in rec['case_embeddings']};support=[];affected=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:
            assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in HALVES]
        if selected:affected.append((p['index'],selected))
    assert len(support)==len(saved) and len(affected)==67
    rows=[];outcomes=Counter();full=[];variants=0;totalmaps=0
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        totalmaps+=len(maps)
        if len(c['marked_component'])==15:
            proof=full_case(c,core);full.append(index);variants+=len(proof['variants'])
            assert len(maps)==8
            outcomes['whole_affine_half_coupon_maps']+=len(maps)
            proof['maps']=maps
        else:
            proof=proper_case(c,core,maps,target)
            for q in proof['maps']:outcomes[q['outcome']]+=1
        rows.append(dict(index=index,core_key=c['core_key'],method=c['method'],proof=proof))
    assert full==FULL and variants==20 and totalmaps==408
    assert outcomes==Counter(projection_exclusion=232,two_entry_exclusion_or_shortcut=64,
        forced_lorentz_index_exclusion=16,whole_affine_half_coupon_maps=96)
    costs=dict(data[C15]['simultaneous_order15_costs']);assert costs==dict(A=75,B=94,C=91,D=96)
    costs['D']=94;assert max(costs.values())==94
    paths=[Path(__file__),Path(D.__file__),Path(X.__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER15_D94_EXISTING_HALF_COUPON_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='D14',diagram_order=15,target=target,automorphisms=autos,old_cost=96,new_cost='94',
        half_weight_pairs=sorted(HALVES),weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
        basic_pairs_retained=sorted(basic),short8_pairs_retained=sorted(short8),triangle_pairs_retained=triangle['required_pairs'],
        complete_eligible_patterns=len(patterns),reconstructed_high_support=support,
        affected_contexts=len(rows),affected_maps=totalmaps,whole_affine_contexts=full,whole_affine_variants=variants,
        outcomes=dict(outcomes),cases=rows,simultaneous_order15_costs=costs,new_order15_envelope='94/15',
        finite_weights_unchanged=True,added_finite_coupons=[],used_finite_same_component_units=False,
        exact_consumed_finite_interfaces='Existing half coupon: whole A4 or A5, any marked root, versus whole D11 marked at a short fork leaf.',
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete selected-map domain and all supporting-root controls. Proper factor completions excluded; already-whole factors have20 sealed variants using exactly one existing A4/A5-D11 half plus one ideal half. No finite D16/D17 unit consumed; basic/short/triangle and inherited nonselected/later interfaces retained. All-label core census inherited. No new ambient bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['affected_maps'],out['outcomes'],out['whole_affine_variants'],flush=True)
