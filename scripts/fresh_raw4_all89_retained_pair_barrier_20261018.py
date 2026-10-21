#!/usr/bin/env python3
"""One retained raw4 pair closes the all89 fixed-mass pair-gate opening.

No covering LP, new clique census, global optimizer, or physical-realization claim.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_new5121_raw4_all89_cover_barrier_20261018 as B

Q=B.Q;Z=Q.Z;R=Q.R;N=Q.A;ROOT=Q.ROOT
TAG='fresh_raw4_all89_retained_pair_barrier_20261018'
VERTEX=495;FACE=14368;SELECTED=(0,9)


def build():
    barrier=Z.M.read(B.TAG);rawdata=Z.M.read(Q.TAG);base=Z.M.read(Z.TAG)
    assert barrier['status']=='PASS_ALL89_RAW4_CORNER_MASS_AND_PAIR_ONLY_CARDINALITY_BARRIER'
    raw=next(r for r in rawdata['forest_records'] if r['vertex']==VERTEX)
    assert VERTEX not in {r['vertex'] for r in base['forest_records']}
    p=Z.prepare();assert p['base']['face_multipliers']==base['face_multipliers']
    p['mu']={v:F(x) for v,x in base['vertex_multipliers']}
    u,g,owner,pairs,oldedges,demands=Z.merge(raw,p)
    inventory=next(r for r in rawdata['cover_inventory']['forests'] if r['vertex']==VERTEX)
    assert u['raw4_to_merged_star_indices']==inventory['raw4_to_merged_star_indices']
    assert u['added_occurrences']==inventory['actual_occurrences']
    assert u['baseline']==inventory['baseline']
    assert u['old_literal_stars']==306 and u['new_literal_stars']==16
    assert u['new_pair_positions']==5016 and len(u['corner_orbits'])==15
    assert all(o['origin']['face']==FACE for o in u['added_occurrences'])
    assert len(u['added_occurrences'])==2
    assert [o['orbit'] for o in u['added_occurrences']]==[13,14]
    pair=tuple(u['raw4_to_merged_star_indices'][i] for i in SELECTED)
    assert pair==(306,315) and [owner[i] for i in pair]==[13,14]
    assert min(pair)>=u['old_literal_stars'] and pair in pairs
    alpha={f:F(x) for f,x in base['face_multipliers']};mass=alpha[FACE]
    assert all(u['corner_orbits'][i]['positive_face_coefficients']==[[FACE,1]] for i in (13,14))
    assert demands[13]==demands[14]==mass
    assert not set(map(tuple,u['literal_stars'][pair[0]]['occupied_pairs']))&set(map(tuple,u['literal_stars'][pair[1]]['occupied_pairs']))
    # Only this one final pair is tested. The exploratory eight-pair pilot is not
    # imported and none of its other outcomes is installed in the current model.
    physical=Z.physical(u,g,[pair])
    assert physical['newly_owned_pair_domain']==[[306,315,None]]
    assert physical['new_exclusions']==0 and physical['new_unknown_pairs']==1
    assert {tuple(sorted(e)) for e in g.edges() if max(e)<306}==oldedges
    selected=[u['literal_stars'][i] for i in pair]
    cov=[N.native(s) for s in selected]
    assert cov[0][0]==cov[1][0] and len(set(cov[0]+cov[1]))==3
    # Save every distinct external-wall product, including the two inherited
    # within-star products, not just the first possible rejection.
    target=R.P.forest_graph(raw['forest'],dict(R.C.models()))
    assert Z.normal(target)==raw['target']
    walls=list(dict.fromkeys(cov[0]+cov[1]));J,active,blocks=R.projection_data(target)
    assert all(i in active for wall in walls for i,edge in wall)
    heights=[]
    for c in walls:
        field=N.pair_field((c,));z=R.projection(c,c,J,field)
        assert not any(z[1:]);h=z[0]-2;assert h.sign()>0;heights.append(h)
    products=[]
    for i,j in combinations(range(3),2):
        field=N.pair_field((walls[i],walls[j]));z=R.projection(walls[i],walls[j],J,field);h=heights[i]*heights[j]
        sg=field.signs_difference(z,h);sg1=field.signs_difference(R.qaddk(z,1),h)
        code,reason=field.arithmetic_wall_gate(tuple((a.a,a.b) for a in z),h.a,h.b)
        assert sg<=0 and not(sg<0 and sg1>0) and code is None
        products.append(dict(walls=[i,j],projection=R.qout(z),height_product=h.out(),
            pair_field_radicals=[field.tau.out(),field.sigma.out()],sign=sg,sign_plus_one=sg1,
            arithmetic_gate_code=code,arithmetic_gate_reason=reason))
    K=Q.K;tau=K(F(7,2),F(1,2));S=K(F(13,2),F(5,2))
    assert heights==[K(-F(1,4),F(1,4)),K(F(3,2),1),K(1,1)]
    assert heights[0]*heights[1]==tau/4 and heights[0]*heights[2]==K(1)
    assert heights[1]*heights[2]==S and (S-4).sign()>0 and Q.M.integral(S)
    # This additionally exhibits the only nontrivial native squareclass cycle.
    assert tau*S==K(3,2)*K(3,2)
    num,den=F(barrier['objective_numerator']),F(barrier['current_normalizer'])
    records=[];gain=F();changed=[]
    for row in barrier['all89']['forests']:
        r=dict(row);v=r['vertex'];t=F(r['baseline_local_cover']);mu=F(r['baseline_vertex_mass']);e=F(r['actual_corner_demand'])
        saving=mass if v==VERTEX else F()
        assert e>=2*saving
        upper=max(mu,t+e-saving);inc=upper-mu;old=F(r['maximum_possible_mu_increase'])
        assert 0<=inc<=old;gain+=inc
        r.update(pair_cover_mass=str(saving),possible_vertex_mass_upper=str(upper),maximum_possible_mu_increase=str(inc),
                 saving_after_multiplier_slack=str(old-inc))
        records.append(r)
        if saving:changed.append(r)
    assert len(changed)==1 and F(changed[0]['saving_after_multiplier_slack'])==mass
    assert F(barrier['all89']['maximum_possible_mu_increase'])-gain==mass
    assert num>12*(den+gain)
    # Retain every actual current multiplier; this is an upper repair covering
    # all89, not a newly optimized or installed complete local cone.
    mu={v:F(x) for v,x in base['vertex_multipliers']};upper=dict(mu)
    for r in records:upper[r['vertex']]=F(r['possible_vertex_mass_upper'])
    assert all(upper.get(v,F())>=x for v,x in mu.items()) and sum(upper.values(),F())==den+gain
    paths=set(p['sources'])
    for tag in (B.TAG,Q.TAG,Z.TAG):
        d=Z.M.read(tag)
        for path,h in d['sources'].items():assert Z.M.digest(path)==h,(tag,path)
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):paths.add(f'{folder}/{tag}.{ext}')
    for mod in (B,Q,Z,N,R,R.T):paths.add(str(Path(mod.__file__).relative_to(ROOT)))
    paths.add(f'scripts/{TAG}.py')
    out=dict(status='PASS_ALL89_CURRENT_PAIR_GATES_RETAIN_ABOVE12_WITNESS',
        sources={path:Z.M.digest(path) for path in sorted(paths)},
        selected_vertex=VERTEX,selected_face=FACE,raw_star_indices=list(SELECTED),merged_star_indices=list(pair),
        forest=raw['forest'],target=raw['target'],baseline=u['baseline'],complete_selected_occurrences=u['added_occurrences'],
        complete_selected_orbits=[u['corner_orbits'][i] for i in (13,14)],selected_native_stars=selected,
        old_graph_edges_sha256=Z.L.jhash(sorted(oldedges)),new_native_star_count=16,old_native_star_count=306,
        pair_outcome=physical,distinct_external_contacts=walls,squared_external_heights=[h.out() for h in heights],
        complete_distinct_wall_products=products,projection_components=blocks,
        explicit_products=dict(AB='0',AC='-1',BC_negative_square=S.out(),tau=tau.out(),cycle_positive_K_root=K(3,2).out()),
        face_multipliers=base['face_multipliers'],initial_vertex_multipliers=base['vertex_multipliers'],
        all89_cover_upper_records=records,upper_vertex_multipliers=[[v,str(x)] for v,x in sorted(upper.items()) if x],
        pair_cover_mass=str(mass),saving_after_multiplier_slack=str(mass),objective_numerator=str(num),
        initial_normalizer=str(den),maximum_possible_mu_increase=str(gain),normalizer_upper=str(den+gain),
        ratio_lower_bound=str(num/(den+gain)),strict_numerator_gap_over12=str(num-12*(den+gain)),
        original_all89_ratio_lower_bound=barrier['all89']['ratio_lower_bound'],
        pair_tests=1,new_clique_enumerations=0,covering_LPs=0,global_model_optimizations=0,
        authoritative_classes=324,authoritative_forests=600,new_dimension_bound=None,
        scope='A retained pair under the current exact port/angle/two-cycle gates saves one full face mass in the complete all89 singleton-cover upper repair. All current old outcomes and unknowns remain fixed. The above12 witness survives every pair-only extension using those gates and retaining this pair; subset extensions follow by projection. Unknown is not physical realizability or immunity to a future stronger geometric or higher-incidence exclusion. No complete new local graph/optimum is claimed or installed.')
    print(out['status'],'pair',pair,'saving',mass,'ratio_lower',float(num/(den+gain)),flush=True)
    return Z.normal(out)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/f'results/{TAG}.json'
    d=build()
    if a.verify:assert d==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(d,separators=(',',':'))+'\n')
