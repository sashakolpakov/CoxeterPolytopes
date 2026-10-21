#!/usr/bin/env python3
"""Order17 D120 to119, with a typed quarter coupon on whole E6 and D11."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_d16_partial_weight_transfer_20261017 as D
import fresh_d17_partial_weight_transfer_20261017 as OLD

ROOT=D.ROOT;L=D.L;Q=D.Q;P=D.P
TAG='fresh_order17_affine_d_three_quarter_transfer_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
CHANGE={(3,15),(3,16),(9,13),(9,14)}
EXPECTED=[506,509,512,515,1037,1051,1827,1830,1907,1981,2220,2223,2300,2374,2549,2551]
FULL=[1981,2374,2549,2551]
STEMS=(BASE,OLD.TAG,D.TAG,D.SMALL,D.OVER,OLD.OLD,
       'fresh_finite8_ideal5_combined_payments_20261017',
       'fresh_distance9_supporting_extensions_review_20261017')


def proper(c,core,maps,target):
    g,control=D.core_constraints(c,core);n=len(g);base=c['marked_component']
    assert c['method']=='lorentz_residual'
    assert control['kind']=='at_most_one_core_contact_root'
    J=Q.inverse(g);rows=[]
    for mp in maps:
        assert all((int(2-g[a][b]) if a!=b and g[a][b] else 0)==target[mp[i]][mp[j]]
                   for i,a in enumerate(base) for j,b in enumerate(base))
        image={t:base[i] for i,t in enumerate(mp)};missing=sorted(set(range(17))-set(mp));entries=[]
        for new in missing:
            ns=[image[t] for t in image if target[new][t]]
            assert len(ns)<=1
            if ns:
                i=ns[0];assert i in c['possible_contacts']
                entries.append(dict(new_target_node=new,core_node=i,native_label=3,
                    projection_square=str(J[i][i])))
        assert entries
        row=dict(map=mp,missing_target_nodes=missing,forced_entries=entries)
        if len(entries)>1:
            row['outcome']='excluded_multiple_core_contact_roots';rows.append(row);continue
        entry=entries[0]
        if F(entry['projection_square'])>2:
            row.update(outcome='excluded_projection_above_norm',projection_square=entry['projection_square']);rows.append(row);continue
        assert entry['core_node'] in control['candidate_nodes']
        image.update({t:n+i for i,t in enumerate(missing)});nn=n+len(missing)
        h=[r+[F(0)]*len(missing) for r in g]+[[F(0)]*nn for _ in missing]
        for a in missing:
            h[image[a]][image[a]]=F(2)
            for b in range(17):
                if target[a][b]:h[image[a]][image[b]]=h[image[b]][image[a]]=F(-1)
        ids=[image[t] for t in range(17)]
        assert Q.inertia(L.sub(h,ids))==[16,0,1]
        sig=Q.inertia(h);assert sig[1]>1
        row.update(outcome='excluded_forced_lorentz_index',target_to_literal=ids,gram=L.encode(h),signature=sig)
        rows.append(row)
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
        complete_residual_contact_control=control,maps=rows)


def e6_long_leaves(g,block):
    ids=block['nodes'];assert block['type']=='E6' and len(ids)==6
    adj={i:[j for j in ids if i!=j and g[i][j]] for i in ids}
    forks=[i for i in ids if len(adj[i])==3];assert len(forks)==1
    center=forks[0];dist=Q.graph_dist(L.sub(g,ids),[ids.index(center)])
    leaves=[i for i in ids if len(adj[i])==1]
    assert sorted(dist[ids.index(i)] for i in leaves)==[1,2,2]
    long=[i for i in leaves if dist[ids.index(i)]==2];assert len(long)==2
    labels=[[3 if i!=j and g[i][j] else 0 for j in ids] for i in ids]
    autos=P.embeddings(labels,labels)
    assert {ids[a[ids.index(long[0])]] for a in autos}==set(long)
    return dict(fork_center=center,long_arm_endpoint_leaves=long,automorphisms=autos)


def sealed(c,core,maps):
    g=[[F(x) for x in row] for row in core['gram']];n=len(g);assert n==19
    assert core['gauge']==[1]*n and Q.inertia(g)==[17,1,1]
    assert c['method']=='radical' and c['core_key']=='(25, 25, ((9, 9),))'
    z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(n))-black);assert len(black)==4
    assert Q.inertia(L.sub(g,white))==[15,0,0]
    for b in black:assert Q.inertia(L.sub(g,sorted(set(white)|{b})))==[16,0,0]
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    affine=c['marked_component'];assert affine==sorted(set(range(n))-{3,13})
    assert Q.inertia(L.sub(g,affine))==[16,0,1]
    ab=next(b for b in c['blocks'] if b['nodes']==affine)
    kernel=list(map(F,ab['positive_kernel']));assert all(x>0 for x in kernel)
    assert all(sum(g[i][j]*kernel[t] for t,j in enumerate(affine))==0 for i in affine)
    assert c['possible_contacts']==[3,13] and z[3]*z[13]<0
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v)
    shortcut=min(du[3]+dv[13],du[13]+dv[3])+2;assert shortcut==2
    corners=[];total=F()
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks));sig=Q.inertia(L.sub(g,ids))
        if marks==c['preferred_corner']:
            assert sig==[16,0,1]
            corners.append(dict(marks=marks,signature=sig,kind='affineD16_order17',
                used_ideal_payment='3/4',used_new_finite_coupon='0'));continue
        assert sig==[17,0,0]
        blocks=D.finite_blocks(g,ids);selected=[next(b for b in blocks if m in b['nodes']) for m in marks]
        coupons=[]
        for a,b,am,bm in ((selected[0],selected[1],marks[0],marks[1]),
                          (selected[1],selected[0],marks[1],marks[0])):
            if a['type']=='E6' and b['type']=='D11' and b['roles'].get(str(bm))=='short_fork_leaf':
                roles=e6_long_leaves(g,a);assert am in roles['long_arm_endpoint_leaves']
                short=[int(i) for i,role in b['roles'].items() if role=='short_fork_leaf'];assert len(short)==2
                coupons.append(dict(small_type='E6',small_nodes=a['nodes'],small_mark=am,
                    small_role='endpoint_of_an_arm_of_length_two',small_role_check=roles,
                    large_type='D11',large_nodes=b['nodes'],large_mark=bm,
                    large_role='short_fork_leaf',all_large_short_leaves=short,weight='1/4',whole_components=True,
                    added_total_per_component_pair='1'))
        assert len(coupons)<=1
        value=F(1,4) if coupons else F();total+=value
        unused=[]
        if selected[0] is selected[1]:
            unused.append(dict(type=selected[0]['type'],nodes=selected[0]['nodes'],marks=marks,
                marked_distance=Q.graph_dist(g,[marks[0]])[marks[1]],used_in_this_proof=False))
        corners.append(dict(marks=marks,signature=sig,kind='finite',components=blocks,
            marked_component_indices=[blocks.index(b) for b in selected],used_ideal_payment='0',
            used_new_finite_coupon=str(value),new_coupon_roles=coupons,unused_same_component_payments=unused))
    assert total==F(1,4) and sum(F(q['used_ideal_payment']) for q in corners)==F(3,4)
    return dict(core_gram=core['gram'],signature=[17,1,1],radical=list(map(str,z)),white=white,
        affine_nodes=affine,positive_affine_kernel=list(map(str,kernel)),allowed_contacts=[3,13],
        opposite_radical_coefficients=[str(z[3]),str(z[13])],both_contacts_shortcut=shortcut,
        whole_component_sealing='Every further white is orthogonal to the entire19-root core unless the length2 path pays.',
        maps=maps,corners=corners,used_new_finite_coupon_total='1/4',used_ideal_total='3/4',used_total='1',
        used_same_component_finite_units=False)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    small=data[D.SMALL];over=data[D.OVER];prior=data[OLD.TAG]
    rec=next(r for r in small['records'] if r['order']==17 and r['kind']=='D')
    target=rec['target'];assert target==prior['target']
    weights={tuple(r['pair']):F(r['weight']) for r in prior['weighted_pairs']}
    assert sum(weights.values(),F())==120 and all(weights[p]==1 for p in CHANGE)
    for p in CHANGE:weights[p]=F(3,4)
    assert sum(weights.values(),F())==119
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all(dist[i][j]==10 for i,j in CHANGE)
    assert all({tuple(sorted(a[i] for i in p)) for p in CHANGE}==CHANGE for a in autos)
    basic=L.S.basic_pairs(target);assert all(weights[p]==1 for p in basic)
    triangles=prior['triangle_pairs_retained_at_weight_one'];assert all(weights[tuple(p)]==1 for p in triangles)
    assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=8)
    patterns=[p for p in small['patterns'] if p['order_bound']>=17 and len(p['nodes'])<=17]
    assert len(patterns)==140
    saved={r['case']:r for r in rec['case_embeddings']};support=[];affected=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:
            assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in CHANGE]
        if selected:affected.append((p['index'],selected))
    assert len(support)==len(saved) and [i for i,mp in affected]==EXPECTED
    rows=[];outcomes=Counter();full=[];totalmaps=0
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']];totalmaps+=len(maps)
        if len(c['marked_component'])==17:
            proof=sealed(c,core,maps);full.append(index);assert len(maps)==8
            outcomes['sealed_E6_D11_quarter_coupon_maps']+=len(maps)
        else:
            proof=proper(c,core,maps,target)
            for q in proof['maps']:outcomes[q['outcome']]+=1
        rows.append(dict(index=index,core_key=c['core_key'],method=c['method'],proof=proof))
    assert full==FULL and totalmaps==96
    assert outcomes==Counter(excluded_projection_above_norm=32,excluded_multiple_core_contact_roots=16,
        excluded_forced_lorentz_index=16,sealed_E6_D11_quarter_coupon_maps=32)
    costs=dict(data[BASE]['small_order_costs']['17']);assert costs==dict(A=85,B=111,C=108,D=120)
    costs['D']=119;assert max(costs.values())==119
    paths=[Path(__file__),Path(D.__file__),Path(OLD.__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER17_D119_TYPED_QUARTER_COUPON_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='D16',diagram_order=17,target=target,automorphisms=autos,old_cost='120',new_cost='119',
        changed_pairs=sorted(CHANGE),old_changed_weight='1',new_changed_weight='3/4',
        weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
        basic_pairs_retained=sorted(basic),triangle_pairs_retained=triangles,
        complete_eligible_patterns=len(patterns),reconstructed_high_support=support,
        affected_contexts=len(rows),affected_maps=totalmaps,whole_affine_contexts=full,outcomes=dict(outcomes),cases=rows,
        new_order17_costs=costs,new_order17_envelope='7',
        added_finite_coupon=dict(weight='1/4',small_type='whole E6',
            small_mark='one of its two endpoint leaves at distance2 from its unique degree3 node',
            large_type='whole D11',large_mark='one of its two short fork leaves',
            number_of_marked_pairs_per_component_pair=4,added_total_per_component_pair='1'),
        finite_budget_proved_here=False,used_finite_same_component_units=False,
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete140-pattern selected support and all96 affected maps. Proper factor completions excluded; four actual whole-factor cores seal the exact E6/D11 coupon. New finite cost1 per E6-D11 component pair must be budgeted separately. Existing finite D16/D17 same-component units unused. All-label core/no-bridge census inherited; no ambient theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['affected_maps'],out['outcomes'],flush=True)
