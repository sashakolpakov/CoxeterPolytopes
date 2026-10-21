#!/usr/bin/env python3
"""Two affine-B singleton payments vanish by impossible actual completions.

The proof uses norm-two roots throughout. No unspecified contacts are set to
zero, no finite coupon is added, and no ambient dimension bound is asserted.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_order13_affine_b_singleton_transfer_20261017 as L

ROOT=L.ROOT; P=L.P; Q=L.Q
TAG='fresh_b14_b15_singleton_completion_transfer_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
OVER='fresh_distance8_high_overlap_20261017'
PRIOR_B15='fresh_order15_affine_b_pair_reduction_20261017'
COMBINED='fresh_finite8_ideal5_combined_payments_20261017'
DROPS={14:(0,10),15:(1,11)}
EXPECTED={
    14:[366,380,392,531,555,587,614,617,704,1115,1159,1241,1246,1285,1503,1505,1508,1541,1544,1550,1619,1624,1642],
    15:[366,611,986,1031,1159,1503,1614,1619,1642,1683,2502]}


def constraints(c,core):
    g=[[F(x) for x in row] for row in core['gram']]; n=len(g)
    assert core['gauge']==[1]*n
    assert all(g[i][i]==2 for i in range(n))
    assert set(x for row in g for x in row)<={F(-1),F(0),F(2)}
    assert Q.inertia(g)==[n-1,1,0]
    inv=Q.inverse(g); assert L.encode(inv)==core['inverse']
    u=c['left'][1]; v=[c['right_copy_map'][i] for i in c['right'][1]]
    corner=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key,value in corner.items(): assert c[key]==value,(c['index'],key)
    du=Q.graph_dist(g,u); dv=Q.graph_dist(g,v)
    ends=c['end_sets']; allowed=c['possible_contacts']
    assert set(ends[0]).isdisjoint(ends[1]) and set(ends[0]+ends[1])==set(allowed)
    shortcut=max((du[i] for i in ends[0]),default=0)+max((dv[i] for i in ends[1]),default=0)+2
    assert shortcut<=8
    assert all(inv[i][j]>0 for end in ends for i in end for j in end)
    black=set(u+v); old_corner=(set(range(n))-black)|set(c['preferred_corner'])
    other=sorted(old_corner-set(c['marked_component']))
    assert all(g[i][j]==0 for i in c['marked_component'] for j in other)
    assert sorted(set(range(n))-old_corner)==c['omitted_black']
    end_proof=dict(allowed_contacts=allowed,end_sets=ends,both_ends_shortcut_bound=shortcut,
                   inverse_end_blocks=[L.encode(L.sub(inv,end)) for end in ends],
                   other_old_corner_nodes=other,old_nodes_absent_from_corner=c['omitted_black'])
    return g,inv,du,dv,end_proof


def inspect_map(c,g,inv,du,dv,target,mp):
    m=len(target); n=len(g); base=c['marked_component']
    assert len(set(mp))==len(mp)==len(base)
    assert all((int(2-g[a][b]) if a!=b and g[a][b] else 0)==target[mp[i]][mp[j]]
               for i,a in enumerate(base) for j,b in enumerate(base))
    image={mp[t]:i for t,i in enumerate(base)}
    missing=sorted(set(range(m))-set(image))
    # These are only forced target incidences, not a completed Gram with
    # unspecified outside contacts artificially set to zero.
    entries=[]
    for new in missing:
        for old in sorted(image):
            bond=target[new][old]
            if bond:
                assert bond in (3,4)
                node=image[old]; q2=F(1 if bond==3 else 2)
                assert node in c['possible_contacts']
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,
                                    coxeter_label=bond,norm_two_contact_square=str(q2),
                                    full_core_inverse_diagonal=str(inv[node][node]),
                                    forced_projection_lower_bound=str(q2*inv[node][node])))
    assert entries and len({r['new_target_node'] for r in entries})==len(entries)
    row=dict(map=mp,missing_target_nodes=missing,forced_entries=entries)
    # Additional contacts in the other end already give a retained short path.
    # Within one end the inverse is entrywise positive, so this is a lower bound.
    bad=next((r for r in entries if F(r['forced_projection_lower_bound'])>2),None)
    assert bad is not None
    assert bad['new_target_node']==m-3 and bad['old_target_node']==m-4
    assert bad['coxeter_label']==4 and bad['norm_two_contact_square']=='2'
    assert bad['core_node'] in c['preferred_corner']
    row.update(outcome='projection_exclusion',witness=bad)
    return row


def run():
    stems=[SMALL,OVER,PRIOR_B15,COMBINED]
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,j in inputs.items():
        assert j['status'].startswith('PASS'),stem
        for name,h in j['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    records=[]; over=inputs[OVER]
    for m,drop in DROPS.items():
        rec=next(r for r in inputs[SMALL]['records'] if r['order']==m and r['kind']=='B')
        target=rec['target']; autos=P.embeddings(target,target); _,ds,features=P.distance_data(target)
        assert all(tuple(sorted(a[i] for i in drop))==drop for a in autos)
        assert ds[drop[0]][drop[1]]==10
        paid=set(map(tuple,rec['retained_paid_pairs'] if m==14 else inputs[PRIOR_B15]['new_paid_pairs']))
        assert len(paid)==(87 if m==14 else 97) and drop in paid
        affected=[r for r in rec['case_embeddings'] if list(drop) in r['marked_images']]
        assert [r['case'] for r in affected]==EXPECTED[m]
        rows=[]; outcomes=Counter(); count=0
        for r in affected:
            c=over['isolation']['cases'][r['case']]; core=over['isolation']['cores'][c['core_key']]
            g,inv,du,dv,end_proof=constraints(c,core)
            maps=[]
            for mp in r['embeddings']:
                marked=tuple(sorted(mp[c['marked_component'].index(i)] for i in c['preferred_corner']))
                if marked!=drop: continue
                row=inspect_map(c,g,inv,du,dv,target,mp)
                maps.append(row); outcomes[row['outcome']]+=1; count+=1
            assert maps
            rows.append(dict(index=c['index'],core_key=c['core_key'],method=c['method'],
                             left=c['left'],right=c['right'],match=c['match'],core_gram=core['gram'],core_inverse=core['inverse'],
                             core_norms=[2]*len(g),base_nodes=c['marked_component'],preferred_corner=c['preferred_corner'],
                             contact_domain=c['possible_contacts'],end_sets=c['end_sets'],end_projection_proof=end_proof,maps=maps))
        reduced=paid-{drop}
        triangle=next(r for r in inputs[COMBINED]['triangle_pair_retention'] if r['order']==m and r['kind']=='B')
        triangle_pairs=set(map(tuple,triangle['required_pairs'])); assert triangle_pairs<=reduced
        short8={p for p in combinations(range(m),2) if ds[p[0]][p[1]]<=8}; assert short8<=reduced
        basic=L.S.basic_pairs(target); assert basic<=reduced
        assert len(reduced)==(86 if m==14 else 96)
        assert count==(46 if m==14 else 18) and outcomes=={'projection_exclusion':count}
        records.append(dict(order=m,affine_type='B'+str(m-1),target=target,automorphisms=autos,
                            inherited_cost=len(paid),new_cost=len(reduced),removed_pair=list(drop),
                            new_paid_pairs=[list(p) for p in sorted(reduced)],affected_cases=EXPECTED[m],
                            affected_map_count=count,map_outcomes=dict(outcomes),cases=rows,
                            triangle_pairs_retained=[list(p) for p in sorted(triangle_pairs)],
                            short8_pairs_retained=[list(p) for p in sorted(short8)],basic_pairs_retained=True))
    paths=[Path(__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_B14_B15_TWO_SINGLETON_COMPLETION_EXCLUSIONS',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},records=records,
                root_normalization='Every old core root and every actual extra root has square2; a target label4 contact has square2, not1.',
                added_finite_coupons=[],triangle_minimum=2,quadrilateral_minimum=1,
                scope='Complete support and exact necessary-contact exclusions for two invariant distance10 singletons. No full Gram is completed by assigning unspecified contacts zero. No new finite budget or global dimension bound is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--verify',action='store_true'); a=p.parse_args()
    out=json.loads(json.dumps(run())); path=ROOT/f'results/{TAG}.json'
    if a.verify: assert out==json.loads(path.read_text())
    else: path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(r['order'],r['affected_map_count'],r['map_outcomes'],r['new_cost']) for r in out['records']],flush=True)
