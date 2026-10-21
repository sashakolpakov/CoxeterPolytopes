#!/usr/bin/env python3
"""Read saved D15 support and inspect eleven narrow Lorentz cores.

This is a diagnostic, not a transfer theorem. No native-map census or DP.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_finite_rank15_existing_half_transfer_20261018 as R

ROOT=R.ROOT
TAG='fresh_finite_d15_residual_orbit_pilot_20261018'
SCOUT=R.R.TAG
RULE='fresh_rank15_sparse_payment_exchange_20261018'
CURRENT='fresh_order13_sparse_ambient60_frontier_20261018'
TARGETS={(0,11),(0,12),(1,12)}


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def key(r):return r['source'],json.dumps(r['case'],sort_keys=True)


def run():
    stems=(SCOUT,RULE,CURRENT,R.A.OVER)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    scout=docs[SCOUT];rec=next(z for z in scout['records'] if z['kind']=='D')
    rule=next(z for z in docs[RULE]['rules'] if z['kind']=='D')
    weights={tuple(z['pair']):F(z['weight']) for z in rule['weighted_pairs']}
    assert sum(weights.values())==102 and rec['complete_eligible_patterns']==536
    current={key(z):z for z in docs[CURRENT]['combined_literal_interface']['rows']}
    assert len(current)==216==len(scout['literal_interfaces'])
    orbits=[]
    for orbit in rec['long_pair_orbits']:
        pairs=set(map(tuple,orbit['pairs']));old=weights[next(iter(pairs))]
        assert all(weights[p]==old for p in pairs)
        trials=[]
        for new in sorted({F(),old/2,F(3,4)}&{F(i,8) for i in range(9)}):
            if new>=old:continue
            hits=[];minimum=None
            for row in scout['literal_interfaces']:
                now=current[key(row)];assert len(row['corners'])==len(now['corners'])
                debit=F();indices=[]
                for i,(before,after) in enumerate(zip(row['corners'],now['corners'])):
                    assert before['marks']==after['marks']
                    images=set(map(tuple,before['whole_rank15_images'].get('D',[])))
                    if images&pairs:
                        assert images<=pairs
                        debit+=old-new;indices.append(i)
                value=F(now['total_lower_bound'])-debit
                minimum=value if minimum is None else min(minimum,value)
                if indices:hits.append(dict(source=row['source'],case=row['case'],corners=indices,
                    current_total=now['total_lower_bound'],debit=debit,hypothetical_total=value))
            trials.append(dict(new_weight=new,factor_saving=len(pairs)*(old-new),minimum_total=minimum,
                affected=hits,violations=[h for h in hits if h['hypothetical_total']<1]))
        orbits.append(dict(pairs=sorted(pairs),distance=orbit['distance'],current_weight=old,
            selected_obligations=orbit['complete_selected_obligations'],trials=trials))
    selected={r['index']:r for r in rec['selected_patterns']}
    ids=sorted({r['index'] for o in orbits if set(o['pairs'])&TARGETS for r in o['selected_obligations']})
    assert ids==[1030,1147,1269,1283,1284,1613,1622,1623,1641,1682,2501]
    over=docs[R.A.OVER];proofs=[]
    for i in ids:
        c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        g,J,du,dv,end=R.X.constraints(c,core)
        base=c['marked_component'];assert all(g[j][j]==2 for j in range(len(g)))
        left=c['left'][1];right=[c['right_copy_map'][j] for j in c['right'][1]]
        white=sorted(set(range(len(g)))-set(left+right))
        inertia=R.Q.inertia(g);assert inertia==[len(g)-1,1,0]
        maps=[]
        for m in selected[i]['maps']:
            if tuple(m['marked_images']) not in TARGETS:continue
            image={j:k for k,j in zip(base,m['images'])}
            missing=sorted(set(range(15))-set(image));forced=[]
            for t in missing:
                contacts=[(j,image[j],rec['target'][t][j]) for j in image if rec['target'][t][j]]
                assert all(label==3 for _,_,label in contacts)
                q=sum((J[a][b] for _,a,_ in contacts for _,b,_ in contacts),F())
                forced.append(dict(target_node=t,old_contacts=contacts,forced_projection=q,
                    contact_inverse_diagonals=[[a,J[a][a]] for _,a,_ in contacts],
                    old_contact_status=[a in c['possible_contacts'] for _,a,_ in contacts]))
            maps.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,forced=forced))
        assert maps
        remaining=sorted(set(c['possible_contacts'])-set(base))
        corners=[]
        for co in c['all_corners']:
            marks=co['preferred_corner'];corner=sorted(set(white)|set(marks))
            sig=R.Q.inertia(R.A.sub(g,corner));blocks=[]
            if sig==[len(corner),0,0]:
                for b in R.D.finite_blocks(g,corner):
                    mm=[m for m in marks if m in b['nodes']]
                    if not mm:continue
                    typ=b['type'];images=[]
                    if typ[0] in 'ABD' and typ[1:].isdigit():
                        images=R.T.mark_pairs(g,b['nodes'],mm,int(typ[1:]),typ[0])
                    blocks.append(dict(type=typ,nodes=b['nodes'],marks=mm,canonical_images=images))
            corners.append(dict(marks=marks,inertia=sig,marked_blocks=blocks))
        proofs.append(dict(index=i,core_key=c['core_key'],core_order=len(g),inertia=inertia,
            displayed_order=len(base),whole_selected=len(base)==15,white=white,
            allowed_contacts=c['possible_contacts'],remaining_contacts=remaining,
            remaining_inverse_diagonals=[[j,J[j][j]] for j in remaining],end_projection=end,
            required_completion_maps=maps,corners=corners))
    paths={Path(__file__),Path(R.__file__),Path(R.X.__file__),Path(R.T.__file__),Path(R.A.__file__),Path(R.Q.__file__),Path(R.D.__file__)}
    paths.update(ROOT/f'results/{s}.json' for s in stems)
    return enc(dict(status='PASS_D15_RESIDUAL_ORBIT_PILOT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        current_D15_cost='102',current_literal_records=216,orbits=orbits,
        inspected_pairs=sorted(TARGETS),exact_core_diagnostics=proofs,
        new_payment_theorem=False,new_forest_bound=None,
        scope='Saved complete536-pattern support inventory and current216 literal debit trials, plus eleven exact Lorentz-core diagnostics. Required-contact projections are not arbitrary-contact exclusions by themselves; displayed alternative factors are not whole without sealing. No map/census/DP or ambient theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'orbits',len(result['orbits']),'cores',len(result['exact_core_diagnostics']),flush=True)
