#!/usr/bin/env python3
"""Use actual contact exclusions directly, without inverse end-set gates.

Radical sign constraints prune contact nodes; a protected leaf or an isolated
marked component then persists under arbitrarily many supporting roots.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_distance7_high_supporting_frontier_20261017 as S

ROOT,H=S.ROOT,S.H
JOINT='fresh_distance7_joint_contact_domains_20261017'
TAG='fresh_distance7_direct_contact_features_20261017'


def prune(allowed,zs,sc,du,dv):
    allowed=set(allowed);rounds=[]
    while True:
        witnesses=[]
        for i in sorted(allowed):
            for zi,z in enumerate(zs):
                if not z[i]:continue
                opposite=[k for k in sorted(allowed)if z[i]*z[k]<0]
                tests=[dict(node=k,shortcut=min(du[i]+dv[k],dv[i]+du[k])+2)for k in opposite]
                if all(t['shortcut']<=S.THRESHOLD for t in tests):
                    witnesses.append(dict(node=i,radical_index=zi,reason='opposite_sign_shortcuts',opposite_sign_tests=tests));break
                # Every core wall shares an actual edge with the new white.
                # In physical norm2 gauge, each nonzero product magnitude is
                # >=1 and <2. Equality of the bounds is already impossible.
                upper=2*sum((abs(z[k])for k in opposite),F())
                if all(sc[k]==1 for k in allowed)and abs(z[i])>=upper:
                    assert opposite
                    witnesses.append(dict(node=i,radical_index=zi,reason='strict_finite_contact_amplitudes',
                        opposite_sign_nodes=opposite,forced_same_sign_minimum=str(abs(z[i])),
                        strict_opposite_sign_upper_bound=str(upper)));break
        if not witnesses:return sorted(allowed),rounds
        rounds.append(dict(domain_before=sorted(allowed),excluded_contacts=witnesses))
        allowed-={w['node']for w in witnesses}


def direct_features(g,sc,corner,allowed):
    ns=corner['marked_component']
    if ns is None:return []
    allowed=set(allowed);forbidden=set(ns)-allowed;out=[]
    if not(allowed&set(ns)):
        out.append(dict(**corner,method='marked_component_isolated_by_actual_contact_domain',
            final_contact_domain=sorted(allowed),component_order_bound=len(ns)))
    adj={i:[k for k in ns if k!=i and g[i][k]]for i in ns}
    leaves=[i for i in ns if len(adj[i])==1 and i in forbidden]
    simple=[i for i in leaves if adj[i][0]in forbidden and len(adj[adj[i][0]])==2
        and g[i][adj[i][0]]**2==sc[i]*sc[adj[i][0]]]
    marked=[i for i in corner['preferred_corner']if i in leaves]
    distances=[H.graph_dist(H.extract(g,ns),[i])for i in range(len(ns))]
    a,b=corner['preferred_corner'];distance=distances[ns.index(a)][ns.index(b)]
    features=[i for i in ns if len(adj[i])>=3 or any(g[i][k]**2>sc[i]*sc[k]for k in adj[i])]
    near=min((distances[ns.index(x)][ns.index(f)]for x in(a,b)for f in features),default=999)
    forced=bool(simple and len(ns)>=10)
    if distance<=9 and marked and (forced or near<=5):
        out.append(dict(**corner,method='protected_leaf_by_actual_contact_domain',
            final_contact_domain=sorted(allowed),forbidden_component_contacts=sorted(forbidden),
            protected_marked_leaves=marked,protected_simple_leaves=simple,
            protected_simple_leaf_neighbors=[adj[i][0]for i in simple],base_order=len(ns),distance=distance,
            feature_distance=near,forced_finite=forced,ideal_basic_pair=near<=5))
    return out


def run():
    parent=json.loads((ROOT/f'results/{S.TAG}.json').read_text())
    joint=json.loads((ROOT/f'results/{JOINT}.json').read_text())
    for name,d in [(S.TAG,parent),(JOINT,joint)]:
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(name,p)
    pending=[r for r in joint['cases']if r['status']=='unresolved_supporting_extensions']
    assert len(pending)==548
    cores={};counts=Counter();bounds=Counter();rows=[]
    for src in pending:
        old=parent['cases'][src['parent_index']];key=src['core_key'];c=parent['cores'][key]
        if key not in cores:
            native=[[F(x)for x in row]for row in c['native_graph']]
            g,sc=H.native_to_rational(native);zs=H.nullspace(g)
            assert H.encode(g)==c['gram']and sc==c['gauge']and H.inertia(g)==c['signature']
            assert c['signature'][1]==1
            assert [[str(x)for x in z]for z in zs]==c['radical']
            cores[key]=(native,g,sc,zs)
        native,g,sc,zs=cores[key]
        du,dv=H.graph_dist(native,src['left_marks']),H.graph_dist(native,src['right_marks'])
        assert min(du[i]for i in src['right_marks'])>S.THRESHOLD
        assert set(src['joint_contact_domain'])<=set().union(*(set(e['nodes'])for e in src['actual_edge_certificates']))
        allowed,rounds=prune(src['joint_contact_domain'],zs,sc,du,dv)
        if rounds:counts['radical_pruned_contexts']+=1
        if any(w['reason']=='strict_finite_contact_amplitudes'for rd in rounds for w in rd['excluded_contacts']):counts['amplitude_pruned_contexts']+=1
        direct=[];opts=[]
        for corner in old['all_corners']:
            a,b=corner['preferred_corner']
            assert H.inspect_corner(g,sc,old['left_marks']+old['right_marks'],a,b)==corner
            direct.extend(direct_features(g,sc,corner,src['joint_contact_domain']))
            opts.extend(direct_features(g,sc,corner,allowed))
        if direct:counts['certified_without_radical_pruning']+=1
        row=dict(parent_index=src['parent_index'],base_context=src['base_context'],core_key=key,
            initial_contact_domain=src['joint_contact_domain'],radical_pruning_rounds=rounds,
            final_contact_domain=allowed,direct_without_radical=bool(direct))
        if opts:
            chosen=min(opts,key=lambda x:(x.get('component_order_bound',100),not x.get('forced_finite',False),
                len(x['final_contact_domain']),x['method'],x['preferred_corner']))
            row.update(status='structural_certificate',chosen=chosen)
            counts['method_'+chosen['method']]+=1
            if 'component_order_bound'in chosen:bounds[chosen['component_order_bound']]+=1
            else:
                counts['protected_forced_finite'if chosen['forced_finite']else'protected_finite_or_ideal_feature']+=1
            if not direct:counts['additional_after_radical_pruning']+=1
        else:row['status']='unresolved_supporting_extensions'
        counts[row['status']]+=1;rows.append(row)
    assert counts['structural_certificate']+counts['unresolved_supporting_extensions']==548
    paths=[Path(__file__)]+[ROOT/f'{d}/{s}.{e}'for s in [S.TAG,JOINT]for d,e in [('scripts','py'),('results','json'),('results','md')]]
    paths += [Path(H.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    return dict(status='PASS_DISTANCE7_DIRECT_ACTUAL_CONTACT_FEATURES',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        input_pending_contexts=548,previous_certified_contexts=6155,counts=dict(counts),
        component_bounds={str(k):v for k,v in sorted(bounds.items())},distinct_pending_core_matrices=len(cores),cases=rows,
        combined_certified_contexts=6155+counts['structural_certificate'],remaining_unresolved_contexts=counts['unresolved_supporting_extensions'],
        scope='All548 joint-domain pending cases. Actual Lorentz-span radical relations, threshold7 paths and strict1<=theta<2 finite-contact bounds in physical gauges prune impossible contacts. A marked component with no possible direct entry is whole. Otherwise forbidden actual contacts preserve marked leaves, simple leaves with degree-two neighbors, and existing branch/heavy features without needing positive inverse end sets. These are structural outputs with no contact realization assumed; finite7 selected payments, low/disjoint coverage, final cost and dimension comparison remain separate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],out['component_bounds'],flush=True)
