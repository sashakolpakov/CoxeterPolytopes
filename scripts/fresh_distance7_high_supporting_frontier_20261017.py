#!/usr/bin/env python3
"""Threshold-seven supporting-root certificates on the complete cross domain.

Unresolved cases are retained explicitly. This does not assert a complete
quadrilateral-payment theorem or choose the final finite/ideal pair weights.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_distance7_high_cross_completions_20261017 as X

H=X.H;ROOT=X.ROOT;TAG='fresh_distance7_high_supporting_frontier_20261017'
THRESHOLD=7


def split_ends(allowed,du,dv,J):
    for k in range(THRESHOLD-1):
        left=[i for i in allowed if du[i]<=k]
        right=[i for i in allowed if dv[i]<=THRESHOLD-2-k]
        if set(left)|set(right)!=set(allowed) or set(left)&set(right):continue
        if all(J[i][j]>0 for block in (left,right) for i in block for j in block):yield [left,right]


def ordinary(sc,J,zs,c,du,dv):
    allowed=c['possible_contacts'];size=len(c['marked_component'])
    if zs:
        for z in zs:
            if any(not z[i] for i in allowed):continue
            pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
            bounds=[max((du[i] for i in pos),default=0)+max((dv[i] for i in neg),default=0)+2,
                    max((dv[i] for i in pos),default=0)+max((du[i] for i in neg),default=0)+2]
            if min(bounds)<=THRESHOLD:
                return dict(**c,method='radical',used_radical=list(map(str,z)),shortcut_bound=min(bounds),component_order_bound=size)
        return None
    ends=next(split_ends(allowed,du,dv,J),None)
    if ends is None:return None
    pp={i:sc[i]*J[i][i] for i in allowed};cand=[i for i in allowed if pp[i]<=2]
    if any(pp[i]<=1 for i in cand):return None
    if any(pp[i]+pp[j]+2*min(sc[i],sc[j])*J[i][j]<=2 for end in ends for i,j in combinations(end,2)):return None
    tests=[]
    for i,j in product(cand,repeat=2):
        a=J[i][j];rho=(2-pp[i])*(2-pp[j]);sq=sc[i]*sc[j]*a*a
        dist=min(du[i]+dv[j],du[j]+dv[i])+3
        if sq<=rho or (a<0 and dist>THRESHOLD):return None
        tests.append(dict(nodes=[i,j],projection_square=str(sq),residual_product=str(rho),sign=1 if a>0 else -1,shortcut=dist if a<0 else None))
    roots=[]
    for i in cand:
        rho=2-pp[i];assert 0<=rho<1
        roots.append(dict(contact_node=i,projection_square=str(pp[i]),residual_square=str(rho),maximum_attached_orthogonal_rank=rho//(1-rho)))
    bound=size+max((1+r['maximum_attached_orthogonal_rank'] for r in roots),default=0)
    return dict(**c,method='lorentz_residual',end_sets=ends,single_contact_candidates=roots,pair_tests=tests,component_order_bound=bound)


def relevant(J,c,du,dv):
    allowed=c['possible_contacts'];ends=next(split_ends(allowed,du,dv,J),None)
    if ends is None:return None
    if any(J[i][i]+J[j][j]+2*J[i][j]<=2 for end in ends for i,j in combinations(end,2)):return None
    cand=[i for i in allowed if J[i][i]<=2]
    if any(J[i][i]<1 for i in cand):return None
    good=[i for i in cand if i in c['marked_component']]
    if any(J[i][i]<=1 for i in good):return None
    records=[]
    for i in good:
        rho=2-J[i][i];outside=[];checks=[]
        for j in cand:
            sig=2-J[j][j];a=J[i][j];dist=min(du[i]+dv[j],du[j]+dv[i])+3
            if a*a>rho*sig and a>0:reason='positive'
            elif a*a>rho*sig and a<0 and dist<=THRESHOLD:reason='short_path'
            elif a>0 and (a+1)**2>rho*sig:
                if j in good:return None
                outside.append(j);reason='forced_orthogonal'
            else:return None
            checks.append(dict(node=j,reason=reason,projection=str(a),residual_product=str(rho*sig),shortcut=dist))
        doubles=[]
        for j in cand:
            if J[j][j]!=1:continue
            a=J[i][j];dist=min(du[i]+dv[j],du[j]+dv[i])+3
            if not(a>0 or (a<0 and dist<=THRESHOLD)):return None
            doubles.append(dict(node=j,sign=1 if a>0 else -1,shortcut=dist))
        if len(outside)>1:return None
        assert 0<=rho<1;joint=None
        if outside:
            j=outside[0];sig=2-J[j][j];a=J[i][j];det=rho*sig-a*a
            assert det>=0 and 0<=sig<=1
            if not det:joint=dict(node=j,determinant='0',kind='positive_residual_relation',tail=0)
            else:
                invdiag=sig/det
                if invdiag<2:return None
                joint=dict(node=j,determinant=str(det),inverseGG=str(invdiag),kind='saturated_single_tail' if invdiag==2 else 'no_tail',tail=int(invdiag==2))
        records.append(dict(node=i,projection=str(J[i][i]),rho=str(rho),alone_tail_bound=rho//(1-rho),pair_checks=checks,double_checks=doubles,joint=joint))
    bound=len(c['marked_component'])+max((1+max(r['alone_tail_bound'],(r['joint'] or {}).get('tail',0)) for r in records),default=0)
    return dict(**c,method='relevant_component',relevant_nodes=good,contact_domains=cand,end_sets=ends,extension_cases=records,component_order_bound=bound)


def isolated_or_leaf(g,sc,J,c,du,dv):
    ends=next(split_ends(c['possible_contacts'],du,dv,J),None)
    if ends is None:return None
    ns=c['marked_component'];allowed=set(c['possible_contacts'])
    if not any(i in allowed and sc[i]*J[i][i]<=2 for i in ns):
        return dict(**c,method='marked_component_isolated_by_end_inverses',end_sets=ends,component_order_bound=len(ns))
    forbidden={i for i in ns if i not in allowed or sc[i]*J[i][i]>2}
    adj={i:[j for j in ns if j!=i and g[i][j]] for i in ns}
    leaves=[i for i in ns if len(adj[i])==1 and i in forbidden]
    simple=[i for i in leaves if adj[i][0] in forbidden and len(adj[adj[i][0]])==2 and g[i][adj[i][0]]**2==sc[i]*sc[adj[i][0]]]
    marked=[i for i in c['preferred_corner'] if i in leaves]
    p,q=c['preferred_corner'];distances=[H.graph_dist(H.extract(g,ns),[i]) for i in range(len(ns))]
    distance=distances[ns.index(p)][ns.index(q)]
    features=[i for i in ns if len(adj[i])>=3 or any(g[i][j]**2>sc[i]*sc[j] for j in adj[i])]
    near=min((distances[ns.index(x)][ns.index(f)] for x in (p,q) for f in features),default=999)
    if not(distance<=9 and marked and ((simple and len(ns)>=10) or near<=5)):return None
    return dict(**c,method='protected_leaf_or_ideal_feature',end_sets=ends,protected_marked_leaves=marked,
        protected_simple_leaves=simple,base_order=len(ns),distance=distance,feature_distance=near,
        forced_finite=bool(simple and len(ns)>=10),ideal_basic_pair=near<=5)


def run():
    stems=(X.C.TAG,X.TAG,H.CENSUS,H.TAG,'fresh_distance9_supporting_extensions_review_20261017')
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for s,d in data.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    diagrams={d['index']:d for d in data[H.CENSUS]['diagrams']};graphs={i:H.bond_graph(d) for i,d in diagrams.items()}
    inputs=[]
    for r in data[X.C.TAG]['cases']:
        b=r['base_context'];di,u=b['left'];dj,v=b['right'];native,mp=H.unite(graphs[di],graphs[dj],b['match'])
        inputs.append(dict(base_context=r['index'],extra_edges=[],native=native,left_marks=u,right_marks=[mp[i] for i in v]))
    for r in data[X.TAG]['nonempty_cross_assignments']:
        if r['status']=='survives_native_core_tests':
            inputs.append(dict(base_context=r['base_context'],extra_edges=r['edges'],native=[[F(x) for x in a] for a in r['native_graph']],left_marks=r['left_marks'],right_marks=r['right_marks']))
    assert len(inputs)==6703
    cores={};matrices={};rows=[];counts=Counter();bounds=Counter()
    @lru_cache(None)
    def corner(key,black,p,q):
        g,sc=matrices[key]
        return H.inspect_corner(g,sc,list(black),p,q)
    for source in inputs:
        native=source['native'];u=source['left_marks'];v=source['right_marks'];black=tuple(u+v)
        key=sha256(json.dumps(H.encode(native),separators=(',',':')).encode()).hexdigest()
        row={k:x for k,x in source.items() if k!='native'};row['core_key']=key
        if key not in cores:
            try:g,sc=H.native_to_rational(native)
            except AssertionError:
                proof=X.native_inertia(tuple(map(tuple,native)))
                cores[key]=dict(native_graph=H.encode(native),native_inertia=proof,method='surd_core')
            else:
                sig=H.inertia(g);zs=H.nullspace(g);iv=None if zs else H.inverse(g)
                cores[key]=dict(native_graph=H.encode(native),gram=H.encode(g),gauge=sc,signature=sig,
                    radical=[list(map(str,z)) for z in zs],inverse=H.encode(iv) if iv else None,method='rational_core')
                matrices[key]=(g,sc)
        core=cores[key]
        sig=core.get('signature',core.get('native_inertia',{}).get('signature'));assert sig[1]>=1
        if sig[1]>1:row['status']='excluded_multiple_negative_directions'
        elif core['method']=='surd_core':row['status']='unresolved_surd_supporting_extensions'
        else:
            g,sc=matrices[key];zs=[[F(x) for x in z] for z in core['radical']]
            iv=[[F(x) for x in r] for r in core['inverse']] if core['inverse'] else None
            du,dv=H.graph_dist(native,u),H.graph_dist(native,v);assert min(du[i] for i in v)>THRESHOLD
            white=sorted(set(range(len(g)))-set(black));bad=None
            for p in black:
                ids=sorted(white+[p]);s=H.inertia(H.extract(g,ids))
                if s!=[len(ids),0,0]:bad=dict(kind='actual_edge',nodes=ids,signature=s);break
            if bad:row.update(status='excluded_actual_principal',obstruction=bad)
            else:
                cs=[]
                for p,q in product(u,v):
                    ids=sorted(white+[p,q]);s=H.inertia(H.extract(g,ids))
                    if s[1]:bad=dict(kind='actual_corner',nodes=ids,signature=s);break
                    cs.append(corner(key,black,p,q))
                if bad:row.update(status='excluded_actual_principal',obstruction=bad)
                else:
                    opts=[]
                    for c in cs:
                        if c['marked_component'] is None:continue
                        a=ordinary(sc,iv,zs,c,du,dv)
                        if a:opts.append(a)
                        if iv and set(sc)=={1}:
                            a=relevant(iv,c,du,dv)
                            if a:opts.append(a)
                        if iv:
                            a=isolated_or_leaf(g,sc,iv,c,du,dv)
                            if a:opts.append(a)
                    row['all_corners']=cs
                    if opts:
                        chosen=min(opts,key=lambda c:(c.get('component_order_bound',100),len(c['possible_contacts']),c['method'],c['preferred_corner']))
                        row.update(status='structural_certificate',chosen=chosen)
                        counts['method_'+chosen['method']]+=1
                        if 'component_order_bound' in chosen:bounds[chosen['component_order_bound']]+=1
                    else:row['status']='unresolved_supporting_extensions'
        counts[row['status']]+=1;rows.append(row)
        if len(rows)%250==0:print('SUPPORT',len(rows),'/6703',dict(counts),flush=True)
    paths=[Path(__file__),Path(X.__file__),Path(H.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths=list(dict.fromkeys(paths))
    return dict(status='PASS_COMPLETE_DISTANCE7_HIGH_SUPPORTING_FRONTIER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        shortcut_threshold=THRESHOLD,total_contexts=len(rows),distinct_literal_cores=len(cores),counts=dict(counts),
        component_bounds={str(k):v for k,v in sorted(bounds.items())},cores=cores,cases=rows,
        scope='Every6691 zero-cross base context and12 surviving fully specified cross assignments tested. Exact full/principal inertia exclusions and supporting-root certificates use threshold7 explicitly, with no cutoff imposed on indirect tails. Unresolved cases remain listed. Structural bounds/protected features do not by themselves prove a retained finite/ideal payment; selected-map and all other branches remain separate. No complete finite7 quadrilateral theorem or new dimension bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['distinct_literal_cores'],out['counts'],out['component_bounds'],flush=True)
