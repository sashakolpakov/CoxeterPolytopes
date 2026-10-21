#!/usr/bin/env python3
"""Complete high-core overlaps under absence of a marked path of length <=8.

New exact overlap and supporting-extension certificate. No dimension theorem.
All disjoint and low-core domains are separate frozen inputs.
"""
import argparse,json,hashlib
from collections import Counter
from functools import lru_cache
from itertools import combinations,product,permutations
from math import comb
from pathlib import Path
import fresh_selective_distance9_dimension108_20261017_review as R
from fresh_selective_distance10_dimension116_20261017_review import (F, extract, inverse, nullspace, psd, graph_dist, bond_graph, all_contexts, unite, obstruction, inertia, native_to_rational, combined, inspect_corner, encode)
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_high_overlap_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
NINE='fresh_distance8_adaptive_leaf_overlap_20261017'
def edge_bad_native(native,black):
    white=[i for i in range(len(native))if i not in black]
    for mark in black:
        ids=sorted(white+[mark]);principal=extract(native,ids)
        try:a,gauge=native_to_rational(principal);mode='positive_diagonal_congruence'
        except AssertionError:
            a=[[F(2 if i==j else -int(bool(principal[i][j])))for j in range(len(ids))]for i in range(len(ids))]
            gauge=[1]*len(ids);mode='unit_bond_relaxation'
        o=obstruction(a)
        if o:return dict(mark=mark,nodes=ids,gauge_squared=gauge,comparison=mode,**o)
    return None


def classify(diagrams):
    ds={d['index']:d for d in diagrams};graphs={i:bond_graph(d)for i,d in ds.items()};contexts=all_contexts(diagrams)
    assert max(max(r['distances'])for r in contexts)==7
    @lru_cache(None)
    def diagonals(di,black,node):
        a=graphs[di];out=[]
        for b in black:
            ids=[i for i in range(len(a))if i!=b];h,scale=native_to_rational(extract(a,ids));psd(h,True)
            j=inverse(h);k=ids.index(node);out.append(scale[k]*j[k][k])
        return tuple(out)
    @lru_cache(None)
    def white_diagonal(di,black,node):
        a=graphs[di];ids=[i for i in range(len(a))if i not in black]
        h,scale=native_to_rational(extract(a,ids));psd(h,True);iv=inverse(h);k=ids.index(node)
        return scale[k]*iv[k][k]
    high=[r for r in contexts if max(r['distances'])>=2]
    assert all(sum(x>=2 for x in r['distances'])<=7 for r in high)
    records=[];candidate_count=0
    for left,right in product(high,repeat=2):
        ga,gb=graphs[left['diagram']],graphs[right['diagram']]
        choices=[(i,j)for i,d in enumerate(left['distances'])for j,e in enumerate(right['distances'])if d+e>=9]
        size_bound=min(sum(x>=2 for x in left['distances']),sum(x>=2 for x in right['distances']))
        candidate_count+=sum(comb(len(choices),k)for k in range(1,min(size_bound,len(choices))+1))
        ls=sorted({i for i,j in choices});rs=sorted({j for i,j in choices})
        for size in range(1,min(size_bound,len(ls),len(rs))+1):
            for aa,bb in product(combinations(ls,size),combinations(rs,size)):
                for pp in permutations(bb):
                    match=tuple(zip(aa,pp))
                    if any(pair not in choices for pair in match):continue
                    if any(ga[i][k]!=gb[j][l]for i,j in match for k,l in match):continue
                    row=dict(left=[left['diagram'],list(left['black'])],right=[right['diagram'],list(right['black'])],match=[list(x)for x in match])
                    rejects=[]
                    for i,j in match:
                        for side,x,y,rr,other,shared in [('left',i,j,left,gb,set(pp)),('right',j,i,right,ga,set(aa))]:
                            ext=[k for k in range(len(other))if k!=y and other[y][k]and k not in shared]
                            if ext and max(diagonals(rr['diagram'],rr['black'],x))>=2:
                                assert all(k not in(right['black']if side=='left'else left['black'])for k in ext)
                                rejects.append([side,x,ext])
                    if rejects:row.update(kind='inverse_rejection',inverse_rejections=rejects)
                    else:
                        native,mp=unite(ga,gb,match)
                        u=list(left['black']);v=[mp[i]for i in right['black']];bad=edge_bad_native(native,u+v)
                        if bad:row.update(kind='edge_rejection',edge_obstruction=bad)
                        else:row['kind']=str((left['diagram'],right['diagram'],match))
                    records.append(row)
    counts=Counter(r['kind'] for r in records)
    assert candidate_count==6390952 and len(records)==82197
    assert counts['inverse_rejection']==75947 and counts['edge_rejection']==3636
    assert sum(v for k,v in counts.items() if k not in ('inverse_rejection','edge_rejection'))==2614
    print('OVERLAP_DOMAIN_PASS',len(records),candidate_count,flush=True)
    return dict(marked_contexts=len(contexts),candidate_matchings=candidate_count,
      consistent_matchings=len(records),overlap_counts=dict(counts),overlap_records=records)

def split_ends(allowed,du,dv,iv):
    # Enumerate every distance split rather than assume symmetric end radii.
    for k in range(7):
        left=[i for i in allowed if du[i]<=k];right=[i for i in allowed if dv[i]<=6-k]
        if set(left)|set(right)!=set(allowed)or set(left)&set(right):continue
        if all(iv[i][j]>0 for block in(left,right)for i in block for j in block):yield [left,right]


def ordinary(gauge,iv,radical,c,du,dv):
    allowed=c['possible_contacts'];size=len(c['marked_component'])
    if radical:
        z=radical[0]
        if any(not z[i]for i in allowed):return None
        pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
        ds=[max((du[i]for i in pos),default=0)+max((dv[i]for i in neg),default=0)+2,
            max((dv[i]for i in pos),default=0)+max((du[i]for i in neg),default=0)+2]
        if min(ds)>8:return None
        return dict(**c,method='radical',shortcut_bound=min(ds),component_order_bound=size)
    ends=next(split_ends(allowed,du,dv,iv),None)
    if ends is None:return None
    projection={i:gauge[i]*iv[i][i]for i in allowed};candidates=[i for i in allowed if projection[i]<=2]
    if any(projection[i]<=1 for i in candidates):return None
    # min(gauge_i,gauge_j) underestimates sqrt(gauge_i*gauge_j), safely.
    if any(projection[i]+projection[j]+2*min(gauge[i],gauge[j])*iv[i][j]<=2 for end in ends for i,j in combinations(end,2)):return None
    checks=[]
    for i,j in product(candidates,repeat=2):
        a=iv[i][j];rho=(2-projection[i])*(2-projection[j]);sq=gauge[i]*gauge[j]*a*a
        ds=min(du[i]+dv[j],du[j]+dv[i])+3
        if sq<=rho or(a<0 and ds>8):return None
        checks.append(dict(nodes=[i,j],projection_square=str(sq),residual_product=str(rho),sign=1 if a>0 else-1,shortcut=ds if a<0 else None))
    roots=[]
    for i in candidates:
        rho=2-projection[i];assert 0<=rho<1
        tail=rho//(1-rho)
        roots.append(dict(contact_node=i,projection_square=str(projection[i]),residual_square=str(rho),maximum_attached_orthogonal_rank=tail))
    bound=size+max((1+r['maximum_attached_orthogonal_rank']for r in roots),default=0)
    if bound>17:return None
    return dict(**c,method='lorentz_residual',end_sets=ends,single_contact_candidates=roots,pair_tests=checks,component_order_bound=bound)


def relevant(iv,c,du,dv):
    allowed=c['possible_contacts'];ends=next(split_ends(allowed,du,dv,iv),None)
    if ends is None:return None
    if any(iv[i][i]+iv[j][j]+2*iv[i][j]<=2 for end in ends for i,j in combinations(end,2)):return None
    candidates=[i for i in allowed if iv[i][i]<=2]
    assert all(iv[i][i]>=1 for i in candidates)
    good=[i for i in candidates if i in c['marked_component']]
    if any(iv[i][i]<=1 for i in good):return None
    records=[]
    for i in good:
        rho=2-iv[i][i];outside=[];checks=[]
        for j in candidates:
            sig=2-iv[j][j];a=iv[i][j];ds=min(du[i]+dv[j],du[j]+dv[i])+3
            if a*a>rho*sig and a>0:reason='positive'
            elif a*a>rho*sig and a<0 and ds<=8:reason='short_path'
            elif a>0 and (a+1)**2>rho*sig:
                if j in good:return None
                outside.append(j);reason='forced_orthogonal'
            else:return None
            checks.append(dict(node=j,reason=reason,projection=str(a),residual_product=str(rho*sig),shortcut=ds))
        doubles=[]
        for j in candidates:
            if iv[j][j]!=1:continue
            a=iv[i][j];ds=min(du[i]+dv[j],du[j]+dv[i])+3
            if not(a>0 or(a<0 and ds<=8)):return None
            doubles.append(dict(node=j,sign=1 if a>0 else-1,shortcut=ds))
        if len(outside)>1:return None
        assert 0<=rho<1;alone=rho//(1-rho);joint=None
        if outside:
            j=outside[0];sig=2-iv[j][j];a=iv[i][j];det=rho*sig-a*a
            assert det>=0 and 0<=sig<=1
            if not det:joint=dict(node=j,determinant='0',kind='positive_residual_relation',tail=0)
            else:
                invdiag=sig/det
                if invdiag<2:return None
                joint=dict(node=j,determinant=str(det),inverseGG=str(invdiag),kind='saturated_single_tail'if invdiag==2 else'no_tail',tail=int(invdiag==2))
        records.append(dict(node=i,projection=str(iv[i][i]),rho=str(rho),alone_tail_bound=alone,pair_checks=checks,double_checks=doubles,joint=joint))
    bound=len(c['marked_component'])+max((1+max(r['alone_tail_bound'],(r['joint']or{}).get('tail',0))for r in records),default=0)
    if bound>17:return None
    return dict(**c,method='relevant_component',relevant_nodes=good,contact_domains=candidates,end_sets=ends,extension_cases=records,component_order_bound=bound)


def isolated_component(gauge,iv,c,du,dv):
    ends=next(split_ends(c['possible_contacts'],du,dv,iv),None)
    if ends is None:return None
    relevant=[i for i in c['possible_contacts']if i in c['marked_component']and gauge[i]*iv[i][i]<=2]
    if relevant:return None
    return dict(**c,method='marked_component_isolated_by_end_inverses',end_sets=ends,
                component_order_bound=len(c['marked_component']))


def protected_leaf(g,gauge,iv,c,du,dv):
    ends=next(split_ends(c['possible_contacts'],du,dv,iv),None)
    if ends is None:return None
    ns=c['marked_component'];allowed=set(c['possible_contacts'])
    forbidden={i for i in ns if i not in allowed or gauge[i]*iv[i][i]>2}
    graph={i:[j for j in ns if j!=i and g[i][j]]for i in ns}
    leaves=[i for i in ns if len(graph[i])==1 and i in forbidden]
    simple=[i for i in leaves if graph[i][0]in forbidden and len(graph[graph[i][0]])==2
            and g[i][graph[i][0]]**2==gauge[i]*gauge[graph[i][0]]]
    marked=[i for i in c['preferred_corner']if i in leaves]
    p,q=c['preferred_corner'];ds=[graph_dist(extract(g,ns),[i])for i in range(len(ns))]
    distance=ds[ns.index(p)][ns.index(q)]
    features=[i for i in ns if len(graph[i])>=3 or any(g[i][j]**2>gauge[i]*gauge[j]for j in graph[i])]
    near=min((ds[ns.index(x)][ns.index(f)]for x in(p,q)for f in features),default=999)
    if not(distance<=9 and marked and ((simple and len(ns)>=10)or near<=5)):return None
    return dict(**c,method='finite_leaf9_or_ideal_feature5',end_sets=ends,protected_marked_leaves=marked,
        protected_simple_leaves=simple,base_order=len(ns),distance=distance,feature_distance=near,
        forced_finite=bool(simple and len(ns)>=10),ideal_basic_pair=near<=5,
        conclusion='finite marked leaf at distance9, or ideal distance9 with a persistent feature within5')


def isolation(diagrams,data):
    ds={d['index']:d for d in diagrams};nine=json.loads((ROOT/'results'/f'{NINE}.json').read_text())
    def k(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['match'])))
    nines={k(r):r for r in nine['cases']};cores={};rows=[];counts=Counter();short_total=pd_total=0
    for row in data['overlap_records']:
        if row['kind']in('inverse_rejection','edge_rejection'):continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']));key=str((di,dj,match))
        native,g,sc,mp=combined(ds,di,dj,match);v=[mp[i]for i in v0]
        du,dv=graph_dist(native,u),graph_dist(native,v);assert min(du[x]for x in v)>8
        if key not in cores:
            zs=nullspace(g);iv=None if zs else inverse(g);sig=inertia(g);assert sig[1]==1
            cores[key]=dict(gram=encode(g),gauge=sc,signature=sig,radical=[[str(x)for x in z]for z in zs],inverse=encode(iv)if iv else None)
        core=cores[key];zs=[[F(x)for x in z]for z in core['radical']];iv=[[F(x)for x in a]for a in core['inverse']]if core['inverse']else None
        corners=[inspect_corner(g,sc,u+v,p,q)for p,q in product(u,v)];opts=[]
        for c in corners:
            if c['marked_component'] is None:continue
            x=ordinary(sc,iv,zs,c,du,dv)
            if x:opts.append(x)
        if not opts and not zs and set(sc)=={1}:
            for c in corners:
                if c['marked_component'] is None:continue
                try:x=relevant(iv,c,du,dv)
                except AssertionError:x=None
                if x:opts.append(x)
        if not opts and iv:
            for c in corners:
                if c['marked_component'] is None:continue
                x=isolated_component(sc,iv,c,du,dv)
                if x:opts.append(x)
        if not opts and iv:
            for c in corners:
                if c['marked_component'] is None:continue
                x=protected_leaf(g,sc,iv,c,du,dv)
                if x:opts.append(x)
        if opts:
            chosen=min(opts,key=lambda c:(c.get('component_order_bound',100),len(c['possible_contacts']),c['preferred_corner']))
        else:
            assert k(row)in nines,(row['left'],row['right'],match)
            n=nines[k(row)];assert n['gram']==encode(g)and n['inverse']==encode(iv)and n['all_corners']==corners
            chosen=dict(method='nine_context_adaptive_finite_payment',local_case=n['index'],payment=n['payment'])
        # A minimal simple extra bond is enough: its nonnegative Z witness survives every larger native bond.
        left=[i for i in range(ds[di]['order'])if i not in{a for a,b in match}]
        right=[mp[i]for i in range(ds[dj]['order'])if i not in{b for a,b in match}]
        shorts=0;obs=[]
        for i,j in product(left,right):
            if native[i][j]:continue
            h=[a.copy()for a in native];h[i][j]=h[j][i]=-1
            if min(graph_dist(h,u)[x]for x in v)<=8:shorts+=1;continue
            proof=edge_bad_native(h,u+v);assert proof,(di,u,dj,v0,match,i,j)
            obs.append(dict(edge=[i,j],**proof))
        counts[chosen['method']]+=1;short_total+=shorts;pd_total+=len(obs)
        rows.append(dict(index=len(rows),left=row['left'],right=row['right'],match=row['match'],right_copy_map=mp,
            core_key=key,all_corners=corners,extra_cross_shortcuts=shorts,extra_cross_pd_rejections=obs,**chosen))
        if len(rows)%100==0:print('OVERLAP_EXTENSIONS',len(rows),dict(counts),flush=True)
    assert len(rows)==2614 and len(cores)==154
    assert counts=={'radical':1329,'lorentz_residual':823,'relevant_component':259,
                   'marked_component_isolated_by_end_inverses':120,'finite_leaf9_or_ideal_feature5':74,
                   'nine_context_adaptive_finite_payment':9},counts
    assert short_total==154502 and pd_total==5528
    return dict(cores=cores,cases=rows,method_counts=dict(counts),extra_cross_shortcuts=short_total,
                extra_cross_pd_rejections=pd_total,
                scope='All2614 contexts have a finite8/leaf9/whole17 payment or an ideal whole-component bound<=17 / basic ideal5-feature5 payment. Exact ideal pair compression is a separate census.')


def run():
    diagrams=json.loads((ROOT/'results'/f'{CENSUS}.json').read_text())['diagrams']
    domain=classify(diagrams);extensions=isolation(diagrams,domain)
    stems=[CENSUS,NINE,'fresh_distance8_high_disjoint_20261017','fresh_distance8_low_core_20261017',
           'fresh_distance9_supporting_extensions_review_20261017','fresh_finite8_a17_d11_coupon_budget_20261017']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    paths += [ROOT/'scripts/fresh_selective_distance9_dimension108_20261017_review.py',ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    return dict(status='PASS_COMPLETE_HIGH_DISTANCE8_OVERLAP_DOMAIN_AND_EXTENSIONS',census=domain,isolation=extensions,
       sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
       scope='No dimension bound. Ideal components of order14..17 require a new exact selected-pair census before importing a compressed ideal budget.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],flush=True)
