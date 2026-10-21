#!/usr/bin/env python3
"""Independent Fraction reconstruction of distance9 high cores and budgets.

Uses only the preceding independent Fraction elementary linear algebra;
imports no producer arithmetic or symbolic/numerical algebra package.
The complete low-order domain is bound as a separately replayed certificate.
"""
from fresh_selective_distance10_dimension116_20261017_review import (F, matrix, extract, inverse, nullspace, psd, graph_dist, rational_gram, bond_graph, all_contexts, overlap_key, unite, positive_kernel, obstruction, inertia, joined, native_to_rational, combined, pieces, inspect_corner, polynomial, sha, encode, square)
import argparse,json
from collections import Counter
from itertools import combinations,product,permutations
from functools import lru_cache
from math import comb,prod
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_selective_distance9_dimension108_20261017_review'
BASE='fresh_selective_distance9_dimension108_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
LOW='fresh_distance9_low_core_20261017'
LOCAL='fresh_distance9_supporting_extensions_review_20261017'


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


def overlap_type(di,dj,match):return str((di,dj,match))

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
    ports=[]
    for r in contexts:
        for node,distance in enumerate(r['distances']):
            if distance<2:continue
            vals=diagonals(r['diagram'],r['black'],node)
            ports.append(dict(diagram=r['diagram'],black=list(r['black']),node=node,distance=distance,
                              inverse_diagonals=list(map(str,vals)),allowed=max(vals)<2))
    allowed=[p for p in ports if p['allowed']];assert len(ports)==1262 and len(allowed)==378
    joins=[];good={};types={}
    for a,b in product(allowed,repeat=2):
        if a['distance']+b['distance']<9:continue
        di,dj=a['diagram'],b['diagram'];u,v=a['black'],b['black'];p,q=a['node'],b['node']
        native,mp=unite(graphs[di],graphs[dj],());native[p][mp[q]]=native[mp[q]][p]=-1;black=u+[mp[x]for x in v]
        white=[i for i in range(len(native))if i not in black];wg,wsc=native_to_rational(extract(native,white));bad=obstruction(wg)
        row=dict(left=[di,u,p],right=[dj,v,q],white=white)
        product_inverse=white_diagonal(di,tuple(u),p)*white_diagonal(dj,tuple(v),q)
        assert bool(bad)==(product_inverse>=1)
        if bad:row['white_obstruction']=bad
        else:
            key=(di,tuple(u),dj,tuple(v));good.setdefault(key,[]).append((p,q))
            maxima=[max(map(F,x['inverse_diagonals']))for x in(a,b)]
            assert max(maxima)>=1 or 2*product_inverse>=1
            row['larger_bond_obstruction']=dict(endpoint_maxima=list(map(str,maxima)),twice_white_inverse_product=str(2*product_inverse))
            typ=tuple(sorted(((di,p),(dj,q))))
            if typ not in types:
                (aa,bb),(cc,dd)=typ;h,gauge=joined(ds[aa],ds[cc],bb,dd)
                sig=inertia(h);assert sig[1]==2
                types[typ]=dict(ports=list(map(list,typ)),inertia=sig,gram=encode(h),gauge=gauge)
        joins.append(row)
    multiple=[]
    for (di,u,dj,v),edges in sorted(good.items()):
        native,mp=unite(graphs[di],graphs[dj],());black=list(u)+[mp[i]for i in v]
        for size in range(2,len(edges)+1):
            for chosen in combinations(edges,size):
                h=[r.copy()for r in native]
                for i,j in chosen:h[i][mp[j]]=h[mp[j]][i]=-1
                ob=edge_bad_native(h,black);assert ob
                multiple.append(dict(left=[di,list(u)],right=[dj,list(v)],edges=list(map(list,chosen)),**ob))
    assert len(joins)==4442 and sum('white_obstruction'in r for r in joins)==1103
    assert Counter(map(len,good.values()))=={1:2007,2:348,3:180,4:24}
    assert len(types)==85 and len(multiple)==1332
    print('HIGH_DISJOINT_PASS',len(joins),len(multiple),flush=True)
    high=[r for r in contexts if max(r['distances'])>=3]
    assert all(sum(x>=3 for x in r['distances'])<=5 for r in high)
    records=[];candidate_count=0
    for left,right in product(high,repeat=2):
        ga,gb=graphs[left['diagram']],graphs[right['diagram']]
        choices=[(i,j)for i,d in enumerate(left['distances'])for j,e in enumerate(right['distances'])if d+e>=10]
        size_bound=min(sum(x>=3 for x in left['distances']),sum(x>=3 for x in right['distances']))
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
                        else:row['kind']=overlap_type(left['diagram'],right['diagram'],match)
                    records.append(row)
    counts=Counter(r['kind']for r in records)
    assert candidate_count==312316 and len(records)==20718,(candidate_count,len(records),counts)
    assert counts['inverse_rejection']==18721 and counts['edge_rejection']==1048
    assert sum(v for k,v in counts.items()if k not in('inverse_rejection','edge_rejection'))==949
    print('OVERLAP_PASS',dict(counts),flush=True)
    return dict(marked_contexts=384,ports=ports,allowed_ports=378,disjoint_joins=joins,
        surviving_disjoint_types=list(types.values()),multiple_edge_obstructions=multiple,
        candidate_matchings=candidate_count,consistent_matchings=len(records),overlap_counts=dict(counts),overlap_records=records)


def split_ends(allowed,du,dv,iv):
    # Enumerate every distance split rather than assume symmetric end radii.
    for k in range(8):
        left=[i for i in allowed if du[i]<=k];right=[i for i in allowed if dv[i]<=7-k]
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
        if min(ds)>9:return None
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
        if sq<=rho or(a<0 and ds>9):return None
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
            elif a*a>rho*sig and a<0 and ds<=9:reason='short_path'
            elif a>0 and (a+1)**2>rho*sig:
                if j in good:return None
                outside.append(j);reason='forced_orthogonal'
            else:return None
            checks.append(dict(node=j,reason=reason,projection=str(a),residual_product=str(rho*sig),shortcut=ds))
        doubles=[]
        for j in candidates:
            if iv[j][j]!=1:continue
            a=iv[i][j];ds=min(du[i]+dv[j],du[j]+dv[i])+3
            if not(a>0 or(a<0 and ds<=9)):return None
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


def isolation_checks(diagrams,data):
    ds={d['index']:d for d in diagrams};out=[];cores={};total_short=total_pd=0
    for nr,row in enumerate(data['overlap_records']):
        if row['kind']in('inverse_rejection','edge_rejection'):continue
        di,u=row['left'];dj,v0=row['right'];match=tuple(map(tuple,row['match']))
        native,g,sc,mp=combined(ds,di,dj,match);v=[mp[i]for i in v0]
        du,dv=graph_dist(native,u),graph_dist(native,v);assert min(du[i]for i in v)>9
        key=str((di,dj,match))
        if key not in cores:
            ks=nullspace(g);iv=None if ks else inverse(g);sig=inertia(g);assert sig[1]==1
            cores[key]=dict(gram=encode(g),gauge=sc,signature=sig,radical=[list(map(str,z))for z in ks],inverse=encode(iv)if iv else None)
        saved=cores[key];ks=[[F(x)for x in z]for z in saved['radical']];iv=[[F(x)for x in r]for r in saved['inverse']]if saved['inverse']else None
        choices=[];allcorners=[]
        for p,q in product(u,v):
            c=inspect_corner(g,sc,u+v,p,q);allcorners.append(c)
            if c['marked_component']is None:continue
            option=ordinary(sc,iv,ks,c,du,dv)
            if option:choices.append(option)
        if not choices:
            assert not ks and set(sc)=={1}
            for c in allcorners:
                if c['marked_component']is None:continue
                option=relevant(iv,c,du,dv)
                if option:choices.append(option)
        assert choices,(di,u,dj,v0,match)
        if choices[0]['method']=='relevant_component':chosen=min(choices,key=lambda c:(c['component_order_bound'],len(c['extension_cases']),c['preferred_corner']))
        else:chosen=min(choices,key=lambda c:(c['component_order_bound'],len(c['possible_contacts'])))
        left=[i for i in range(ds[di]['order'])if i not in{x for x,y in match}]
        right=[mp[i]for i in range(ds[dj]['order'])if i not in{y for x,y in match}]
        shorts=0;obs=[]
        for i,j in product(left,right):
            if native[i][j]:continue
            h=[r.copy()for r in native];h[i][j]=h[j][i]=-1
            if min(graph_dist(h,u)[k]for k in v)<=9:shorts+=1;continue
            proof=edge_bad_native(h,u+v);assert proof,(di,dj,match,i,j)
            obs.append(dict(edge=[i,j],**proof))
        total_short+=shorts;total_pd+=len(obs)
        out.append(dict(type=row['kind'],left=[di,u],right=[dj,v0],match=row['match'],right_copy_map=mp,
            all_corners=allcorners,extra_cross_shortcuts=shorts,extra_cross_pd_rejections=obs,**chosen))
        if len(out)%100==0:print('INDEPENDENT_CORNERS',len(out),flush=True)
    assert len(out)==949 and len(cores)==54 and total_short==60831 and total_pd==888
    assert Counter(r['method']for r in out)['relevant_component']==69
    assert max(r['component_order_bound']for r in out)==17
    print('ISOLATION_PASS',949,54,total_short,total_pd,flush=True)
    return dict(cores=cores,preferred_corners=out,extra_cross_shortcuts=total_short,extra_cross_pd_rejections=total_pd)


def profile(ranks):
    d=sum(ranks);h=polynomial(ranks);v=prod(r+1 for r in ranks)
    assert h==h[::-1]and sum(h)==v
    g=sum(x*min(i,d-i)for i,x in enumerate(h));curvature=F(5*g-d*v,2*v)
    moment=0;s=len(ranks)
    for k in range(s+1):
        for subset in combinations(ranks,k):
            shift=sum(x+1 for x in subset);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**k*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curvature==F(3*d,4)-F(5*moment,2*v)
    cost=sum((min(F(9),F(r,2))for r in ranks),F())
    return dict(ranks=list(ranks),normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))


def budgets():
    rows=[];mins=[]
    for d in range(108,116):
        shapes=set()
        for k in range(d//18+1):
            for low in range(18):
                high=d-18*k-low
                if high==0 or high>=19:shapes.add(tuple(sorted([18]*k+([low]if low else[])+([high]if high else[]))))
        group=[profile(rs)for rs in sorted(shapes)];assert all(F(r['margin'])>0 for r in group)
        rows+=group;n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8);assert beta>9*n
        mins.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=str(min(F(r['margin'])for r in group)),finite_slack=str(beta-9*n)))
    assert len(rows)==756 and min(F(r['margin'])for r in rows)==F(45873,2476099)
    boundary=profile((18,)*5+(17,));assert F(boundary['margin'])==F(-4659,260642)
    coeff=[]
    for start in(79,80):
        f=lambda t:((start+2*t-1)*(start+2*t-5)if start%2 else(start+2*t)*(start+2*t-6))-72*(start+2*t+1)
        vals=[f(i)for i in range(3)];a=(vals[2]-2*vals[1]+vals[0])//2;b=vals[1]-vals[0]-a
        coeff.append([vals[0],b,a]);assert all(x>0 for x in coeff[-1])
    assert coeff==[[12,160,4],[88,164,4]]
    return dict(profiles=rows,minima=mins,boundary=boundary,compact_facet_slack_coefficients=coeff)

def check_witness(native,record):
    ids=record['nodes'];h=extract(native,ids)
    if record.get('comparison')=='unit_bond_relaxation':
        a=[[F(2 if i==j else-int(bool(h[i][j])))for j in range(len(h))]for i in range(len(h))];sc=[1]*len(h)
    else:a,sc=native_to_rational(h)
    assert list(map(F,record['gauge_squared']))==list(map(F,sc))
    vec=list(map(F,record['vector']));assert len(vec)==len(a)and any(vec)and all(x>=0 for x in vec)
    q=square(a,vec);assert q<=0 and q==F(record['square'])
    return 1


def compare_saved(diagrams,data,iso,bs,p):
    ds={d['index']:d for d in diagrams};graphs={i:bond_graph(d)for i,d in ds.items()}
    assert data['ports']==p['census']['ports']
    ours={overlap_key(r):r for r in data['overlap_records']};theirs={overlap_key(r):r for r in p['census']['overlap_records']}
    assert set(ours)==set(theirs);witnesses=0
    for key,row in theirs.items():
        kind='inverse_rejection'if'inverse_rejection'in row else'edge_rejection'if'edge_pd_rejection'in row else row['surviving_type']
        assert ours[key]['kind']==kind
        if 'edge_pd_rejection'in row:
            native,_=unite(graphs[key[0]],graphs[key[2]],key[4]);witnesses+=check_witness(native,row['edge_pd_rejection'])
    def join_key(r):return(r['left'][0],tuple(r['left'][1]),r['left'][2],r['right'][0],tuple(r['right'][1]),r['right'][2])
    ours_join={join_key(r):r for r in data['disjoint_joins']};theirs_join={join_key(r):r for r in p['disjoint_obstructions']['rows']}
    assert set(ours_join)==set(theirs_join)
    for key,row in theirs_join.items():
        assert ('white_obstruction'in ours_join[key])==('white_rejection'in row)
        if 'white_rejection'in row:
            native,mp=unite(graphs[key[0]],graphs[key[3]],());native[key[2]][mp[key[5]]]=native[mp[key[5]]][key[2]]=-1
            witnesses+=check_witness(native,row['white_rejection'])
        else:assert ours_join[key]['larger_bond_obstruction']==row['larger_bond_obstruction']
    def multi_key(r):return(r['left'][0],tuple(r['left'][1]),r['right'][0],tuple(r['right'][1]),tuple(map(tuple,r['edges'])))
    ownmulti={multi_key(r)for r in data['multiple_edge_obstructions']};assert len(ownmulti)==1332
    assert ownmulti=={multi_key(r)for r in p['disjoint_obstructions']['multiple_edge_obstructions']}
    for row in p['disjoint_obstructions']['multiple_edge_obstructions']:
        di,u=row['left'];dj,v=row['right'];native,mp=unite(graphs[di],graphs[dj],())
        for i,j in row['edges']:native[i][mp[j]]=native[mp[j]][i]=-1
        witnesses+=check_witness(native,row['edge_pd_rejection'])
    types={tuple(map(tuple,r['ports'])):r for r in data['surviving_disjoint_types']}
    for row in p['disjoint_obstructions']['underlying_inertia_types']:
        own=types[tuple(map(tuple,row['ports']))];assert own['inertia']==row['signature']
        assert own['gram']==row['gram']and own['gauge']==row['gauge_squared']
    owniso={overlap_key(r):r for r in iso['preferred_corners']};assert len(owniso)==949
    for row in p['isolation']['cases']:
        key=overlap_key(row);own=owniso[key]
        for field in('preferred_corner','omitted_black','blocks','possible_contacts','marked_component','method','component_order_bound'):
            assert own[field]==row[field],(key,field,own[field],row[field])
        if row['method']=='lorentz_residual':
            for a,b in zip(own['single_contact_candidates'],row['single_contact_candidates']):
                assert all(a[f]==b[f]for f in a)
        if row['method']=='relevant_component':
            assert own['relevant_nodes']==row['relevant_nodes']
            for a,b in zip(own['extension_cases'],row['extension_cases']):
                assert all(a[f]==b[f]for f in('node','projection','rho','alone_tail_bound'))
                aa=a['joint'];bb=b['coexisting_outside'];assert bool(aa)==bool(bb)
                if aa:
                    assert aa['node']==bb['outside_node']
                    for f in('determinant','kind','tail'):assert aa[f]==bb[f]
                    if 'inverseGG'in aa:assert aa['inverseGG']==bb['inverseGG']
        native,mp=unite(graphs[key[0]],graphs[key[2]],key[4]);cross=row['extra_cross_bond_checks']
        assert own['extra_cross_shortcuts']==cross['shortcut_rejections']
        assert {tuple(r['edge'])for r in own['extra_cross_pd_rejections']}=={tuple(r['edge'])for r in cross['edge_pd_z_rejections']}
        for certificate in cross['edge_pd_z_rejections']:
            h=[r.copy()for r in native];i,j=certificate['edge'];h[i][j]=h[j][i]=-1
            witnesses+=check_witness(h,certificate)
    for key,row in p['isolation']['cores'].items():
        own=iso['cores'][key];assert own['gram']==row['gram']and own['gauge']==row['gauge_squared']and own['signature']==row['signature']
        if 'inverse'in row:assert own['inverse']==row['inverse']
        if 'radical'in row:
            g=[[F(x)for x in r]for r in row['gram']];z=list(map(F,row['radical']));assert any(z)and all(sum(a*b for a,b in zip(r,z))==0 for r in g)
    ownbudget={tuple(sorted(r['ranks'])):r for r in bs['profiles']}
    for row in p['reduced_profiles']:
        own=ownbudget[tuple(sorted(row['ranks']))]
        assert all(own[k]==row[k]for k in('normalized_curvature','normalized_cost','margin'))
    assert bs['minima']==p['dimension_minima']and bs['compact_facet_slack_coefficients']==p['compact_facet_parity_slack_times8']
    assert witnesses==4371,witnesses
    return dict(ports=1262,joins=4442,multiple_edge_subsets=1332,consistent_overlaps=20718,
        complete_isolation_cases=949,critical_native_cores=54,saved_nonnegative_witnesses=witnesses,reduced_profiles=756)


def compare_local(iso,local):
    rows={overlap_key(r):r for r in iso['preferred_corners']if r['method']=='relevant_component'}
    other={overlap_key(r):r for r in local['contexts']};assert set(rows)==set(other)and len(rows)==69
    for key,a in rows.items():
        b=other[key]
        for f in('preferred_corner','omitted_black','possible_contacts','marked_component','end_sets','contact_domains','relevant_nodes'):assert a[f]==b[f],(key,f)
        assert a['component_order_bound']==b['component_bound']
        assert len(a['blocks'])==len(b['blocks'])
        for x,y in zip(a['blocks'],b['blocks']):
            assert x['nodes']==y['nodes']and x['kind']==y['kind']
            if x['kind']=='finite':assert x['normalized_inverse_diagonals']==y['inverse_diagonals']
        assert len(a['extension_cases'])==len(b['extensions'])
        for x,y in zip(a['extension_cases'],b['extensions']):
            for f in('node','projection','rho','alone_tail_bound'):assert x[f]==y[f]
            j=x['joint'];k=y['coexisting_outside'];assert bool(j)==bool(k)
            if j:
                assert j['node']==k['outside_node']and j['determinant']==k['determinant']and j['tail']==k['tail_bound']and j['kind']==k['kind']
                if 'inverseGG'in j:assert j['inverseGG']==k['inverse'][0][0]
    return dict(contexts=69,relevant_root_cases=sum(len(r['extension_cases'])for r in rows.values()),max_component=max(r['component_order_bound']for r in rows.values()))


def run(prepare=False):
    diagrams=json.loads((ROOT/f'results/{CENSUS}.json').read_text())['diagrams']
    data=classify(diagrams);iso=isolation_checks(diagrams,data);bs=budgets()
    local=json.loads((ROOT/f'results/{LOCAL}.json').read_text());low=json.loads((ROOT/f'results/{LOW}.json').read_text())
    for cert in(local,low):
        assert cert['status'].startswith('PASS')
        for path,h in cert['sources'].items():assert sha(ROOT/path)==h,(path,h)
    local_comparison=compare_local(iso,local)
    result=dict(status='PASS_INDEPENDENT_FRACTION_DISTANCE9_PREPARED',classification=data,isolation=iso,budgets=bs,local_extension_comparison=local_comparison)
    if prepare:return result
    producer=json.loads((ROOT/f'results/{BASE}.json').read_text());comparison=compare_saved(diagrams,data,iso,bs,producer)
    for path,h in producer['sources'].items():assert sha(ROOT/path)==h,path
    paths=[f'scripts/{TAG}.py','scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    for stem in(BASE,CENSUS,LOCAL,LOW):paths += [f'scripts/{stem}.py',f'results/{stem}.json',f'results/{stem}.md']
    paths += ['results/fresh_cusp_link_shape_curvature_20261017.md','results/fresh_selective_distance10_dimension116_20261017.md']
    result.update(status='PASS_INDEPENDENT_FRACTION_SELECTIVE_DISTANCE9_DIMENSION108_REVIEW',sources={path:sha(ROOT/path)for path in paths},
        producer_source_bindings=len(producer['sources']),saved_comparison=comparison,dimension_bound=108,compact_facet_ambient_bound=79,
        low_core_certificate=LOW,supporting_extension_certificate=LOCAL,
        scope='Independent Fraction high-catalogue, joins, all overlap and extension-scalar reconstruction, saved Z-witness checks and curvature arithmetic. Separately frozen low-core and69-context extension proofs are bound; their full computations are not rerun here. Mathematical face/cusp/classification/convex-order inputs remain explicit dependencies. No external refereeing, formalization, priority or construction claim.')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--prepare',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run(args.prepare)));path=Path('/private/tmp/distance9_review_prepared.json')if args.prepare else ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'1262ports/20718overlaps/949corners/756profiles',flush=True)
