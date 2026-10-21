#!/usr/bin/env python3
"""Saved support diagnostic and three redundant order15 affine B payments.

The local proof uses the complete frozen all-label supporting-white domains.
It does not rerun the overlap/map census or assert a new dimension bound.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal5_small_factor_payments_20261017 as S
import fresh_selective_distance10_dimension116_20261017_review as Q

P=S.P;ROOT=P.ROOT
TAG='fresh_order15_affine_b_pair_reduction_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
OVERLAP='fresh_distance8_high_overlap_20261017'
DROP={(0,11),(11,13),(11,14)}

def extract(g,ids):return [[g[i][j]for j in ids]for i in ids]
def encode(g):return [[str(x)for x in r]for r in g]
def component(g,ids,start):
    out={start};todo=[start]
    while todo:
        i=todo.pop()
        for j in ids:
            if j not in out and g[i][j]:out.add(j);todo.append(j)
    return sorted(out)

def finite_type(g,ids):
    assert Q.inertia(extract(g,ids))==[len(ids),0,0]
    adj={i:{j for j in ids if j!=i and g[i][j]}for i in ids}
    assert all(g[i][j]==-1 for i in ids for j in adj[i])
    assert sum(map(len,adj.values()))==2*(len(ids)-1)
    forks=[i for i in ids if len(adj[i])==3]
    if not forks:
        assert max(map(len,adj.values()))<=2
        return 'A'+str(len(ids))
    assert len(forks)==1 and all(len(a)<=3 for a in adj.values())
    assert len([i for i in adj[forks[0]]if len(adj[i])==1])==2
    return 'D'+str(len(ids))

def support_diagnostic(new,old,over):
    records=[]
    for m in(14,15,16):
        for kind in 'BCD':
            rec=next(r for r in new['records']if r['order']==m and r['kind']==kind)
            prev=next(r for r in old['records']if r['order']==m and r['kind']==kind)
            target=rec['target'];autos=P.embeddings(target,target)
            adj,ds,features=P.distance_data(target)
            paid=set(map(tuple,rec['retained_paid_pairs']));basic=set(map(tuple,prev['basic_pairs']))
            todo=paid-basic;orbits=[]
            while todo:
                p=min(todo);orbit={tuple(sorted(a[i]for i in p))for a in autos}
                assert orbit<=todo;todo-=orbit;cases=[]
                for r in rec['case_embeddings']:
                    if set(map(tuple,r['marked_images']))&orbit:
                        c=over['isolation']['cases'][r['case']]
                        cases.append({k:c[k]for k in('index','left','right','match','preferred_corner','component_order_bound','method','core_key')})
                orbits.append(dict(representative=list(p),distance=ds[p[0]][p[1]],
                    pairs=[list(p)for p in sorted(orbit)],contexts=cases))
            records.append(dict(order=m,kind=kind,automorphism_count=len(autos),
                paid_cost=len(paid),basic_pair_count=len(basic),extra_pair_orbits=orbits))
    b=next(r for r in records if r['order']==15 and r['kind']=='B')
    assert b['paid_cost']==100 and b['basic_pair_count']==91 and b['automorphism_count']==2
    assert [len(o['pairs'])for o in b['extra_pair_orbits']]==[1,1,1,2,2,2]
    assert [len(o['contexts'])for o in b['extra_pair_orbits']]==[13,3,11,53,17,4]
    return records

def run():
    stems=[SMALL,OVERLAP,S.TAG,'fresh_finite8_ideal5_combined_payments_20261017',
        'fresh_ideal5_piecewise_dimension75_20261017']
    inputs={s:S.read(s)for s in stems}
    for stem,data in inputs.items():
        assert data['status'].startswith('PASS'),stem
        for name,h in data['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    new=inputs[SMALL];old=inputs[S.TAG];over=inputs[OVERLAP]
    diagnostics=support_diagnostic(new,old,over)
    rec=next(r for r in new['records']if r['order']==15 and r['kind']=='B')
    target=rec['target'];autos=P.embeddings(target,target);adj,ds,features=P.distance_data(target)
    oldset=set(map(tuple,rec['retained_paid_pairs']));reduced=oldset-DROP
    assert len(reduced)==97 and all(ds[a][b]>9 for a,b in DROP)
    assert all({tuple(sorted(a[i]for i in p))for p in DROP}==DROP for a in autos)
    basic=S.basic_pairs(target);assert basic<=reduced
    affected=[r['case']for r in rec['case_embeddings']if set(map(tuple,r['marked_images']))&DROP]
    assert affected==[1284,1623,1641,1646,1666,1933,2326]
    rows=[]
    for index in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        assert c['method']=='lorentz_residual' and c['component_order_bound']==15
        g=[[F(v)for v in row]for row in core['gram']];n=len(g)
        assert all(F(x)==1 for x in core['gauge'])
        assert Q.inertia(g)==[n-1,1,0]
        inv=Q.inverse(g);assert encode(inv)==core['inverse']
        nodes=c['marked_component'];marks=c['preferred_corner']
        assert all(m in nodes for m in marks)
        black=set(c['left'][1])|{c['right_copy_map'][i]for i in c['right'][1]}
        assert len(black)==4
        white=sorted(set(range(n))-black);corner=sorted(set(white)|set(marks))
        assert Q.inertia(extract(g,white))==[n-4,0,0]
        assert Q.inertia(extract(g,corner))==[n-2,0,0]
        assert component(g,corner,marks[0])==nodes
        assert set(corner)-set(nodes)<=set(white)
        candidates={r['contact_node']:r for r in c['single_contact_candidates']}
        for node,r in candidates.items():
            assert inv[node][node]==F(r['projection_square'])
            assert 2-inv[node][node]==F(r['residual_square'])
            # All nonzero norm2 Coxeter contacts are at least1; the next
            # possibility is sqrt2, already above each norm bound here.
            assert 2*inv[node][node]>2
        if index!=1641:
            assert n==17 and set(candidates)=={3,13}
            assert not(set(nodes)&set(candidates))
            assert finite_type(g,nodes)=='D14'
            assert all(F(r['projection_square'])==2 and F(r['residual_square'])==0 for r in candidates.values())
            outcome=dict(method='no_direct_contact_to_marked_component',whole_finite_type='D14',maximum_order=14)
        else:
            assert n==16 and finite_type(g,nodes)=='A12'
            assert set(candidates)=={1,3,11,13} and set(nodes)&set(candidates)=={1,11}
            assert all(sum(bool(g[i][j])for j in nodes if i!=j)==1 for i in(1,11))
            pair_rows=[];extensions=[]
            for i in(1,11):
                assert inv[i][i]==2
                for j in sorted(candidates):
                    value=inv[i][j]
                    assert value>0 or value<=-2
                    pair_rows.append(dict(first=i,second=j,projection_product=str(value),
                        exclusion='acuteness'if value>0 else'strict_positive_definiteness_of_two_supporting_whites'))
                extended=[r+[F(-int(k==i))]for k,r in enumerate(g)]
                extended.append([F(-int(k==i))for k in range(n)]+[F(2)])
                ids=nodes+[n];assert finite_type(extended,ids)=='A13'
                extensions.append(dict(contact_node=i,whole_finite_type='A13',nodes=ids))
            outcome=dict(method='zero_residual_endpoint_extension',base_finite_type='A12',maximum_order=13,
                exclusions=pair_rows,extensions=extensions,
                orthogonal_whites='The residual-zero added root lies in the core span, hence is orthogonal to every core-orthogonal white.')
        rows.append(dict(index=index,left=c['left'],right=c['right'],match=c['match'],
            core_key=c['core_key'],gram=encode(g),signature=Q.inertia(g),marks=marks,
            actual_white=white,actual_corner=corner,marked_component=nodes,
            whole_component_complement_is_white=True,allowed_single_contacts=c['single_contact_candidates'],
            inherited_end_sets=c['end_sets'],outcome=outcome,finite_payment=1))
    assert len(rows)==7
    boundary=F(-16873,171366)
    gain_two=boundary+F(2,15);gain_three=boundary+F(3,15)
    assert boundary+F(1,15)<0<gain_two<gain_three
    paths=[Path(__file__),Path(Q.__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER15_B_THREE_PAIR_LOCAL_REDUCTION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        support_diagnostic=diagnostics,affine_type='B14',diagram_order=15,
        old_cost=100,new_cost=97,removed_pairs=[list(p)for p in sorted(DROP)],
        new_paid_pairs=[list(p)for p in sorted(reduced)],affected_cases=affected,cases=rows,
        triangle_and_short8_pairs_retained=True,new_finite_coupons=[],
        quantitative_boundary=dict(ambient_dimension=75,ranks=[12,12,12,12,12,14],old_margin=str(boundary),
            one_pair_margin=str(boundary+F(1,15)),two_pair_margin=str(gain_two),three_pair_margin=str(gain_three),
            warning='Only this one formal profile is compared. All other profiles and concentration moves require a separate new comparison.'),
        scope='Complete saved support diagnostic and seven all-support local finite alternatives, using the frozen complete all-label/no-extra-cross domain. No new finite budget or ambient dimension bound is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'7 finite alternatives, B14 order15 cost100->97',flush=True)
