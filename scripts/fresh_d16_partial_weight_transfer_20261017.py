#!/usr/bin/env python3
"""Ten affine-D15 pair weights fall to one half after actual-factor completion.

The frozen complete overlap and induced-map censuses are inputs. This checker
reconstructs every affected completion, its residual obstruction or sealing
argument, and its whole finite-component payments with rational arithmetic.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_order13_affine_b_singleton_transfer_20261017 as L

ROOT=L.ROOT; P=L.P; Q=L.Q
TAG='fresh_d16_partial_weight_transfer_20261017'
OVER='fresh_distance8_high_overlap_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
HALF_PAIRS={(0,10),(1,11),(1,14),(1,15),(10,12),(10,13),
            (2,14),(2,15),(9,12),(9,13)}
EXPECTED=[369,373,504,509,515,662,766,989,997,1019,1031,1036,1037,1050,1051,
          1164,1177,1266,1285,1312,1537,1614,1624,1647,1648,1667,1668,1796,
          1830,1906,1907,1912,1934,1935,1947,1975,1979,2035,2044,2189,2223,
          2299,2300,2305,2327,2328,2340,2368,2372,2495,2505,2512,2536,2540]


def components(g,ids):
    todo=set(ids); out=[]
    while todo:
        ns={min(todo)}; stack=list(ns)
        while stack:
            i=stack.pop()
            for j in sorted(todo-ns):
                if g[i][j]:
                    ns.add(j); stack.append(j)
        todo-=ns; ns=sorted(ns)
        out.append(ns)
    return out


def finite_blocks(g,ids):
    out=[]
    for ns in components(g,ids):
        assert Q.inertia(L.sub(g,ns))==[len(ns),0,0]
        adj={i:[j for j in ns if i!=j and g[i][j]] for i in ns}
        assert sum(map(len,adj.values()))==2*(len(ns)-1)
        assert all(g[i][j]==-1 for i in ns for j in adj[i])
        forks=[i for i in ns if len(adj[i])==3]
        if forks:
            assert len(forks)==1 and max(map(len,adj.values()))==3
            center=forks[0]; leaves=[i for i in adj[center] if len(adj[i])==1]
            if len(leaves)>=2:
                kind='D'+str(len(ns))
                roles={str(center):'fork_center',**{str(i):'short_fork_leaf' for i in leaves}}
                assert (len(ns)==4 and len(leaves)==3) or len(leaves)==2
            else:
                assert len(ns) in (6,7,8)
                kind='E'+str(len(ns)); roles={}
        else:
            assert max(map(len,adj.values()))<=2
            kind='A'+str(len(ns)); roles={str(i):'endpoint' for i in ns if len(adj[i])<=1}
        out.append(dict(nodes=ns,type=kind,roles=roles))
    return out


def core_constraints(c,core):
    """Recheck all numeric hypotheses of the inherited residual lemmas."""
    g=[[F(x) for x in row] for row in core['gram']]; n=len(g)
    assert core['gauge']==[1]*n and set(x for row in g for x in row)<={F(-1),F(0),F(2)}
    assert Q.inertia(g)==[n-1,1,0]
    J=Q.inverse(g); assert L.encode(J)==core['inverse']
    u=c['left'][1]; v=[c['right_copy_map'][i] for i in c['right'][1]]
    corner=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert corner[key]==c[key],(c['index'],key)
    du=Q.graph_dist(g,u); dv=Q.graph_dist(g,v)
    allowed=c['possible_contacts']; ends=c['end_sets']
    assert set(ends[0]).isdisjoint(ends[1]) and set(ends[0]+ends[1])==set(allowed)
    assert max((du[i] for i in ends[0]),default=0)+max((dv[i] for i in ends[1]),default=0)+2<=8
    assert all(J[i][j]>0 for end in ends for i in end for j in end)
    assert all(J[i][i]+J[j][j]+2*J[i][j]>2 for end in ends for i,j in combinations(end,2))
    candidates=[i for i in allowed if J[i][i]<=2]
    assert all(J[i][i]>=1 for i in candidates)
    result=dict(allowed_contacts=allowed,end_sets=ends,
                candidate_nodes=candidates,projection_diagonals={str(i):str(J[i][i]) for i in allowed})
    if c['method']=='lorentz_residual':
        assert all(J[i][i]>1 for i in candidates)
        assert candidates==[r['contact_node'] for r in c['single_contact_candidates']]
        checks=[]
        for i,j in product(candidates,repeat=2):
            a=J[i][j]; residual=(2-J[i][i])*(2-J[j][j])
            shortcut=min(du[i]+dv[j],du[j]+dv[i])+3
            assert a*a>residual and (a>0 or shortcut<=8)
            checks.append(dict(nodes=[i,j],projection_square=str(a*a),residual_product=str(residual),
                               sign=1 if a>0 else -1,shortcut=shortcut if a<0 else None))
        assert checks==c['pair_tests']
        for r in c['single_contact_candidates']:
            p=J[r['contact_node']][r['contact_node']]; rho=2-p
            assert r['projection_square']==str(p) and r['residual_square']==str(rho)
            assert r['maximum_attached_orthogonal_rank']==rho//(1-rho)
        result.update(kind='at_most_one_core_contact_root',pair_tests=checks)
    else:
        assert c['method']=='relevant_component'
        assert candidates==c['contact_domains'] and len(c['relevant_nodes'])==1
        good=c['relevant_nodes'][0]
        assert [i for i in candidates if i in c['marked_component']]==[good]
        assert 1<J[good][good]<2 and len(c['extension_cases'])==1
        r=c['extension_cases'][0]; assert r['node']==good
        rho=2-J[good][good]
        assert r['projection']==str(J[good][good]) and r['rho']==str(rho)
        assert r['alone_tail_bound']==rho//(1-rho)
        checks=[]; outside=[]
        for j in candidates:
            a=J[good][j]; residual=rho*(2-J[j][j]); shortcut=min(du[good]+dv[j],du[j]+dv[good])+3
            if a*a>residual and a>0: reason='positive'
            elif a*a>residual and a<0 and shortcut<=8: reason='short_path'
            else:
                assert a>0 and (a+1)**2>residual and j not in c['marked_component']
                reason='forced_orthogonal'; outside.append(j)
            checks.append(dict(node=j,reason=reason,projection=str(a),residual_product=str(residual),shortcut=shortcut))
        assert checks==r['pair_checks'] and len(outside)==1
        doubles=[]
        for j in candidates:
            if J[j][j]!=1: continue
            a=J[good][j]; shortcut=min(du[good]+dv[j],du[j]+dv[good])+3
            assert a>0 or (a<0 and shortcut<=8)
            doubles.append(dict(node=j,sign=1 if a>0 else -1,shortcut=shortcut))
        assert doubles==r['double_checks']
        j=outside[0]; S=[[rho,-J[good][j]],[-J[good][j],2-J[j][j]]]
        det=S[0][0]*S[1][1]-S[0][1]**2
        assert r['joint']['node']==j and r['joint']['determinant']==str(det)
        if det==0:
            assert S[0][0]>0 and S[1][1]>0 and S[0][1]<0
            assert r['joint']['kind']=='positive_residual_relation' and r['joint']['tail']==0
        else:
            assert S[1][1]/det==2 and r['joint']['inverseGG']=='2'
            assert r['joint']['kind']=='saturated_single_tail' and r['joint']['tail']==1
        result.update(kind='one_relevant_root_with_optional_exterior_root',relevant_node=good,
                      exterior_node=j,residual_gram=L.encode(S),joint_tail_bound=r['joint']['tail'],
                      joint_component_order_bound=len(c['marked_component'])+1+r['joint']['tail'],
                      pair_checks=checks,double_checks=doubles)
        assert result['joint_component_order_bound']<16
    return g,result


def complete(c,core,g0,constraints,target,mp):
    n=len(g0); base=c['marked_component']; image={mp[t]:i for t,i in enumerate(base)}
    missing=sorted(set(range(16))-set(image)); image.update({i:n+t for t,i in enumerate(missing)})
    nn=n+len(missing); g=[row+[F(0)]*len(missing) for row in g0]+[[F(0)]*nn for _ in missing]
    for i in missing:
        g[image[i]][image[i]]=2
        for j in range(16):
            if target[i][j]: g[image[i]][image[j]]=g[image[j]][image[i]]=F(-1)
    required=[dict(new_node=i,core_nodes=[j for j in range(n) if g[i][j]])
              for i in range(n,nn) if any(g[i][j] for j in range(n))]
    row=dict(target_to_literal=[image[i] for i in range(16)],new_roots=list(range(n,nn)),required_core_contacts=required)
    assert all(len(r['core_nodes'])==1 for r in required)
    if len(required)>1:
        assert c['method']=='lorentz_residual'
        row.update(outcome='excluded_multiple_core_contact_roots')
        return row,g
    assert len(required)==int(bool(missing))
    if required:
        node=required[0]['core_nodes'][0]
        permitted=(constraints['candidate_nodes'] if c['method']=='lorentz_residual'
                   else [constraints['relevant_node']])
        if node not in permitted:
            assert c['method']=='lorentz_residual'
            assert node in constraints['allowed_contacts']
            projection=F(constraints['projection_diagonals'][str(node)])
            assert projection>2
            row.update(outcome='excluded_projection_above_norm',projection_square=str(projection))
            return row,g
        if c['method']=='relevant_component':
            assert len(missing)>1+constraints['joint_tail_bound']
            row['exterior_root_exclusion']=dict(required_factor_order=16,
                maximum_order_if_exterior_root_present=constraints['joint_component_order_bound'])
        sealing='Forced attaching root exhausts the permitted core contacts; the positive affine kernel seals its entire completion.'
    else:
        assert c['method']=='lorentz_residual' and constraints['candidate_nodes']==[]
        sealing='No core-contact root is possible; the displayed affine factor is already whole.'
    affine=[image[i] for i in range(16)]
    assert Q.inertia(L.sub(g,affine))==[15,0,1]
    kernel=Q.nullspace(L.sub(g,affine)); assert len(kernel)==1
    z=kernel[0]; z=[x if z[0]>0 else -x for x in z]
    assert all(x>0 for x in z)
    signature=Q.inertia(g)
    if signature[1]>1:
        row.update(outcome='excluded_lorentz_index',gram=L.encode(g),signature=signature,
                   forced_completion_reason=sealing)
        return row,g
    assert signature[1]==1
    black=set(c['left'][1])|{c['right_copy_map'][i] for i in c['right'][1]}
    white=sorted(set(range(nn))-black)
    assert len(black)==4 and Q.inertia(L.sub(g,white))==[len(white),0,0]
    for b in black:
        ids=sorted(set(white)|{b}); assert Q.inertia(L.sub(g,ids))==[len(ids),0,0]
    corners=[]
    for corner in c['all_corners']:
        marks=corner['preferred_corner']; ids=sorted(set(white)|set(marks)); sig=Q.inertia(L.sub(g,ids))
        if marks==c['preferred_corner']:
            assert sig==[len(ids)-1,0,1] and sorted(affine) in components(g,ids)
            corners.append(dict(marks=marks,signature=sig,kind='affineD15_order16',ideal_weight='1/2',finite_unit=0,finite_half='0'))
            continue
        assert sig==[len(ids),0,0]
        cc=finite_blocks(g,ids); selected=[next(r for r in cc if m in r['nodes']) for m in marks]
        unit=int(selected[0] is selected[1] and len(selected[0]['nodes'])<=17)
        coupon=F(0); roles=[]
        for i,j in ((0,1),(1,0)):
            a,b=selected[i],selected[j]; role=b['roles'].get(str(marks[j]))
            if len(a['nodes'])<=5 and b['type'] in {'D'+str(t) for t in range(11,17)} and role in ('fork_center','short_fork_leaf'):
                coupon=F(1,2)
                roles.append(dict(small_type=a['type'],small_mark=marks[i],large_type=b['type'],large_mark=marks[j],large_role=role))
        corners.append(dict(marks=marks,signature=sig,kind='finite',components=cc,
                            marked_component_indices=[cc.index(r) for r in selected],ideal_weight='0',
                            finite_unit=unit,finite_half=str(coupon),coupon_roles=roles))
    finite=sum(r['finite_unit']+F(r['finite_half']) for r in corners)
    assert finite>=F(1,2)
    row.update(outcome='sealed_completion',sealing=sealing,gram=L.encode(g),signature=signature,white=white,
               affine_nodes=affine,positive_affine_kernel=list(map(str,z)),corners=corners,
               ideal_total='1/2',finite_payment_total=str(finite),total_payment=str(finite+F(1,2)))
    return row,g


def run():
    stems=[OVER,SMALL,'fresh_finite8_ideal5_combined_payments_20261017',
           'fresh_finite8_a17_d11_coupon_budget_20261017',
           'fresh_distance9_supporting_extensions_review_20261017',
           'fresh_d17_partial_weight_transfer_20261017']
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,j in inputs.items():
        assert j['status'].startswith('PASS'),stem
        for name,h in j['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    over=inputs[OVER]
    rec=next(r for r in inputs[SMALL]['records'] if r['order']==16 and r['kind']=='D')
    target=rec['target']; autos=P.embeddings(target,target); _,ds,features=P.distance_data(target)
    paid=set(map(tuple,rec['retained_paid_pairs'])); assert len(paid)==111
    assert HALF_PAIRS=={p for p in paid if ds[p[0]][p[1]]>9}
    assert all({tuple(sorted(a[i] for i in p)) for p in HALF_PAIRS}==HALF_PAIRS for a in autos)
    affected=[r for r in rec['case_embeddings'] if set(map(tuple,r['marked_images']))&HALF_PAIRS]
    assert [r['case'] for r in affected]==EXPECTED
    rows=[]; outcomes=Counter(); map_count=0
    for r in affected:
        c=over['isolation']['cases'][r['case']]; core=over['isolation']['cores'][c['core_key']]
        g0,constraints=core_constraints(c,core)
        configurations=[]; keys={}
        for mp in r['embeddings']:
            marked=tuple(sorted(mp[c['marked_component'].index(i)] for i in c['preferred_corner']))
            if marked not in HALF_PAIRS: continue
            assert len(set(mp))==len(mp)
            assert all((int(2-g0[a][b]) if a!=b and g0[a][b] else 0)==target[mp[i]][mp[j]]
                       for i,a in enumerate(c['marked_component']) for j,b in enumerate(c['marked_component']))
            row,g=complete(c,core,g0,constraints,target,mp)
            key=str(g)
            if key not in keys:
                keys[key]=len(configurations); row['embeddings']=[]; configurations.append(row)
                outcomes[row['outcome']]+=1
            configurations[keys[key]]['embeddings'].append(dict(map=mp,marked_image=list(marked)))
            map_count+=1
        assert configurations
        rows.append(dict(index=c['index'],core_key=c['core_key'],method=c['method'],left=c['left'],right=c['right'],
                         match=c['match'],preferred_corner=c['preferred_corner'],base_nodes=c['marked_component'],
                         core_gram=core['gram'],core_inverse=core['inverse'],constraints=constraints,completions=configurations))
    weights={p:(F(1,2) if p in HALF_PAIRS else F(1) if p in paid else F(0)) for p in combinations(range(16),2)}
    assert sum(weights.values())==106
    triangles={p for p in combinations(range(16),2) if ds[p[0]][p[1]]<=5 or
               (ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
    short8={p for p in combinations(range(16),2) if ds[p[0]][p[1]]<=8}
    assert all(weights[p]==1 for p in triangles|short8)
    assert all(weights[tuple(p)]==0 for p in rec['dropped_pairs'])
    orbits=[]; todo=set(HALF_PAIRS)
    while todo:
        p=min(todo); orbit={tuple(sorted(a[i] for i in p)) for a in autos}; assert orbit<=todo; todo-=orbit
        cases=[r['case'] for r in affected if set(map(tuple,r['marked_images']))&orbit]
        orbits.append(dict(representative=list(p),distance=ds[p[0]][p[1]],pairs=[list(q) for q in sorted(orbit)],cases=cases))
    assert [len(r['pairs']) for r in orbits]==[2,4,4]
    assert map_count==328 and outcomes=={
        'excluded_lorentz_index':22,'excluded_multiple_core_contact_roots':20,
        'excluded_projection_above_norm':14,'sealed_completion':10}
    paths=[Path(__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_AFFINE_D15_TEN_HALF_WEIGHT_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='D15',diagram_order=16,target=target,automorphisms=autos,
                inherited_cost=111,new_cost='106',half_weight_pairs=[list(p) for p in sorted(HALF_PAIRS)],
                weighted_pairs=[dict(pair=list(p),weight=str(w)) for p,w in sorted(weights.items())],
                half_weight_orbits=orbits,affected_cases=EXPECTED,affected_map_count=map_count,
                completion_outcomes=dict(outcomes),cases=rows,
                triangle_pairs_retained_at_weight_one=[list(p) for p in sorted(triangles)],
                short8_pairs_retained_at_weight_one=[list(p) for p in sorted(short8)],
                inherited_zero_pairs=rec['dropped_pairs'],added_finite_coupons=[],
                scope='Complete saved-map support, exact actual-factor completions, residual exclusions and sealed whole finite payments. Ten ideal pairs have weight one half; the previous nine zeros are unchanged. No new ambient bound or scalar comparison is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--verify',action='store_true'); a=p.parse_args()
    out=json.loads(json.dumps(run())); path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_map_count'],'maps;',out['completion_outcomes'],'; cost106',flush=True)
