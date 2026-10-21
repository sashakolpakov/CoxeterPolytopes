#!/usr/bin/env python3
"""Finite A16/D16/D17 long-pair transfers with sparse whole-factor coupons.

No ambient bound: the new coupons require their own complete finite budget.
The mathematical proof supplies arbitrary-root and universal-interface scope.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d16_d17_two_unit_exchange_20261018 as E
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=E.ROOT; T=E.T; A=T.A; P=T.P; Q=T.Q; D=T.D
TAG='fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
CURRENT='fresh_order11_order17_dimension63_20261018'
CHANGES={('A',16):{(0,12),(3,15)},('D',16):{(0,10),(0,13)},
         ('D',17):{(1,12),(1,13)}}
COUPONS=[
    dict(name='A4end_A12end',small='A4',small_nodes=[0,3],large='A12',large_nodes=[0,11],weight='1',pair_cost='4'),
    dict(name='A4end_A13second',small='A4',small_nodes=[0,3],large='A13',large_nodes=[1,11],weight='1',pair_cost='4'),
    dict(name='D5long_A12second',small='D5',small_nodes=[0],large='A12',large_nodes=[1,10],weight='1',pair_cost='2'),
    dict(name='A2_A13end',small='A2',small_nodes=[0,1],large='A13',large_nodes=[0,12],weight='1/2',pair_cost='2'),
    dict(name='A2_A14second',small='A2',small_nodes=[0,1],large='A14',large_nodes=[1,12],weight='1/2',pair_cost='2'),
    dict(name='E6long_A11second',small='E6',small_role='length-two-arm endpoint',small_role_count=2,
         large='A11',large_nodes=[1,9],weight='1',pair_cost='4')]
STEMS=tuple(dict.fromkeys((*E.STEMS,E.TAG,CURRENT,'fresh_e6_b10_half_coupon_forest_budget_20261018')))


def e6_roles(g,ids):
    adj={i:[j for j in ids if j!=i and g[i][j]] for i in ids}
    centers=[i for i in ids if len(adj[i])==3]
    assert len(ids)==6 and len(centers)==1
    center=centers[0]; ds=Q.graph_dist(A.sub(g,ids),[ids.index(center)])
    leaves=[i for i in ids if len(adj[i])==1]
    assert sorted(ds[ids.index(i)] for i in leaves)==[1,2,2]
    return sorted(i for i in leaves if ds[ids.index(i)]==2)


def marked_role(g,block,mark):
    kind=block['type'][0]; rank=len(block['nodes'])
    if block['type']=='E6':
        return dict(type='E6',nodes=block['nodes'],mark=mark,
                    long_endpoint=mark in e6_roles(g,block['nodes']),canonical_images=[])
    images=T.mark_pairs(g,block['nodes'],[mark],rank,kind) if kind in 'ABD' else []
    return dict(type=block['type'],nodes=block['nodes'],mark=mark,
                long_endpoint=False,canonical_images=sorted(p[0] for p in images))


def coupon_matches(g,ids,marks):
    blocks=D.finite_blocks(g,ids)
    a=next(b for b in blocks if marks[0] in b['nodes'])
    b=next(b for b in blocks if marks[1] in b['nodes'])
    if a is b:return []
    pa,pb=marked_role(g,a,marks[0]),marked_role(g,b,marks[1]); found=[]
    for x,y in ((pa,pb),(pb,pa)):
        for spec in COUPONS:
            if x['type']!=spec['small'] or y['type']!=spec['large']:continue
            good=x['long_endpoint'] if x['type']=='E6' else x['canonical_images']==spec['small_nodes']
            if good and y['canonical_images']==spec['large_nodes']:
                found.append(dict(name=spec['name'],weight=spec['weight'],small=x,large=y))
    assert len(found)<=1
    return found


def radical_transfer(c,core):
    g=A.matrix(core['gram']); n=len(g)
    assert c['method']=='radical' and c['core_key']=='(25, 25, ((9, 9),))'
    assert n==19 and core['gauge']==[1]*n and Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1]; v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(n))-black)
    assert len(white)==15 and Q.inertia(A.sub(g,white))==[15,0,0]
    for b in black:assert Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0]
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);shortcuts=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8;shortcuts.append(dict(nodes=[i,j],length=length))
    assert shortcuts
    corners=[];total=F()
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[17,0,0]
        coupons=coupon_matches(g,ids,marks);value=sum((F(t['weight']) for t in coupons),F())
        total+=value
        corners.append(dict(marks=marks,components=D.finite_blocks(g,ids),new_coupons=coupons,
                            used_new_payment=str(value)))
    assert total>=1,(c['index'],total)
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=list(map(str,z)),
                white=white,allowed_contacts=allowed,mixed_sign_shortcuts=shortcuts,
                whole_core_sealed_unless_short8=True,corners=corners,coupon_total=str(total))


def proper_transfer(c,core,maps,target):
    g,inv,du,dv,endproof=X.constraints(c,core);base=c['marked_component'];rows=[]
    for saved in maps:
        if tuple(saved['marked_images'])!=(0,10):continue
        mp=saved['images'];image={mp[t]:i for t,i in enumerate(base)}
        missing=sorted(set(range(16))-set(image));entries=[]
        for new in missing:
            contacts=[(old,image[old]) for old in image if target[new][old]]
            for old,node in contacts:
                assert target[new][old]==3 and node in c['possible_contacts']
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,
                                    inverse_diagonal=str(inv[node][node])))
        assert entries
        record=dict(map=mp,missing_target_nodes=missing,forced_entries=entries)
        if any(F(e['inverse_diagonal'])>2 for e in entries):
            assert len(entries)==1 and F(entries[0]['inverse_diagonal'])==6
            record.update(outcome='projection_exceeds_two',lower_bound='6')
        elif len(entries)==2:
            assert {e['new_target_node'] for e in entries}=={14,15}
            assert len({e['core_node'] for e in entries})==1
            k=entries[0]['core_node'];assert inv[k][k]==F(3,2) and target[14][15]==0
            # Both residual squares <=1/2 but mutual product has magnitude >=3/2.
            record.update(outcome='two_orthogonal_fork_leaves_cauchy',
                          residual_square_upper='1/2',cross_magnitude_lower='3/2',
                          residual_determinant_upper='-2')
        else:
            assert len(entries)==1 and entries[0]['new_target_node']==13
            k=entries[0]['core_node'];assert inv[k][k]==F(4,3)
            assert missing==[13,14,15] and target[13][14]==target[13][15]==3 and target[14][15]==0
            same=next(end for end in c['end_sets'] if k in end)
            assert all(inv[k][k]+inv[j][j]+2*inv[k][j]>2 for j in same if j!=k)
            assert 2*inv[k][k]>2
            outside=sorted(set(c['possible_contacts'])-set(base));tests=[]
            for j in outside:
                if inv[j][j]>2:reason='projection_exceeds_two'
                elif inv[j][j]==2 and inv[k][j]>0:reason='zero_residual_positive_actual_product'
                else:
                    length=min(du[k]+dv[j],dv[k]+du[j])+3
                    assert length<=8;reason='required_bond_shortcut'
                tests.append(dict(node=j,inverse_diagonal=str(inv[j][j]),
                                  center_projection_product=str(inv[k][j]),reason=reason,
                                  required_bond_shortcut=min(du[k]+dv[j],dv[k]+du[j])+3))
            assert all(inv[j][j]>=F(4,3) for j in outside)
            for end in c['end_sets']:
                assert all(inv[a][a]+inv[b][b]+2*inv[a][b]>2
                           for a,b in combinations(set(end)&set(outside),2))
            # The center's sole core contact is k. Old corner roots outside the
            # marked block are therefore orthogonal to it and cannot be leaves.
            black=set(c['left'][1]+[c['right_copy_map'][i] for i in c['right'][1]])
            old_corner=(set(range(len(g)))-black)|set(c['preferred_corner'])
            old_other=sorted(old_corner-set(base));assert k in base
            residual=[[F(2,3),F(-1),F(-1)],[F(-1),F(2),F(0)],[F(-1),F(0),F(2)]]
            assert Q.inertia(residual)==[2,1,0]
            record.update(outcome='core_orthogonal_two_leaf_fork',center_contact_node=k,
                          center_residual_square='2/3',outside_base_contact_checks=tests,
                          old_other_corner_nodes=old_other,old_omitted_nodes=c['omitted_black'],
                          residual_gram=[[str(v) for v in row] for row in residual],
                          residual_signature=[2,1,0])
        rows.append(record)
    assert rows
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],inverse=core['inverse'],
                base_nodes=base,preferred_corner=c['preferred_corner'],end_projection_proof=endproof,maps=rows)


def literal_audit(data):
    old={(r['source'],E.normalize_case(r['case'])):r for r in data[E.TAG]['literal_interfaces']['rows']}
    assert len(old)==216;rows=[]
    for stem,index,g,white,cs in T.literal_inputs(data):
        saved=old[(stem,E.normalize_case(index))];corners=[]
        for j,(marks,unit,_) in enumerate(cs):
            assert marks==saved['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));finite=Q.inertia(A.sub(g,ids))==[len(ids),0,0]
            blocks=D.finite_blocks(g,ids) if finite else []
            shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
            debit=F();images=[]
            if shared:
                key=(shared['type'][0],len(shared['nodes']))
                if key in CHANGES:
                    images=T.mark_pairs(g,shared['nodes'],marks,key[1],key[0])
                    if set(images)&CHANGES[key]:
                        assert set(images)<=CHANGES[key];debit=F(1)
            coupons=coupon_matches(g,ids,marks) if finite else []
            previous=F(saved['corners'][j]['retained_weight'])
            value=previous-debit+sum((F(c['weight']) for c in coupons),F())
            assert value>=0
            corners.append(dict(marks=marks,old_weight=str(previous),debit=str(debit),
                                marked_images=images,new_coupons=coupons,retained_weight=str(value)))
        total=sum(F(c['retained_weight']) for c in corners);assert total>=1,(stem,index,total)
        rows.append(dict(source=stem,case=index,corners=corners,total=str(total)))
    affected=[r for r in rows if any(F(c['debit']) or c['new_coupons'] for c in r['corners'])]
    assert len(rows)==216 and len(affected)==4
    assert sorted(r['case'] for r in affected if r['source']=='fresh_order16_finite_cross_payments_20261017')==[394,648]
    return dict(complete_configurations=216,affected_configurations=4,rows=rows,
                minimum_total=str(min(F(r['total']) for r in rows)))


def run():
    data={stem:json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem,doc in data.items():
        assert doc['status'].startswith('PASS'),stem
        for path,digest in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,(stem,path)
    assert data[CURRENT]['dimension_bound']==63
    over=data[A.OVER];rules=[];supports=[];proper=[];radical_ids=set(E.EXPECTED)
    for (kind,n),drop in CHANGES.items():
        support=T.selected(over,n,kind);target=T.target(n,kind);adj,dist,_=P.distance_data(target)
        affected=[r for r in support['rows'] if set(map(tuple,r['marked_images']))&drop]
        expected={('A',16):[1651,2026],('D',16):[501,504,1031,1282,1534,1537,1614,1621,1688,2547],('D',17):[]}[kind,n]
        assert [r['index'] for r in affected]==expected
        for r in affected:
            c=over['isolation']['cases'][r['index']];core=over['isolation']['cores'][c['core_key']]
            if c['method']=='radical':radical_ids.add(c['index'])
            else:
                assert (kind,n)==('D',16)
                proper.append(proper_transfer(c,core,r['maps'],target))
        if kind=='A':
            old=next(r for r in data[T.TAG]['A16_B16_records'] if r['kind']=='A')
            weights={p:F(int(list(p) in old['retained_pairs'])) for p in combinations(range(n),2)}
        else:
            old=next(r for r in data[E.TAG]['rules'] if r['rank']==n)
            weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
        before=sum(weights.values());assert before==(127 if n==17 else 116)
        assert all(weights[p]==1 for p in drop)
        for p in drop:weights[p]=F()
        autos=P.embeddings(target,target)
        assert all({tuple(sorted(m[i] for i in p)) for p in drop}==drop for m in autos)
        basic={p for p in weights if dist[p[0]][p[1]]<=8 or
               (dist[p[0]][p[1]]==9 and any(len(adj[i])==1 for i in p))}
        assert not basic&drop and all(weights[p]==1 for p in basic)
        rules.append(dict(kind=kind,rank=n,target=target,automorphisms=autos,
                          deleted_pairs=sorted(drop),old_cost=str(before),new_cost=str(sum(weights.values())),
                          weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
                          retained_basic_pairs=sorted(basic)))
        supports.append(dict(kind=kind,rank=n,eligible_patterns=support['complete_eligible_patterns'],
                             complete_maps=support['complete_maps'],affected=affected))
    assert sorted(radical_ids)==[1282,1621,1651,1681,1685,1688,2026,2500,2532,2547]
    radicals=[radical_transfer(over['isolation']['cases'][i],
              over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in sorted(radical_ids)]
    literals=literal_audit(data)
    counts=Counter(m['outcome'] for c in proper for m in c['maps'])
    paths=[Path(__file__),Path(E.__file__),Path(T.__file__),Path(X.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_A16_D16_D17_SPARSE_COUPON_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                rules=rules,coupons=COUPONS,selected_supports=supports,
                proper_completion_cases=proper,proper_outcomes=dict(counts),
                radical_cases=radicals,literal_interfaces=literals,
                baseline_ideal_rule=CURRENT,ideal_weights_unchanged=True,
                new_finite_component_costs=dict(A16=114,B16=104,D16=114,D17=125),
                current_coupons_retained=True,new_dimension_bound=None,
                scope='Complete induced selected maps, six proper-completion exclusions, ten radical transfers and all216 literal interfaces. Written proof handles arbitrary contacts and baseline63 universal short8/deep11 interfaces. Sparse new whole-component coupons require a separate finite budget; no ambient theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_finite_component_costs'],out['proper_outcomes'],flush=True)
