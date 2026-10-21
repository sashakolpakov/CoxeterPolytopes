#!/usr/bin/env python3
"""Complete one-touch irreducible finite-component rank completion.

All-size A/B/D families are symbolic. Rank29 instances are a logically
separate finite export using the published compact dimension upper bound.
No sealing against more external attachments is asserted.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015 as C
import fresh_sqrt5_quad_filters_20261015 as A
import fresh_sqrt5_single_chain_completion_20261015 as S
import fresh_sqrt5_singular_attachment_closure_20261015 as ATT
import fresh_sqrt5_sealed_weights_20261015_geometry as GEO

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rooted_component_completion_20261015';K=A.K


def diagram(kind,n):
    g={p:2 for p in combinations(range(n),2)}
    if kind in ('A','B','F','H','I'):
        for i in range(n-1):g[i,i+1]=3
        if kind=='B':g[n-2,n-1]=4
        if kind=='F':assert n==4;g[1,2]=4
        if kind=='H':assert n in (3,4);g[0,1]=5
        if kind=='I':assert n==2;g[0,1]=5
    elif kind=='D':
        assert n>=4
        for i in range(n-3):g[i,i+1]=3
        g[n-3,n-2]=g[n-3,n-1]=3
    elif kind=='E':
        assert n in (6,7,8);start=1
        for length in (1,2,n-4):
            prev=0
            for i in range(start,start+length):g[tuple(sorted((prev,i)))]=3;prev=i
            start+=length
        assert start==n
    else:raise ValueError(kind)
    return dict(nodes=n,labels=[g[p] for p in combinations(range(n),2)])


def tree_key(row,root):
    g=C.edges(row);adj={i:{} for i in range(row['nodes'])}
    for (i,j),m in g.items():
        if m!=2:adj[i][j]=adj[j][i]=m
    def rec(i,p):return tuple(sorted((m,rec(j,i)) for j,m in adj[i].items() if j!=p))
    return rec(root,None)


FIXED=[('A',5,2),('A',6,2),('A',7,2),('B',3,2),
       ('D',5,4),('D',6,5),('D',7,6),('E',6,3),('E',7,6),('H',3,2),('I',2,0)]


def allowed_templates(n):
    out=[('A',n,0)]
    if n>=3:out.append(('A',n,1))
    if n>=2:out.append(('B',n,0))
    if n>=4:out.append(('D',n,0))
    out.extend(t for t in FIXED if t[1]==n)
    return out


def leaf_test(row,root):
    n=row['nodes'];g=C.edges(row)
    g.update({(i,n):3 if i==root else 2 for i in range(n)})
    return C.E.finite_word(n+1,tuple(g[p] for p in combinations(range(n+1),2)))


def classification_checks():
    # A/B/D proofs below are all-size; this finite replay checks every root
    # through rank35 and every exceptional type against the rooted list.
    counts=Counter()
    for n in range(1,36):
        expect={tree_key(diagram(k,r),i) for k,r,i in allowed_templates(n)}
        types=['A']+(['B'] if n>=2 else [])+(['D'] if n>=4 else [])
        if n in (6,7,8):types+=['E']
        if n==4:types+=['F']
        if n in (3,4):types+=['H']
        if n==2:types+=['I']
        actual=set()
        for kind in types:
            row=diagram(kind,n)
            assert C.E.finite_word(n,tuple(row['labels']))
            for root in range(n):
                ok=leaf_test(row,root);counts['root_tests']+=1
                if ok:actual.add(tree_key(row,root));counts['passing_root_tests']+=1
        assert actual==expect,(n,actual-expect,expect-actual)
    for label in (6,10):
        assert not C.E.finite_word(3,(label,3,2));counts['isolated_high_dihedral_failures']+=1
    return dict(counts)


def graft(core,contact,spec):
    kind,n,root=spec;ext=diagram(kind,n);old=core['nodes'];total=old+n
    g={p:2 for p in combinations(range(total),2)};g.update(C.edges(core))
    for (i,j),m in C.edges(ext).items():g[old+i,old+j]=m
    for i,m in enumerate(contact):g[i,old+root]=m
    return dict(nodes=total,white=list(range(4,total)),labels=[g[p] for p in combinations(range(total),2)])


def inverse_cost(spec):
    kind,n,root=spec;row=diagram(kind,n);co=A.coloring(row);assert co['passes']
    G=A.scaled_matrix(row,co['colors']);Gi=S.inverse(G)
    return G[root][root]*Gi[root][root]


def tail_corner_witness(core,contact):
    g=C.edges(core);out=[]
    for a,b in product((0,2),(1,3)):
        ids=core['white']+[a,b];touched=[i for i in ids if contact[i]!=2]
        if not touched:out.append(dict(marked=[a,b],component_nodes=[],path_length=0));continue
        assert len(touched)==1 and contact[touched[0]]==3
        root=touched[0];adj={i:[j for j in ids if i!=j and g[tuple(sorted((i,j)))]!=2] for i in ids}
        seen={root};todo=[root]
        while todo:
            for j in adj[todo.pop()]:
                if j not in seen:seen.add(j);todo.append(j)
        assert len(adj[root])<=1 and all(len(adj[i])<=2 for i in seen)
        assert all(g[tuple(sorted((i,j)))]==3 for i in seen for j in adj[i])
        out.append(dict(marked=[a,b],attachment=root,component_nodes=sorted(seen),path_length=len(seen)))
    return out


def full_test(row):
    n=row['nodes'];g=C.edges(row)
    if not C.E.corners_finite(row['white'],g):return dict(status='REJECT_NONFINITE_CORNER')
    co=A.coloring(row);assert co['passes'];G=A.scaled_matrix(row,co['colors'])
    pi=A.inertia(G);ci=A.inertia([[x.conj() for x in r] for r in G])
    out=dict(physical_inertia=pi,conjugate_inertia=ci)
    if pi!=[n-2,1,1] or ci!=[n-1,0,1]:out['status']='REJECT_FULL_INERTIA';return out
    if (w:=C.affine_witness(row)):out.update(status='REJECT_AFFINE',affine_witness=w)
    elif not C.E.no_short(row['white'],g):out['status']='PAID_SHORT'
    elif C.h4_payments(row):out['status']='PAID_H4'
    else:out['status']='UNPAID_PRINCIPAL_STAR'
    return out


def run():
    load=lambda tag:json.loads((ROOT/f'results/{tag}.json').read_text())
    single=load(S.TAG);old=load('fresh_sqrt5_quad_contacts_20261015')
    union=load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    lookup={r['class_index']:r for r in old['records']}
    original=load('fresh_sqrt5_forest_sensitive_weights_20261015');sealed=load('fresh_sqrt5_sealed_weights_20261015')
    extended=load('fresh_sqrt5_extended_forest_weights_20261015')
    catalogue=classification_checks();fixed_costs={t:inverse_cost(t) for t in FIXED}
    # Norm-two leaf elimination gives reduced diagonal1 at the distant
    # weighted end/fork, and1 is a fixed point along every simple bond.
    assert K(2)-K(2)/2==1 and K(2)-K(1)/2-K(1)/2==1 and K(2)-1/K(1)==1
    counts=Counter();records=[];families=[];instances=[];seen_instances=set()
    def save_instance(profile_index,kind,rank,root,row,origin):
        # White paths/components are literal; no cross-source isomorphism
        # quotient is claimed for this conservative weight input.
        if row['nodes']>31:return
        key=(profile_index,kind,rank,root)
        if key in seen_instances:return
        seen_instances.add(key)
        assert C.E.corners_finite(row['white'],C.edges(row))
        pay=ATT.payment(row,original,sealed,extended)
        gap=31-row['nodes']
        instances.append(dict(index=len(instances),source_profile=profile_index,
            external_type=kind,external_rank=rank,attachment_root=root,row=row,origin=origin,
            corner_forests_and_orbits=pay['corners'],current_sums=pay['sums_without_with_external_a1'],
            rank29_orthogonal_spectator_rank=gap,
            possible_spectator_A1_counts=[j for j in range(gap+1) if gap-j!=1]))
    for ri,saved in enumerate(single['records']):
        core=union[lookup[saved['source_class_index']]['union_index']];contact=saved['contacts']
        d=K(*saved['first_norm']);c=K(*saved['projection_cost']);q=c/d
        out=dict(source_profile=ri,source_class_index=saved['source_class_index'],trial_index=saved['trial_index'],
                 contacts=contact,cost_ratio=q.enc(),finite_candidates=[])
        specs=[]
        for branch,den,num,minimum,root in (
                ('A_END',2*q-1,K(1),1,0),('A_NEXT',4*q-1,4*q+1,3,1)):
            n=None if not den else num/den
            out[branch]=dict(rank_in_field=None if n is None else n.enc())
            if n is not None and not n.b and n.a.denominator==1 and n.a>=minimum:
                specs.append(('A',int(n.a),root))
        specs.extend(t for t,z in fixed_costs.items() if q*z==1)
        assert len(specs)==len(set(specs))
        for spec in specs:
            assert q*inverse_cost(spec)==1
            row=graft(core,contact,spec);test=full_test(row);counts[test['status']]+=1
            item=dict(external_type=spec[0],rank=spec[1],root=spec[2],**test)
            if test['status']=='UNPAID_PRINCIPAL_STAR':
                item['row']=row;save_instance(ri,*spec,row,origin='FINITE_EXACT_RANK_SOLUTION')
            out['finite_candidates'].append(item)
        if q==F(1,2):
            witnesses=tail_corner_witness(core,contact)
            assert sum(m!=2 for m in contact)==1 and 3 in contact
            # The effective old-core Schur matrix is independent of tail rank
            # and type: root inverse*root norm is exactly2.
            gb=graft(core,contact,('B',2,0));gd=graft(core,contact,('D',4,0))
            tb=full_test(gb);td=full_test(gd)
            assert tb['status']==td['status']=='UNPAID_PRINCIPAL_STAR'
            geometry=GEO.geometry(gb);other=GEO.geometry(gd)
            assert geometry['projected_gram']==other['projected_gram']
            fam=dict(index=len(families),source_profile=ri,source_class_index=saved['source_class_index'],
                     trial_index=saved['trial_index'],contacts=contact,corner_attachment_paths=witnesses,
                     root_inverse_times_norm='2',representative_B2=gb,representative_D4=gd,
                     B2_checks=tb,D4_checks=td,common_projected_gram=geometry['projected_gram'],
                     actual_quadrilateral_geometry=geometry)
            families.append(fam);out['tail_family_index']=fam['index']
            for kind,start in (('B',2),('D',4)):
                for n in range(start,32-core['nodes']):
                    row=graft(core,contact,(kind,n,0))
                    save_instance(ri,kind,n,0,row,origin='SYMBOLIC_ALL_RANK_B_D_TAIL')
        records.append(out)
        if (ri+1)%50==0:print('profiles',ri+1,'finite',dict(counts),'families',len(families),'instances',len(instances),flush=True)
    assert len(records)==458 and len(families)==37
    assert counts==Counter(UNPAID_PRINCIPAL_STAR=209,REJECT_NONFINITE_CORNER=25)
    assert len(instances)==1715
    paths=[f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',f'results/{S.TAG}.json',f'results/{S.TAG}_review.json',
           'scripts/fresh_sqrt5_quad_contacts_20261015.py','scripts/fresh_sqrt5_quad_filters_20261015.py',
           f'scripts/{ATT.TAG}.py',f'results/{ATT.TAG}.json',
           'scripts/fresh_sqrt5_sealed_weights_20261015_geometry.py',
           'results/fresh_sqrt5_sealed_weights_20261015_geometry.json',
           'results/fresh_sqrt5_sealed_weights_20261015_review.json']
    for tag in('fresh_sqrt5_quad_contacts_20261015','fresh_sqrt5_quad_exhaustion_20261015',
               'fresh_sqrt5_forest_sensitive_weights_20261015','fresh_sqrt5_sealed_weights_20261015',
               'fresh_sqrt5_extended_forest_weights_20261015'):
        paths.append(f'results/{tag}.json')
    ATT.Q.manifest_check(paths)
    return dict(status='PASS_ALL_RANK_ONE_TOUCH_FINITE_COMPONENT_NECESSITY',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        rooted_classification_checks=catalogue,fixed_rooted_types=[dict(type=k,rank=n,root=r,
            inverse_times_norm=fixed_costs[k,n,r].enc()) for k,n,r in FIXED],
        source_profiles=458,finite_candidate_counts=dict(counts),tail_profile_families=37,
        all_size_B_families=37,all_size_D_families=37,records=records,families=families,
        tail_elimination_identities=['B weighted leaf:2-2/2=1','D fork:2-1/2-1/2=1',
                                     'every preceding simple bond:2-1/1=1'],
        rank29_relevant_instances=len(instances),rank29_instances=instances,
        scope='Complete rank necessity for exactly one connected finite external white component touching one of the130 old nonsingular finite unpaid cores at exactly one external root, with a nonzero accepted first-contact profile. Finite-type classification and exact equations treat all ranks; B/D endpoint families are unbounded. Other external components may be orthogonal spectators only. No sealing against additional touching components, branching attachments with multiple touching nodes, dotted cores or arbitrary global completion is claimed. The separate finite instance export uses full-star nodes<=31 for compact dimension29, includes complete corner forests/orbits, and is only conservative necessary input for weights.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite_candidate_counts'],out['tail_profile_families'],out['rank29_relevant_instances'])
