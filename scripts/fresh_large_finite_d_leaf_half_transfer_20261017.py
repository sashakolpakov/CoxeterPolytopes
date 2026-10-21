#!/usr/bin/env python3
"""Large finite D short-leaf-nine halves and a small-component/D10 coupon.

Exact local interfaces only. The inherited all-label geometric reductions and
residual inverse-rank lemma are stated inputs, not replaced by this census.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_finite_d17_long_pair_support_20261017 as T
import fresh_b14_b15_singleton_completion_transfer_20261017 as B
import fresh_distance8_high_overlap_20261017 as H

ROOT=T.ROOT
TAG='fresh_large_finite_d_leaf_half_transfer_20261017'
OVER='fresh_distance8_high_overlap_20261017'
ADAPT='fresh_distance8_adaptive_leaf_overlap_20261017'
STEMS=(OVER,ADAPT,T.TAG,'fresh_finite8_cycle_replacement_20261017',
       'fresh_finite8_ideal5_local_replacement_20261017',
       'fresh_distance9_supporting_extensions_review_20261017',
       'fresh_ideal5_weighted_dimension66_20261017')


def target(n):
    return T.P.labelled(n,[(i,i+1,3) for i in range(n-3)]+
                       [(n-3,n-2,3),(n-3,n-1,3)])


def degree(g,i):
    return sum(bool(g[i][j]) for j in range(len(g)) if j!=i)


def protected_options(c,g,inv,du,dv):
    out=[]
    for z in c['all_corners']:
        if not z['marked_component']:
            continue
        ans=H.protected_leaf(g,[1]*len(g),inv,z,du,dv)
        if ans:
            out.append(ans)
    return out


def small_coupon(c,g,inv,end_proof):
    """The sixteen one-protected-corner cases, conditional on whole finite D>=20."""
    ns=c['marked_component']; marks=c['preferred_corner']
    leaves=c['protected_marked_leaves']
    assert len(leaves)==1
    b=leaves[0]; a=next(i for i in marks if i!=b)
    # The alternate black is an old full-core leaf attached only to a.
    ap=next(i for i in c['omitted_black'] if g[i][a])
    assert [i for i in range(len(g)) if i!=ap and g[ap][i]]==[a]
    assert g[ap][a]==-1
    bp=next(i for i in c['omitted_black'] if i!=ap)
    assert set(ns)==set(range(len(g)))-{ap,bp}
    alternate=next(z for z in c['all_corners'] if set(z['preferred_corner'])=={ap,b})
    assert alternate['marked_component'] is None
    # The base is already a complete finite D_s with its fork present.
    s=len(ns); assert 13<=s<=17
    maps=T.mapped(g,ns,marks,target(s)); assert len(maps)==2
    long_end=next(i for i in ns if sum(bool(g[i][j]) for j in ns if j!=i)==1 and i not in
                  [ns[z['images'].index(s-2)] for z in maps]+[ns[z['images'].index(s-1)] for z in maps])
    before=next(j for j in ns if j!=long_end and g[long_end][j])
    assert set(c['possible_contacts'])=={ap,bp,long_end,before}
    r=inv[ap][ap]
    assert r in (F(5,4),F(3,2),F(7,4),F(9,4))
    assert inv[long_end][long_end]==1 and inv[before][before]==2
    assert inv[ap][long_end]==F(1,2) and inv[ap][before]==1
    # Every allowed end has positive inverse entries; all multi-supports fail
    # at native lower magnitude one. Contacts meeting both ends were resolved.
    supports=[]
    for end in c['end_sets']:
        for k in range(2,len(end)+1):
            for sub in combinations(end,k):
                q=sum(inv[i][j] for i in sub for j in sub)
                assert q>2
                supports.append(dict(nodes=list(sub),unit_projection_lower_bound=str(q)))
    assert all(inv[i][i]>=1 for i in c['possible_contacts'])
    profiles=[(i,q) for i in c['possible_contacts'] for q in (1,2) if q*inv[i][i]<=2]
    # Label>=5 has squared magnitude >(5/2)>2; diagonals>=1 exclude it.
    if r>2:
        assert all(i!=ap for i,q in profiles)
        component_bound=1; interaction=[]; residual=None
    else:
        assert [(i,q) for i,q in profiles if i==ap]==[(ap,1)]
        rho=2-r; interaction=[]
        for i,q in profiles:
            projection=inv[ap][i]
            residual_product=rho*(2-q*inv[i][i])
            square=q*projection**2
            if i==long_end and q==1:
                assert projection==F(1,2)
                # A native negative product would have magnitude at least1.
                assert (1+projection)**2>residual_product
                outcome='only_native_product_zero'
            elif projection>0:
                assert square>residual_product
                outcome='positive_actual_product_impossible'
            else:
                assert i==bp and projection<0 and square>residual_product
                assert ap in c['omitted_black'] and bp in c['omitted_black']
                outcome='forced_bond_gives_marked_path_length3'
            interaction.append(dict(core_node=i,contact_square=q,
                                    projection_coefficient=str(projection),
                                    projection_square=str(square),
                                    residual_product=str(residual_product),outcome=outcome))
        # G is unique. If a unit tail-end B coexists, two B residuals would
        # be opposites, incompatible with their common -1/2 product with G.
        assert r*r>rho*rho
        gram=[[rho,F(-1,2)],[F(-1,2),F(1)]]
        det=rho-F(1,4); assert det>=0
        if det:
            ri=B.Q.inverse(gram); assert ri[0][0]>=2
            assert all(x>0 for row in ri for x in row)
            if ri[0][0]==2:
                assert ri==[[F(2),F(1)],[F(1),F(3,2)]]
                # At equality E=-2*g-b. Thus a further root touching B
                # but right to G has positive product with E, impossible.
                assert [sum(gram[i][j]*v for j,v in enumerate((-2,-1))) for i in range(2)]==[-1,0]
            residual=dict(gram=B.L.encode(gram),inverse=B.L.encode(ri),
                          conclusion='No indirect B-tail bridge reaches G; at most one terminal orthogonal-white root touches G.')
        else:
            assert [sum(gram[i][j]*v for j,v in enumerate((2,1))) for i in range(2)]==[0,0]
            residual=dict(gram=B.L.encode(gram),positive_kernel=[2,1],
                          conclusion='Every core-orthogonal white is orthogonal to G and B.')
        # If B is absent, the endpoint-inverse rank lemma controls the ENTIRE
        # attached forest, not only immediate neighbours. Two nonempty parts
        # cost at least1, and one of rank4 costs at least4/5>rho.
        assert rho<=F(3,4)<F(4,5) and rho<1
        component_bound=5
    # A whole D_n containing this D_s extends ONLY its existing long arm.
    # Check the canonical D20 representatives, and expose the same fixed D10
    # after removing the nonleaf mark. The report proves the all-rank step.
    exposed=[]
    for row in T.mapped(g,ns,marks,target(20)):
        image=dict(zip(row['images'],ns)); mp={v:k for k,v in image.items()}
        assert mp[b] in (18,19) and mp[a]==9
        part=set(range(10,20)); assert part<=set(image)
        d10=sorted(image[i] for i in part)
        assert len(T.mapped(g,d10,[b, d10[0]],target(10)))==2
        assert all(g[ap][i]==0 for i in d10)
        exposed.append(dict(map=row['images'],whole_D10_nodes=d10,
                            D10_mark=b,small_mark=ap,small_component_order_bound=component_bound))
    assert len(exposed)==2
    return dict(index=c['index'],core_key=c['core_key'],core_gram=B.L.encode(g),
                core_inverse=B.L.encode(inv),selected_marks=marks,selected_base_nodes=ns,
                old_opposite_leaf=ap,other_omitted_black=bp,long_end=long_end,
                preceding_long_node=before,old_leaf_inverse_diagonal=str(r),
                end_projection_proof=end_proof,multisupport_exclusions=supports,
                complete_single_profiles=[dict(node=i,squared_magnitude=q) for i,q in profiles],
                interactions_with_G=interaction,two_residual_certificate=residual,
                small_component_order_bound=component_bound,exposed_D10=exposed,
                retained_total='1',payments=['large_D_short_leaf9_half','small_component_D10_short_leaf_half'])


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(s,name)
    over=data[OVER]; rows=[]; coupons=[]; counts=Counter()
    cases=over['isolation']['cases']; assert len(cases)==2614
    for c in cases:
        if c.get('component_order_bound') is not None:
            assert c['component_order_bound']<=17
            counts['whole_component_at_most17']+=1;continue
        if c['method']=='nine_context_adaptive_finite_payment':
            counts['adaptive_simple_leaf_neighbor_degree2']+=1;continue
        assert c['method']=='finite_leaf9_or_ideal_feature5'
        core=over['isolation']['cores'][c['core_key']]
        g,inv,du,dv,end_proof=B.constraints(c,core)
        options=protected_options(c,g,inv,du,dv)
        assert options
        if c['forced_finite']:
            assert c['protected_simple_leaves']
            forbidden={i for i in c['marked_component'] if i not in c['possible_contacts'] or inv[i][i]>2}
            for leaf in c['protected_simple_leaves']:
                adjacent=[j for j in c['marked_component'] if j!=leaf and g[leaf][j]]
                assert len(adjacent)==1 and leaf in forbidden and adjacent[0] in forbidden
                assert sum(bool(g[adjacent[0]][j]) for j in c['marked_component'] if j!=adjacent[0])==2
            method='protected_simple_leaf_neighbor_degree2'
        elif len(options)>=2:
            assert len({tuple(z['preferred_corner']) for z in options})>=2
            method='two_protected_leaf9_corners'
        else:
            coupons.append(small_coupon(c,g,inv,end_proof));method='new_small_component_D10_half'
        counts[method]+=1
        rows.append(dict(index=c['index'],method=method,protected_options=options))
    assert counts==Counter(whole_component_at_most17=2531,adaptive_simple_leaf_neighbor_degree2=9,
                          protected_simple_leaf_neighbor_degree2=18,two_protected_leaf9_corners=40,
                          new_small_component_D10_half=16)
    adaptive=[]
    for c in data[ADAPT]['cases']:
        p=c['payment']
        leaves=[p] if 'leaf' in p else p['leaf_alternatives'] if 'leaf_alternatives' in p else [p['if_far_present']]
        assert all(z['distance']==9 and z['base_order']>=10 for z in leaves)
        if 'if_far_present' in p:assert p['whole_component_bound']==16
        adaptive.append(dict(index=c['index'],method=c['method'],leaves=leaves))
    # The cycle theorem's D_m single profile has distance L-1. To need9,
    # L=10 is necessarily the longer arc: shorter P is6..9, hence m<=19.
    cycle=[dict(shorter_arc=p,longer_arc=10,finite_D_order=p+10) for p in (6,7,8,9)]
    assert max(r['finite_D_order'] for r in cycle)==19
    paths=[Path(__file__),Path(T.__file__),Path(B.__file__),Path(H.__file__),Path(B.Q.__file__),Path(T.P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_LARGE_FINITE_D_SHORT_LEAF9_HALF_TRANSFER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                changed_rule='Whole finite D_n with n>=20: short-fork-leaf pairs at distance9 have weight1/2; all other old payments unchanged.',
                added_coupon='Weight1/2 between any mark of a whole finite component of order<=5 and a short-fork leaf of a separate whole finite D10.',
                new_large_D_cost='8*n-33',routing=dict(counts),protected_cases=rows,
                coupon_cases=coupons,adaptive_interfaces=adaptive,cycle_D_leaf9_orders=cycle,
                triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
                scope='Exact affected high-core interfaces and residual inequalities, with an all-rank written whole-component proof. Inherited all-label core classifications, universal transfers and connected inverse-rank lemma are not rerun. No global finite maximum or new dimension theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['routing'],flush=True)
