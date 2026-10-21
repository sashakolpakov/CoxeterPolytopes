#!/usr/bin/env python3
"""Complete the35 remaining selected D12 contexts, in the no-short8 branch.

Nonselected and earlier-priority interfaces remain separate obligations.
No universal face rule, forest count or new dimension bound is asserted.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_d12_radical_saturated_half_replacements_20261018 as D

ROOT=D.ROOT; A=D.A; Q=D.Q; T=D.T
TAG='fresh_d12_remaining_selected_half_replacements_20261018'
STEMS=(D.TAG,D.S.TAG,D.R.TAG,A.OVER)


def validate(c,core):
    g,J,du,dv,end=D.R.X.constraints(c,core)
    left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(len(g)))-set(left+right))
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in left+right)
    return g,J,du,dv,end,white,left,right


def extend(g,white,left,right,contacts):
    n=len(g);k=len(contacts)
    h=[row+[F(-int(i==q))for q in contacts]for i,row in enumerate(g)]
    for j,q in enumerate(contacts):
        h.append([F(-int(i==q))for i in range(n)]+[F(2*int(i==j))for i in range(k)])
    ww=white+list(range(n,n+k))
    assert Q.inertia(A.sub(h,ww))==[len(ww),0,0]
    assert all(Q.inertia(A.sub(h,ww+[b]))==[len(ww)+1,0,0]for b in left+right)
    return h,ww


def used(row):return set().union(*(set(r['nodes'])for r in row['roles']))


def whole(c,core,weights):
    g,J,du,dv,end,white,left,right=validate(c,core)
    assert len(c['marked_component'])==12
    remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert all(J[i][i]>=1 for i in remaining)
    possible={i for i in remaining if J[i][i]<=2}
    # In the no-short8 branch a new white contacts only one positive end
    # block. Two contacts already have projection square strictly above2.
    two=[]
    for ee in c['end_sets']:
        for i,j in combinations(sorted(set(ee)&set(remaining)),2):
            lower=J[i][i]+2*J[i][j]+J[j][j]
            assert lower>2
            two.append(dict(nodes=[i,j],projection_lower=lower))
    pay=D.corners(c,g,white,weights)
    alternatives=[p for p in pay if p['marks']!=c['preferred_corner']and p['payment']>=F(1,2)]
    common=dict(index=c['index'],core_key=c['core_key'],gram=g,white=white,
        selected_marks=c['preferred_corner'],selected_factor=c['marked_component'],
        end_projection=end,remaining_inverse_diagonals=[[i,J[i][i]]for i in remaining],
        possible_single_contacts=sorted(possible),two_contact_lower_bounds=two)
    direct=[p for p in alternatives if not used(p)&possible]
    if direct:
        return common|dict(method='protected_distinct_corner',paying_corner=direct[0],
            all_possible_contact_positions_disjoint_from_paying_factors=True)
    for p in alternatives:
        touched=used(p)&possible
        if len(touched)!=1:continue
        q=next(iter(touched))
        if J[q][q]!=2:continue
        h,ww=extend(g,white,left,right,[q])
        after=next(z for z in D.corners(c,h,ww,weights)if z['marks']==p['marks'])
        if after['payment']<F(1,2):continue
        pairs=[]
        for i in sorted(possible-{q}):
            length=min(du[q]+dv[i],dv[q]+du[i])+3
            assert J[q][i]>=0 or length<=8
            pairs.append(dict(other_contact=i,projection_product=J[q][i],
                shortcut_length=length if J[q][i]<0 else None,
                outcome='short8_escape'if J[q][i]<0 else'positive_product_exclusion'if J[q][i]>0 else'orthogonal'))
        return common|dict(method='one_saturated_contact_at_paying_factor',contact=q,
            before= p,after=after,extended_gram=h,extended_white=ww,
            new_root_projection_square='2',new_root_residual_square='0',
            all_other_contact_interactions=pairs,
            core_orthogonal_whites_are_orthogonal_to_new_root=True)
    # The four cases with a unit inverse diagonal at each end cannot be
    # sealed separately. If both ends are active, no-short8 forces their
    # extra roots to be orthogonal. Residual Cauchy forces a shared unit
    # residual and exactly one root per end.
    assert c['index']in(1648,1668,1935,2328)
    ends=[sorted(set(ee)&possible)for ee in c['end_sets']]
    one_sided=[]
    for active in range(2):
        choices=[p for p in alternatives if not used(p)&set(ends[active])]
        assert choices
        one_sided.append(dict(active_end=active,paying_corner=choices[0]))
    mixed=[];survivors=[]
    for i,j in product(*ends):
        length=min(du[i]+dv[j],du[j]+dv[i])+3;assert length<=8 and J[i][j]<0
        rhs=(2-J[i][i])*(2-J[j][j]);gap=J[i][j]**2-rhs;assert gap>=0
        if gap==0:
            assert J[i][i]==J[j][j]==1 and J[i][j]==-1
            survivors.append((i,j))
        mixed.append(dict(nodes=[i,j],two_root_shortcut=length,inverse_diagonals=[J[i][i],J[j][j]],
            inverse_cross=J[i][j],unit_residual_Cauchy_gap=gap))
    assert len(survivors)==1
    i,j=survivors[0];h,ww=extend(g,white,left,right,[i,j])
    assert Q.inertia(h)==[len(g),1,1]
    shared_length=min(du[i]+dv[j],du[j]+dv[i])+4;assert shared_length<=8
    newpay=D.corners(c,h,ww,weights)
    return common|dict(method='both_ends_shared_unit_residual',one_sided_branches=one_sided,
        cross_end_constraints=mixed,unique_possible_contact_pair=[i,j],
        both_native_contacts='1',residual_Gram=[['1','1'],['1','1']],
        same_end_duplicate_root_product='2',common_tail_shortcut_length=shared_length,
        extended_gram=h,extended_white=ww,extended_signature=Q.inertia(h),
        extended_corners=newpay,total=sum(p['payment']for p in newpay),
        all_other_whites_orthogonal_to_extended_core_in_no_short8_branch=True)


def proper(c,core,maps,target):
    g,J,du,dv,end,white,left,right=validate(c,core)
    assert len(c['marked_component'])==11
    records=[]
    for mp in maps:
        image={t:i for t,i in zip(mp['images'],c['marked_component'])}
        missing=sorted(set(range(12))-set(image));assert len(missing)==1
        forced=[(new,i,target[new][old])for new in missing for old,i in image.items()if target[new][old]]
        assert len(forced)==1
        new,q,label=forced[0]
        assert label==3 and q in c['possible_contacts']and J[q][q]>2
        records.append(dict(images=mp['images'],marked_images=mp['marked_images'],missing_target_node=new,
            forced_core_contact=q,required_label=label,projection_lower_bound=J[q][q],
            outcome='proper_completion_exceeds_norm2'))
    return dict(index=c['index'],core_key=c['core_key'],method='proper_completion_projection_exclusion',
        gram=g,white=white,end_projection=end,maps=records,
        unspecified_contacts_not_set_to_zero=True)


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,doc in docs.items():
        assert doc['status'].startswith('PASS')
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old,scout,rank14,over=[docs[t]for t in STEMS]
    weights={tuple(x['pair']):F(x['weight'])for rule in rank14['rules']if rule['kind']=='D'for x in rule['weighted_pairs']}
    patterns=[c for c in over['isolation']['cases']if(c.get('component_order_bound')or 0)>=12 and len(c['marked_component'])<=12]
    assert len(patterns)==len(scout['selected_patterns'])==758
    hits=[];mapcount=0;maprows=[]
    for c,saved in zip(patterns,scout['selected_patterns']):
        assert c['index']==saved['index']
        core=over['isolation']['cores'][c['core_key']]
        maps=T.T.mapped(A.matrix(core['gram']),c['marked_component'],c['preferred_corner'],scout['target'])
        assert maps==saved['maps'];mapcount+=len(maps)
        affected=[m for m in maps if tuple(m['marked_images'])in D.S.PAIRS]
        assert affected==saved['affected_maps']
        maprows.append(dict(index=c['index'],maps=maps,affected_maps=affected))
        if affected:hits.append(c['index'])
    assert mapcount==938 and len(hits)==137
    done={p['index']for p in old['proofs']};remaining={r['index']:r for r in old['remaining_selected_contexts']}
    assert len(done)==102 and len(remaining)==35 and not done&remaining.keys()and done|remaining.keys()==set(hits)
    proofs=[]
    for ix,row in sorted(remaining.items()):
        c=over['isolation']['cases'][ix];core=over['isolation']['cores'][c['core_key']]
        proofs.append(whole(c,core,weights)if row['displayed_order']==12 else proper(c,core,row['affected_maps'],scout['target']))
    counts=Counter(p['method']for p in proofs)
    assert counts==dict(protected_distinct_corner=16,one_saturated_contact_at_paying_factor=8,
        both_ends_shared_unit_residual=4,proper_completion_projection_exclusion=7)
    assert sum(len(p.get('maps',[]))for p in proofs)==28
    paths={Path(__file__),Path(D.__file__),Path(D.R.X.__file__),Path(A.__file__),Path(Q.__file__),Path(T.__file__),Path(T.T.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return D.enc(dict(status='PASS_ALL137_SELECTED_D12_CONTEXTS_IN_NO_SHORT8_BRANCH',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        complete_eligible_patterns=758,complete_native_maps=mapcount,selected_patterns=maprows,
        affected_contexts=137,affected_maps=sum(len(r['affected_maps'])for r in maprows),
        inherited_covered_contexts=102,newly_covered_contexts=35,new_method_counts=counts,new_proofs=proofs,
        selected_contexts_remaining=0,proposed_D12_pairs=sorted(D.S.PAIRS),proposed_weight='1/2',
        proposed_component_cost='65',new_half_coupons=old['new_half_coupons'],new_coupon_costs_unchanged=True,
        independent_short8_and_nonselected_interfaces_open=True,universal_transfer_proved=False,
        new_forest_count=False,new_dimension_bound=None,
        scope='All137 affected selected native contexts, including every proper completion, in the branch with no independent adjacent-black white path of length<=8. Earlier-priority and nonselected interfaces remain open; no universal payment rule.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_method_counts'],flush=True)
