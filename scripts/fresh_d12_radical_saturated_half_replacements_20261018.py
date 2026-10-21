#!/usr/bin/env python3
"""Exact replacements in102 selected D12 contexts, not a universal transfer.

The35 remaining selected contexts and all nonselected interfaces stay open.
The report supplies the arbitrary-support reasoning for the covered contexts.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_d12_long_leaf_support_scout_20261018 as S
import fresh_finite_rank14_existing_payment_transfer_20261018 as R

ROOT=S.ROOT; A=S.A; Q=S.Q; T=S.T; P=S.P; C=S.C
TAG='fresh_d12_radical_saturated_half_replacements_20261018'
EXTRA=('fresh_order17_e6_unit_coupon_zero_transfer_20261018',
       'fresh_a11_cycle_coupon_domination_20261017',
       'fresh_large_finite_d_leaf_half_transfer_20261017',
       'fresh_depth6_finite8_zero_transfer_20261018')
STEMS=tuple(dict.fromkeys((S.TAG,R.TAG,C.TAG,A.OVER,*EXTRA)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x


def basic(c,core):
    g=A.matrix(core['gram']);n=len(g);sc=core['gauge']
    assert all(g[i][i]==2*sc[i]for i in range(n))
    left=c['left'][1];right=[c['right_copy_map'][i]for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right))
    assert Q.inertia(g)==core['signature']
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in left+right)
    actual=Q.inspect_corner(g,sc,left+right,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[k]==c[k]
    assert len(c['marked_component'])==12
    return g,white,left,right


def block(g,nodes,mark):
    n=len(nodes);types=[]
    for kind in 'ABD':
        if kind=='D'and n<4 or kind=='B'and n==1:continue
        im=T.mark_pairs(g,nodes,[mark],n,kind)
        if im:types.append(dict(type=kind+str(n),nodes=nodes,mark=mark,images=im))
    assert len(types)<=1
    if types:return types[0]
    # The only exceptional paying role used here is E6's length-two arm end.
    if n==6 and all(g[i][i]==2 for i in nodes)and all(g[i][j]in(-1,0)for i,j in combinations(nodes,2)):
        adj={i:{j for j in nodes if i!=j and g[i][j]}for i in nodes}
        centers=[i for i in nodes if len(adj[i])==3]
        if len(centers)==1:
            center=centers[0];dist=Q.graph_dist(A.sub(g,nodes),[nodes.index(center)])
            ends=[i for i in nodes if len(adj[i])==1]
            if sorted(dist[nodes.index(i)]for i in ends)==[1,2,2]:
                return dict(type='E6',nodes=nodes,mark=mark,
                    long_endpoint=mark in ends and dist[nodes.index(mark)]==2)
    return dict(type='unclassified',nodes=nodes,mark=mark)


def corners(c,g,white,d14):
    rows=[]
    for marks in [r['preferred_corner']for r in c['all_corners']]:
        ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
        value=F();reason='unused';roles=[];images=[]
        if sig==[len(ids),0,0]:
            parts=A.components(g,ids);aa,bb=[next(p for p in parts if m in p)for m in marks]
            if aa is bb:
                b=block(g,aa,marks[0]);typ=b['type'];roles=[b]
                if typ[0]in'ABD'and typ[1:].isdigit():
                    n=int(typ[1:]);images=T.mark_pairs(g,aa,marks,n,typ[0]);assert images
                    if n<=13:
                        assert C.actual_base((typ[0],n))==16*n*(n-1)//2
                        if typ=='D12'and set(images)&S.PAIRS:
                            assert set(images)<=S.PAIRS;value=F(1,2);reason='proposed_D12_half'
                        else:value=F(1);reason='retained_whole_small_finite_unit'
                    elif typ=='D14':
                        values={d14[p]for p in images};assert len(values)==1
                        value=min(values);reason='retained_whole_D14_weight'
            else:
                roles=[block(g,aa,marks[0]),block(g,bb,marks[1])]
                for sm,la in (roles,roles[::-1]):
                    if sm['type']in('A1','A2')and la['type']=='A11'and set(la['images'])=={(0,),(10,)}:
                        value=F(1,2);reason='new_small_A11_endpoint_half'
                    if la['type']=='D11'and set(la['images'])=={(9,),(10,)}:
                        if sm['type']in('A1','A2','A3','A4','A5'):
                            value=max(value,F(1,2));reason='retained_broad_small_D11_leaf_half'
                        if sm['type']=='E6'and sm['long_endpoint']:
                            value=F(1);reason='retained_E6_long_end_D11_leaf_unit'
        assert sig[1]==0
        rows.append(dict(marks=marks,signature=sig,roles=roles,images=images,payment=value,reason=reason))
    selected=next(r for r in rows if r['marks']==c['preferred_corner'])
    assert selected['reason']=='proposed_D12_half'and selected['payment']==F(1,2)
    assert sum(r['payment']for r in rows)>=1,(c['index'],enc(rows))
    assert any(r['payment']>=F(1,2)for r in rows if r is not selected)
    return rows


def radical(c,core,d14):
    g,white,left,right=basic(c,core);n=len(g)
    assert core['signature']==[n-2,1,1]and len(core['radical'])==1
    z=list(map(F,core['radical'][0]));assert all(sum(g[i][j]*z[j]for j in range(n))==0 for i in range(n))
    # The actual whole-D12 premise forbids every new white contact to it.
    allowed=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert all(z[i]for i in allowed)
    du,dv=Q.graph_dist(g,left),Q.graph_dist(g,right);escapes=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2;assert length<=8
            escapes.append(dict(nodes=[i,j],length=length))
    rows=corners(c,g,white,d14)
    return dict(index=c['index'],method='radical_annihilation',core_key=c['core_key'],
        gram=core['gram'],gauge=core['gauge'],white=white,radical=z,
        remaining_contact_positions=allowed,mixed_sign_shortcuts=escapes,
        corners=rows,total=sum(r['payment']for r in rows))


def saturated(c,core,d14):
    g,white,left,right=basic(c,core);n=len(g)
    assert core['gauge']==[1]*n and core['signature']==[n-1,1,0]
    gg,J,du,dv,ends=R.X.constraints(c,core);assert gg==g
    assert ends['both_ends_shortcut_bound']<=5
    remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert all(J[i][i]>=2 for i in remaining)
    equal=[i for i in remaining if J[i][i]==2]
    assert len(equal)in(0,2)
    assert all(J[i][j]==-2 for i,j in combinations(equal,2))
    branches=[]
    for i in [None,*equal]:
        h=[r.copy()for r in g];ww=white.copy()
        if i is not None:
            h=[r+[F(-int(j==i))]for j,r in enumerate(h)]
            h.append([F(-int(j==i))for j in range(n)]+[F(2)]);ww.append(n)
            assert Q.inertia(A.sub(h,ww))==[len(ww),0,0]
            assert all(Q.inertia(A.sub(h,ww+[b]))==[len(ww)+1,0,0]for b in left+right)
        rows=corners(c,h,ww,d14)
        branches.append(dict(extra_contact=i,extra_norm2_root_count=int(i is not None),
            gram=h,white=ww,corners=rows,total=sum(r['payment']for r in rows)))
    return dict(index=c['index'],method='saturated_lorentz_extensions',core_key=c['core_key'],
        gram=core['gram'],white=white,end_projection=ends,remaining_contacts=remaining,
        inverse_diagonals=[[i,J[i][i]]for i in remaining],saturated_positions=equal,
        two_saturated_products=[[i,j,J[i][j]]for i,j in combinations(equal,2)],
        orthogonal_whites_cannot_attach_saturated_root=True,branches=branches)


def actual_cross(a,b):
    pair=frozenset((a,b));extra=16 if pair==frozenset((('A',1),('A',11)))else 32 if pair==frozenset((('A',2),('A',11)))else 0
    return C.actual_cross(a,b)+extra


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    scout=docs[S.TAG];over=docs[A.OVER];d14={tuple(r['pair']):F(r['weight'])for rule in docs[R.TAG]['rules']if rule['kind']=='D'for r in rule['weighted_pairs']}
    proofs=[];remaining=[];reconstructed=0
    for row in scout['affected_contexts']:
        c=over['isolation']['cases'][row['index']];core=over['isolation']['cores'][c['core_key']]
        maps=T.T.mapped(A.matrix(core['gram']),c['marked_component'],c['preferred_corner'],scout['target'])
        assert [m for m in maps if tuple(m['marked_images'])in S.PAIRS]==row['affected_maps'];reconstructed+=len(row['affected_maps'])
        if len(c['marked_component'])==12 and c['method']=='radical':proofs.append(radical(c,core,d14))
        elif len(c['marked_component'])==12 and c['method']=='lorentz_residual':proofs.append(saturated(c,core,d14))
        else:remaining.append(dict(index=c['index'],method=c['method'],displayed_order=len(c['marked_component']),affected_maps=row['affected_maps']))
    assert len(proofs)==102 and len(remaining)==35 and reconstructed==288
    assert Counter(p['method']for p in proofs)=={'radical_annihilation':94,'saturated_lorentz_extensions':8}
    engine=C.C.C.P.Z.B.B;types=engine.types();rep=engine.representative
    for a in types:
        for b in types:assert C.upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)
    costs=[]
    for r in(1,2):
        a,b=('A',r),('A',11)
        costs.append(dict(small_rank=r,small_roles='all vertices',large_rank=11,large_roles='two endpoints',
            weight='1/2',increment=str(F(2*r,2)),old_actual=str(F(C.actual_cross(a,b),16)),
            new_actual=str(F(actual_cross(a,b),16)),unchanged_upper=str(F(C.upper_cross(rep(a),rep(b)),16))))
    assert [x['new_actual']for x in costs]==['2','3']
    paths={Path(__file__),Path(S.__file__),Path(R.__file__),Path(R.X.__file__),Path(A.__file__),Path(Q.__file__),Path(T.__file__),Path(T.T.__file__),Path(P.__file__),Path(C.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_102_SELECTED_D12_HALF_REPLACEMENTS_WITH_35_SELECTED_CONTEXTS_OPEN',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        selected_affected_contexts=137,selected_affected_maps=reconstructed,covered_contexts=102,
        radical_contexts=94,lorentz_contexts=8,saturated_extra_branches=sum(len(p.get('branches',[]))-1 for p in proofs if 'branches'in p),
        proofs=proofs,remaining_selected_contexts=remaining,proposed_D12_pairs=sorted(S.PAIRS),proposed_weight='1/2',
        proposed_component_cost='65',new_half_coupons=costs,representative_increments='0',simultaneous_type_pairs=len(types)**2,
        geometric_scope='Only the102 explicitly covered selected whole-D12 contexts; no universal face-payment rule.',
        remaining_nonselected_and_later_interfaces_open=True,all_216_literals_inherited_only=True,
        universal_transfer_proved=False,new_forest_count=False,new_dimension_bound=None))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['covered_contexts'],'remaining',len(out['remaining_selected_contexts']),flush=True)
