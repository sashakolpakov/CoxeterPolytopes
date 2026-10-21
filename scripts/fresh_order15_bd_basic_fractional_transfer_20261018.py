#!/usr/bin/env python3
"""A B15 basic deletion and a universal economical D15 fractional transfer.

The new D4 end profile, complete native maps and whole-component roles are
checked exactly. One new A1/D14 role costs1/16; forest counting is separate.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order14_d_near_fork_half_transfer_20261018 as N
import fresh_order15_bd_terminal_coupon_transfer_20261018 as M
import fresh_rank13_contextual_coupon_rebalance_20261018 as R

ROOT=N.ROOT;L=N.L;D=N.D;T=N.T;P=N.P;A=N.A;Q=N.Q;X=M.X
TAG='fresh_order15_bd_basic_fractional_transfer_20261018'
FRONT='fresh_order14_contextual_ambient59_frontier_20261018'
F15='fresh_finite_d15_three_singleton_half_transfer_20261018'
BC=N.BC;LEAF=N.LEAF;BELOW=N.BELOW
DELTA=F(1,16);BDROP={(3,12)};DCHANGES={(2,13),(2,14),(8,11),(8,12)}
STEMS=tuple(dict.fromkeys((FRONT,N.TAG,M.TAG,R.TAG,F15,BC,LEAF,BELOW,
    X.SMALL,X.OVER,X.COMBINED,L.TAG,N.C.TAG,
    'fresh_finite7_ideal5_triangle_payments_20261017',
    'fresh_depth6_finite7_adaptive_transfer_20261018')))
enc=N.enc


def affine():
    g=L.gram(15,[(i,i+1,1) for i in range(10)]+[(0,11,1),(0,12,1),(10,13,1),(10,14,1)])
    nu=s.Matrix([1]*11+[s.Rational(1,2)]*4)
    assert g*nu==s.zeros(15,1) and g.rank()==14
    return g,nu


def profiles():
    g,_=affine();K=g.extract([1,0,11,12],[1,0,11,12]).inv()
    assert K==s.Matrix([[1,1,s.Rational(1,2),s.Rational(1,2)],
        [1,2,1,1],[s.Rational(1,2),1,1,s.Rational(1,2)],
        [s.Rational(1,2),1,s.Rational(1,2),1]])
    t,u,v,w,z=s.symbols('t u v w z',nonnegative=True)
    q=s.Matrix([u,v,w,z]);e=(q.T*K*q)[0];h=t+u+v+(w+z)/2
    near=[dict(name='outer',contacts={2:t,1:u,0:v,11:w,12:z},height=h,c0=2-e,moment=0),
        dict(name='inner3',contacts={2:t,3:1,1:u,0:v,11:w,12:z},height=h+1,c0=1-e,moment=1),
        dict(name='inner4',contacts={2:t,3:s.sqrt(2)},height=t+s.sqrt(2),c0=0,moment=s.sqrt(2)),
        dict(name='next3',contacts={2:t,4:1},height=t+1,c0=0,moment=2)]
    far=[dict(name='isolated',contacts={13:t},height=t/2,c0=2,moment=0),
        dict(name='mate3',contacts={13:t,14:1},height=(t+1)/2,c0=1,moment=1),
        dict(name='mate4',contacts={13:t,14:s.sqrt(2)},height=(t+s.sqrt(2))/2,c0=0,moment=s.sqrt(2)),
        dict(name='center3',contacts={13:t,10:1},height=t/2+1,c0=0,moment=2)]
    return near,far,(t,u,v,w,z),K


def formulas():
    g,nu=affine();W=[i for i in range(15) if i not in (2,13)]
    J=g.extract(W,W).inv();sc=s.diag(nu[2],nu[13]);V=-g.extract(W,[2,13])*sc
    assert sc*g.extract([2,13],[2,13])*sc-V.T*J*V==s.Matrix([[1,-1],[-1,1]])/10
    near,far,variables,K=profiles();inverses=[];rows=[];cross=[]
    for ids in [list(range(3,11))+[13,14],[14]+list(reversed(range(11)))+[11,12]]:
        diag=list(g.extract(ids,ids).inv().diagonal())
        assert diag[:2]==[1,2] and all(x>2 for x in diag[2:])
        inverses.append(dict(nodes=ids,diagonal=diag))
    for end,ps,other in [('near',near,V[:,1]),('far',far,V[:,0])]:
        for p in ps:
            x=s.Matrix([p['contacts'].get(i,0) for i in W])
            assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-p['height'])==0
            assert s.simplify(10*(other.T*J*x)[0]-p['moment'])==0
            assert s.simplify(2-(x.T*J*x)[0]-p['c0']-p['moment']**2/10)==0
            rows.append(dict(end=end,**p,cost=s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2)))
    for p,q in product(near,far):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([q['contacts'].get(i,0) for i in W])
        value=s.simplify((x.T*J*y)[0]);assert s.simplify(value-p['moment']*q['moment']/10)==0
        cross.append(dict(near=p['name'],far=q['name'],projection=value))
    h,k,l,m,c0,d0,p,q,r=s.symbols('h k l m c0 d0 p q r',nonzero=True);alpha=s.Rational(1,10)
    G=s.Matrix([[alpha,-alpha,-h+l*alpha,-m*alpha],[-alpha,alpha,-l*alpha,-k+m*alpha],
        [-h+l*alpha,-l*alpha,c0+l*l*alpha-p,-l*m*alpha-q],[-m*alpha,-k+m*alpha,-l*m*alpha-q,d0+m*m*alpha-r]])
    identity=1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/10+(p/h**2-2*q/(h*k)+r/k**2)/10
    assert s.factor(G.det()/(h*h*k*k)-identity)==0
    mult=[]
    for i,j in combinations(range(4),2):
        lower=K[i,i]+K[j,j]+2*K[i,j];assert lower>=3
        mult.append(dict(nodes=[i,j],projection_lower=lower))
    t,u,v,w,z=variables;cost=lambda p:s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2)
    assert cost(near[0]).subs({t:0,u:0,v:0,w:1,z:0})==4
    assert cost(near[0]).subs({t:0,u:1,v:0,w:0,z:0})==1
    assert cost(near[1]).subs({t:0,u:0,v:0,w:0,z:0})==3
    assert cost(far[0]).subs(t,1)==cost(far[1]).subs(t,0)==8
    assert cost(far[0]).subs(t,s.sqrt(2))==cost(far[2]).subs(t,0)==cost(far[3]).subs(t,0)==4
    return dict(resistance=10,near_D4_inverse=K.tolist(),large_inverse_diagonals=inverses,
        profiles=rows,cross_identities=cross,determinant_identity=identity,multiple_outer_exclusions=mult,
        near_extreme_cost=4,near_other_maximum=3,far_extreme_cost=8,far_other_maximum=4,
        far_nonextreme_strict_gap='1/5',nonextreme_near_weighted_energy_maximum=1,
        nonzero_far_spectator_lower_bound=2,both_extreme_energy='1/2')


def immediate(docs):
    f,_=affine();rows=[]
    for name,node,mag,kind,mapping in [('inner4',3,s.sqrt(2),'B',list(range(10,2,-1))+[15,13,14]),
        ('next3',4,1,'D',list(range(10,3,-1))+[13,14,3,15])]:
        g=2*s.eye(16);g[:15,:15]=f;g[15,node]=g[node,15]=-mag
        canonical,_,_=L.affine(kind);gc=g.extract(mapping,mapping);assert gc==canonical
        ker=gc.nullspace()[0];ker=ker if ker[0]>0 else -ker;assert all(x>0 for x in ker)
        pair=tuple(sorted((mapping.index(13),mapping.index(15))))
        rs=docs[BC]['rules'] if kind=='B' else docs[L.TAG]['factor_rules']
        weights={tuple(r['pair']):F(r['weight']) for r in next(r for r in rs if r['kind']==kind)['weighted_pairs']}
        assert weights[pair]==(0 if kind=='B' else F(1,2))
        rows.append(dict(profile=name,type='affine_'+kind+'10',nodes=mapping,marks=[13,15],
            gram=gc.tolist(),positive_kernel=list(ker),canonical_pair=pair,weight=weights[pair],
            payment='independent_BC11_zero_theorem' if kind=='B' else 'retained_half_plus_selected_15/16'))
    return rows


def whole(docs):
    f,_=affine();rows=[]
    for name,omit,attach,marks in [('near_empty',2,11,[13,15]),('far_isolated_empty',13,13,[2,15]),('far_mate_empty',13,14,[2,15])]:
        g=2*s.eye(16);g[:15,:15]=f;g[15,attach]=g[attach,15]=-1
        ids=[i for i in range(16) if i!=omit];assert L.inertia(g.extract(ids,ids))==[15,0,0]
        gg=[[F(x) for x in r] for r in g.tolist()];bs=D.finite_blocks(gg,ids)
        aa,bb=[next(bl for bl in bs if m in bl['nodes']) for m in marks];roles=[]
        if aa is bb:
            assert aa['type']=='D15';images=T.mark_pairs(gg,aa['nodes'],marks,15,'D');assert images==[(0,10)]
            assert next(r['weight'] for r in docs[F15]['weighted_pairs'] if r['pair']==[0,10])=='1'
            roles=[dict(type='D15',pair_images=images)];pay=F(1)
        else:
            for bl in (aa,bb):
                mark=next(m for m in marks if m in bl['nodes'])
                roles.append(dict(type=bl['type'],nodes=bl['nodes'],mark=mark,
                    images=T.mark_pairs(gg,bl['nodes'],[mark],len(bl['nodes']),bl['type'][0])))
            by={r['type']:r for r in roles}
            if name=='near_empty':
                assert set(by)=={'D5','D10'} and by['D5']['images']==[(0,)] and by['D10']['images']==[(8,),(9,)]
                pay=F(1,2)
            else:
                assert set(by)=={'A1','D14'} and by['D14']['images']==[(9,)]
                pay=DELTA
        rows.append(dict(case=name,omitted_mark=omit,marks=marks,gram=g.tolist(),components=bs,roles=roles,
            distinct_payment=pay,selected_payment=1-DELTA,total=1-DELTA+pay))
    return rows


def b_support(docs,target):
    small=docs[X.SMALL];over=docs[X.OVER]
    patterns=[p for p in small['patterns'] if p['order_bound']>=15 and len(p['nodes'])<=15]
    assert len(patterns)==536
    saved={r['index']:r for r in docs[M.TAG]['B_reconstructed_support']};success=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if maps:
            assert maps==saved[p['index']]['embeddings'] and pairs==list(map(tuple,saved[p['index']]['marked_images']))
            assert not set(pairs)&BDROP
            success.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        else:assert p['index'] not in saved
    assert len(success)==len(saved)==215 and sum(len(r['embeddings']) for r in success)==548
    protected=[]
    for c in over['isolation']['cases']:
        if c['method']!='finite_leaf9_or_ideal_feature5':continue
        core=over['isolation']['cores'][c['core_key']];g=A.matrix(core['gram'])
        labels=T.T.labels(g,c['marked_component']);assert labels is not None
        marks=[c['marked_component'].index(i) for i in c['preferred_corner']]
        maps=P.embeddings(labels,target);pairs=sorted({tuple(sorted(mp[i] for i in marks)) for mp in maps})
        assert not set(pairs)&BDROP
        protected.append(dict(index=c['index'],base_nodes=c['marked_component'],labels=labels,marks=marks,maps=maps,marked_images=pairs))
    assert len(protected)==74
    assert docs[M.TAG]['complete_high_overlap_routing']==dict(low=1691,protected=83,too_large=304,eligible=536)
    return dict(eligible_patterns=536,successful_patterns=215,successful_maps=548,
        selected_support=success,changed_selected_maps=0,complete_protected_bases=74,
        protected_support=protected,changed_protected_maps=0,whole_high_overlap_routing=docs[M.TAG]['complete_high_overlap_routing'])


def coupon(g,ids,marks):
    bs=D.finite_blocks(g,ids);aa,bb=[next(b for b in bs if m in b['nodes']) for m in marks]
    if aa is bb:return []
    rows=[]
    for sm,la,m,n in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='D14' and T.mark_pairs(g,la['nodes'],[n],14,'D')==[(9,)]:
            rows.append(dict(small=sm,small_mark=m,large=la,large_mark=n,large_images=[(9,)],increment=DELTA))
    assert len(rows)<=1;return rows


def literal(docs,targets):
    tags=N.S.B.C.R.B.C.L.STEMS;all_docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in tags}
    prior={R.key(r['source'],r['case']):r for r in docs[R.TAG]['combined_literal_interface']['rows']}
    assert docs[N.TAG]['current_literal_interface']['changed']==[]
    rows=[];changes=[]
    for stem,index,g,white,cs in T.literal_inputs(all_docs):
        old=prior[R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
            cp=[];images={};debit=F()
            if sig==[len(ids),0,0]:cp=coupon(g,ids,marks)
            elif sig[1]==0:
                block=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if block and len(block)==15:
                    labels=T.T.labels(g,block)
                    if labels is not None:
                        for kind,target in targets.items():
                            ps=sorted({tuple(sorted(mp[block.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                            if ps:
                                images[kind]=ps;changed=BDROP if kind=='B' else DCHANGES
                                if set(ps)&changed:
                                    assert set(ps)<=changed;debit=F(1) if kind=='B' else DELTA
            inc=DELTA*len(cp);value=F(old['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,images=images,debit=debit,coupons=cp,increment=inc,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['debit'] or c['increment'] for c in corners):changes.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changes,rows=rows),tags


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert docs[FRONT]['current_order14_family_costs']==dict(A=70,B=80,C=80,D=80)
    assert docs[M.TAG]['new_family_costs']==dict(A=75,B=92,C=91,D=92)
    rules=[];targets={}
    for old in docs[M.TAG]['rules']:
        kind=old['kind'];target=old['target'];targets[kind]=target
        weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']};changes=BDROP if kind=='B' else DCHANGES
        assert all(weights[p]==1 for p in changes)
        for p in changes:weights[p]=F() if kind=='B' else 1-DELTA
        assert sum(weights.values())==(91 if kind=='B' else F(367,4))
        _,dist,_=P.distance_data(target);autos=P.embeddings(target,target)
        assert all({tuple(sorted(mp[i] for i in p)) for p in changes}==changes for mp in autos)
        assert all(dist[i][j]==9 for i,j in changes)
        assert all(weights[tuple(p)]==1 for p in old['retained_triangle']+old['retained_short8'])
        rules.append(dict(kind=kind,target=target,changed_pairs=sorted(changes),old_weight='1',
            new_weight='0' if kind=='B' else '15/16',new_cost=str(sum(weights.values())),
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
            retained_triangle=old['retained_triangle'],retained_short8=old['retained_short8']))
    g,_=affine();assert targets['D']==[[3 if i!=j and g[i,j] else 0 for j in range(15)] for i in range(15)]
    assert docs[LEAF]['added_coupon']=='Weight1/2 between any mark of a whole finite component of order<=5 and a short-fork leaf of a separate whole finite D10.'
    assert N.S.actual_cross(('D',5),('D',10))==40 and N.S.actual_cross(('A',1),('D',14))==16
    support=b_support(docs,targets['B']);lit,tags=literal(docs,targets)
    paths={Path(__file__),Path(N.__file__),Path(M.__file__),Path(R.__file__),Path(L.__file__),Path(T.__file__),Path(T.T.__file__),Path(P.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/f'results/{tag}.json' for tag in tags)
    return enc(dict(status='PASS_ORDER15_BD_BASIC_FRACTIONAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        rules=rules,B_complete_support=support,D_formulas=formulas(),D_immediate_corners=immediate(docs),
        D_whole_finite_corners=whole(docs),current_literal_interface=lit,
        new_family_costs=dict(A='75',B='91',C='91',D='367/4'),new_order15_maximum='367/4',new_order15_envelope='367/60',
        new_coupon=dict(small='whole_A1',large='whole_D14',large_role='long_arm_vertex_two_edges_from_fork_center',
            canonical_large_node=9,weight=str(DELTA),eligible_marked_pairs=1,component_pair_increment=str(DELTA),
            old_component_pair_cost='2',new_component_pair_cost='33/16'),
        valid_fractional_family=dict(delta_range='0<delta<=1/2',selected_weight='1-delta',
            new_coupon='delta',ideal_D_cost='92-4*delta',finite_component_increment='delta',chosen_delta=str(DELTA)),
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,finite_budget_recomputed=False,
        scope='Complete B selected/protected native domains plus a universal D4-end profile/spectator proof. One explicit A1/D14 increment1/16. Historical actual-face routing inherited,216 literal interfaces reconstructed; no finite forest or new ambient bound asserted.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_family_costs'],'changed literals',len(out['current_literal_interface']['changed']),flush=True)
