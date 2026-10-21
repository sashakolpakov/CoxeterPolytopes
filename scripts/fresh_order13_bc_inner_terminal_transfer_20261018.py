#!/usr/bin/env python3
"""Universal order13 B71.5/C71 with one absorbed A1/B12 half.

Complete physical one-root profiles, whole-factor alternatives and current
literal/cost interfaces. No selected finite census substitutes for support.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order12_bc_near_terminal_half_transfer_20261018 as H
import fresh_order13_bc_sparse_terminal_transfer_20261018 as O
import fresh_order15_d16_d17_ambient59_frontier_20261018 as C

ROOT=H.ROOT;L=H.L;K=H.K;D=H.D;P=H.P;T=C.T;A=C.A;Q=C.Q
TAG='fresh_order13_bc_inner_terminal_transfer_20261018'
BC11='fresh_order11_bc_zero_transfer_20261018'
D11='fresh_order11_bcd_half_pair_transfer_20261018'
CHANGES={'B':{(8,11),(8,12)},'C':{(2,10)}}
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,H.TAG,O.TAG,BC11,D11,
    'fresh_order13_d_center_leaf_half_transfer_20261018',
    'fresh_finite7_ideal5_triangle_payments_20261017')))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,s.Basic):assert not x.has(s.Float);return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def affine(kind):
    if kind=='B':
        g=L.gram(13,[(i,i+1,s.sqrt(2) if i==9 else 1) for i in range(10)]+[(0,11,1),(0,12,1)])
        nu=s.Matrix([1]*10+[1/s.sqrt(2),s.Rational(1,2),s.Rational(1,2)])
    else:
        g=L.gram(13,[(i,i+1,s.sqrt(2) if i in (0,11) else 1) for i in range(12)])
        nu=s.Matrix([1/s.sqrt(2)]+[1]*11+[1/s.sqrt(2)])
    assert g*nu==s.zeros(13,1) and g.rank()==12
    return g,nu


def near(mark,outer_near,outer_terminal,first,second):
    t=s.symbols('t',nonnegative=True)
    rows=[dict(name='isolated',contacts={mark:t},height=t,moment=0,c0=s.Integer(2),maximum=2,lower_t=1)]
    for name,node,height,maximum in [('outer_near',outer_near,1,1),('outer_terminal',outer_terminal,1/s.sqrt(2),2)]:
        rows.append(dict(name=name,contacts={mark:t,node:1},height=t+height,moment=0,c0=s.Integer(1),maximum=maximum,lower_t=0))
    for name,node,height in [('end3',None,0),('end3_outer_near',outer_near,1),('end3_outer_terminal',outer_terminal,1/s.sqrt(2))]:
        contacts={mark:t,first:1}
        if node is not None:contacts[node]=1
        rows.append(dict(name=name,contacts=contacts,height=t+1+height,moment=1,c0=s.Integer(node is None),maximum=3,lower_t=0))
    rows.extend([dict(name='end4',contacts={mark:t,first:s.sqrt(2)},height=t+s.sqrt(2),moment=s.sqrt(2),c0=s.Integer(0),maximum=2,lower_t=0),
        dict(name='next3',contacts={mark:t,second:1},height=t+1,moment=2,c0=s.Integer(0),maximum=4,lower_t=0)])
    return rows,t


def formulas(kind):
    f,nu=affine(kind);marks=(8,11) if kind=='B' else (2,10);R=10 if kind=='B' else 8
    W=[i for i in range(13) if i not in marks];J=f.extract(W,W).inv().applyfunc(s.simplify)
    va,vb=[-nu[i]*f.extract([i],W).T for i in marks]
    old=s.diag(*[nu[i] for i in marks])*f.extract(marks,marks)*s.diag(*[nu[i] for i in marks])
    assert (old-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)).applyfunc(s.simplify)==s.Matrix([[1,-1],[-1,1]])/R
    left,t=near(8,9,10,7,6) if kind=='B' else near(2,1,0,3,4)
    right,_=H.leaf(11,12,0) if kind=='B' else near(10,11,12,9,8)
    identities=[];bounds=[];u=s.symbols('u',nonnegative=True)
    for side,ps,far in [('left',left,vb),('right',right,va)]:
        for p in ps:
            x=s.Matrix([p['contacts'].get(i,0) for i in W]);h=p['height'];m=s.sympify(p['moment']);c0=p['c0']
            assert s.simplify(sum(nu[i]*z for i,z in p['contacts'].items())-h)==0
            assert s.simplify(R*(far.T*J*x)[0]-m)==0
            assert s.simplify(2-(x.T*J*x)[0]-c0-m*m/R)==0
            certificate=K.positive_polynomial(s.expand((p['maximum']*h*h-2*m*h-c0).subs(t,u+p['lower_t'])),(u,))
            identities.append(dict(side=side,**p,cost=s.factor(2*m/h+c0/h**2)))
            bounds.append(dict(side=side,name=p['name'],certificate=certificate))
    crosses=[]
    for p,q in product(left,right):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([q['contacts'].get(i,0) for i in W])
        assert s.simplify((x.T*J*y)[0]-s.sympify(p['moment'])*q['moment']/R)==0
        crosses.append(dict(left=p['name'],right=q['name'],product=s.simplify((x.T*J*y)[0])))
    domains=[('D10',list(range(8))+[11,12],7,6),('B12',[i for i in range(13) if i!=11],12,0)] if kind=='B' else [
        ('B10',list(range(3,13)),3,4),('B10',list(range(10)),9,8)]
    principal=[]
    for name,ids,end,nxt in domains:
        inv=f.extract(ids,ids).inv().applyfunc(s.simplify);diag=dict(zip(ids,inv.diagonal()))
        assert diag[end]==1 and diag[nxt]==2 and all(v>2 for i,v in diag.items() if i not in (end,nxt))
        principal.append(dict(type=name,nodes=ids,inverse_diagonal=diag,endpoint=end,next_node=nxt))
    outer=[9,10] if kind=='B' else [0,1]
    outer_inv=f.extract(outer,outer).inv().applyfunc(s.simplify)
    assert list(outer_inv.diagonal())==[1,1]
    assert (s.ones(2,1).T*outer_inv*s.ones(2,1))[0]>2
    if kind=='B':
        costs={p['name']:s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2) for p in left}
        assert costs['isolated'].subs(t,1)==2 and costs['outer_terminal'].subs(t,0)==2 and costs['end3'].subs(t,0)==3
        for name in ['outer_near','end3_outer_near','end3_outer_terminal']:
            assert (costs[name].subs(t,0)-2).is_negative
        for name in ['isolated','outer_terminal','end3']:
            lower=s.sqrt(2) if name=='isolated' else 1
            assert (costs[name].subs(t,lower)-2).is_negative
        for p in left:
            if p['name'] in ('isolated','end4','next3'):continue
            assert s.diff(2*p['moment']/p['height']+p['c0']/p['height']**2,t).is_nonpositive
        for p in right:
            lower=s.sqrt(2) if p['name']=='isolated' else 1 if p['name']=='mate3' else 0
            maximum=3 if p['name']=='mate3' else 4
            K.positive_polynomial(s.expand((maximum*p['height']**2-2*p['moment']*p['height']-p['c0']).subs(t,u+lower)),(u,))
        consequences=dict(nonimmediate_near_maximum=3,far_nonextreme_maximum=4,far_extreme_height=s.Rational(1,2),
            far_extreme_cost=8,near_cost_at_least_two=[dict(name='isolated',t=1,height=1,cost=2),
                dict(name='outer_terminal',t=0,height=1/s.sqrt(2),cost=2),dict(name='end3',t=0,height=1,cost=3)],
            spectator_energy_maximum=1,single_support_lower=2,two_support_strict_lower=s.Rational(9,8))
    else:
        assert max(p['maximum'] for p in left if p['name']!='next3')==3
        consequences=dict(no_next3_strict_gap=s.Rational(1,4),next3_payment='whole_affine_B11_retained_distance7_unit',
            triangle_allowed_contacts_left=[0,1,2,3,4],triangle_allowed_contacts_right=[8,9,10,11,12])
        assert not set(consequences['triangle_allowed_contacts_left'])&set(consequences['triangle_allowed_contacts_right'])
    return dict(gram=f.tolist(),kernel=list(nu),marks=marks,white=W,resistance=R,principal_inverse_checks=principal,
        strict_outer_B2_nodes=outer,outer_inverse=outer_inv.tolist(),end_identities=identities,bounds=bounds,cross_identities=crosses,
        determinant_identity=K.determinant(R),consequences=consequences)


def immediate(docs):
    rows=[]
    for kind,name,omitted,node,mag,mark,target_kind,weight in [('B','end4',8,7,s.sqrt(2),11,'B',F()),
        ('B','next3',8,6,1,11,'D',F(1,2)),('C','left_next3',2,4,1,10,'B',F(1)),('C','right_next3',10,8,1,2,'B',F(1))]:
        f,_=affine(kind);g=2*s.eye(14);g[:13,:13]=f;g[13,node]=g[node,13]=-mag
        ids=[i for i in range(14) if i!=omitted];ns=next(b for b in L.components(g,ids) if 13 in b)
        assert len(ns)==11 and mark in ns;gram=g.extract(ns,ns);z=gram.nullspace()[0]
        if z[0]<0:z=-z
        assert all(v>0 for v in z)
        target,_,_=L.affine(target_kind);maps=P.embeddings(D.labels(gram),D.labels(target));assert maps
        pairs=sorted({tuple(sorted(mp[ns.index(i)] for i in (13,mark))) for mp in maps})
        current=next(r for r in docs[BC11]['rules'] if r['kind']=='B') if target_kind=='B' else next(r for r in docs[D11]['factor_rules'] if r['kind']=='D')
        ws={tuple(r['pair']):F(r['weight']) for r in current['weighted_pairs']}
        assert all(ws[p]==weight for p in pairs)
        distances=[L.distances(L.graph(target),i)[j] for i,j in pairs]
        if kind=='C':assert set(distances)=={7}
        rows.append(dict(kind=kind,profile=name,gram=g.tolist(),component=ns,kernel=list(z),marked_pairs=pairs,
            distances=distances,current_weight=str(weight),payment='independent_BC11_zero_first' if not weight else 'distinct_retained_whole_factor'))
    return rows


def finite_corners():
    f,_=affine('B');rows=[]
    for name,node in [('isolated',11),('mate3',12)]:
        g=2*s.eye(14);g[:13,:13]=f;g[13,node]=g[node,13]=-1
        ids=[i for i in range(14) if i!=11];parts=L.components(g,ids);L.pd(g.extract(ids,ids))
        if name=='isolated':
            assert sorted(map(len,parts))==[1,12] and [13] in parts
            b=next(p for p in parts if 8 in p)
            # Intrinsic third-from-doubled-terminal role, independently of orientation.
            adj=L.graph(g.extract(b,b));j=b.index(8);heavy=[i for i in adj if len(adj[i])==1 and any(s.simplify(g[b[i],b[k]]**2)==2 for k in adj[i])]
            assert len(heavy)==1 and L.distances(adj,j)[heavy[0]]==2
            pay=F(1,2);reason='whole_A1/B12_third_from_doubled_terminal_half'
        else:
            assert len(parts)==1 and len(parts[0])==13
            adj=L.graph(g.extract(ids,ids));assert max(map(len,adj.values()))==2
            assert L.distances(adj,ids.index(8))[ids.index(13)]==10
            assert C.actual_base(('B',13))==16*78
            pay=F(1);reason='whole_finite_B13_distance10_unit'
        rows.append(dict(profile=name,gram=g.tolist(),paying_corner=[8,13],components=parts,payment=str(pay),reason=reason,
            selected_half='1/2',total=str(pay+F(1,2))))
    return rows


def coupon(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes']) for m in marks];out=[]
    if a is b:return out
    for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='B12':
            images=T.mark_pairs(g,la['nodes'],[n],12,'B')
            if images==[(9,)]:out.append(dict(small=sm,large=la,marks=[m,n],images=images,weight='1/2'))
    assert len(out)<=1;return out


def literal(docs,targets):
    old={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']};rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        before=old[C.R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids));debit=F();cp=[];images={}
            if sig==[len(ids),0,0]:cp=coupon(g,C.D.finite_blocks(g,ids),marks)
            elif sig[1]==0:
                b=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if b and len(b)==13:
                    labels=T.T.labels(g,b)
                    if labels is not None:
                        for kind,target in targets.items():
                            ps=sorted({tuple(sorted(mp[b.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
                            if ps:
                                images[kind]=ps
                                if set(ps)&CHANGES[kind]:assert set(ps)<=CHANGES[kind];debit=F(1,2) if kind=='B' else F(1)
            inc=F(len(cp),2);pay=F(before['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert pay>=0
            corners.append(dict(marks=marks,images=images,debit=debit,coupons=cp,increment=inc,retained_payment_lower_bound=pay))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['debit'] or c['increment'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in docs.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert docs[O.TAG]['simultaneous_order13_family_costs']==dict(A='65',B='145/2',C='72',D='71')
    assert docs['fresh_order13_d_center_leaf_half_transfer_20261018']['new_cost']=='71'
    assert F(docs[C.TAG]['ideal']['current_envelopes']['12'])==F(145,26)
    fs={kind:formulas(kind) for kind in ('B','C')};rules=[];targets={}
    for kind in ('B','C'):
        old=next(r for r in docs[O.TAG]['records'] if r['kind']==kind);ws={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
        f,_=affine(kind);target=D.labels(f);targets[kind]=target;adj=L.graph(f)
        assert target==old['target'] and all(ws[p]==1 for p in CHANGES[kind])
        new={p:(F(1,2) if kind=='B' else F()) if p in CHANGES[kind] else w for p,w in ws.items()}
        autos=P.embeddings(target,target);assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in autos)
        assert all(L.distances(adj,i)[j]==(9 if kind=='B' else 8) for i,j in CHANGES[kind])
        assert all(new[p]==ws[p] for p in ws if L.distances(adj,*p[:1])[p[1]]<=7)
        assert sum(new.values())==(F(143,2) if kind=='B' else F(71))
        rules.append(dict(kind=kind,target=target,changed_pairs=sorted(CHANGES[kind]),old_cost=str(sum(ws.values())),new_cost=str(sum(new.values())),
            weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(new.items())],automorphisms=autos))
    rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types() if t[1]<=59]
    for a in types:
        for b in types:
            inc=8 if {a,b}=={('A',1),('B',12)} else 0
            assert C.actual_cross(a,b)+inc<=C.upper_cross(rep(a),rep(b))
    assert C.actual_cross(('A',1),('B',12))==0 and C.upper_cross(('A',1),('D',12))==40
    assert docs[C.TAG]['finite']['exact_maximum']=='793/2'
    assert C.value([('A',16),('A',16),('D',27)])==16*F(793,2)
    ideal=[]
    for r in docs[C.TAG]['ideal']['rows']:
        saving=F(r['ranks'].count(12),13);margin=F(r['normalized_margin'])+saving
        ideal.append(dict(ranks=r['ranks'],link_vertices=r['link_vertices'],old_normalized_margin=r['normalized_margin'],
            normalized_saving=saving,normalized_margin=margin,absolute_margin=margin*r['link_vertices']))
    assert len(ideal)==1890
    assert sum(r['normalized_margin']<0 for r in ideal)==3 and all(r['normalized_margin']!=0 for r in ideal)
    paths={Path(__file__),Path(H.__file__),Path(O.__file__),Path(C.__file__),Path(L.__file__),Path(K.__file__),Path(D.__file__),Path(P.__file__),Path(T.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER13_B143_OVER2_C71_INNER_TERMINAL_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        formulas=fs,immediate_affine_roles=immediate(docs),whole_finite_corners=finite_corners(),rules=rules,
        current_literal_interface=literal(docs,targets),new_order13_family_costs=dict(A='65',B='143/2',C='71',D='71'),new_order13_maximum='143/2',
        new_coupon=dict(small='whole_A1',large='whole_B12',large_role='third node from doubled terminal',canonical_large_node=9,
            weight='1/2',eligible_marked_pairs=1,component_pair_increment='1/2',old_actual_total='0',new_actual_total='1/2'),
        finite_comparison=dict(simultaneous_type_pairs=len(types)**2,unchanged_representative_allowance='5/2',new_residual_weights=False,
            inherited_finite_maximum='793/2',unchanged_actual_attainer=['A16','A16','D27'],forest_DP_repeated=False),
        ideal_frontier_update=dict(rows=ideal,complete_saved_domain=1890,new_convolution=False,previous_negative_count=6,
            negative_rows=[r for r in ideal if r['normalized_margin']<0],zero_count=sum(r['normalized_margin']==0 for r in ideal),
            normalized_minimum=min(r['normalized_margin'] for r in ideal),negative_absolute_margins_not_extended=True),
        triangle_minimum=2,quadrilateral_minimum=1,no_selected_census_used=True,new_dimension_bound=None,
        scope='Universal native B2/large-component profiles, arbitrary spectator Cauchy, whole B11/D11 and finite alternatives; independent C-first priority and C triangle exclusion. New actual A1/B12 half absorbed by unchanged representative allowance. No global bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_order13_family_costs'],'literal',len(out['current_literal_interface']['changed']),flush=True)
