#!/usr/bin/env python3
"""Universal affine D13 near-fork halves, with one explicit A1/D13 coupon.

All-native/arbitrary-spectator coverage is proved in the companion report.
No finite forest maximum or ambient dimension theorem is asserted.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order13_d_center_pair_half_transfer_20261018 as H
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as S
import fresh_order14_c_basic_distance9_transfer_20261018 as C

ROOT=H.ROOT;L=H.L;T=S.T;A=S.A;Q=S.Q;D=S.D;P=C.P
TAG='fresh_order14_d_near_fork_half_transfer_20261018'
CURRENT='fresh_depth6_dimension59_20261018'
TAIL='fresh_deep_finite8_dimension64_20261018'
BC='fresh_order11_bc_zero_transfer_20261018'
FIN14='fresh_finite_rank14_existing_payment_transfer_20261018'
BELOW='fresh_spectator_energy_below_one_20261017'
LEAF='fresh_large_finite_d_leaf_half_transfer_20261017'
CHANGED={(1,12),(1,13),(8,10),(8,11)}
STEMS=(CURRENT,TAIL,BC,FIN14,BELOW,LEAF,S.TAG,H.TAG,L.TAG,C.TAG,
       'fresh_order14_b_basic_distance9_transfer_20261018',
       'fresh_finite7_ideal5_triangle_payments_20261017',
       'fresh_depth6_finite7_adaptive_transfer_20261018')


def affine():
    g=L.gram(14,[(i,i+1,1) for i in range(9)]+[(0,10,1),(0,11,1),(9,12,1),(9,13,1)])
    nu=s.Matrix([1]*10+[s.Rational(1,2)]*4)
    assert g*nu==s.zeros(14,1) and g.rank()==13
    return g,nu


def profiles():
    t,u,v,w=s.symbols('t u v w',nonnegative=True)
    K=s.Matrix([[s.Rational(3,4),s.Rational(1,2),s.Rational(1,4)],
                [s.Rational(1,2),1,s.Rational(1,2)],
                [s.Rational(1,4),s.Rational(1,2),s.Rational(3,4)]])
    z=s.Matrix([u,v,w]);energy=(z.T*K*z)[0];height=t+(u+w)/2+v
    near=[dict(name='outer',contacts={1:t,10:u,0:v,11:w},height=height,c0=2-energy,moment=0),
          dict(name='inner3',contacts={1:t,2:1,10:u,0:v,11:w},height=height+1,c0=1-energy,moment=1),
          dict(name='inner4',contacts={1:t,2:s.sqrt(2)},height=t+s.sqrt(2),c0=0,moment=s.sqrt(2)),
          dict(name='next3',contacts={1:t,3:1},height=t+1,c0=0,moment=2)]
    far=[dict(name='isolated',contacts={12:t},height=t/2,c0=2,moment=0),
         dict(name='mate3',contacts={12:t,13:1},height=(t+1)/2,c0=1,moment=1),
         dict(name='mate4',contacts={12:t,13:s.sqrt(2)},height=(t+s.sqrt(2))/2,c0=0,moment=s.sqrt(2)),
         dict(name='center3',contacts={12:t,9:1},height=t/2+1,c0=0,moment=2)]
    return near,far,(t,u,v,w),K


def formulas():
    g,nu=affine();W=[i for i in range(14) if i not in (1,12)];J=g.extract(W,W).inv()
    scale=s.diag(nu[1],nu[12]);V=-g.extract(W,[1,12])*scale
    old=scale*g.extract([1,12],[1,12])*scale-V.T*J*V
    assert old==s.Matrix([[1,-1],[-1,1]])/10
    near,far,variables,K=profiles();assert g.extract([10,0,11],[10,0,11]).inv()==K
    inverse=[]
    for ids in [list(range(2,10))+[12,13],[13]+list(reversed(range(10)))+[10,11]]:
        inv=g.extract(ids,ids).inv();diag=list(inv.diagonal())
        assert diag[:2]==[1,2] and all(x>2 for x in diag[2:])
        inverse.append(dict(nodes=ids,diagonal=list(map(str,diag))))
    rows=[]
    for end,ps,farvec in [('near',near,V[:,1]),('far',far,V[:,0])]:
        for p in ps:
            x=s.Matrix([p['contacts'].get(i,0) for i in W])
            assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-p['height'])==0
            assert s.simplify(10*(farvec.T*J*x)[0]-p['moment'])==0
            assert s.simplify(2-(x.T*J*x)[0]-p['c0']-p['moment']**2/10)==0
            rows.append(dict(end=end,profile=p['name'],height=p['height'],c0=p['c0'],moment=p['moment'],
                cost=s.factor(2*p['moment']/p['height']+p['c0']/p['height']**2)))
    cross=[]
    for p,q in product(near,far):
        x=s.Matrix([p['contacts'].get(i,0) for i in W]);y=s.Matrix([q['contacts'].get(i,0) for i in W])
        value=s.simplify((x.T*J*y)[0]);assert s.simplify(value-p['moment']*q['moment']/10)==0
        cross.append(dict(near=p['name'],far=q['name'],product=value))
    h,k,l,m,c0,d0,p,q,r=s.symbols('h k l m c0 d0 p q r',nonzero=True);a=s.Rational(1,10)
    G=s.Matrix([[a,-a,-h+l*a,-m*a],[-a,a,-l*a,-k+m*a],
        [-h+l*a,-l*a,c0+l*l*a-p,-l*m*a-q],[-m*a,-k+m*a,-l*m*a-q,d0+m*m*a-r]])
    expected=1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/10+(p/h**2-2*q/(h*k)+r/k**2)/10
    assert s.factor(G.det()/(h*h*k*k)-expected)==0
    t,u,v,w=variables;cost=lambda z:s.factor(2*z['moment']/z['height']+z['c0']/z['height']**2)
    assert cost(near[0]).subs({t:0,u:1,v:0,w:0})==5
    assert cost(near[0]).subs({t:0,u:s.sqrt(2),v:0,w:0})==1
    assert cost(near[0]).subs({t:0,u:0,v:1,w:0})==1
    assert cost(near[1]).subs({t:0,u:0,v:0,w:0})==3
    two=[]
    for i,j in combinations(range(3),2):
        lower=K[i,i]+K[j,j]+2*K[i,j];assert lower>=2;two.append(dict(nodes=[i,j],projection_lower=lower))
    assert cost(far[0]).subs(t,1)==8 and cost(far[0]).subs(t,s.sqrt(2))==4
    assert cost(far[1]).subs(t,0)==8 and cost(far[1]).subs(t,1)==3
    assert cost(far[2]).subs(t,0)==cost(far[3]).subs(t,0)==4
    assert s.Rational(1,4)*(5+8-10)==s.Rational(3,4)
    return dict(resistance=10,near_small_inverse=K.tolist(),large_inverse_diagonals=inverse,
        end_profiles=rows,cross_identities=cross,determinant_identity=expected,
        near_multi_contact_strict_edge_exclusions=two,near_extreme_cost=5,near_other_maximum=3,
        far_extreme_cost=8,far_other_maximum=4,far_nonextreme_strict_gap='1/10',
        one_far_extreme_weighted_energy_maximum=1,nonzero_far_support_lower_bound=2,
        both_extreme_unscaled_energy='3/4')


def immediate(data):
    f,_=affine();out=[]
    for name,node,mag,kind,mapping in [
        ('inner4',2,s.sqrt(2),'B',list(range(9,1,-1))+[14,12,13]),
        ('next3',3,1,'D',list(range(9,2,-1))+[12,13,2,14])]:
        g=2*s.eye(15);g[:14,:14]=f;g[14,node]=g[node,14]=-mag
        canonical,_,_=L.affine(kind);Hc=g.extract(mapping,mapping);assert Hc==canonical
        z=Hc.nullspace()[0]
        if z[0]<0:z=-z
        assert all(x>0 for x in z)
        pair=tuple(sorted((mapping.index(12),mapping.index(14))))
        rules=data[BC]['rules'] if kind=='B' else data[L.TAG]['factor_rules']
        rule=next(r for r in rules if r['kind']==kind);weights={tuple(r['pair']):F(r['weight']) for r in rule['weighted_pairs']}
        assert weights[pair]==(0 if kind=='B' else F(1,2))
        out.append(dict(profile=name,nodes=mapping,gram=Hc.tolist(),positive_kernel=list(z),marks=[12,14],
            type='affine_'+kind+'10',canonical_pair=pair,current_weight=str(weights[pair]),
            payment='independent_BC11_zero_theorem' if kind=='B' else 'distinct_half_plus_selected_half'))
    return out


def finite_literals(data):
    f,_=affine();out=[]
    for name,omit,attach,marks,expected in [
        ('near_extreme_empty',1,10,[12,14],('A4','D10')),
        ('far_isolated_empty',12,12,[1,14],('A1','D13')),
        ('far_mate_empty',12,13,[1,14],('D14',))]:
        g=2*s.eye(15);g[:14,:14]=f;g[14,attach]=g[attach,14]=-1
        ids=[i for i in range(15) if i!=omit];assert L.inertia(g.extract(ids,ids))==[len(ids),0,0]
        gg=[[F(x) for x in row] for row in g.tolist()];blocks=D.finite_blocks(gg,ids)
        aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks];roles=[]
        if aa is bb:
            assert expected==('D14',) and aa['type']=='D14'
            images=T.mark_pairs(gg,aa['nodes'],marks,14,'D');assert images==[(0,10)]
            current=next(r for r in data[FIN14]['rules'] if r['kind']=='D')
            assert next(r['weight'] for r in current['weighted_pairs'] if r['pair']==[0,10])=='1'
            pay='1';roles=[dict(type='D14',pair_images=images)]
        else:
            small,large=(aa,bb) if len(aa['nodes'])<len(bb['nodes']) else (bb,aa)
            sm=next(m for m in marks if m in small['nodes']);lm=next(m for m in marks if m in large['nodes'])
            assert (small['type'],large['type'])==expected
            si=T.mark_pairs(gg,small['nodes'],[sm],len(small['nodes']),'A')
            li=T.mark_pairs(gg,large['nodes'],[lm],len(large['nodes']),'D')
            expected_roles=([(0,),(3,)],[(8,),(9,)]) if name=='near_extreme_empty' else ([(0,)],[(9,)])
            assert (si,li)==expected_roles
            pay='1/2';roles=[dict(type=small['type'],images=si),dict(type=large['type'],images=li)]
        out.append(dict(case=name,omitted_mark=omit,marks=marks,components=blocks,roles=roles,
            alternative_payment=pay,selected_half='1/2',total=str(F(pay)+F(1,2))))
    return out


def coupon(g,ids,marks):
    blocks=D.finite_blocks(g,ids);aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if aa is bb:return []
    out=[]
    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='D13':
            images=T.mark_pairs(g,la['nodes'],[lm],13,'D')
            if images==[(9,)]:out.append(dict(small_nodes=sm['nodes'],small_mark=m,large_nodes=la['nodes'],large_mark=lm,large_images=images,increment='1/2'))
    assert len(out)<=1
    return out


def literal(data,target):
    tags=S.B.C.R.B.C.L.STEMS;docs={t:json.loads((ROOT/f'results/{t}.json').read_text()) for t in tags}
    norm=S.B.C.R.B.C.L.N.normalize
    prior={(r['source'],norm(r['case'])):r for r in data[FIN14]['combined_literal_interface']['rows']}
    rows=[];changed=[]
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=prior[stem,norm(index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids));images=[];cp=[]
            if sig==[len(ids),0,0]:cp=coupon(g,ids,marks)
            elif sig[1]==0:
                block=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if block and len(block)==14:
                    labels=T.T.labels(g,block)
                    if labels is not None:images=sorted({tuple(sorted(mp[block.index(m)] for m in marks)) for mp in P.embeddings(labels,target)})
            debit=F(1,2) if set(images)&CHANGED else F();inc=F(len(cp),2)
            if debit:assert set(images)<=CHANGED
            value=F(old['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert value>=0
            corners.append(dict(marks=marks,images=images,debit=debit,coupons=cp,increment=inc,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['debit'] or c['increment'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows),tags


def enc(x):
    if isinstance(x,(F,s.Basic)):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,d in data.items():
        assert d['status'].startswith('PASS'),tag
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    old=next(r for r in data[TAIL]['short7_small_factor_interfaces'] if(r['order'],r['kind'])==(14,'D'))
    b_current=data['fresh_order14_b_basic_distance9_transfer_20261018']
    assert b_current['new_cost']=='80' and data[C.TAG]['new_cost']=='80'
    assert b_current['current_order14_family_costs']==dict(A=70,B=80,C=80,D=82)
    target=old['target'];weights={p:F() for p in combinations(range(14),2)}
    weights.update({tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']})
    assert sum(weights.values())==82 and weights[0,9]==0 and all(weights[p]==1 for p in CHANGED)
    autos=P.embeddings(target,target);assert {tuple(sorted(mp[i] for i in (1,12))) for mp in autos}==CHANGED
    f,_=affine();assert target==[[3 if i!=j and f[i,j] else 0 for j in range(14)] for i in range(14)]
    for p in CHANGED:assert L.distances(L.graph(f),p[0])[p[1]]==9;weights[p]=F(1,2)
    assert sum(weights.values())==80
    assert all(weights[p]==1 for p in weights if L.distances(L.graph(f),p[0])[p[1]]<=8)
    assert data[LEAF]['added_coupon']=='Weight1/2 between any mark of a whole finite component of order<=5 and a short-fork leaf of a separate whole finite D10.'
    assert S.actual_cross(('A',4),('D',10))==32 and S.actual_cross(('A',1),('D',13))==12
    tri=next(r for r in json.loads((ROOT/f'results/{C.X.COMBINED}.json').read_text())['triangle_pair_retention'] if(r['order'],r['kind'])==(14,'D'))
    assert all(weights[tuple(p)]==1 for p in tri['required_pairs'])
    lit,tags=literal(data,target)
    paths={Path(__file__),Path(H.__file__),Path(L.__file__),Path(S.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    paths.update(ROOT/f'results/{tag}.json' for tag in tags);paths.add(ROOT/f'results/{C.X.COMBINED}.json')
    return enc(dict(status='PASS_ORDER14_D_NEAR_FORK_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,changed_pairs=sorted(CHANGED),old_weight='1',new_weight='1/2',old_cost='82',new_cost='80',
        current_order14_family_costs=dict(A=70,B=80,C=80,D=80),current_order14_maximum='80',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        formulas=formulas(),immediate_corners=immediate(data),finite_whole_component_corners=finite_literals(data),
        new_coupon=dict(small='whole_A1',large='whole_D13',large_role='unique_long_arm_neighbor_of_fork_center',
            canonical_large_node=9,weight='1/2',eligible_marked_pairs=1,component_pair_increment='1/2',
            old_actual_component_pair_cost='3/2',new_actual_component_pair_cost='2'),
        current_literal_interface=lit,triangle_required_pairs=tri['required_pairs'],triangle_minimum=2,quadrilateral_minimum=1,
        finite_budget_recomputed=False,new_dimension_bound=None,
        scope='Universal all-native near-fork/far-leaf profile proof with arbitrary spectators and inherited below-one classification. Four ideal halves save2; one explicit A1/D13-neighbor half costs1/2. Current216 literal interfaces reconstructed. Finite/global counting separate.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_cost'],'changed_literals',len(result['current_literal_interface']['changed']),flush=True)
