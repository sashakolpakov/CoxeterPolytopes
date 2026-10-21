#!/usr/bin/env python3
"""Repair two universal ideal theorems after finite D12 long-leaf halving.

No selected-core census, finite forest optimization or cusp convolution.
The complete new argument is the accompanying report, including arbitrary
spectator support. Literal checks corroborate its forced whole components.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s
import fresh_affine_d11_leaf_pair_transfer_20261017 as L
import fresh_affine_d13_long_pair_transfer_20261017 as H
import fresh_order11_bcd_half_pair_transfer_20261018 as K

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_d12_universal_ideal_output_repairs_20261018'
STEMS=(L.TAG,H.TAG,K.TAG,'fresh_d12_radical_saturated_half_replacements_20261018',
 'fresh_order17_e6_unit_coupon_zero_transfer_20261018',
 'fresh_order15_bd_terminal_coupon_transfer_20261018',
 'fresh_spectator_energy_below_one_20261017',
 'fresh_spectator_five_quarters_transfer_20261017',
 'fresh_order10_low8_ambient58_frontier_20261018')

def enc(x):
    if isinstance(x,(F,s.Rational)):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v)for v in x]
    return x

def gram(n,edges):return K.gram(n,[(i,j,q)for i,j,q in edges])
def put(g,i,j,q=1):g[i,j]=g[j,i]=-q
def inertia(g):
    # Physical norm-two graphs are retained for role checks. A diagonal
    # congruence only makes doubled bonds rational for the inertia helper.
    parity={}
    for start in range(g.rows):
        if start in parity:continue
        parity[start]=0;queue=[start]
        for i in queue:
            for j in range(g.rows):
                if i==j or not g[i,j]:continue
                assert g[i,j]in(-1,-s.sqrt(2))
                p=parity[i]^int(g[i,j]==-s.sqrt(2))
                if j in parity:assert parity[j]==p
                else:parity[j]=p;queue.append(j)
    D=s.diag(*[s.sqrt(2)**parity[i]for i in range(g.rows)])
    q=(D*g*D).applyfunc(s.simplify);assert all(x.is_Rational for x in q)
    return K.inertia(q)

def parts(g,ids,marks):
    out=[]
    for ns in K.components(g,ids):
        adj={i:{j for j in ns if j!=i and g[i,j]!=0}for i in ns}
        leaves=[i for i in ns if len(adj[i])<=1]
        branch=[i for i in ns if len(adj[i])==3]
        assert all(g[i,j]in(0,-1)for i in ns for j in ns if i!=j)
        typ=None;short=[];long=[]
        if not branch:typ='A'+str(len(ns))
        elif len(branch)==1:
            center=branch[0];dist=K.distances(adj,center)
            arms=sorted(dist[i]for i in leaves)
            if arms[:2]==[1,1]:typ='D'+str(len(ns));short=[i for i in leaves if dist[i]==1]
            elif arms==[1,2,2]:typ='E6';long=[i for i in leaves if dist[i]==2]
        assert typ is not None,(ns,branch)
        out.append(dict(type=typ,nodes=ns,marks=[m for m in marks if m in ns],leaves=leaves,
            short_leaves=short,long_E6_ends=long,adj=adj))
    return out

def corner(g,W,marks):
    ids=sorted(set(W)|set(marks));K.pd(g.extract(ids,ids))
    ps=parts(g,ids,marks);a,b=[next(p for p in ps if m in p['nodes'])for m in marks]
    value=F();reason='unused'
    if a is b:
        d=K.distances(a['adj'],marks[0])[marks[1]]
        assert len(a['nodes'])<=13
        if a['type']=='D12'and d==10 and any(m in a['short_leaves']for m in marks):
            value=F(1,2);reason='new_D12_long_leaf_half'
        else:value=F(1);reason='unchanged_whole_small_unit'
    else:
        for small,large in((a,b),(b,a)):
            if small['type']in('A1','A2')and large['type']=='A11'and large['marks'][0]in large['leaves']:
                value=F(1,2);reason='new_A1_A2_A11_endpoint_half'
            if large['type']=='D11'and large['marks'][0]in large['short_leaves']:
                if small['type']=='E6'and small['marks'][0]in small['long_E6_ends']:
                    value=F(1);reason='retained_E6_length2_end_D11_leaf_unit'
                elif small['type']=='A5':
                    dd=K.distances(small['adj'],small['marks'][0])
                    assert sorted(dd[i]for i in small['leaves'])==[1,3]
                    value=F(1);reason='retained_A5_second_D11_leaf_unit'
    for p in ps:p.pop('adj')
    return dict(marks=marks,components=ps,payment=value,reason=reason)

def check_star(g,W,a,b,c,d):
    assert inertia(g)==[g.rows-2,1,1]
    K.pd(g.extract(W,W))
    for r in(a,b,c,d):K.pd(g.extract(W+[r],W+[r]))
    for pair in((a,b),(b,c),(c,d),(d,a)):
        sig=inertia(g.extract(W+list(pair),W+list(pair)))
        assert sig[1]==0

def one_sided_spectrum():
    # All-rank A/B/D reasoning is in the report; these exact endpoints
    # and all exceptional diagrams complete its bounded remainder.
    phi=(1+s.sqrt(5))/2
    assert 2<s.Rational(5,2)<phi**2
    assert s.Rational(16,3).q!=1
    e={}
    for name,arms in(('E6',(1,2,2)),('E7',(1,2,3)),('E8',(1,2,4))):
        edges=[];n=1
        for length in arms:
            prev=0
            for _ in range(length):edges.append((prev,n,1));prev=n;n+=1
        e[name]=gram(n,edges)
    e['F4']=gram(4,[(0,1,1),(1,2,s.sqrt(2)),(2,3,1)])
    e['H3']=gram(3,[(0,1,phi),(1,2,1)])
    e['H4']=gram(4,[(0,1,phi),(1,2,1),(2,3,1)])
    rows=[]
    for name,g in e.items():
        J=g.inv();ds=[s.simplify(J[i,i])for i in range(g.rows)]
        assert all(s.simplify(x-s.Rational(5,4)).is_positive for x in ds)
        rows.append(dict(type=name,diagonals=list(map(str,ds))))
    assert s.simplify(2/(4-phi**2)-s.Rational(5,4)).is_positive
    assert F(1,2)+F(3,4)==F(5,4)
    return dict(energy='5/4',one_component='whole D5, one simple short-leaf contact',
        two_components='whole A1 and A3, simple endpoint contacts on the same side',
        exceptional_inverse_diagonals=rows,I2_4_diagonal='1',I2_5_diagonal=str(s.simplify(2/(4-phi**2))),
        A_second_node_equation_forces_noninteger_N='16/3',A_deeper_node_minimum='3/2',
        nonunit_A1_required_squared_magnitude='5/2',native_gap_squared=['2',str(phi**2)])

def d11_repairs():
    f,nu=L.affine();rows=[];kernel=[]
    # Nonextreme end4/next3 already contains an affine principal block.
    for b in(10,11):
        ids=[i for i in range(12)if i!=b]
        for name in('mate4','center3'):
            p=next(p for p in L.profiles(b)if p['name']==name)
            g=2*s.eye(12);g[:11,:11]=f.extract(ids,ids)
            for i,q in p['contacts'].items():
                if i!=b:put(g,ids.index(i),11,q)
            vs=g.nullspace();assert len(vs)==1
            v=vs[0]if vs[0][-1].is_positive else-vs[0]
            assert all(x.is_positive for x in v)
            kernel.append(dict(mark=b,profile=name,kernel=list(map(str,v))))
    phi=(1+s.sqrt(5))/2
    phienergy=s.simplify(2-3*phi**2/4)
    assert phienergy.is_positive and s.simplify(s.Rational(1,2)-phienergy).is_positive
    assert F(5,4)>F(12,11)
    for(a,b),side,kind in product(L.MARKS,('left','right'),('isolated4','mate3_1','D5','A1_A3')):
        c,d=12,13;quiet,other=(c,d)if side=='left'else(d,c)
        near,far=(a,b)if side=='left'else(b,a)
        mate=lambda x:17-x if x<10 else 21-x
        extra={'isolated4':1,'mate3_1':0,'D5':5,'A1_A3':4}[kind]
        g=2*s.eye(14+extra);g[:12,:12]=f;put(g,quiet,mate(near))
        if kind=='isolated4':put(g,other,far,s.sqrt(2));put(g,other,14)
        elif kind=='mate3_1':put(g,other,far);put(g,other,mate(far))
        elif kind=='D5':
            put(g,other,far)
            for i,j in((14,15),(15,16),(16,17),(16,18),(other,17)):put(g,i,j)
        else:
            put(g,other,far)
            for i,j in((other,14),(other,15),(15,16),(16,17)):put(g,i,j)
        W=sorted(set(range(g.rows))-{a,b,c,d})
        check_star(g,W,a,b,c,d)
        first=corner(g,W,(far,quiet));assert first['reason']=='new_D12_long_leaf_half'
        second=corner(g,W,(quiet,other)if kind in('isolated4','mate3_1')else(near,other))
        assert first['payment']+second['payment']>=1
        rows.append(dict(old_marks=[a,b],spectator_free_extreme_side=side,branch=kind,
            gram=[list(map(str,row))for row in g.tolist()],white=W,corners=[first,second],
            total=first['payment']+second['payment']))
    assert len(rows)==32
    return dict(rows=rows,payment_classes=dict(Counter(r['corners'][1]['reason']for r in rows)),
        nonextreme_affine_kernel_exclusions=kernel,isolated_phi_one_sided_energy=str(phienergy),
        both_extreme_mate_edge_residual='12/11',both_extreme_required_spectator_energy='5/4')

def d13_repairs(old):
    f,nu=H.affine();rows=[]
    for r in old['nonextreme_zero_cases']:
        a,b=r['marks'];side=r['center_inner3_side'];leaf=r['opposite_leaf_attachment']
        g=2*s.eye(16);g[:14,:14]=f
        for i,j in(((14,1),(15,leaf))if side=='left'else((14,leaf),(15,8))):put(g,i,j)
        W=sorted(set(range(14))-{a,b});check_star(g,W,a,b,14,15)
        first=corner(g,W,tuple(r['chosen_unit_corner']))
        if first['payment']==1:second=None;total=F(1)
        else:second=corner(g,W,(14,15));total=first['payment']+second['payment'];assert total>=1
        rows.append(dict(old_marks=[a,b],side=side,attachment=leaf,corners=[first,second],total=total))
    assert len(rows)==12
    assert Counter(r['corners'][0]['payment']for r in rows)==Counter({F(1):4,F(1,2):8})
    assert Counter(r['corners'][1]['reason']for r in rows if r['corners'][1])=={
        'new_A1_A2_A11_endpoint_half':4,'unchanged_whole_small_unit':4}
    for r in old['extreme_equality']['rows']:
        for c in r['corners'][1:]:
            assert not any(p['type']=='D12'and len(p['marks'])==2 for p in c['components'])
    return dict(rows=rows,unchanged_sealed_extreme_cases=72,unchanged_finite_extreme_corners=216)

def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    e=docs['fresh_order17_e6_unit_coupon_zero_transfer_20261018']['finite_coupon']
    assert e['new_weight']=='1'and e['large_mark']=='short fork leaf'
    a=docs['fresh_order15_bd_terminal_coupon_transfer_20261018']['new_cross_component_coupons']
    assert any(x['small_type']=='whole A5'and x['small_role']=='second_from_endpoint'and x['new_marked_pair_weight']=='1'for x in a)
    for r in docs[L.TAG]['two_A1_A3_spectators']['rows']:
        for c in r['corners']:
            assert not any(p['type']=='D12'and len(p['marks'])==2 for p in c['components'])
    paths=[Path(__file__),Path(K.__file__)]+[ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_D12_HALVES_UNIVERSAL_D11_D13_IDEAL_OUTPUT_REPAIRS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        one_sided_five_quarters=one_sided_spectrum(),affine_D11=d11_repairs(),affine_D13=d13_repairs(docs[H.TAG]),
        unaffected_D11_two_sided_cases=32,new_coupons_already_in_proposed_D12_rule=True,
        additional_coupon_beyond_proposal='0',ideal_weights_changed=False,
        scope='Two universal ideal-output repairs only; no claim of complete nonselected D12 coverage, forest budget or new dimension bound.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=enc(run());p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['affine_D11']['rows']),len(out['affine_D13']['rows']),flush=True)
