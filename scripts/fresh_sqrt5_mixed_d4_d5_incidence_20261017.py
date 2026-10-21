#!/usr/bin/env python3
"""Exact set accounting for mixed D4/D5 special-vertex incidence bounds."""
import argparse,hashlib,json,math
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_mixed_d4_d5_incidence_20261017'
E='fresh_sqrt5_mixed_d4_d5_migration_20261017'
P='fresh_sqrt5_mixed_d4_d5_quad_patterns_20261017'
G='fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017'
D='fresh_sqrt5_d4_short_face_incidence_20261017'
R='fresh_sqrt5_d4_short_face_incidence_review_20261017'
H='fresh_sqrt5_discount_structural_floor_20261016'

def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def run():
    e,p,g,d,r,h=map(read,[E,P,G,D,R,H]);targets=[]
    assert d['leaf_pair_arithmetic_obstruction']['possible_negative_norms']==[-1,-4,-16,-64,-256]
    for endpoint,pattern,gate in zip(e['targets'],p['targets'],g['targets']):
        name=endpoint['name'];assert name==pattern['name']==gate['name']
        n=endpoint['dimension'];tri={tuple(x)for x in endpoint['possible_two_special_triangle_corners']}
        d4_leaves={}
        for comp in endpoint['components']:
            if comp['type']!='D4':continue
            center,*leaves=comp['nodes']
            for leaf in leaves:
                attach=endpoint['ports'][leaf]['attachments']
                assert len(attach)==1 and attach[0]['contacts']==[center] and attach[0]['copy']
            for pair in combinations(leaves,2):d4_leaves[pair]=center
        retained=set(gate['retained_indices']);quad=set();excluded=[];retained_maps=[]
        for index,m in enumerate(pattern['maps']):
            if m['target'] not in retained:continue
            pair=tuple(m['ports'])
            if pair in d4_leaves:
                center=d4_leaves[pair]
                assert m['first_copy'] and m['second_copy'] and m['new_new_label']==2
                assert m['first_contacts']==m['second_contacts']==[center]
                excluded.append(dict(map_index=index,ports=list(pair),center=center,
                    reason='Literal copied D4 leaf pair: apply the independently reviewed AB=16 norm/congruence obstruction.'))
            else:quad.add(pair);retained_maps.append(index)
        assert tri<=quad
        # A root pair determines one incident two-face. At a special vertex,
        # it cannot simultaneously be a two-special triangle and a quad.
        corner_charges=[]
        for pair in combinations(range(n),2):
            triangle=F(1,2)if pair in tri else F(0)
            quadrilateral=F(1,4)if pair in quad else F(0)
            corner_charges.append(dict(ports=list(pair),two_special_triangle=str(triangle),
                all_special_quadrilateral=str(quadrilateral),maximum=str(max(triangle,quadrilateral))))
        alpha=sum((F(x['maximum'])for x in corner_charges),F(0))
        assert alpha==F(len(quad)+len(tri),4)
        # Check the face payments themselves, including overpayment at faces
        # with several nonspecial corners.
        triangle_payments=[]
        for special in range(3):
            pay=F(3-special)+(F(1)if special==2 else F(0));assert pay>=2
            triangle_payments.append(dict(special_vertices=special,payment=str(pay)))
        quad_payments=[]
        for special in range(5):
            pay=F(1)if special==4 else F(4-special);assert pay>=1
            quad_payments.append(dict(special_vertices=special,payment=str(pay)))
        corners=math.comb(n,2);threshold=F(5*((n-1)**2//4),2)-corners
        left=threshold-alpha;right=corners-threshold;assert left>0 and right>0
        scale=math.lcm(left.denominator,right.denominator)
        L=int(left*scale);Rr=int(right*scale);div=math.gcd(L,Rr);L//=div;Rr//=div
        targets.append(dict(name=name,dimension=n,forest=endpoint['forest'],triangle_corner_count=len(tri),
            quad_corner_count_before_D4_gate=gate['retained_corner_count'],quad_corner_count=len(quad),
            triangle_corners=[list(x)for x in sorted(tri)],quad_corners=[list(x)for x in sorted(quad)],
            D4_leaf_maps_excluded=excluded,retained_quad_map_indices=retained_maps,corner_charges=corner_charges,
            triangle_payments=triangle_payments,quadrilateral_payments=quad_payments,
            special_vertex_charge=str(alpha),other_vertex_charge=str(corners),strict_h_threshold=str(threshold),
            frequency_coefficients=[L,Rr],frequency_inequality=f'{L}*V_special < {Rr}*V_other'))
    paths=[f'scripts/{TAG}.py']+[f'results/{t}.json'for t in[E,P,G,D,R,H]]
    return dict(status='PASS_EXACT_MIXED_D4_D5_SPECIAL_VERTEX_INCIDENCE',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},targets=targets,
        hypotheses=['P is a compact hyperbolic Coxeter polytope arithmetic of simplest type over Q(sqrt5).',
            'For each separately stated inequality, special means its ENTIRE vertex forest equals the displayed forest.',
            'The frozen complete endpoint/triangle and quadrilateral-pattern classifications, rational factor lemma, and independently reviewed D4 leaf-pair obstruction are used.'],
        proof='No triangle has three special vertices: its whole normalized Gram would be rational and physically Lorentzian, contradicting conjugate positivity. Pay one unit at every other-vertex corner of each triangle and quadrilateral. A triangle with two special vertices additionally receives one-half at each special corner; an all-special quadrilateral receives one-quarter at each corner. This pays at least2 per triangle and1 per quadrilateral. At each special vertex, the literal root-pair sets bound its total by alpha; every other vertex contributes at most binom(d,2). The strict h-vector inequality T(d)*V<2*T+Q then gives (T(d)-alpha)*V_special<(binom(d,2)-T(d))*V_other.',
        scope='Four separate restrictions on special-vertex frequencies, including the already known3D4 calibration. They exclude polytopes all of whose vertices have the indicated forest. They do not bound dimension, exclude the occurrence of these forests, or establish a simultaneous inequality for their union.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],[(x['name'],x['frequency_inequality'])for x in data['targets']],flush=True)
