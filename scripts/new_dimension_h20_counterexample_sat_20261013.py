#!/usr/bin/env python3
"""Minimal exact obstruction to all subsets of the 580-wall H20 dictionary.

This is a replay of two independently discovered deletion escapes, not an
LP search. Integer arithmetic proves each mandatory-wall clause. Boolean
enumeration proves the three-clause contradiction; Z3 also exports its proof.
No imported project verifier or floating arithmetic is used.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from fractions import Fraction
import itertools
import json
from pathlib import Path
import z3

ROOT=Path(__file__).resolve().parents[1]
BOUNDARY=ROOT/'data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json'
GEOMETRY=ROOT/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUTPUT=ROOT/'data/enumerations/new_dimension_h20_counterexample_sat_20261013.json'
CENTER=[61,-136,-92,-182,-270,-220,-168,-114,-58,230,152,306,462,390,320,252,186,-58,-114,-84,-46,122]
ESCAPES={
0:[61,-136,-92,-182,-270,-220,-168,-114,-58,230,152,306,462,390,320,252,186,-58,-114,-84,72499954,122],
74:[998269353,-2453846426,-1661923260,-3323846518,-4907692848,-3970769670,-3028269566,-2085769458,-1087500116,3881538922,2587692612,5175385228,7763077848,6525000780,5392885256,4260769734,3128654218,-909038578,-1773461766,-1277115552,-713846246,1996538706]}

def geometry():
    raw=json.loads(GEOMETRY.read_text())['coxeter_input']
    Q=[[-int(v) for v in row] for row in raw['basis_gram']]
    f=list(map(int,raw['wall_coordinates_as_columns_transposed'][0]))
    boundary=json.loads(BOUNDARY.read_text())['six_tile']
    R=boundary['primitive_outward_normals']
    assert len(R)==580
    def pairing(a,b):
        return sum(a[i]*sum(Q[i][j]*b[j] for j in range(22)) for i in range(22))
    return f,R,boundary,pairing

def verify(data):
    f,R,boundary,pairing=geometry()
    c=data['common_strict_timelike_center']
    assert len(c)==22 and all(type(x) is int for x in c)
    assert pairing(c,f)==0 and pairing(c,c)<0
    assert max(pairing(r,c) for r in R)<0
    ids=[]
    for row in data['mandatory_wall_witnesses']:
        index=row['index_zero_based'];v=row['vector'];ids.append(index)
        assert 0<=index<len(R) and len(v)==22 and all(type(x) is int for x in v)
        assert pairing(v,f)==0 and pairing(v,v)>0
        values=[pairing(r,v) for r in R]
        assert values[index]>0 and max(values[:index]+values[index+1:])<0
        # Thus the necessary escape-hitting clause contains exactly this wall.
        assert [i for i,p in enumerate(values) if p>0]==[index]
        assert row['square']==pairing(v,v)
        assert row['maximum_other_wall_pairing']==max(values[:index]+values[index+1:])
        assert row['own_wall_pairing']==values[index]
    assert len(ids)==2 and len(set(ids))==2
    i,j=sorted(ids)
    G=[[pairing(R[a],R[b]) for b in (i,j)] for a in (i,j)]
    assert G==data['conflicting_pair_gram']
    value=Fraction(G[0][1]**2,G[0][0]*G[1][1])
    assert G[0][1]<0 and 0<value<1
    assert value not in {Fraction(0),Fraction(1,4),Fraction(1,2),Fraction(3,4)}
    # The original genuine ridge survives every deletion retaining these walls.
    actual=[r for r in boundary['all_actual_ridges'] if r[:2]==[i,j]]
    assert len(actual)==1 and actual[0][2:] == data['actual_ridge'][2:]
    assert actual[0][2]=='acute' and Fraction(actual[0][3])==value
    # Independent exhaustive truth-table replay of x_i, x_j, not(x_i and x_j).
    assert not any(a and b and not(a and b) for a,b in itertools.product((False,True),repeat=2))
    return {'walls':len(R),'mandatory_pair':[i,j],'cosine_squared':str(value),
            'exact_escape_checks':'PASS','actual_ridge_check':'PASS',
            'Boolean_truth_table':'UNSAT'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:
        data=json.loads(OUTPUT.read_text())
    else:
        f,R,boundary,pairing=geometry()
        rows=[]
        for index,v in ESCAPES.items():
            values=[pairing(r,v) for r in R]
            rows.append({'index_zero_based':index,'vector':v,'square':pairing(v,v),
                         'maximum_other_wall_pairing':max(values[:index]+values[index+1:]),
                         'own_wall_pairing':values[index]})
        data={'status':'ALL_580_WALL_SUBSETS_EXCLUDED',
              'source_geometry':str(GEOMETRY.relative_to(ROOT)),
              'source_boundary':str(BOUNDARY.relative_to(ROOT)),
              'witness_provenance':'Independent all-single-facet-deletion escape search by construction_audit; this file independently replays the minimal pair.',
              'common_strict_timelike_center':CENTER,
              'mandatory_wall_witnesses':rows,
              'conflicting_pair_gram':[[pairing(R[a],R[b]) for b in (0,74)] for a in (0,74)],
              'actual_ridge':next(r for r in boundary['all_actual_ridges'] if r[:2]==[0,74]),
              'clauses':['x_0','x_74','not x_0 or not x_74'],
              'scope':'No subset of these580 outward halfspaces is both finite-volume and Coxeter. Adding genuinely new planes or different gluing constructions is not excluded.'}
    result=verify(data)
    bad=deepcopy(data);bad['mandatory_wall_witnesses'][0]['vector']=data['common_strict_timelike_center']
    try:verify(bad)
    except AssertionError:result['forged_spacelike_certificate']='REJECTED'
    else:raise AssertionError('Forged timelike witness was accepted')
    z3.set_param(proof=True)
    solver=z3.Solver();a,b=z3.Bools('x_0 x_74')
    solver.add(a,b,z3.Or(z3.Not(a),z3.Not(b)))
    assert solver.check()==z3.unsat
    result['Z3_status']='UNSAT'
    if not args.verify:
        data['verification']=result
        data['z3_smt2']=solver.to_smt2();data['z3_proof']=str(solver.proof())
        OUTPUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':main()
