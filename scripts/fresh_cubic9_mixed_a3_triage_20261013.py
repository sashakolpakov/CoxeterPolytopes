#!/usr/bin/env python3
"""Archive exact field-space equivalence triage for the rejected A3 pivot.

Requires PARI/GP for local Hilbert symbols. This does not identify the mixed
lattice or its full integral reflection group with an earlier B^4+1 lattice.
"""
import ast
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import shutil
import subprocess
import fresh_cubic9_dihedral_independent_audit_20261013 as f
ROOT=Path(__file__).resolve().parents[1]

def run():
    gp=shutil.which('gp') or '/opt/homebrew/bin/gp'
    command=[gp,'-q',str(ROOT/'scripts/fresh_cubic9_mixed_a3_triage_20261013.gp')]
    proc=subprocess.run(command,capture_output=True,text=True,check=True)
    rows=[ast.literal_eval(line) for line in proc.stdout.splitlines() if line.startswith('[')]
    assert rows==[[0,2,1,1,1,1],[0,3,1,1,1,1],[1,2,1,-1,1,-1],[1,3,1,-1,1,-1]]
    A3=tuple(tuple(f.scale(f.ONE,2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(3)) for i in range(3))
    P=tuple(tuple(f.scale(f.ONE,x) for x in row) for row in ((0,1,0),(1,-1,1),(-1,0,1)))
    assert tuple(tuple(f.sum_field(f.mul(P[k][i],P[k][j]) for k in range(3)) for j in range(3)) for i in range(3))==A3
    assert f.determinant(P) in (f.TWO,f.neg(f.TWO))
    records=[]
    for k in (0,1):
        unit=f.mul(f.THETA,f.BETA) if k else f.THETA
        residues=[]
        for raw in product(range(2),repeat=3):
            a=tuple(map(F,raw));q=f.scale(f.mul(unit,f.mul(a,a)),F(-1,2))
            if f.integral(q):residues.append(list(raw))
        assert residues==[[0,0,0]]
        # Overclosed M determinant / matched old determinant, before scaling.
        ratio=f.ONE if k else f.inverse(f.mul(f.BETA,f.BETA))
        records.append({'new_alpha':f.encoded(f.scale(unit,2)),
            'old_alpha':f.encoded(f.scale(f.THETA,2) if k else f.scale(f.mul(f.THETA,f.BETA),2)),
            'old_form_scale':f.encoded(f.BETA if k else f.ONE),
            'unscaled_overclosed_determinant_ratio':f.encoded(ratio),
            'time_residues_allowing_integral_order_two_overclosure':residues,
            'field_space_verdict':'SIMILAR_BY_HASSE_MINKOWSKI' if k else 'ISOMETRIC_BY_HASSE_MINKOWSKI',
            'integral_identification_with_old_B4_host':'UNRESOLVED'})
    out={'status':'PASS','field':'Q(T), T^3-3T-1=0',
        'local_hasse_columns':['case','rational_prime','prime_index','mixed_overclosed','old','scaled_old'],
        'local_hasse_records':rows,'a3_embedding_in_I3':[f.vecencode(r) for r in P],
        'a3_overclosure_index_ideal':'(2)','a3_overclosure_index_cardinality':8,
        'cases':records,
        'argument':'At all finite primes outside 2 and 3 every diagonal entry is a unit, so local Hasse symbols are +1. The determinant ratios are squares, and the real signatures agree; beta is totally positive. The displayed local symbols therefore establish the stated field isometry/similarity. The unique integral index-(2) A3 overclosure gives O(<-alpha>+B^3+A3)=O(<-alpha>+B^3+I3), since the smaller even lattice is characterized by q in2O. This does not identify either group with the older B^4+1 full integral group.',
        'scope':'Triage of two proposed hosts only; no integral old-host duplication or nonreflectivity claim.'}
    path=ROOT/'results/fresh_cubic9_mixed_a3_triage_20261013.json';path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','cases':2,'old_integral_group_identification':'UNRESOLVED'}))
    return out
if __name__=='__main__':run()
