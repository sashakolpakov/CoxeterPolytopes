#!/usr/bin/env python3
"""Exact transport to the earlier saturated silver K3 carrier.

The integral lattice transport has index sqrt(2), not index one.
Equality of full integral orthogonal groups follows separately from the
intrinsic mod2 descriptions and unique integral overclosures in the report.
"""
from itertools import permutations
import json
from pathlib import Path

from scripts.fresh_silver_f4_h9_20261013 import S,Z,O,SQRT2,UNIT,F4,Host,pair,mv,inverse,enc,ray_key
from coxeter_polytopes.exact.silver import silver_determinant
from coxeter_polytopes.geometric.silver_k3_carrier import audit_silver_c5_saturation
from coxeter_polytopes.search.silver_k3_vinberg import silver_k3_alexeev_field_gram

ROOT=Path(__file__).resolve().parents[1]
def matmul(a,b): return tuple(tuple(sum((x*y for x,y in zip(row,col)),Z) for col in zip(*b)) for row in a)
def transpose(a):return tuple(zip(*a))
def decode(rows):return tuple(tuple(S.from_json(x) for x in row) for row in rows)

def main():
    calibration=json.loads((ROOT/'results/fresh_silver_f4_h4_calibration_20261013.json').read_text())
    B=decode(calibration['gram']); W=decode(calibration['all_walls'])
    G=tuple(tuple(pair(r,B,s) for s in W) for r in W)
    old=tuple(tuple(S.coerce(x) for x in row) for row in audit_silver_c5_saturation().field_gram)
    permutation=next(p for p in permutations(range(5)) if all(G[p[i]][p[j]]==old[i][j] for i in range(5) for j in range(5)))
    A=tuple(tuple(W[permutation[j]][i] for j in range(5)) for i in range(5))
    saturation=tuple(tuple(O if i==j and j<4 else SQRT2/2 if j==4 and i in (1,2,4)
                          else Z for j in range(5)) for i in range(5))
    T5=matmul(A,saturation)
    assert all(x.is_integral for row in T5 for x in row)
    assert all(x.is_integral for row in inverse(T5) for x in row)
    carrier=silver_k3_alexeev_field_gram(); new=Host(UNIT,2,True).gram
    T=tuple(tuple(T5[i][j] if i<5 and j<5 else SQRT2 if i==j==9 else O if i==j and 5<=i<9 else Z
                  for j in range(10)) for i in range(10))
    assert matmul(matmul(transpose(T),new),T)==carrier
    assert silver_determinant(T)**2==2
    assert all(x.is_integral for row in T for x in row)
    assert all((SQRT2*x).is_integral for row in inverse(T) for x in row)
    # The two new height30 searches are literally identical oriented walls
    # after the time-coordinate inclusion and positive parallel rescaling.
    unit=json.loads((ROOT/'results/fresh_silver_f4_h9_unit_h30_20261013.json').read_text())
    doubled=json.loads((ROOT/'results/fresh_silver_f4_h9_2unit_h30_20261013.json').read_text())
    unit_walls=decode(unit['all_walls']); doubled_walls=decode(doubled['all_walls'])
    transported=[(SQRT2*r[0],*r[1:]) for r in doubled_walls]
    assert {ray_key(r) for r in transported}=={ray_key(r) for r in unit_walls}
    # Explicit field isometry to the already scanned diagonal quadratic
    # space. It is not an integral lattice isometry: A4 has denominatorsqrt2.
    A4=((O,Z,Z,-O/SQRT2),(-O,O,Z,-O/SQRT2),(Z,-O,SQRT2,-O/SQRT2),(Z,Z,Z,-O/SQRT2))
    assert matmul(transpose(A4),A4)==F4
    assert silver_determinant(A4)**2==1
    out={'status':'EXACT_CARRIER_TRANSPORT_AND_DUPLICATE_HEIGHT30_WALLS',
         'old_carrier_source':'silver_k3_alexeev_field_gram',
         'calibration_to_unsaturated_C5_permutation':list(permutation),
         'saturated_C5_to_time_plus_F4':[enc(row) for row in T5],
         'old_carrier_to_F4_squared_scalar1':[enc(row) for row in T],
         'transport_determinant':silver_determinant(T).to_json(),
         'old_carrier_gram':[enc(row) for row in carrier],
         'target_gram':[enc(row) for row in new],
         'F4_to_diagonal_field_isometry':[enc(row) for row in A4],
         'diagonal_space_note':'For every alpha, this transports F4²+1 to I9 over Q(sqrt2). Denominator sqrt2 prevents an integral identification; full group equivalence with the diagonal lattice remains unproved.',
         'equal_height30_oriented_walls':len(unit_walls),
         'scope':'Index2 as an O-module cardinality, not integral isometry. Full-group equality uses the separate mod2 and unique-overclosure proof.'}
    path=ROOT/'results/fresh_silver_f4_identifications_20261013.json'
    path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'],len(unit_walls),'walls; determinant',silver_determinant(T))

if __name__=='__main__':main()
