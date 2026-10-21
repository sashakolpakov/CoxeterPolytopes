# Joint A4/scalar completion and three targeted relative-norm orbits

No finite-volume example was obtained. The prime-five cubic-nine H9 host now has a **203-wall exact Coxeter cone**, including six explicit primitive norm-five mirrors found by targeted relative norm equations. Its 197-wall H5 face and the full H9 cone retain an exact spacelike escape. This is constructive progress in the existing changed host, not a new polytope or a proof that the full host is nonreflective.

The final artifact is [fresh_h9_joint_arithmetic_20261013_pell3.json](fresh_h9_joint_arithmetic_20261013_pell3.json). Its [independent Fraction audit](fresh_h9_joint_arithmetic_20261013_pell3_independent.json) checks all **20,503 H9 wall pairs** and **19,306 H5 wall pairs**, lattice integrality, every finite Coxeter angle, strict timelike interiors, and strict spacelike perturbations of the saved escape.

## Host and scope

Put `T=2cos(pi/9)`, `T³−3T−1=0`, `beta=4−T²`, and `B=[[2,−T],[-T,2]]`. The host remains

`L=<-T> + B² + A4 + <2>`.

The lattice, full zero-height census, norm-class proof, and determinant comparisons were established in [the preceding positive arithmetic report](fresh_cubic9_positive_20261013.md). The present round does not relabel it as a new host. In particular, the earlier global nonreflectivity theorem for a conductor-nine `B⁴+<1>` host does not apply to this A4 host.

The old 22-wall escape lies in the time+A4+A1 complement, with both B coordinates zero. This complement has full primitive norm pool `{1,2,5,10}`. The full H9 pool additionally allows these classes times beta; beta classes cannot occur in a root supported only on the complement because reflection integrality at beta would make every coordinate divisible by beta, contradicting its odd norm valuation. The dyadic imprimitive-vector argument excludes classes divisible by 4. Thus freeing all four A4 covectors and the scalar is a complete search of this complement at the stated height, not a smaller hand-picked root-norm search.

## Completing the missing directions at the same bound

[Producer](../scripts/fresh_h9_joint_arithmetic_20261013.py), [artifact](fresh_h9_joint_arithmetic_20261013.json), [independent census](../scripts/fresh_h9_joint_arithmetic_20261013_independent.py), [audit](fresh_h9_joint_arithmetic_20261013_independent.json).

The search retains the previous height bound 1000 and frees **all four A4 Cartan covectors plus the scalar**. It conditions the dictionary only on the original 13-wall chamber, allowing chronological replacement rather than assuming that all previous high walls must remain. In the resulting selection all 22 old walls do, in fact, remain.

There are 53 exact time/norm controls. They yield 190 roots after the simple chamber bounds, 186 roots compatible with the original 13 walls, and 184 chronologically selected additions. Fourteen dictionary roots cut the old 22-wall escape. The resulting cones have 197 H9 walls and 191 H5 walls. The dictionary contains 2 norm-one and 184 norm-two roots; the complete norm-five and norm-ten controls produce no dictionary root at this bound.

The producer uses the inverse-A4 identity, for covectors `(a,b,c,d)`,

`E=(4a²+6ab+6b²)/5 + (4d²+6dc+6c²)/5 + 2(a+2b)(d+2c)/5`.

Integrality is equivalent to `a+2b = d+2c (mod 5)`. A residue-indexed pair join therefore replaces a large four-variable loop. Each adjacent pair also satisfies its exact Schur bound `2(a²+ab+b²)/3 <= total_energy` at all three embeddings. The inherited cap inequalities are `a+b+c+d <= T*x` and `2z <= T*x`.

The independent audit uses standard-library Fraction arithmetic, a trace-ellipsoid field-integer alphabet, sequential A2/A3 Schur bounds, and a final Cartesian mod-five join. It reproduces all 53 controls, all 190 raw counts, the 186-root dictionary, the selected wall set, and the old-escape metadata. Its enumeration algorithm is different from the producer's pair join. This is a complete conditional complement census through 1000; **it is not a complete H9 root census**.

## Three constructive norm-equation steps

All six mirrors have norm 5. Since the ideal (5) is squarefree in this field, an integral vector of norm 5 is primitive. Every reflection coefficient `2Gr/5` is checked integral. PARI supplies proposals; the saved exact root identities and all geometric assertions are independently replayed in Python and do not depend on a PARI completeness claim.

1. **Middle fundamental-weight gap.** The 191-wall H5 face has escape `(1,(2T/3)omega3,0)`. The relative norm equation

   `30t²−Tx²=5`

   has the explicit integral solution

   `x=8715+30785T+16380T²`, `t=2181+7705T+4100T²`.

   The normal `(x,t(2,4,6,3),0)` and its A4 reversal are compatible with every old wall. The first cuts the escape by pairing `20+69T+35T²>0`. These give 199 H9 and 193 H5 walls. The enormous exact height is archived; there is no chronological-minimality claim. [Producer](../scripts/fresh_h9_joint_arithmetic_20261013_pell.py), [PARI proposal](../scripts/fresh_h9_joint_arithmetic_20261013_pell.gp), [artifact](fresh_h9_joint_arithmetic_20261013_pell.json). The full enumeration audit above also replays this pair.

2. **End fundamental-weight plus scalar gap.** The persistent H9 escape suggests the primitive spatial direction

   `d=((4,3,2,1),10+5T−5T²)`, `d²=120−50T`.

   An integral solution of `(120−50T)t²−Tx²=5` produces a new mirror and its reversal, both passing all 199 old wall tests. The first cuts that escape. The exact coefficients, of roughly 65 digits, are saved without rounding. This gives 201 H9 and 195 H5 walls. [Producer](../scripts/fresh_h9_joint_arithmetic_20261013_pell2.py), [PARI proposal](../scripts/fresh_h9_joint_arithmetic_20261013_pell2.gp), [artifact](fresh_h9_joint_arithmetic_20261013_pell2.json), [independent audit](fresh_h9_joint_arithmetic_20261013_pell2_independent.json).

3. **Three A4 covectors plus scalar.** Another simple H5 escape remains after both preceding orbits. Its exact spatial ray has norm `4700−1530T−440T²`. Direct PARI Hilbert tests did not propose a solution in the four allowed norm classes; this observation is not used as an independently certified nonrepresentation theorem. Instead the scout tests 40 explicitly prescribed perturbations of that ray by `±D^j*5omega_i` or `±D^j*5e_scalar`, for `D=2+T−T²`, `j=1,…,4`. Thirty-six pass exact continuous cutting/compatibility inequalities.

   Eleven of these directions receive positive norm-five Hilbert tests and are passed to the relative norm solver. Rational solutions with denominators are not accepted as lattice roots. Direction 33 gives an integral norm-five solution with 926-digit coefficients. Its spatial vector has norm

   `delta=5820−3760T+430T²`.

   The equation `delta*t²−Tx²=5` yields two reversal-related mirrors. Both pass every previous wall test, and the first cuts the three-covector H5 escape. The final cones have 203 H9 and 197 H5 walls. [Direction scout](../scripts/fresh_h9_joint_arithmetic_20261013_pell_scout.py), [40-direction artifact](fresh_h9_joint_arithmetic_20261013_pell_scout.json), [exploratory PARI log](fresh_h9_joint_arithmetic_20261013_pell_scout_pari.json), [finite PARI norm-solver replay](../scripts/fresh_h9_joint_arithmetic_20261013_pell_scout_norms.gp), [selected proposal](fresh_h9_joint_arithmetic_20261013_pell3_proposal.json), [producer](../scripts/fresh_h9_joint_arithmetic_20261013_pell3.py), [PARI proposal](../scripts/fresh_h9_joint_arithmetic_20261013_pell3.gp).

The final independent audit also rebuilds all 40 direction tests and verifies the 36 exact passes. It treats the Hilbert-symbol log as an external proposal record, not as a negative theorem. The finite-angle gate throughout allows only `4cos²(pi/m)` for `m=2,3,4,6,9,18` when the mirrors meet; a nonpositive pairing alone is never accepted as proof of a Coxeter angle.

## Exact remaining obstruction

In full H9 coordinates the rational point

`v=(1,0,0,0,0,499/500,1497/2000,499/1000,499/2000,287/500)`

satisfies all 203 wall inequalities and has

`v²=1903957/1000000−T > 0`.

It belongs to the same time+A4+A1 complement, so its projection also escapes the 197-wall H5 face. The independent checker moves it toward the archived strict timelike center and verifies a **strict** spacelike escape, which proves infinite volume for these particular cones. The simple rational witness replaces a much larger algebraic vertex certificate.

The three targets have distinct spatial supports: a middle A4 weight, an end weight with scalar, and three A4 covectors with scalar. The last surviving direction is again in an end-weight/scalar plane, now limited by one of the newly inserted mirrors. No integral isometry identifying a periodic full-chamber gap has been found or certified. Consequently the full H5 and H9 reflection chambers remain open; these results neither establish global nonreflectivity nor justify declaring a new compact polytope.

## Reproduction

```sh
PYTHONPATH=src:.:scripts python scripts/fresh_h9_joint_arithmetic_20261013.py
PYTHONPATH=src:.:scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell.py
PYTHONPATH=scripts python scripts/fresh_h9_joint_arithmetic_20261013_independent.py
PYTHONPATH=src:.:scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell2.py
PYTHONPATH=scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell_independent.py --round 2
PYTHONPATH=src:.:scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell_scout.py
PYTHONPATH=src:.:scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell3.py
PYTHONPATH=scripts python scripts/fresh_h9_joint_arithmetic_20261013_pell_independent.py --round 3
```

The full independent complement census took approximately 40 seconds; each later complete root/geometry replay took approximately 2 seconds on this run. PARI is optional for verification: the three saved positive proposals are exact, and their independent checks use no PARI or project field-arithmetic implementation.
