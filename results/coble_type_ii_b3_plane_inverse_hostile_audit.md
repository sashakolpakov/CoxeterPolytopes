# Hostile audit: raw H10/30 to the Choi–Lam/B3 strong-Coble surface

Status: **GREEN, with one characteristic-3 correction** (2026-08-09).

The independent replay confirms the project-new claim in its proper scope:
the already known abstract noncompact finite-volume H10/30 chamber admits a
strict, equation-free-at-input transition from its raw lattice and walls to a
7-proper-plus-3-near blowup of the plane.  The forced antibicanonical curve is
the known Choi–Lam sextic.  The replay then returns exactly the same marked
thirty-wall chamber.

This is a new realization/cross-identification in this project.  It is not a
new abstract Coxeter polytope, not a new dimension, not compact, and not yet a
general necessary-and-sufficient polytope-to-surface theorem.

## Independent-input rule

The hostile program does not import
`certify_coble_type_ii_b3_plane_inverse.py`, its `build_certificate`, or its
helpers.  Before the final comparison with the frozen target it reads only:

- the rank-11 lattice Gram matrix and thirty wall vectors in
  `coble_one_boundary_full_ns_h10.json`;
- independent exact-cone, Vinberg-face, and P-first recognizer archives;
- the explicit algebra it reconstructs from those raw data.

The target JSON is read only after the reconstruction, to compare headline
claims.  Both branches transport all thirty wall vectors exactly, so this is
not merely a comparison of hashes or summary counts.

## Raw chamber reconstruction

Exact arithmetic recovers an odd unimodular lattice of signature `(1,10)` and
rank-11 wall span.  Among the eighteen norm-minus-one walls, wall 0 is the
unique characteristic vector.  Calling it `b=-K`, the other walls split
intrinsically as

- twelve norm-minus-two roots orthogonal to `b`;
- seventeen norm-minus-one vectors `E` with `b.E=1`, for which `f=E+b` is a
  primitive isotropic vector orthogonal to `b`.

The sequential contraction search excludes `b` and adds a wall only when its
class in the current quotient has square `-1`.  Deduplication is by the exact
wall subset.  The complete state counts by depths zero through ten are

```text
1, 17, 87, 329, 948, 2094, 3472, 4417, 4455, 3291, 573.
```

Thus all 573 numerical P2 endpoints, not a preselected contraction, are
tested.  They have 28 degree signatures.  Exactly two have boundary degree 3,
root degrees `0^3 1^9`, and contracted original norms `(-2)^3(-1)^7`:

```text
{1,5,9,13,15,17,19,22,26,28}
{3,7,11,13,15,17,19,22,26,28}.
```

Their lexicographically first valid contraction orders are respectively

```text
13,5,15,17,19,9,22,26,28,1
13,11,15,17,19,3,22,26,28,7.
```

The orthogonal plane classes in the raw basis are

```text
(13,9,19,27,23,20,15,12,7,3,3)
(12,9,18,26,22,19,14,11,7,3,3).
```

Appending the ten sequential total exceptional transforms gives Gram
`diag(1,-1^10)` and determinants `+1,-1`.  Hence both are integral
unimodular plane bases.  In each basis `K=(-3,1^10)`, `b=(3,-1^10)`, and every
one of the thirty raw walls transports integrally and reconstructs exactly.

## The forced B3 incidence object

In each branch, three roots have degree zero and become
`E_parent-E_near`; nine have degree one and become `H-E_i-E_j-E_k`.
The resulting colored incidence graph has exactly twelve isomorphisms to the
same 7+3 design.

Normalize three proper points as

```text
A=(1:0:0), D=(0:1:0), G=(0:0:1), B=(1:1:1).
```

Solving five incidences gives

```text
C=(1:1:c), E=(1/c:1:1), F=(c:1:c).
```

The remaining determinant is `(1-c^2)/c`.  The branch `c=1` collides `C`
with `B`; the unique nondegenerate solution is `c=-1` (so this incidence
argument excludes characteristic 2).  The nine root lines are precisely the
projectivized B3 reflection arrangement.  Thirty of its 36 line pairs meet at
the seven blown proper points; the other six meet transversely at distinct
unblown points.

## Exact surface checks over Q and C

The ten cubic conditions have rank 10, while the thirty sextic conditions
have rank 27 in a 28-dimensional space.  Therefore

```text
h0(-K)=0,   h0(-2K)=1.
```

The unique sextic is

```text
F = x^4 y^2 + y^4 z^2 + z^4 x^2 - 3 x^2 y^2 z^2.
```

Singular exact factorization, absolute factorization, saturated Jacobian, and
radical computations give:

- one absolute component, with exponent one;
- singular scheme length 13 and reduced support degree 7;
- reduced support exactly `A,D,G,B,C,E,F`;
- three A3 tacnodes at `A,D,G`, whose declared first-near blowups leave
  ordinary nodes;
- four A1 nodes at `B,C,E,F`.

The seven proper and three declared near blowups resolve the curve.  Its
strict transform is smooth, geometrically integral, rational, has class
`-2K`, and square `-4`.  The twelve root curves are nine line transforms and
three exceptional-chain transforms, hence smooth geometrically integral
rational minus-two curves.  Their 66 pairings are `15` ones and `51` zeros.

Every root has intersection zero with the boundary.  Absolute irreducibility
excludes a root line as a boundary component, and the strict transform has no
exceptional component.  Distinct effective irreducible curves on a smooth
surface have nonnegative local intersection, so total intersection zero
forces actual disjointness.

## Finite volume and exact same-P return

The hostile program independently feeds the raw 30-by-11 inequality matrix
to GMP cddlib.  It obtains

```text
1540 extreme rays = 1333 finite + 207 ideal + 0 hyperideal.
```

The active-facet histogram is

```text
finite_10: 1333
ideal_10: 48
ideal_11: 98
ideal_12: 52
ideal_13: 9.
```

All thirty walls are irredundant facets.  The sorted primitive-ray digest is

```text
f88e2b4397b6e6c27cee784cd1ecad4f534dd39decf219774e3f3f12b83e94f9.
```

This independently confirms finite volume and noncompactness; absence of
hyperideal extreme rays is the relevant complete double-description test.
The separate Vinberg archive agrees: all 10,030 elliptic rank-9 edges have
two endpoints, with `ff=4423`, `fi=4484`, `ii=1123`, and all 436 connected
affine components pass their full parabolic-neighborhood test.

Finally, the raw inventory returns as `1+12+17=30`: boundary, explicit roots,
and the seventeen `E=f-b` cusp lifts.  The independent P-first recognizer
finds the unique Type-II candidate, the same complete seventeen cusps, root
span index two in E10, and the unique E10 saturation up to diagram
automorphism.  Coupled with the exact thirty-wall transports, this is the
same marked H10/30 chamber, not merely an isomorphic root subdiagram.

## Characteristic 3: the corrected rejection

Do not extend the phrase “three A3 plus four A1” to every characteristic other
than 2.  In characteristic 3 the four discriminants `-48` vanish and the
tangent cones at `B,C,E,F` become repeated.  The saturated singular scheme
has length 21 rather than 13.

An earlier scratch computation incorrectly selected only the first generator
of Singular's `sat(...)` ideal and thereby suggested a positive-dimensional
singular locus.  That is false.  Correct saturation gives affine-cone
dimension 1 (hence zero-dimensional projective support), radical degree 7,
and the same seven reduced points.  A first blowup makes the strict curve
smooth at each of the four repeated-tangent points, but they are not A1 nodes.
The audit makes no claim here about absolute irreducibility over
`F3bar`.  The proved surface statement remains over Q/C.

## What this contributes to a P-to-surface theorem

This example supports a reusable sufficient gate, formulated entirely from a
marked chamber before the final AG checks:

1. Find a primitive characteristic `b^2=-1` wall and the intrinsic
   `b + roots + cusp-lifts` partition.
2. Enumerate sequential rank-ten contraction endpoints and require a
   unimodular plane basis.
3. Turn transported root classes into an incidence realization scheme and
   find a nondegenerate point of that scheme.
4. Check the anticanonical/antibicanonical Hilbert functions, boundary
   irreducibility, and resolution conditions.
5. Apply a complete curve/wall exhaustion theorem and an exact cone test to
   prove that the surface returns all—and only—the walls of P.

For this marked H10/30 chamber those conditions are sufficient, and the
lattice partition, plane marking, effective root incidence, Hilbert
functions, and full-wall return are necessary for any carrier of this same
marked rational type.  This is genuine theoretical leverage beyond “a
surface happened to give a polytope”: P algorithmically forces a finite
incidence scheme, and AG decides its realization and boundary.

What remains for an L4 theorem is substantial and explicit: prove that an
allowed list of carrier types and contraction signatures is exhaustive for a
variable class of arithmetic chambers, then classify the resulting incidence
schemes and effectivity strata.  Empty schemes/strata would give
nonexistence and potentially dimension bounds; rational points would give new
examples.  The current H10/30 result alone says neither that every compact
arithmetic polytope has a surface carrier nor that a higher-dimensional
bound follows.

## Lineage

The sextic is the cyclic ternary sextic associated with Man-Duen Choi and
Tsit-Yuen Lam, *Extremal Positive Semidefinite Forms*, **Mathematische
Annalen** 231 (1977), 1–18,
[DOI 10.1007/BF01360024](https://doi.org/10.1007/BF01360024).  The equation is
not new.  As a formula-level cross-check, Blekherman–Kozhasov–Reznick print
this exact cyclic sextic and explicitly attribute its study to Choi and Lam in
[Forum of Mathematics, Sigma (2026)](https://doi.org/10.1017/fms.2026.10221).
Nor is the abstract D16-D5 H10/30 chamber.  The project-new content
claimed here is the raw-P derivation of the two intrinsic plane branches, the
forced B3/Choi–Lam identification, the exact strong-Coble verification, and
the same-marked-P return.  No worldwide priority claim is made for every
possible historical Choi–Lam/Coble connection.

## Frozen artifacts and replay

```bash
python3 -m py_compile scripts/audit_coble_type_ii_b3_plane_inverse_hostile.py
python3 scripts/audit_coble_type_ii_b3_plane_inverse_hostile.py \
  --check-archive data/enumerations/coble_type_ii_b3_plane_inverse_hostile_audit.json
python3 -m pytest -q tests/test_coble_type_ii_b3_plane_inverse_hostile_audit.py
```

The live replay requires exact `Singular`, GMP `cddexec_gmp`, SymPy, and
NetworkX.  On the audited machine it takes roughly half a minute.
