# Hostile audit: smooth-toric full-nef Coxeter theorem

**Verdict (2026-08-09): the scoped direct theorem is correct.**  The sharp
ceiling is `H6`, the endpoint is a finite-volume noncompact `H6/9` chamber,
and its exact ray profile is `18 finite + 3 ideal + 0 hyperideal`.  The
abstract polytope is known from Im Hof's Napier-cycle classification.  The
current inverse is an exact finite recognition theorem, but it remains below
the project's strict `L4` polytope-to-surface leverage target.

## 1. Logical directions

Let `D_i` be the invariant boundary curves of a smooth projective toric
surface `X`.  They generate `NE(X)`, and adjunction gives

```text
(-K_X).D_i = 2 + D_i^2.
```

Thus `-K_X` is nef if and only if every `D_i^2>=-2`.  Every negative
irreducible curve is an actual extremal Mori ray: for an ample `A`, a suitable
positive multiple of the curve added to `A` supports precisely that ray (and
torically this is the contraction obtained by deleting its fan ray).  Hence a
boundary curve with square at most `-3` cannot hide as a redundant negative
inequality.

Each `D_i` is primitive because it intersects an adjacent boundary curve
once.  Since `NS(X)` is unimodular, reflection in `D_i^perp` is integral if
and only if `D_i^2` divides `2`, hence exactly for squares `-1` and `-2`.
This proves every direction of

```text
all hyperbolic nef facets integral
  <=> all invariant boundary squares >= -2
  <=> -K_X nef.
```

The pairwise angle assertion is also complete: disjoint negative boundary
curves give `m=2`, adjacent `(-2,-2)` gives `m=3`, adjacent `(-1,-2)` gives
`m=4`, and adjacent `(-1,-1)` gives asymptotic mirrors (`m=infinity`).

One wording restriction is essential.  A square-zero or positive invariant
curve need not be redundant in the **linear** nef-cone presentation.  Its
inequality is automatic only after restricting to the chosen sheet of
hyperbolic space: every nonzero effective causal class lies in the same
future cone as an ample class, and therefore pairs strictly positively with
every point of that hyperbolic sheet.

## 2. Volume and noncompactness

The finite-volume implication does not need to be smuggled in as a Vinberg
assumption.  The complete toric nef cone is rational polyhedral.  Every nef
class has nonnegative square, by approximation with ample classes.  A compact
affine slice of the cone therefore has finitely many vertices in the closed
Klein ball and no hyperideal vertex.  In hyperbolic dimension at least two
(`rho>=3`), its intersection with hyperbolic space has finite volume.

Noncompactness is separate.  Contracting a rational toric surface to a
Hirzebruch surface, or through `F_1` when the minimal endpoint is `P2`, and
pulling back a ruling gives a nonzero nef class of square zero.  A null ray in
a future Lorentz cone is automatically extremal, so it is an ideal vertex.

## 3. Sharp fan bound

For `-K_X` nef, it is also big, and the anticanonical model is a Gorenstein
toric del Pezzo surface.  Define `Delta` to be its **fan polygon in `N`**.
The smooth fan is the crepant regular subdivision, so its rays are exactly
the lattice points on `partial Delta`; with this orientation,

```text
r = B(Delta),   B(Delta)+B(Delta^vee)=12,   B(Delta^vee)>=3.
```

Consequently `r<=9`, `rho=r-2<=7`, and the hyperbolic dimension is at most
six.  One must not write `B(Delta)=r` for an unspecified anticanonical or
polar polygon; the fan-polygon convention (or the reversed convention with
the roles exchanged) is part of the proof.

An independent toric blowup replay starts from `P2,F0,F2`, applies
`(a,b)->(a-1,-1,b-1)` only when `a,b>=-1`, and gives the exact fan-ray
histogram

```text
r=3,4,5,6,7,8,9 : 1,3,2,4,2,3,1.
```

The unique endpoint cycle is `(-2,-1,-2)^3`.

## 4. Endpoint and literal Vinberg cusp test

The independent replay uses the nine displayed divisor classes in
`I_(1,6)` but imports no production-certifier code.  Exhausting all
`binom(9,6)=84` active subsets gives exactly 21 cone rays:

```text
18 timelike, 3 lightlike, 0 spacelike.
```

The three null rays are `H-E_1`, `H-E_2`, and `H-E_3`.  Each has seven active
walls.  Its complete active Gram has

```text
rank 5, inertia (0 positive, 5 negative, 2 zero),
```

and splits into two connected affine components: a four-node component of
rank three (labels `4-3-4`, affine `C_3`) and a three-node component of rank
two (labels `4-4`, affine `C_2`).  Their affine ranks add to `3+2=5=n-1`.
Thus every full zero-neighborhood is a complete Euclidean cusp; no partial
affine subdiagram has been substituted for Vinberg's criterion.

Reproduce this independent check with

```bash
PYTHONPATH=.:src python3 scripts/independent_toric_full_nef_coxeter_hostile_audit.py
```

## 5. Worldwide lineage

The endpoint diagram is a nine-cycle.  Reading clockwise from the top edge
of the dimension-six lower-right diagram in Im Hof's Table III gives

```text
3,4,4,3,4,4,3,4,4,
```

which is cyclically `(4,4,3)^3`.  Im Hof proves that every listed
crystallographic Napier cycle produces a finite-covolume Coxeter group.
Therefore the abstract `H6/9` polytope is known, not a worldwide-new
example.  The defensible novelty is the project-new toric full-nef
realization and the sharp classification of this carrier species.

Primary lineage: H.-C. Im Hof, *Napier cycles and hyperbolic Coxeter groups*,
Bull. Soc. Math. Belg. Ser. A **42** (1990), 523--545, Table III.  The compact
three-page diagram extract maintained by Felikson and Tumarkin is useful for
visual comparison but should not be presented as the primary article itself.

## 6. Inverse status under the project's strict contract

The twelve-row procedure is a genuine terminating `P -> surface` **recognizer**:
for a marked integral lattice polytope it tests a full marked-lattice
isometry, primitive wall lengths, fan recurrence, cyclic closure, winding,
and positive-cone orientation, and acceptance returns the same full nef
chamber.  A normalized Gram match alone is not sufficient.

This merits `L3` (exact same-`P` finite carrier recognition), not yet strict
`L4`.  The accepted class is defined by matching one of the already
classified twelve toric rows, so the subsequent `H6` bound does not shrink an
independently defined variable `P`-only preclass; the pre/post AG delta is
zero at that level.  A strict `L4` promotion would require an independent,
non-singleton `P`-only class--for example a saturated cyclic-completion
condition stated without surface existence or the twelve-row list--together
with an iff carrier theorem and a surface-classification consequence that
removes candidates or proves a new property of those original polytopes.

## 7. Required archival cautions

1. State redundancy of nonnegative boundary curves only on the hyperbolic
   sheet.
2. Define which reflexive polygon has boundary points equal to fan rays.
3. Record the three full affine cusp neighborhoods, not only the ray signs.
4. Orient candidate rays by all wall pairings; using only the first vector
   coordinate is not a formally complete cone-enumeration convention.
5. Describe the inverse as a full marked-lattice test and rate it `L3` until
   a nonzero classwide AG leverage delta is exhibited.

