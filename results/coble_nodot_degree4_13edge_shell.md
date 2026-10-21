# The exact degree-four, thirteen-edge nodot shell is hyperideal

**Status (2026-08-09): proved P-side shell theorem; no new finite-volume
polytope and no algebraic-geometric leverage delta.**

The executable certificate is
[`certify_coble_nodot_degree4_13edge_shell.py`](../scripts/certify_coble_nodot_degree4_13edge_shell.py),
with frozen
[`JSON`](../data/enumerations/coble_nodot_degree4_13edge_shell.json)
and
[`tests`](../tests/test_coble_nodot_degree4_13edge_shell.py).

## The intrinsic shell

Let `S_Delta=4,e=13,kappa<=3` consist of weighted root data with:

- twelve norm-minus-two wall symbols;
- connected ordinary support with thirteen edges and maximum degree exactly
  four;
- off-diagonal root pairings in `{0,1,2}`, with exactly `kappa` pairing-two
  edges for `0 <= kappa <= 3`;
- exact root-Gram inertia `(1,9,2)`.

No determinant, `E10` saturation, boundary, surface, carrier, or known-diagram
condition occurs in membership.  This is the first frozen shell obtained by
relaxing the prior degree-at-most-three restriction while holding the
thirteen-edge and small-pairing coordinates fixed.  It is not claimed to be
all degree-four diagrams at every edge count.

The connected degree-at-most-four census has 16,977 ordinary supports.  The
pinned degree-at-most-three lane has 2,888, leaving exactly 14,089 supports
of maximum degree four.

## Exact census

| pairing-two count `kappa` | weighted trials | Lorentzian marked presentations | weighted isometry classes | finite-volume classes |
|---:|---:|---:|---:|---:|
| 0 | 14,089 | 1 | 1 | 0 |
| 1 | 183,157 | 26 | 24 | 0 |
| 2 | 1,098,942 | 189 | 145 | 0 |
| 3 | 4,029,454 | 422 | 330 | 0 |
| **total** | **5,325,642** | **638** | **500** | **0** |

Every exact Lorentzian presentation is retained, not merely determinant-four
or even-unimodular cases.  The 500 quotient determinants occupy 31 values;
in particular the census includes 45 determinant-minus-four classes and
seven determinant-minus-one classes.  Hence the negative conclusion is not
an artefact of an early Coble-lattice filter.

## Theorem

> **Degree-four shell theorem.**  No weighted diagram in
> `S_Delta=4,e=13,kappa<=3` defines a finite-volume hyperbolic root chamber.
> Each of the 500 exact Lorentzian weighted-diagram isometry classes has an
> exact hyperideal extreme ray in its saturated radical quotient.

For each class the certificate performs the following exact replay.

1. It isolates the saturated rank-ten radical quotient by a unimodular HNF
   transformation and checks quotient inertia `(1,9,0)`.
2. An exact Perron iteration supplies a timelike point with strictly positive
   wall slacks.  Thus the chamber is nonempty in the chosen time cone.
3. GMP `cddlib` enumerates the complete extreme-ray set.  Every input wall is
   irredundant, and all causal rays lie in one time cone.
4. At least one extreme ray has negative Lorentzian square in every class.
   This is a hyperideal endpoint, so Vinberg's finite-volume test stops and
   fails.

Across all classes the complete extreme-ray ledger is

```text
finite       2,588
ideal        1,182
hyperideal  12,666
```

The presence of finite vertices or cusps in an individual cone does not
compensate for a hyperideal ray.  A single hyperideal extreme direction gives
an open end beyond the absolute and therefore infinite hyperbolic volume.
Passing to an integral overlattice cannot repair this: saturation changes the
arithmetic coordinate lattice, not the real wall cone or its ray squares.

Every class record contains its complete negative Cartan matrix, an exact
strict-interior certificate, one hyperideal-ray witness, the exact ray-type
counts, and hashes of both the complete and hyperideal ray lists.

## Before/after ledger and honest leverage score

```text
intrinsic degree-four weighted classes        500
after the exact Vinberg finite-volume gate      0
before algebraic geometry                       0
after algebraic geometry                        0
AG leverage delta                               0
```

The computation stops before carrier enumeration.  It proves a useful search
boundary, but it does **not** prove a `POLYTOPE -> SURFACE -> new theorem on
POLYTOPE` result.  The appropriate score is

```text
SCOPED_P_SIDE_NEGATIVE_L2_L3_NOT_L4.
```

There are no new compact or finite-volume noncompact examples in this shell.
There is also no claim about degree-four supports with other edge counts,
four or more pairing-two edges, dotted edges, or other root norms.

## Smallest credible broadened-class target

There are two statements here which must not be conflated.

### What is already proved

The broader, surface-free `D_CobVin` class already has a finite exhaustive
all-contraction incidence stack and an intrinsic exact-return criterion:
over an algebraically closed field of characteristic different from two, an
input has an exact marked strong-Coble full-nef carrier if and only if its
compiled incidence stack has a geometric point.  Every geometric point
returns literally the original oriented `P`.

The mixed calibration is

```text
raw arithmetic inputs:                 Type I, Type II, classical P10
complex one-boundary Coble carriers:   Type I, Type II
AG-obstructed Coble carrier:           classical P10
```

This is a nonzero delta for the **carrier locus**, not a nonexistence theorem
for polytopes: classical `P10` remains an actual hyperbolic polytope and has a
different rational-surface realization.  It therefore does not cross the
requested flagship gate.

### What remains a proposed target

The smallest rigorous Coble target is a P-computable subdenominator
`C_free` inside `D_CobVin`.  At least one P-derived contraction branch should
satisfy, without solving a surface instance:

1. an acyclic dependency order of free blow-up points and prescribed
   transverse intersections;
2. base-point-freeness, `H1` vanishing, and evaluation-map surjectivity for
   each successive curve class;
3. independent expected codimension for all ordinary and infinitely-near
   conditions;
4. a numerical exact-sequence certificate forcing
   `h0(-K)=0` and `h0(-2K)=1`.

A uniform realization lemma for that predicate would force an incidence
point.  The existing exact-return theorem would then apply, and Kondō's
independent classification would push back

```text
boundary count in {1,2,5,8};
hyperbolic dimension in {10,11,14,17};
no members in dimensions 12,13,15,16, or above 17.
```

For a genuinely nonzero AG leverage delta, a surface-free census of
`C_free` must leave apparent P-side candidates outside those four dimensions
which Kondō, rather than the Coxeter cone test, removes.  Neither the uniform
realization lemma nor that broadened census is presently proved.

Classical `P10` cannot belong to any subclass on which a **complex
one-boundary strong-Coble** carrier is forced: its empty Coble carrier stack
is already certified.  A theorem class required to contain `P10` must enlarge
the carrier species, for example to suitable rational Halphen or
pencil-at-infinity surfaces.  For that broader species, finite hidden-curve
exhaustion and an independent classification strong enough to push back to
the original P class are still missing.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_degree4_13edge_shell.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_degree4_13edge_shell.py
```
