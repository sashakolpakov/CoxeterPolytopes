# TOOT: the `X2 H10/37` chamber has an explicit full-Picard K3 family

## Positive theorem

The already certified `X2-RAT-C2-H10-37` finite-volume noncompact Coxeter
chamber now has an equation-level characteristic-zero carrier whose **full**
Neron--Severi lattice is the chamber lattice.

On a nonempty Zariski-open part of the normalized nine-parameter family

\[
  X:\quad y^2=x^3+A(s^2)x+B(s^2),
\]

where

\[
\begin{aligned}
 A(t)&=-3+a_1t+a_2t^2+a_3t^3+a_4t^4,\\
 B(t)&=2+(1-a_1)t+b_2t^2+b_3t^3+b_4t^4+b_5t^5+b_6t^6,
\end{aligned}
\]

the very general member is a K3 surface with

\[
  NS(X)=F=U+E_8(-2)+\langle-2\rangle .
\]

The exact transport in the certificate puts all 37 primitive facet normals
in this full Picard lattice.  They generate all of `F`, and every one of the
37 reflections is integral on `F`.  Thus the 37-wall chamber is a full-rank
mixed divisor-reflection subsystem on every very-general marked member of
this explicit family.

This is a new project-level equation and full-Picard deformation realization
of the known `H10/37` chamber.  It is not a new abstract Coxeter polytope.

## Why this is more than the original `X_2` occurrence

The original construction lives inside the rank-20 Picard lattice of
Shimada's named singular K3 surface `X_2`: the rank-11 lattice `F` is the
saturated invariant divisor lattice of a rational involution.  That already
gave a genuine algebraic fixed-section construction, but `F` was not the
whole Picard lattice of `X_2`.

The present composition uses the pinned determinant-one isometry

\[
  F_{\rm standard}\longrightarrow NS(X_2)^\iota
\]

to transport every one of the 37 archived walls to the geometric basis of
the quadratic-base-change family.  The independent family theorem proves
that the very-general Picard lattice is exactly `F`, not a proper
sublattice of a larger Neron--Severi group.  This closes the full-Picard gate
for the 37-wall chamber.

The same equation family already carries the distinct 21-wall
`A15_11-C2-H10-21` chamber.  One explicit nine-dimensional complex-K3 family
therefore carries two different certified `H10` mixed divisor-reflection
chambers, with 21 and 37 facets.

## Exact lattice and wall checks

The standard lattice has

```text
rank                 11
signature            (1,10)
determinant          512
Smith invariants     1^2 2^9
Nikulin invariants   (r,a,delta)=(11,9,1)
K3 moduli dimension  9
```

An explicit primitive embedding in
`Lambda_K3=U^3+E8(-1)^2` is also archived:

- embed `U` in the first ambient `U`;
- embed `E8(-2)` diagonally in the two `E8(-1)` factors;
- send the generator of `<-2>` to `e-f` in the second ambient `U`.

The orthogonal complement is `U+E8(-2)+<2>`.  A unit minor proves
primitivity, and the combined embedding/complement index is 512, exactly the
discriminant order.

For the walls, the exact split is

```text
  6  primitive (-2), divisibility 1
 21  primitive (-2), divisibility 2
 10  primitive (-4), divisibility 2
```

The 37-by-11 wall-coordinate matrix has Smith diagonal `1^11`, so its rows
span the full Picard lattice with index one.  The complete transported facet
Gram matrix agrees entry by entry with the original `X_2` Gram matrix.  A
primitive square-12 positive class has all 37 oriented wall pairings
positive, with minimum pairing one; this independently fixes the common
halfspace orientation.

## Finite volume is checked separately

No volume conclusion is inferred merely from the lattice.  The package pins
the literal Vinberg finite-special-subdiagram computation:

```text
finite vertices                  3,236
ideal vertices / cusps             261
elliptic rank-9 edges           24,350
bad edge endpoint audits             0
bad affine-neighbourhood audits      0
```

Every elliptic rank-9 edge has exactly two finite-or-ideal endpoints, and
every irreducible affine zero-neighbourhood is a full rank-9 parabolic cusp
diagram.  The cusps prove noncompactness.  These are precisely the two
Vinberg gates needed for finite volume in this setting.

## Exact boundary

This is a mixed divisor carrier, **not the full nef cone** of the
very-general K3.  In particular, the ten square-minus-four mirrors are
divisor cuts, not K3 negative-curve walls.  After Weyl re-marking,
Riemann--Roch can orient the 27 square-minus-two classes effectively, but no
irreducibility statement is made.

The 37 integral Picard reflections are **not surface automorphisms** and are
not asserted to extend to Hodge isometries.  The original rational
involution on the special surface `X_2`, the deck involution of the explicit
family, and the abstract reflection group are three distinct pieces of
structure.  No worldwide priority claim is made for the abstract polytope,
the lattice, or the equation family.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_x2_h10_37_ue8m2a1_full_picard_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_x2_h10_37_ue8m2a1_full_picard_carrier.py
```

Frozen payload SHA-256:

`50c48a6ef7dcaaa4719b21d91f303a6d1983d5ab15cd998e4832d498a7ae4e08`

