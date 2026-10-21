# TOOT: one explicit K3 family now carries the 13-, 21-, and 37-wall `H10` chambers

## Positive theorem

The known finite-volume noncompact `K3-S3-C2-H10-13` chamber has a second,
equation-level realization.  On the nonempty open part of the normalized
nine-parameter family

\[
 X:\quad y^2=x^3+A(s^2)x+B(s^2),
\]

with

\[
\begin{aligned}
 A(t)&=-3+a_1t+a_2t^2+a_3t^3+a_4t^4,\\
 B(t)&=2+(1-a_1)t+b_2t^2+b_3t^3+b_4t^4+b_5t^5+b_6t^6,
\end{aligned}
\]

a very general member is a K3 surface with full Picard lattice

\[
 NS(X)=F=U+E_8(-2)+\langle-2\rangle .
\]

An explicit determinant-one integral isometry transports the complete
13-wall system from the published rank-19 K3 transposition fold to the
geometric basis of `F`.  The thirteen primitive normals span full Picard and
all thirteen reflections act integrally on it.  Thus every very-general
marked member carries this known chamber as a full-rank mixed
divisor-reflection subsystem.

This is a new project-level equation and very-general full-Picard carrier of
the known chamber, not a new abstract Coxeter polytope.

## The missing integral bridge is now explicit

Let `G` be the archived Gram matrix of the saturated transposition-fixed
lattice.  Its displayed divisibility-one isotropic pair is

```text
v = [0,  0, 0,-1,-1,-2,-2,-3,-2,-1,1]
w = [0, -1, 0,-1,-1,-2,-2,-3,-2,-1,1]
```

with `v^2=w^2=0` and `(v,w)=1`.  Splitting off this integral `U` and solving
the remaining exact positive-definite rank-nine isometry gives the archived
matrix `Q` satisfying

\[
 \det Q=1,\qquad Q^tGQ=U+E_8(-2)+\langle-2\rangle.
\]

Both `Q` and its integral inverse are stored entry by entry.  This is the
machine gate that was absent from the earlier quotient-nef certificate.

The resulting lattice data are

```text
rank                 11
signature            (1,10)
determinant          512
Smith invariants     1^2 2^9
Nikulin invariants   (r,a,delta)=(11,9,1)
K3 moduli dimension  9
```

The package also replays an explicit primitive embedding in
`Lambda_K3=U^3+E8(-1)^2`; its orthogonal complement is
`U+E8(-2)+<2>`.

## Exact wall transport

In the geometric `F` basis the primitive mirror types are

```text
2  of square -2 and divisibility 1
2  of square -2 and divisibility 2
9  of square -4 and divisibility 2
```

The 13-by-11 coordinate matrix has Smith diagonal `1^11`, proving index-one
generation of the full Picard lattice.  The transported facet Gram agrees
entry by entry with the source Gram.  The positive class

```text
[-42,-40,2,-6,-7,-11,-8,-14,3,21,19]
```

has square 358 and all thirteen wall pairings positive, with minimum two.
This fixes the common halfspace orientation without floating point.

## Finite volume was rerun from the Gram

The new certificate does not merely copy the earlier classification.  It
reruns three exact lanes from the transported 13-wall Gram:

1. literal all-subset CoxIter Proposition 1.10 extension enumeration;
2. the independent optimized Vinberg full-neighbourhood test; and
3. exact enumeration of all projective extreme rays.

All three return 34 finite vertices, four cusps, and no hyperideal rays.  The
literal enumeration finds 208 elliptic rank-nine edges, each with exactly
two finite-or-ideal endpoints.  Its endpoint histogram is

```text
(2 finite, 0 ideal)  137
(1 finite, 1 ideal)   66
(0 finite, 2 ideal)    5
```

The optimized lane independently performs seven affine-neighbourhood
audits; all seven are full rank-nine parabolic cusp neighbourhoods.  Hence
the chamber has finite volume and is noncompact.

## Three chambers on one equation family

The same explicit nine-dimensional family already carries the certified
`A15_11-C2-H10-21` and `X2-RAT-C2-H10-37` chambers.  It now carries three
distinct mixed divisor-reflection chambers with 13, 21, and 37 facets.  This
is a useful positive criterion: once an exact chamber lattice is identified
integrally with `F`, the equation family supplies hyperbolicity, primitive K3
embedding, and very-general full Picard automatically; only wall transport,
span, reflection integrality, and Vinberg finite volume remain discrete
checks.

## Exact boundary

Kondo's geometric `S3` action and the quotient-nef interpretation belong to
the special rank-19 source surface.  They are not asserted to persist on a
very-general member of the quadratic-base-change family.

This mixed chamber is **not the full nef cone** of the very-general K3.  The
nine square-minus-four mirrors are integral divisor reflections, not K3
negative-curve walls.  After Weyl re-marking, Riemann--Roch can orient the
four square-minus-two classes effectively, but irreducibility is not
claimed.  The thirteen Picard reflections are **not surface automorphisms**
and are not asserted to extend to Hodge isometries.  No worldwide priority
claim is made.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.py
```

Frozen payload SHA-256:

`c033a684aa04b2c6d0bfe80d1b5393f0c52f5a9c9c1039140ca674ea9cc67a34`
