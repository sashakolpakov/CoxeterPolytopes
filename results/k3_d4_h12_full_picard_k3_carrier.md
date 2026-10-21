# TOOT: the `K3-D4 H12/17` chamber has a full-Picard complex-K3 carrier

## Positive theorem

The known 17-wall finite-volume noncompact chamber obtained from the
order-three `D4` triality fold now has a very-general complex-K3 carrier whose
**full** Neron--Severi lattice is the chamber lattice:

\[
 S=U\oplus E_8(-1)\oplus A_2(-1)\oplus A_1(-1).
\]

There is a seven-dimensional family of projective K3 surfaces with
very-general `NS(X)=S`.  On every very-general marked member, the 17 archived
primitive normals define integral reflections of full `NS(X)` and generate
the whole Picard lattice.  The exact chamber remains the full reflection
chamber of `S`, with 93 finite vertices and three cusps.

This is a new project-level full-Picard divisor-reflection realization of a
known abstract chamber.  It is not a new lattice, reflection group, or
hyperbolic Coxeter polytope.

## The fixed kernel is explicit and saturated

The source K3 Picard lattice is

\[
 U\oplus E_8(-1)\oplus D_4(-1)\oplus A_1(-1).
\]

The order-three divisor-lattice action is `D4` triality.  Its saturated fixed
kernel has rank 13, signature `(1,12)`, determinant 6, and Smith invariants
`1^12,6`.  The pinned determinant-one change of basis identifies it exactly
with

\[
 U\oplus A_2(-1)\oplus A_1(-1)\oplus E_8(-1).
\]

This uses the lattice action only to construct the fixed divisor lattice and
its walls.  Triality is **not a surface automorphism** of the very-general
carrier, and no such assertion is needed for the theorem.

## Explicit primitive K3-lattice embedding

The certificate gives a literal 22-by-13 embedding matrix into

\[
 \Lambda_{K3}=U^3\oplus E_8(-1)^2.
\]

In words:

- the target `U` enters the first ambient `U`;
- the target `E8(-1)` enters the first ambient `E8(-1)`;
- `A2(-1)` uses two adjacent simple-root basis vectors in the second
  `E8(-1)`;
- `A1(-1)` uses a third simple-root basis vector orthogonal to those two.

Those 13 ambient coordinates contain a determinant-one minor, so the
embedding is primitive.  The matrix equality with the target Gram is checked
entry by entry.  Its orthogonal complement has signature `(2,7)`, hence the
corresponding K3 period domain has dimension seven.  Avoiding the countable
Noether--Lefschetz hyperplanes gives `NS(X)=S` on the very general member,
rather than only an embedding in a larger Picard lattice.

## All 17 mirrors live on full Picard

The exact primitive mirror types in the standard `S` basis are

```text
13  square -2, divisibility 1
 2  square -2, divisibility 2
 2  square -6, divisibility 3
```

Every reflection is checked for integrality, determinant `-1`, order two,
and preservation of the full Picard Gram.  The 17-by-13 coordinate matrix
has Smith diagonal `1^13`, so the mirrors span `S` with index one.

The transported facet Gram agrees literally with the original 17-wall orbit
Gram.  A primitive positive class has square 1,110 and wall pairings

```text
2 (thirteen times), 54, 54, 58, 58.
```

It is therefore a strict interior class and fixes all halfspace
orientations.  After Weyl re-marking, K3 Riemann--Roch gives effective signs
for the fifteen square-minus-two classes; irreducibility is not claimed.

## Literal finite-volume theorem

Volume is not inferred from signature or from the K3 embedding.  The new
package pins both independent archived volume lanes.  The literal
CoxIter-Proposition-1.10/Vinberg calculation gives

```text
elliptic rank-11 edges       623
finite vertices               93
ideal vertices / cusps         3
hyperideal extreme rays         0
bad edge extensions             0
bad affine cusp neighbourhoods  0
```

Every edge has exactly two endpoints.  Their type histogram is

```text
496  (two finite, zero ideal)
124  (one finite, one ideal)
  3  (zero finite, two ideal).
```

Every affine zero-neighbourhood is a full rank-11 parabolic cusp diagram.
An optimized Vinberg replay independently returns the same 623 edges, 93
finite vertices, and three cusps.

The separate all-reflective-root height audit proves that these 17 roots are
the complete simple system of the full reflection group.  Thus the subgroup
index is one.  Scharlau--Kirschmer already record the same determinant-six
lattice and chamber as `List12/Gram12` entry 6.

## Exact boundary

This is a mixed divisor-reflection carrier, **not the full nef cone** of the
very-general K3.  In particular, the two square-minus-six mirrors are
integral divisor cuts but cannot be negative-curve walls on a K3 surface.

Neither the order-three triality action nor the 17 Picard reflections are
claimed to extend to Hodge isometries or surface automorphisms of the
very-general carrier.  The result is period-level rather than equation-level;
constructing a seven-moduli Weierstrass family with geometric
`E8+A2+A1` marking remains a separate positive upgrade.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_k3_d4_h12_full_picard_k3_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_k3_d4_h12_full_picard_k3_carrier.py
```

Frozen payload SHA-256:

`d41e18b550d600d42ec75a7e068073e8bbd51366af975c25e6efdc798ef62550`
