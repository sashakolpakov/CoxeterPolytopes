# TOOT: the Hessian `H10/27` chamber has a very-general full-Picard K3 carrier

## Positive theorem

Let `S_H` be the saturated fixed divisor lattice of an `S5` transposition in
the rank-16 Picard lattice of a general Hessian quartic K3.  There is a
**nine-dimensional** family of projective complex K3 surfaces whose very
general member has

\[
 \operatorname{NS}(X)=S_H.
\]

All 27 primitive mirrors of the generic-Hessian transposition fold are
integral reflections of this full Picard lattice, and the mirrors generate
`S_H` with index one.  They cut out the known finite-volume noncompact
`D16-D4-C7-H10/27` chamber.

This is a new project-level very-general full-Picard carrier of a known
chamber.  It strengthens the existing Hessian divisor-stock realization; it
does not add a new slot to the maintained geometric count and is **not a new
abstract polytope**.

## Exact lattice and discriminant form

The certifier reconstructs the fixed kernel from the twenty nodal and line
curves, rather than copying its invariants.  In the saturated source basis its
Gram is even of rank 11 and signature `(1,10)`, with

```text
determinant       256
Smith diagonal    1^4 2^6 4
discriminant      (Z/4) + (Z/2)^6.
```

A determinant-minus-one change of basis splits a literal hyperbolic plane:

\[
 S_H=U\oplus K_H,
\]

where the certificate displays the full negative-definite rank-nine Gram of
`K_H`.  More sharply, it gives seven explicit dual-lattice generators and
puts the finite quadratic form in the orthogonal normal form

\[
 q_{S_H}=\langle 3/4\rangle\oplus u(2)^{\oplus2}\oplus v(2).
\]

Here the cyclic generator has order four and quadratic value `3/4`; each
`u(2)` pair has generator values zero and mutual bilinear value `1/2`; the
`v(2)` pair has generator values one and mutual bilinear value `1/2`.  The
seven displayed generators exhaust all 256 discriminant classes.

## Explicit primitive K3 embedding

The embedding is a literal graph-glue construction, not a stable-range
shortcut.  Write `P_H` for the reconstructed rank-16 Hessian Picard Gram and
set

\[
 T_H=U\oplus U(2)\oplus A_2(-2).
\]

The certificate exhibits generators of orders `6,2,2,2` for both
discriminant groups.  They generate all 48 classes, and their quadratic and
bilinear forms are entrywise opposite.  Gluing their graph produces a
displayed 22-by-22 integral even Gram of determinant minus one and signature
`(3,19)`.  It is therefore the K3 lattice

\[
 \Lambda_{K3}=U^3\oplus E_8(-1)^2.
\]

The certificate then gives the 22-by-16 embedding of `P_H` and the induced
22-by-11 embedding of `S_H`.  Both matrices are integral and have Smith
diagonal consisting entirely of ones; both embeddings are primitive and
their Gram identities hold entry by entry.  The rational orthogonal
complement of `S_H` has signature `(2,9)`.

Consequently the `S_H`-polarized K3 period domain is nonempty and has
dimension nine.  Away from the countable union of Noether--Lefschetz
hyperplanes, the full Neron--Severi lattice is exactly `S_H`, not a larger
ambient Hessian lattice.

## All 27 oriented mirrors span full Picard

The primitive fixed-lattice mirror census is

```text
 8  square -2, divisibility 1
10  square -4, divisibility 2
 9  square -4, divisibility 4.
```

For every mirror the package records its source wall orbit and all eleven
Picard coordinates.  Every reflection is checked for integrality,
determinant minus one, order two, preservation of the full Picard Gram, and
negation of its primitive normal.  The 27-by-11 coordinate matrix has Smith
diagonal `1^11`, proving index-one generation.

The invariant Hessian class

```text
[0,-3,2,2,2,0,0,-3,-3,-3,5]
```

has square 20 and pairs positively with all 27 primitive normals, with
minimum pairing one.  Thus it fixes the orientations without importing a
floating-point chamber choice.

## Two independent finite-volume gates

The literal Vinberg edge-and-parabolic program is rerun on the primitive
27-wall Gram.  It returns

```text
finite vertices                    895
cusps                               67
elliptic rank-nine edges         6,221
bad edge extensions                  0
bad parabolic cusp neighbourhoods    0.
```

Every elliptic edge has exactly two finite or ideal endpoints, and every
maximal parabolic neighbourhood passes the cusp test.

Independently, GMP rational double description enumerates all 962 extreme
rays of the chamber cone:

```text
finite     895
ideal       67
hyperideal   0.
```

The finite rays have ten active facets.  Among the ideal rays, 6, 58, and 3
have respectively 11, 12, and 13 active facets.  Every one of the 27 walls
has exact extreme-ray incidence.  The zero hyperideal count independently
confirms finite volume, while the 67 ideal rays prove noncompactness.

## Exact boundary

The normalized Coxeter relation matrix is entrywise the already frozen
`D16-D4-C7-H10/27` diagram.  The containing known 14-wall reflection chamber
and index `2016` are unchanged.  Thus there is no new abstract polytope,
arithmetic commensurability class, or maximal reflection group here.

This period carrier is not equation-level and the imported chamber is **not
the full nef cone**.  In particular, its nineteen primitive square-minus-four
mirrors are divisor cuts, not negative-curve walls of a K3 surface.  No
effectivity or irreducibility is asserted for them.

The `S5` transposition is used to derive the lattice and wall system.  It is
**not a surface-automorphism** claim for the very-general carrier.  Likewise,
the 27 integral Picard reflections are not asserted to extend to Hodge
isometries or K3 automorphisms.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_ghess_h10_27_full_picard_k3_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_ghess_h10_27_full_picard_k3_carrier.py
```

Frozen payload SHA-256:

`f8e2f2e128be4d212370b8b99dc298e0725da8d0b8e560decbb03db32019d036`
