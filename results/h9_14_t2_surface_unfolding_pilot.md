# Surface-free unfolding theorem for the T2/C01 H9/14 datum

Status: **exact IC0 uniqueness theorem, weak K3 realization, and marked bridge
PASS** (2026-08-08).

This is not a new hyperbolic polytope and not yet a class-wide theorem saying
that every suitable Coxeter polytope comes from an algebraic surface.  It is a
closed pilot for that proposed reverse direction: beginning only with a
normalized Coxeter Gram, a finite compiler recovers the unique compatible
fixed lattice and an explicit K3-lattice unfolding.

## Surface-free input and finite theorem

The logical input is only the 14-by-14 facet Gram of `A18_4_0_C01` in
`alexeev_k3_h9_14_factory.json`.  Its canonical string-JSON SHA-256 is

```text
6979b16de7583c20277077a142df863e9f9a7209c443e9252505234b15d66c56
```

It has rank ten, inertia `(1,9,4)`, and normalized relation census

```text
75*m2 + 8*m3 + 8*m4.
```

Contracting the `m3` components and solving every oriented crystallographic
ratio equation on the eight `m4` edges gives exactly 18 primitive minimal even
normal scalings.  Their six diagram-automorphism orbits have sizes

```text
1, 1, 2, 4, 4, 6,
```

and the full labelled diagram automorphism group has order 48.  Sixteen
scalings contain square-minus-eight normals and cannot be built from either a
fixed `(-2)` root or an orthogonal exchanged pair of `(-2)` roots.  Exactly two
one-/two-sheet lanes remain:

```text
archived orientation   4*(-2) + 10*(-4)
dual orientation      10*(-2) +  4*(-4).
```

The verifier then enumerates every reflection-stable quadratically isotropic
subgroup of the root-span discriminant module.  This is a complete marked
overlattice enumeration, not a sample:

```text
                         archived       dual
stable even candidates       16            5
glue dimensions          0:1,1:7,2:7,3:1  0:1,1:3,2:1.
```

Three logically separate gates are applied to each candidate:

1. distinct long-wall discriminant classes, needed for an acute rootless
   two-sheet lift;
2. an exact embedding of the long-class quadratic space into the 120
   nonsingular classes of `A_E8(-2)` represented by pairs of opposite E8
   roots;
3. discriminant length at most four after minimal rank-18 graph glue, needed
   for a primitive rank-four K3 complement.

The archived root span has a seven-dimensional long-class space with
three-dimensional polar radical.  A seven-subspace of a nondegenerate
eight-dimensional symplectic space has radical dimension at most one, so it
cannot embed in `A_E8(-2)`.  The dual root span does embed, but its minimally
glued rank-18 discriminant length is eight, too large for a rank-four K3
complement.  Repeated-class candidates fail independently.

Exactly one of all 21 marked fixed-lattice candidates survives.  It is the
index-two archived extension with glue vector `0x28a`, determinant `-256`,
Smith form `1^2 2^8`, and invariants

```text
(r,a,delta) = (10,8,0).
```

The indefinite even 2-elementary classification therefore forces

```text
F = U + E8(-2).
```

This is the finite surface-free theorem: under the minimal Nikulin
one-/two-sheet model and the weak K3 length gate, the divisor-lattice candidate
is not guessed from the known surface; it is uniquely reconstructed from the
Coxeter Gram.

## Explicit rank-18 and K3-lattice realization

The ten long classes span dimension six.  The compiler embeds a canonical
basis as E8 discriminant classes

```text
01, 02, 04, 10, 20, 40
```

and gives exact E8 roots in all ten required classes.  Gluing the corresponding
half-sums into `F + E8(-2)` has index 64 and produces an even lattice `N` with

```text
signature        (1,17)
determinant         -16
Smith             1^14 2^4
(r,a,delta)       (18,4,0).
```

The sheet involution is integral.  Its fixed and coinvariant embeddings are
primitive, with fixed lattice `U+E8(-2)` and coinvariant lattice `E8(-2)`.
The compiler constructs 24 specified `(-2)` wall roots—not a claim that these
are every `(-2)` root of `N`—with off-diagonal pairings

```text
244*0 + 32*1.
```

Four roots are fixed and ten pairs are exchanged.  Their orbit sums reproduce
the input 14-by-14 Gram entry for entry.

Finally, an exact discriminant anti-isometry glues `N` to `U(2)+U(2)`.  The
result is even unimodular of signature `(3,19)` and determinant `-1`; hence it
is the K3 lattice.  The involution extends integrally with characteristic
polynomial

```text
(lambda-1)^14 (lambda+1)^8.
```

This closes the weak lattice embedding and involution-extension gates.  It does
not by itself prove effectivity, irreducibility, extremality, or completeness
of the displayed divisor walls on a surface.

## Exact bridge to the known C01 surface datum

Only after the surface-free theorem is complete does the certifier read the
pinned Alexeev 48-wall ancillary matrix and the existing C01 strong-promotion
archive.  The 14 active source-root orbits are

```text
[1] [3,9] [4,8] [5,7] [6] [11,17] [12,16]
[13,15] [14] [18,24] [19,25] [22,32] [23,31] [30].
```

All `2^10` sheet assignments are checked.  Exactly four give integral marked
isometries.  The canonical unswapped witness has determinant `-1` and matrix
SHA-256

```text
15f683beded0f7366fe8be371298e5c17bf42238e36c2bec7833398acf6d6e85.
```

It maps all 24 compiled roots to the archived source roots, gives exact Gram
congruence, and intertwines the compiled involution with the transpose of the
archived row action.  The action hash agrees with the pinned strong promotion,
so there is no remaining coordinate gap to that already certified
datum-specific K3/nef realization.

The top-level factory geometry—not its older pre-volume section status—records
the known folded polytope as finite-volume noncompact.  Literal Vinberg gives
90 finite and 20 ideal endpoints, with zero bad elliptic edges and zero bad
parabolic neighborhoods.

## What this contributes to the surface-correspondence program

The pilot isolates a plausible finite front end for a reverse carrier theorem:

```text
Coxeter Gram
  -> crystallographic normal scalings
  -> stable even divisor-lattice overclosures
  -> rootless coinvariant sheet classes
  -> primitive K3-lattice embedding/action
  -> geometric wall realization
  -> literal Vinberg polytope gate.
```

The first four arrows are exact and finite here.  The existing C01 promotion
supplies the last two for this one marked datum.  A class-wide theorem still
needs algebraic-geometric predicates ensuring effectivity, the correct
irreducible or reducible wall interpretation, a common ample orientation,
complete-root/hidden-wall separation, and exact extraction of the full nef
section.  On the polytope side it also needs full Vinberg edge-endpoint and
parabolic-neighborhood checks after redundancy removal.

Thus the useful conclusion is a criterion template, not universality yet:
surface moduli can replace unrestricted polytope fishing only after the
geometric wall-realization and reverse-extraction arrows are made automatic.
This example shows that the lattice portion of that idea can be both finite and
rigid enough to force the expected surface lattice.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_t2_surface_unfolding_pilot.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_t2_surface_unfolding_pilot.py
```

Artifacts:

- [`certify_h9_14_t2_surface_unfolding_pilot.py`](../scripts/certify_h9_14_t2_surface_unfolding_pilot.py)
- [`h9_14_t2_surface_unfolding_pilot.json`](../data/enumerations/h9_14_t2_surface_unfolding_pilot.json)
- [`test_h9_14_t2_surface_unfolding_pilot.py`](../tests/test_h9_14_t2_surface_unfolding_pilot.py)
