# Hostile audit: compact complete K3 root return

**Audit status (2026-08-09): PASS after certificate and exposition repairs.**

Audited objects:

- [`compact_complete_k3_root_return_classification.md`](compact_complete_k3_root_return_classification.md);
- [`certify_compact_k3_root_return.py`](../scripts/certify_compact_k3_root_return.py); and
- [`test_compact_k3_root_return.py`](../tests/test_compact_k3_root_return.py).

## Quantifier and theorem-chain verdict

The input class is genuinely surface-free: an even lattice `N` of signature
`(1,n)`, `n>=2`, together with a compact chamber of its **complete**
square-minus-two mirror arrangement.  No K3 embedding, surface, Weyl vector,
finite automorphism group, or dimension bound is assumed.

The return chain is valid for exactly that class:

1. A complete root chamber is a fundamental chamber for `W_2(N)`.  If it is
   compact, `W_2(N)` is a cocompact lattice.  Since it is contained in the
   arithmetic lattice `O^+(N)`, it has finite index.
2. Compactness of the resulting arithmetic quotient and Godement's criterion
   make `N_Q` anisotropic.  Meyer then gives `rank(N)<=4`, hence `dim(P)<=3`.
3. The stable primitive-embedding range is therefore automatic: `N` has one
   primitive embedding orbit in the K3 lattice under `O(Lambda_K3)`.  The
   type-IV period domain minus the countable Noether--Lefschetz hyperplanes is
   nonempty and has dimension `20-rank(N)`, giving `NS(X)=N` exactly.
4. For an exact-Picard K3, Riemann--Roch and K3 Weyl theory identify its nef
   cone with a chamber of the complete `(-2)` arrangement.  Simple
   transitivity of `W_2(N)` on chambers gives a unique re-framing, after the
   marking has been oriented by the input positive-cone component.  Each root
   reflection extends integrally to
   `Lambda_K3` by the same formula and fixes `N^perp`, so this re-framing is
   compatible with the fixed period/embedding package.  Thus the extraction
   is literally `P -> X -> same P`, with no hidden root cut.
5. The cited Nikulin rank-three and Vinberg rank-four compact exhaustions,
   restated with full curve Grams by Roulleau, leave exactly
   `S1,...,S6,L24,L27`.  Since all chambers of a fixed full root arrangement
   are `W_2(N)`-conjugate, an exhausted lattice does not split into extra
   marked chamber types.

The bridge to the surface classification is also sound: finite index of
`W_2(N)` makes the integral nef-chamber stabilizer finite, the kernel of the
K3 automorphism action on `NS` is finite, and compactness excludes nef
isotropic rays and hence elliptic fibrations.

## Repairs made by this audit

The original matrix replay compared the determinant of a selected wall basis
with `det(NS)`.  That is not enough to prove integral equality of lattices.
The certifier now freezes and checks, for every row, an explicit unimodular
change of basis `U` satisfying

```text
U^T * Gram(NS) * U = Gram(selected facet walls).
```

It reconstructs every wall integrally in that basis and rechecks the full
wall Gram.  In particular, the `L27` identification is no longer implicit.
It also canonicalizes every wall Gram under simultaneous row/column
permutation and verifies that the eight metric chambers remain pairwise
distinct after forgetting the frozen wall order.

The cone replay was also strengthened.  It now checks that every extreme ray
is timelike in the same positive-cone component as the interior point and
that finite facet pairs agree exactly with vertex incidences.  This directly
certifies the stated polygon/polyhedron `f`-vectors, including
`L24=(8,12,6)` and `L27=(12,18,8)`.

The proof text now states the `O(Lambda_K3)` scope of embedding uniqueness,
the simple-transitivity/unique-reframing argument, and the finite-
automorphism bridge instead of leaving those steps implicit.

## Classification and Weyl-vector checks

The pinned Roulleau arXiv-v5 source contains both literal exhaustion
statements: six compact rank-three types and only two compact rank-four
types.  Its lattice/curve data agree with all eight frozen matrices.  The
source archive and main TeX hashes replay exactly.

The division of labor is explicit: the finite matrix code proves
compactness and irredundancy **after** accepting the facet lists.  Their
full-root completeness is a cited classification input, not something
inferred from the word "finite" or from the Gram alone.  Roulleau's complete
irreducible `(-2)`-curve lists give all nef facets; mirrors of reducible
effective roots may pass through lower-dimensional faces, but their
inequalities are redundant and do not cut the returned chamber.

The rational Weyl-vector test is exact linear algebra:

```text
rank(G) = rank([G | 1]).
```

It gives precisely five positive rows `S1,S2,S3,S5,L24` and three negative
rows `S4,S6,L27`.  For the five solutions, the independently reconstructed
vector squares are respectively `1/2,1/4,1,2,2`, so all five lie in the
positive cone rather than merely solving a formal singular-Gram equation.
The three negative rows are new strata relative to the project's
earlier Weyl-vector theorem, not new worldwide polytopes; the main report
states that lineage correctly.

## Replayed tests

```text
python3 scripts/certify_compact_k3_root_return.py --check-source
  source hashes and both exhaustion statements: PASS

python3 -m pytest -q tests/test_compact_k3_root_return.py
  3 passed
```

## Scope boundary

This audit does **not** turn the result into a theorem about arbitrary
compact arithmetic Coxeter polytopes.  It proves the exhaustive inverse
return only for complete integral norm-two root chambers over `Q`, in
hyperbolic dimension at least two.  Mixed root norms, totally real
eigencones, invariant folds, and Bugaenko's compact `H8` example remain
outside the quantifier.  The Nikulin--Vinberg exhaustion is cited and pinned,
not recomputed by the certifier.
