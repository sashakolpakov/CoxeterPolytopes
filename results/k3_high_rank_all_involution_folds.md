# Exhaustive involution folding of the rank-18 and rank-19 K3 nef chambers

## Result

The complete exact census produces **no new hyperbolic Coxeter polytope**.
It recovers precisely the two already archived K3 quotient-nef folds:

| ambient K3 wall graph | full weighted graph group | nonidentity involutions | involution classes | folded dimension | outcome |
|---|---:|---:|---:|---:|---|
| rank 18, 19 walls | `C2` | 1 | 1 | \(\mathbb H^9\) | exactly the published 10-facet fold |
| rank 19, 24 walls | `S3` | 3 | 1 | \(\mathbb H^{10}\) | the published 13-facet fold and its two `S3` conjugates |

Thus the strict target \(\dim\mathbb H>9\) contains three labeled actions,
but only one geometric conjugacy class.  It is the already published
rank-19 transposition fold.  No `TOOT` is warranted.

This is nevertheless a completion result for the specified search space:
there is no untested weighted-wall-graph involution of either complete
configuration.

## Exhaustive graph census

The search reads the complete integral curve Gram matrices in
`k3_high_rank_nef_chambers.json`.  It enumerates all bijections \(p\) such
that

\[
G_{p(i),p(j)}=G_{i,j}
\]

by exact backtracking on the complete weighted matrix.  Vertex signatures
are used only to prune the search; every pair relation is checked before a
permutation is accepted.  Closure, inverses, and element orders are then
checked exactly.

The rank-18 group has order two and order histogram

\[
1^1,2^1.
\]

Its unique nonidentity action is exactly the archived involution

\[
(1\ 18)(2\ 17)(3\ 15)(4\ 16)(5\ 14)(6\ 13)
(7\ 12)(8\ 11)(10\ 19),
\]

with wall 9 fixed.

The rank-19 group has order six and order histogram

\[
1^1,2^3,3^2,
\]

so it is `S3`.  Its three involutions form a single conjugacy class.  One is
the archived transposition; the other two are obtained by rotating the
threefold-symmetric wall graph.  Exact conjugators and the induced orbit
permutations are frozen in the certificate.  After transport, every
orbit-sum Gram entry agrees—not merely its determinant or spectrum.

## Direct Coxeter and finite-volume audits

Every nonidentity involution was audited directly, even though conjugacy
would suffice for the two extra rank-19 actions.  For every orbit
\(\mathcal O\), the normal is

\[
R_{\mathcal O}=\sum_{C\in\mathcal O}C.
\]

All negative-square two-curve orbits consist of disjoint curves.  Hence the
two ambient root reflections commute and restrict to the integral reflection
in \(R_{\mathcal O}\).  In each rank-19 action, one of the three pairs among
walls 22, 23, and 24 has curve intersection 6 and orbit-normal square 8.  Its
hyperplane misses hyperbolic space, so that inequality is redundant.  The
three choices are exchanged by `S3`.

For every retained orbit wall the computation checks:

1. Lorentzian inertia of the fixed span;
2. an exact Coxeter label for every facet pair;
3. integrality of every folded reflection;
4. irredundancy by spanning the equality face with exact extreme rays; and
5. Vinberg--CoxIter Proposition 1.10 using the complete folded wall list.

The endpoint results are:

| fold | facets | Coxeter edges | finite vertices | ideal vertices | elliptic edge diagrams | endpoint count |
|---|---:|---|---:|---:|---:|---|
| \(\mathbb H^9\) | 10 | \(8\times3,1\times4\) | 8 | 2 | 45 | exactly 2 for all 45 |
| \(\mathbb H^{10}\) | 13 | \(9\times3,2\times4,2\times\infty\) | 34 | 4 | 208 | exactly 2 for all 208 |

All four labeled folds are therefore finite-volume and noncompact.  The
rank-18 fold fails only the requested strict dimension gate because its
dimension is 9, not greater than 9.

## Algebraic realization versus graph symmetry

A graph involution by itself is not an algebraic automorphism of a K3
surface.  The realization gate here is separate:

- the unique rank-18 graph involution is exactly the surface involution
  realized by Kondō;
- Kondō realizes the full rank-19 `S3` graph action, so all three
  transpositions are surface automorphisms, not merely combinatorial
  symmetries.

Consequently the fixed nef sections really are quotient-surface divisor
cones.  But this same realized action also settles novelty: the three
rank-19 folds are algebraically conjugate and define isometric copies of the
same published Coxeter polytope.

## Search boundary

The two selected complete K3 configurations are now exhausted under
involution folding.  A new example must change at least one input:

- use a different complete K3 curve configuration;
- use a finite subgroup or correspondence not conjugate to these
  involutions; or
- move beyond ordinary fixed-space orbit sums while retaining an independent
  algebraic-realization theorem.

The reproducible certifier is
`scripts/certify_k3_high_rank_all_involution_folds.py`; its frozen output is
`data/candidates/k3_high_rank_all_involution_folds.json`, with regression
tests in `tests/test_k3_high_rank_all_involution_folds.py`.
