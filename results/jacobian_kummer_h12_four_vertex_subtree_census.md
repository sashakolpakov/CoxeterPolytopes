# TOOT: 2,644 four-chamber Jacobian--Kummer polytopes in `H12`

Status: **exact positive census; production and independent hostile
certificates green** (2026-08-08).  Worldwide priority remains open pending a
publication-complete lineage audit.

## Complete positive theorem

Let `Q=JK-H12/172`, and use its certified pairwise disjoint doubling-wall
stock `W61,...,W172`.  Up to tree automorphism and the exact 384-element seed
action, every four-vertex subtree belongs to exactly one of

```text
1,964 colored paths P4,
  680 colored stars K1,3.
```

Exact reconstruction proves that the resulting **2,644 finite-volume
noncompact hyperbolic Coxeter polytopes in `H12` are pairwise
nonisomorphic**, even at the level of ordinary Coxeter matrices where
parallel and ultraparallel pairs are both assigned the single entry
`infinity`.

The facet theorem is unusually simple.  Split the selected generators into

```text
A = W61,...,W76      (16 colors),
B = W77,...,W172     (96 colors).
```

If `B(T)` counts the three edge-color occurrences of `T` that lie in the
second class, with multiplicity, then

\[
             f(T)=607+9B(T).
\]

Consequently the complete type census is:

| facets | `P4` types | `K1,3` types | all types |
|---:|---:|---:|---:|
| 607 | 17 | 6 | **23** |
| 616 | 127 | 43 | **170** |
| 625 | 626 | 218 | **844** |
| 634 | 1,194 | 413 | **1,607** |
| total | **1,964** | **680** | **2,644** |

The four facet counts `607,616,625,634` are new to the project ledger.
Every polytope consists of four seed chambers, so its reflection index and
volume multiplier are exactly four.

## Correct coloring denominator

A proper edge coloring of a subtree of the free-product Cayley tree has
distinct colors on adjacent edges at a vertex.

For a path, write the colors in order as `(a,b,c)`.  The exact conditions are

```text
a != b, b != c;     a = c is allowed;
(a,b,c) ~ (c,b,a)  by path reversal.
```

There are 1,379,952 admissible ordered words and 696,192 colored paths after
reversal.  The diagonal seed action together with reversal has 1,964 orbits;
47 have repeated outer color `a=c`.  Treating the word as reversible is
geometrically exact: left translation by the appropriate terminal word
reverses the four chamber path.

For a star, the three colors incident to the central chamber must be
distinct and are unordered.  Thus there are

\[
 \binom{112}{3}=227,920
\]

colored stars, split into 680 seed-action orbits.  In total the computation
covers 924,112 unrooted colored four-vertex subtrees before quotienting by
the seed action.  No adjacent repeated color is admitted: that would
backtrack across the same wall and fail to give four distinct vertices.

Weighted by these raw unrooted colorings, the facet histogram is

| facets | raw colored trees |
|---:|---:|
| 607 | 2,480 |
| 616 | 47,616 |
| 625 | 293,376 |
| 634 | 580,640 |

## Ordinary-Coxeter nonisomorphism theorem

The certifier forms the full ordinary Coxeter matrix of every representative,
using the five symbols

```text
m2, m3, m4, m6, infinity.
```

It first canonically sorts the multiset of the five exact degree counts at
all facets.  This elementary invariant already gives 2,535 classes.  Only
102 buckets, containing 211 representatives and of size at most four,
remain.

Those 211 matrices receive canonical one-dimensional color refinement.  In
195 cases the coloring becomes discrete, yielding a literal canonical
adjacency word.  The remaining 16 have nontrivial final cells, but their
canonical stable quotient invariants are already different from every other
member of their initial degree bucket.  Hence every residual pair is also
nonisomorphic.  Combining the exact initial profiles with the refinement
words gives 2,644 distinct invariants.

Here “canonical” includes the numeric color IDs, not merely the unnumbered
partition.  The proof is a short induction.  Initially, the five-component
label-degree vector of a vertex is invariant under simultaneous row/column
permutation, and lexicographically sorting the distinct vectors assigns the
same IDs in every presentation.  Inductively, if the old IDs are canonical,
then the signature consisting of the old ID and all `(edge label, neighbor
color)` counts is also permutation-invariant; lexicographically sorting the
distinct signatures canonically assigns the next IDs.  At stability, sorting
vertices by IDs therefore gives a canonical adjacency word in the discrete
case.  In the nondiscrete case, the sorted inventory of `(edge label, smaller
cell ID, larger cell ID, multiplicity)`, together with canonical cell sizes
and round counts, is an isomorphism invariant (not asserted to be a complete
canonical form).  The executable certificate permutes one actual discrete
and one actual nondiscrete residual matrix in two nontrivial ways and checks
both transported color IDs and the exact final byte word; the unit tests do
the same on small independent fixtures.

This proof is stronger than a metric dotted-edge comparison: it forgets all
ultraparallel distances and even the distinction between parallel and
ultraparallel walls.  It also separates every path type from every star type.
The script compares the underlying exact byte words before freezing their
SHA-256 digests.

For reference, the complete ordinary relation inventories depend only on
the facet count:

| facets | `m2` | `m3` | `m4` | `m6` | `infinity` |
|---:|---:|---:|---:|---:|---:|
| 607 | 10,467 | 492 | 3,062 | 0 | 169,900 |
| 616 | 10,609 | 508 | 3,093 | 0 | 175,210 |
| 625 | 10,751 | 524 | 3,124 | 0 | 180,601 |
| 634 | 10,893 | 540 | 3,155 | 0 | 186,073 |

These coarse global inventories do not prove the type theorem—the refined
vertex invariants do.

## Boundary, volume, and carrier gates

The exact chamber sets are

```text
P4:    {1, R_a, R_a R_b, R_a R_b R_c},
K1,3:  {1, R_a, R_b, R_c}.
```

For each of all 2,644 representatives the verifier constructs the integral
reflection matrices, deletes exactly three internal mirrors, transports all
172 seed roots, and deduplicates primitive hyperplanes.  It then checks:

1. the summed chamber interior is timelike and strictly satisfies every
   boundary inequality;
2. every boundary root is primitive, negative, and integrally reflective;
3. every boundary pair has one of the five ordinary Coxeter relations; and
4. the unchanged exterior core `W1,...,W60` remains on the boundary and has
   Smith diagonal `1^13`.

The largest coordinate encountered is 128,830 and the largest absolute
boundary pairing is 143,888.  Before boundary reconstruction, let `M=336`
be the largest reflection-matrix entry, `X=30` the largest seed-root entry,
and `I=17` the largest seed-interior entry.  For rank `r=13`, the verifier
uses the elementary entry bound

\[
 |(A_1\cdots A_k)_{ij}|\le r^{k-1}M^k.
\]

In particular it checks, before the corresponding products, the bounds

```text
three-reflection chamber transform       r^2 M^3       =       6,410,686,464
four-reflection internal-mirror product  r^3 M^4       =  28,001,878,474,752
transformed seed root                    r^3 M^3 X     =   2,500,167,720,960
transformed seed interior                r^3 M^3 I     =   1,416,761,708,544
sum of four transformed interiors        4 r^3 M^3 I   =   5,667,046,834,176
```

All are strictly below `2^63-1 = 9,223,372,036,854,775,807`.  Subsequent
boundary-pairing and Gram products retain their existing instance-specific
pre-multiplication bounds.  Thus neither the four-reflection mirror test nor
the transformed-interior sum relies on unchecked `int64` arithmetic.

Allcock's finite-subtree theorem now gives a Coxeter fundamental polytope of
reflection index and volume multiplier four.  Since the seed lattice is
nonuniform, every finite-index descendant is nonuniform and therefore
noncompact.  The fixed Smith core proves that every boundary still spans

\[
 U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle
\]

with index one.  Thus all 2,644 types retain the same seven-dimensional
very-general full-Picard K3 family as a **weak divisor-reflection carrier**.
No nef-cone, effectivity, Hodge, or surface-automorphism claim is made.

## Lineage boundary

Allcock's theorem and exponential `H12` abundance are classical.  The
project-level gain is the complete order-four classification for this exact
surface-derived 112-wall factory, four new explicit facet counts, 2,644
ordinary-Coxeter type certificates, and their common full-Picard carrier.

The pinned CoxIter lanes contain 1,214 ordinary files and 2,959 index-two
files; their `H12` headers are respectively `13,14,16,18` and
`18,20,25,35`.  None of `607,616,625,634` occurs.  This is bounded
same-header evidence only, not a publication-complete novelty theorem.
Worldwide priority for these exact diagrams remains open.

The independent hostile replay reconstructs the raw Dolgachev--Shimada seed,
the 384-element action, both Burnside/BFS orbit counts, every boundary, and
all carrier and Allcock gates.  Degree profiles already distinguish 2,535
types; for the remaining 102 collision groups, the multiset of labelled
two-walk entries separates all 211 matrices and survives explicit vertex
permutations.  The archive replay passes and the focused hostile suite reports
`4 passed in 368.99s`; see the
[`hostile audit`](jacobian_kummer_h12_four_vertex_subtree_hostile_audit.md).

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_four_vertex_subtree_census.py
```

Artifacts:

- [`jacobian_kummer_h12_four_vertex_subtree_census.json`](../data/enumerations/jacobian_kummer_h12_four_vertex_subtree_census.json)
- [`certify_jacobian_kummer_h12_four_vertex_subtree_census.py`](../scripts/certify_jacobian_kummer_h12_four_vertex_subtree_census.py)
- [`test_jacobian_kummer_h12_four_vertex_subtree_census.py`](../tests/test_jacobian_kummer_h12_four_vertex_subtree_census.py)

Frozen JSON SHA-256:

`2591d5f078421b5b4902045d37fa57e33548e2450f090f5ea8513c564c1a4b24`

Hostile archive SHA-256:

`83efbd60b5ded8886d1cca201a323f6a3ae4eacdd66fe110f41adb028796f7d8`
