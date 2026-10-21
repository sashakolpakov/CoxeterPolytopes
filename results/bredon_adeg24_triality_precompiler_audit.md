# Audit of the 24-row Bredon ADEG triality precompiler

## Status

The census in
[`certify_bredon_adeg24_triality_census.py`](../scripts/certify_bredon_adeg24_triality_census.py)
is sound for its declared purpose.  It is an exact **integral-triality and
K3-lattice precompiler**, not a theorem that any unfolded diagram is the nef
chamber of a K3 surface.  In particular, it certifies **zero** same-polytope
or full-nef returns.

The source is Naomi Bredon, *On ADEG-polyhedra in hyperbolic spaces*,
[arXiv:2507.05153v1](https://arxiv.org/abs/2507.05153), Theorem 4.2 and
Table 2.  The arXiv record still has only version 1 as of the audit date.

## 1. Source and normalized-Gram audit

Table 2 was checked diagram-by-diagram against the 24 in-code records.  The
dimension distribution is

| dimension | 2 | 3 | 5 | 6 | 7 | 9 | 11 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| diagrams | 4 | 7 | 4 | 1 | 5 | 2 | 1 |

and the facet counts and placements of every label-six edge agree with the
source, including both 5-dimensional octagons and the 14-facet polyhedron
`P_star` in dimension nine.

For every source diagram `Sigma`, the script independently forms twice the
normalized Coxeter Gram matrix

```text
G_ii = 2,    G_ij = -1 on a 3-edge,
             G_ij = -sqrt(3) on a 6-edge,
             G_ij = 0 otherwise.
```

Exact congruence elimination in `Q(sqrt(3))` gives, for all 24 rows,

```text
rank(G) = n + 1,
inertia(G) = (n positive, 1 negative, facets - n - 1 zero).
```

Thus the table transcription passes the expected hyperbolic rank and inertia
test in every row.  This audit uses Bredon's finite-volume classification as
input; it does not redo Vinberg's finite-volume proof for the 24 polytopes.

## 2. Exact integral-triality criterion

Give each folded wall square either `2` or `6`.  The Coxeter inner product is

```text
-sqrt(q_i q_j) cos(pi/m_ij).
```

It is integral on a 3-edge exactly when the two squares agree, and it is
integral on a 6-edge exactly when they differ.  Consequently an assignment
exists if and only if, after contracting all 3-edge components, the 6-edge
quotient is loop-free and bipartite.  This proves both necessity and
sufficiency of the implemented combinatorial gate.

An independent union-find and bipartite-colouring replay gives:

| outcome | diagrams | orientations |
|---|---:|---:|
| passes | **21** | **42** |
| fails | **3** | 0 |

The failures are:

- `ADEG-N2-02-triangle-633`: its 6-edge becomes a loop after 3-edge
  contraction;
- `ADEG-N2-04-triangle-666`: its quotient is an odd 6-edge cycle;
- `ADEG-N3-06-square-one-6`: its 6-edge becomes a loop after contraction.

Each passing source graph is connected, so its connected quotient has exactly
the two complementary long/short orientations.  This accounts for all
`21 * 2 = 42` orientation packages, without quotient automorphisms being used
to identify the two choices.

## 3. Integral folding, unfolding, and quotient lattices

For every one of the 42 orientations, exact replay verifies the integral
folded Gram:

- diagonal `2` on short walls and `6` on long walls;
- inner product `-1` on short-short 3-edges;
- inner product `-3` on long-long 3-edges and long-short 6-edges.

Equivalently, if `D_ii` is `1` on a short wall and `sqrt(3)` on a long wall,
then the integral matrix is exactly `D G D`.

Every square-six wall is then replaced by three mutually orthogonal
square-two roots.  A long-long 3-edge becomes three phase-matched 3-edges,
and a long-short 6-edge joins the short root to all three roots in the long
orbit.  In all 42 cases the orbit sums recover the integral folded Gram.

The certifier uses HNF coordinates to compute the integral lattice generated
by the unfolded roots modulo the integral radical.  As an independent audit,
the nonzero Smith factors of the full unfolded Gram were recomputed directly
for all 42 packages.  They agree exactly with every archived Smith diagonal,
discriminant, and discriminant length.

Terminology must remain precise here.  This is the **root-generated quotient
lattice**.  The internal helper name `saturated_wall_lattice` does not mean
that a primitive closure inside an already specified K3 Neron--Severi lattice
has been constructed.  The archive's public name `root_generated_lattice` is
the correct interpretation.

## 4. The 31 stable-range K3 carrier candidates

The unfolded roots have square `+2`.  For the K3 test one applies Nikulin's
primitive-embedding criterion to the **sign-reversed** even lattice.  If the
unfolded quotient has rank `r`, Lorentzian inertia `(r-1,1)`, and discriminant
length `ell`, the implemented sufficient gate is

```text
r <= 19,    22 - r >= ell + 2.
```

The exact outcome is:

| source dimension | 2 | 3 | 5 | 6 | 7 | 9 | 11 | total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| stable packages | 4 | 12 | 7 | 1 | 5 | 1 | 1 | **31** |

Of the other eleven orientation packages, nine unfolded forms are not
Lorentzian.  The remaining two are merely outside this easy theorem: one has
rank 20, while one `P_star` orientation has rank 18 and discriminant length
four, so `22-r = 4 < 6 = ell+2`.  The precompiler does not claim that either
of these two lattices is non-embeddable.

For a stable package the lattice theorem supplies a primitive embedding of
the sign-reversed root lattice into the K3 lattice.  It does **not** by itself
prove any of the following:

1. that the marked roots are simultaneously effective irreducible curves;
2. that their half-spaces contain an ample class;
3. that there are no additional effective `(-2)` roots or hidden nef walls;
4. that folding the resulting divisor chamber returns Bredon's original
   polytope with the same facets; or
5. that the displayed chamber is the complete nef cone.

Accordingly the archive correctly records `31` stable-range carrier
**candidates** and `0` strong same-`P` K3 returns.

## 5. Reproducibility

- archive:
  [`bredon_adeg24_triality_census.json`](../data/enumerations/bredon_adeg24_triality_census.json)
- canonical payload SHA-256:
  `00ce6afa559a4fd84217b6e38cfa57867747da19964bbe0c1e03e94b1fc14722`
- complete live `--check-archive` replay: passed;
- focused independent regression:
  [`test_bredon_adeg24_triality_census_audit.py`](../tests/test_bredon_adeg24_triality_census_audit.py),
  **5 passed**.

The regression freezes the source transcription, repeats all 24 exact
`Q(sqrt(3))` rank/inertia checks, independently recomputes the 21/24 quotient
gate and all 42 orientations, reconstructs every folded and unfolded Gram,
and obtains the quotient discriminants through Smith form rather than the
certifier's HNF path.
