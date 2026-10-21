# TOOT: complete Jacobian–Kummer H12 five-vertex star census

Status: **PASS/FREEZE; all exact shape-census checks pass** (2026-08-08).

## Positive result

The complete properly colored `K1,4` finite-subtree construction over the
112-wall doubling stock of the 172-facet Jacobian–Kummer seed gives

```text
16,872 pairwise nonisomorphic ordinary Coxeter polytopes in H^12.
```

Every polytope is finite-volume noncompact, has reflection index five, and
has volume exactly five times the seed volume.  The complete histogram is:

| B-colored incident edges | facets `752+9B` | raw unrooted colorings | ordinary Coxeter types |
|---:|---:|---:|---:|
| 0 | 752 | 1,820 | 19 |
| 1 | 761 | 53,760 | 169 |
| 2 | 770 | 547,200 | 1,582 |
| 3 | 779 | 2,286,080 | 6,117 |
| 4 | 788 | 3,321,960 | 8,985 |
| **total** |  | **6,210,820** | **16,872** |

Here `A` means a color in `W61,...,W76`, and `B` means a color in
`W77,...,W172`.  The formula

```text
f(K1,4) = 172 + 145*4 + 9B = 752 + 9B
```

was checked on every reconstructed boundary; it was not extrapolated from
the five-path pilot.

## Exact coloring denominator

A proper `K1,4` coloring is an unordered four-subset of distinct incident
edge colors, so there are `C(112,4)=6,210,820` raw colorings.  Two independent
quotients agree:

- explicit orbit cover under the exact 384-element seed action: 16,872
  orbits with orbit-size sum 6,210,820;
- Burnside from the ten exact action cycle types: fixed-four-subset sum
  6,478,848, hence `6,478,848/384=16,872`.

The orbit-size histogram is

```text
12^1, 16^3, 24^9, 48^23, 64^10, 96^125,
128^12, 192^1111, 384^15578.
```

The separate all-shape denominator closes the combinatorial size of the
remaining five-vertex search:

| shape | seed-and-tree color orbits |
|---|---:|
| `P5` | 201,747 |
| `K1,4` | 16,872 |
| degree-3 fork | 199,362 |
| **total** | **417,981** |

This denominator is not itself a polytope census: different color orbits
must still be compared as ordinary Coxeter matrices.  The present result
closes precisely the whole `K1,4` shape.

## Direct Coxeter and lattice gates

For every orbit representative, the certifier starts from the exact
rank-13 Picard Gram matrix and uses the five chamber transformations

```text
1, r_a, r_b, r_c, r_d.
```

It verifies that multiplication by the selected reflections induces exactly
the marked star, with four internal mirrors and no extra adjacency.  All 860
transported local walls are deduplicated as primitive hyperplanes and
oriented by the sum of the five transported strict interiors.

Every one of the 16,872 rebuilt boundaries passes all of the following:

- exactly four internal mirrors;
- a strict timelike interior, with minimum wall pairing at least 5;
- every exterior normal primitive, spacelike, and reflective on the full
  integral lattice;
- all off-diagonal pairs in the ordinary label set
  `m2,m3,m4,m6,infinity`, where infinity combines parallel and
  ultraparallel pairs;
- the fixed boundary core `W1,...,W60`, whose Smith diagonal is `1^13` and
  therefore spans the full Picard lattice.

The ordinary relation inventory depends only on the facet count:

| facets | m2 | m3 | m4 | m6 | infinity |
|---:|---:|---:|---:|---:|---:|
| 752 | 13,016 | 608 | 3,816 | 0 | 264,936 |
| 761 | 13,158 | 624 | 3,847 | 0 | 271,551 |
| 770 | 13,300 | 640 | 3,878 | 0 | 278,247 |
| 779 | 13,442 | 656 | 3,909 | 0 | 285,024 |
| 788 | 13,584 | 672 | 3,940 | 0 | 291,882 |

Every row sums to `f(f-1)/2`.

The production engine uses signed 64-bit integer arrays only after explicit
pre-multiplication bounds below `2^63`.  The largest boundary coordinate
actually observed is 502, and the largest boundary pairing is 8,716.

## Pairwise nonisomorphism

Seed-color orbits were not assumed to be Coxeter types.  The certifier first
formed the exact multiset, over all facets, of the five ordinary-label degree
counts.  This gave 16,355 initial invariant classes.  The residue consisted
of 498 collision groups, 1,015 members in total, and maximum group size four.

For every residual member the complete ordinary Coxeter matrix was rebuilt
and passed through canonically numbered one-dimensional color refinement.
The stable invariant separated every member inside its initial exact class:

```text
16,872 color orbits = 16,872 pairwise nonisomorphic ordinary matrices.
```

Of the 1,015 refined matrices, 963 obtained a discrete canonical coloring.
For the remaining 52, distinct stable quotient inventories suffice to prove
nonisomorphism within their respective initial classes; no claim that the
nondiscrete quotient is a full canonical form is needed.

## Why finite volume and noncompactness are valid

The frozen seed certificate directly establishes finite-volume
noncompactness by both the literal Vinberg criterion and an independent
exact cone computation.  It also directly checks that all 112 selected walls
are doubling walls and that all `C(112,2)=6,216` pairs are parallel or
ultraparallel in Allcock's sense.

A proper finite subtree therefore gives a Coxeter fundamental polytope for
a reflection subgroup, and the number of seed chambers is its reflection
index and volume multiplier.  Every star here has five vertices, hence index
and volume multiplier five.  Finite index in the nonuniform seed lattice
preserves nonuniformity, so these are finite-volume **noncompact**, not
infinite-volume, examples.

## Divisor-lattice carrier

Every example lives in the same full Néron–Severi lattice

```text
U(2) + D5(-1)^2 + <-4>
```

of the pinned seven-dimensional K3 family.  The Smith-`1^13` core rules out a
proper-sublattice artifact.  This is a weak full-Picard divisor-reflection
carrier: it does not assert that each glued chamber is a full nef cone, that
all desired classes are effective, or that a surface automorphism realizes
the reflection subgroup.

Thus the result is directly relevant to the inverse-carrier program while
respecting the current theorem boundary: hyperbolicity and the integral
reflection arrangement are carried by the surface lattice; finite covolume
still enters through the separately certified seed plus Allcock's theorem.

## Relation to the five P5 pilots

The five earlier `P5` pilots have the same five facet counts and relation
histograms.  A separate audit rebuilds those paths from the raw rank-13 seed
and compares SHA-256 digests of their exact ordinary degree-profile
multisets against all 16,872 archived star profiles.  No collision occurs.
This gives a certificate-grade project count of 16,877 pairwise distinct
five-chamber examples.  The pure exact theorem of this report remains the
16,872 mutually nonisomorphic star matrices even without relying on hash
collision resistance.

## Claim boundary and lineage

This is a complete `K1,4` census, not a complete census of all five-vertex
trees.  The `P5` and degree-3 fork shapes remain open at 201,747 and 199,362
color orbits respectively.  Worldwide novelty remains open: Allcock already
proved exponential abundance in dimension 12, and the present diagrams have
not received a publication-complete literature denominator.  The verified
positive increment is the exact surface-derived star layer and its complete
ordinary-matrix classification.

## Reproduction

The focused tests validate every archived row structurally, rebuild one
representative for each `B=0,...,4`, replay the exact all-shape Burnside
calculation, and rebuild all five cross-shape pilots:

```bash
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_five_vertex_tree_denominators.py \
  tests/test_jacobian_kummer_h12_five_vertex_star_census.py \
  tests/test_jacobian_kummer_h12_five_vertex_star_vs_path_pilot_audit.py
```

Fresh result: `9 passed in 14.37s`.

The publication-strength archive replay recomputes all 16,872 boundaries
and all 1,015 refinement cases and therefore takes roughly 30 minutes on the
current machine:

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py --check-archive
```

Artifacts:

- [`certify_jacobian_kummer_h12_five_vertex_star_census.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_star_census.py)
- [`jacobian_kummer_h12_five_vertex_star_census.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_star_census.json)
- [`test_jacobian_kummer_h12_five_vertex_star_census.py`](../tests/test_jacobian_kummer_h12_five_vertex_star_census.py)
- [`certify_jacobian_kummer_h12_five_vertex_tree_denominators.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_tree_denominators.py)
- [`jacobian_kummer_h12_five_vertex_tree_denominators.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json)
- [`audit_jacobian_kummer_h12_five_vertex_star_vs_path_pilot.py`](../scripts/audit_jacobian_kummer_h12_five_vertex_star_vs_path_pilot.py)
- [`jacobian_kummer_h12_five_vertex_star_vs_path_pilot_audit.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_star_vs_path_pilot_audit.json)
