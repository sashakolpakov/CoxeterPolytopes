# TOOT: Jacobian–Kummer H12 five-chamber positive pilot

Status: **PASS/FREEZE; all exact checks pass** (2026-08-08).

## Positive result

The 112-wall doubling stock of the 172-facet Jacobian–Kummer seed produces
the following five explicit finite-volume noncompact Coxeter polytopes in
`H^12`:

| B edges | properly colored P5 edge walls | facets | volume |
|---:|---|---:|---|
| 0 | `W61,W62,W63,W64` | 752 | `5 Vol(P_seed)` |
| 1 | `W61,W62,W63,W77` | 761 | `5 Vol(P_seed)` |
| 2 | `W61,W62,W77,W78` | 770 | `5 Vol(P_seed)` |
| 3 | `W61,W77,W78,W79` | 779 | `5 Vol(P_seed)` |
| 4 | `W77,W78,W79,W80` | 788 | `5 Vol(P_seed)` |

Thus the proposed pilot formula is exact without correction:

```text
f(P5) = 172 + 145(5-1) + 9B = 752 + 9B
```

on these five paths.  Since the ordinary Coxeter matrices have five distinct
orders, they are pairwise nonisomorphic.

This is intentionally a minimal positive pilot.  It is not a census of all
properly colored five-vertex trees and does not assert that the formula holds
for every five-vertex coloring or tree shape.

## Exact construction

The certificate starts from the exact rank-13 seed lattice and its 172
oriented roots.  It does not read a four-tree boundary or extrapolate one.
For edge-color word `(a,b,c,d)`, it constructs the five chamber transforms

```text
1, r_a, r_a r_b, r_a r_b r_c, r_a r_b r_c r_d.
```

All four colors in each displayed path are distinct.  The transforms are
five distinct integral lattice isometries, and multiplication by the four
selected reflections induces exactly the marked `P5`, with no extra chamber
adjacency or backtracking.

The verifier transports all 172 seed walls through all five transforms.
It identifies the four internal mirrors in two independent ways:

- adjacency in the marked chamber tree;
- cancellation of opposite local wall orientations.

The two four-element sets agree in every row.  All other hyperplanes are
deduplicated and oriented by the sum of the five transported strict
interiors.  The resulting direct facet counts are exactly those in the
table.

## Coxeter and lattice gates

Every exterior normal is primitive, spacelike, and defines an integral
reflection of the full lattice.  The summed interior is strictly timelike
and has positive pairing with every exterior wall; the minimum pairing is
`5` in all five examples.

Every boundary contains the fixed seed walls `W1,...,W60`.  Their root
matrix has Smith diagonal

```text
1^13,
```

so each boundary spans the entire fixed Picard lattice rather than a proper
sublattice.

Every off-diagonal pair was classified directly.  The ordinary relation
inventories are:

| facets | m2 | m3 | m4 | m6 | infinity |
|---:|---:|---:|---:|---:|---:|
| 752 | 13,016 | 608 | 3,816 | 0 | 264,936 |
| 761 | 13,158 | 624 | 3,847 | 0 | 271,551 |
| 770 | 13,300 | 640 | 3,878 | 0 | 278,247 |
| 779 | 13,442 | 656 | 3,909 | 0 | 285,024 |
| 788 | 13,584 | 672 | 3,940 | 0 | 291,882 |

Here `infinity` combines parallel and ultraparallel pairs, as in the
ordinary Coxeter matrix.  Each row sums to `f(f-1)/2`.

All arithmetic uses Python arbitrary-precision integers and exact SymPy
rationals.  The largest boundary coordinate appearing in the five pilots is
`115,772`; there is no fixed-width overflow assumption.

## Why finite volume and why index five

The certificate rechecks the two relevant stock properties directly from
the seed Gram matrix:

- all 112 walls `W61,...,W172` are doubling walls;
- all `C(112,2)=6,216` pairs are disjoint in Allcock's sense—parallel or
  ultraparallel.

A properly edge-colored finite subtree in the corresponding chamber tree
therefore satisfies Allcock's finite-subtree construction.  Each displayed
`P5` is a fundamental polytope for a reflection subgroup of index five and
is the union of five seed chambers with disjoint interiors.  Consequently
its volume is exactly five times the seed volume.

Finite volume is not inferred from the divisor lattice.  The seed is
separately pinned as finite-volume noncompact by both the literal Vinberg
certificate and the independent exact GMP cone certificate.  A finite-index
subgroup of this nonuniform lattice is again nonuniform, so every pilot is
finite-volume noncompact.

## Algebro-geometric scope

The fixed Smith core shows that all five reflection arrangements live in the
same full Néron–Severi lattice

```text
U(2) + D5(-1)^2 + <-4>
```

carried by the pinned seven-dimensional K3 family.  This is a weak
full-Picard divisor-reflection carrier.  No full-nef realization, effectivity
of a prescribed chamber, or surface-automorphism realization is claimed.

This distinction is useful for the broader converse program: these examples
certify the soft lattice-carrier level.  Turning such a chamber into the nef
cone of a surface remains the separate no-extra-effective-wall problem.

## Bounded lineage

The pinned project predecessors have the JK-H12 facet counts

```text
172;
317,326;
462,471,480;
607,616,625,634.
```

The pinned project-lineage snapshot has H12 headers
`14,17,18,19,23,172`.  None contains `752,761,770,779,788`.  The existing
path-cube lane begins at `941` facets through the formula `886N+55`.

At CoxIter commit `b800be48d8240fbabd2d1c9542e61054584ad7fe`, the pinned
1,214-file ordinary lane has H12 headers `13,14,16,18`, and the pinned
2,959-file index-two lane has H12 headers `18,20,25,35`.  The five pilot
headers are absent from both.

This is bounded project/header absence only.  It is not a
publication-complete denominator, Allcock already proved exponential H12
abundance, and worldwide priority remains open.

## Reproduction

```bash
python3 scripts/certify_jacobian_kummer_h12_five_chamber_pilot.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_five_chamber_pilot.py
```

Machine-readable archive:
[`jacobian_kummer_h12_five_chamber_pilot.json`](../data/enumerations/jacobian_kummer_h12_five_chamber_pilot.json).

Certifier:
[`certify_jacobian_kummer_h12_five_chamber_pilot.py`](../scripts/certify_jacobian_kummer_h12_five_chamber_pilot.py).

Current SHA-256 values:

```text
script  d1e678dc2eea7c67998fff1d5fcef34a8ad45281fac87424f74f3876ea4c9f6d
JSON    809302f4aa10fcf723c68061c34a4d857c9d5c64a4c03e04227c4c54c917ebec
test    4a8589decfc63dfce057a2d98f7fb764d59009904bb995d08b200f30f0cf9d6b
```

Fresh focused test result: `4 passed in 6.99s`.

Fresh exact archive replay: `archive_check=PASS`, payload SHA-256
`809302f4aa10fcf723c68061c34a4d857c9d5c64a4c03e04227c4c54c917ebec`.
