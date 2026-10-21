# Hostile replay and lineage of the five parked rank-ten `H9` sections

## Decision

An independent replay confirms all five archived sections as exact
finite-volume, noncompact hyperbolic Coxeter polytopes in `H9`. It does not
import or call the batch builder. Instead, it reconstructs the oriented
mirrors from the immutable batch, reruns GMP double description, and checks
every ray incidence, time-sheet sign, elliptic edge endpoint, and affine cusp
neighborhood.

The lineage result is deliberately narrower than a novelty claim:

- `A14_8_0_C01` is exactly Im Hof's classical doubly truncated orthoscheme
  `9-imhof85_03`;
- the other four labelled diagrams have no exact match in the pinned ordinary
  and index-two CoxIter lanes, but worldwide priority remains open;
- every section lies in a known arithmetic class after passing to its unique
  maximal reflection-stable even overlattice, so there are zero new arithmetic
  classes.

The executable audit is
[`audit_alexeev_parked_rank10_h9_examples_lineage.py`](../scripts/audit_alexeev_parked_rank10_h9_examples_lineage.py),
the frozen output is
[`alexeev_parked_rank10_h9_examples_hostile_lineage.json`](../data/enumerations/alexeev_parked_rank10_h9_examples_hostile_lineage.json),
and the regression gate is
[`test_alexeev_parked_rank10_h9_examples_hostile_lineage.py`](../tests/test_alexeev_parked_rank10_h9_examples_hostile_lineage.py).

## Independent geometry replay

The input batch is pinned by SHA-256
`e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430`.
For every case, the audit reconstructs the saturated rank-ten Gram matrix and
all active primitive mirrors, verifies every reflection is integral and
unimodular, and independently checks one strict timelike witness for each
facet. A fresh `cddexec_gmp --rep` run reproduces the archived primitive ray
set exactly.

| ID | facets | fresh rays | finite + ideal + hyperideal | elliptic edges | affine components |
|---|---:|---:|---:|---:|---:|
| `A14_8_0_C01` | 12 | 49 | `44+5+0` | 260 | 10 |
| `A14_8_1_C02` | 16 | 188 | `165+23+0` | 1,112 | 55 |
| `A14_10_1_C01` | 13 | 77 | `70+7+0` | 393 | 14 |
| `A14_10_1_C02` | 15 | 138 | `123+15+0` | 778 | 35 |
| `A16_14_1_C01` | 15 | 109 | `100+9+0` | 588 | 23 |

Every extreme ray is in the chosen closed time cone, with strictly positive
pairing against the archived interior vector. There are isotropic rays and
no hyperideal rays. Independently, every elliptic rank-eight edge has its
two required finite-or-ideal endpoints, while every connected affine
component has a complete zero-neighborhood equal to an ideal vertex. Thus
both literal Vinberg failure counts are zero for every case.

## Arithmetic overclosure

For each displayed reflection group, the audit exhausts the integral
isotropic subgroups of the fixed lattice's discriminant group that are stable
under all facet reflections. Each case has one maximal stable even
overlattice. Exact signature, Smith group, discriminant quadratic form, and
the uniqueness theorem for indefinite even 2-elementary lattices identify
the normalized lattice with the pinned Scharlau--Kirschmer `List9` entry shown
below.

| ID | stable even subgroup orders | maximal index | maximal determinant | `(r,a,delta)` | normalization | `List9` entry |
|---|---|---:|---:|---|---|---:|
| `A14_8_0_C01` | `1,4` | 4 | -64 | `(10,6,0)` | scaled dual `M*(2)` | 23 |
| `A14_8_1_C02` | `1,2,4,8` | 8 | -16 | `(10,4,0)` | sign reversal | 23 |
| `A14_10_1_C01` | `1,2,4` | 4 | -256 | `(10,8,0)` | scaled dual `M*(2)` | 7 |
| `A14_10_1_C02` | `1,2,4,8` | 8 | -64 | `(10,6,0)` | scaled dual `M*(2)` | 23 |
| `A16_14_1_C01` | `1,16,32` | 32 | -16 | `(10,4,0)` | sign reversal | 23 |

The Scharlau--Kirschmer snapshots are pinned by hashes
`b7f6982a5c03b1ff07a5a9f383517438a926a499d27fb33ea89a8a69d7fb66c3`
for `List9` and
`772187d28a067c2d10f2fe00d7cee07a063a6c3d79be5cc8400c2c0ca8c96c36`
for `Gram9`. Entry 7 has a full 10-facet chamber with eight finite vertices
and two cusps; entry 23 has a full 11-facet chamber with 24 finite vertices
and two cusps. The five sections may still be interesting reflection
subgroups and divisor-lattice carriers, but they do not define new maximal
arithmetic reflection classes.

## Exact corpus collision and bounded exclusions

The pinned ordinary CoxIter archive contains 1,214 records and has `H9`
header counts `38,10,13` at 10, 11, and 12 facets. The pinned index-two lane
contains 2,959 records and has eight `H9` records at each of 12, 13, 15, 19,
27, and 43 facets.

There is exactly one labelled-graph collision among the five candidates:

```text
A14_8_0_C01
  = CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_03.coxiter
candidate -> Im Hof walls: 6,11,10,9,8,1,2,3,4,5,7,12
```

The raw Im Hof record has SHA-256
`8f61b20c7e8ccc658073a14d020cec94f67fcb6a8b62ffa7a7dc836cbf5b7171`.
The `_02` and `_04` records have the same global relation histogram but are
not isomorphic to the candidate. The other four candidates have zero exact
matches. All five are pairwise nonisomorphic and have zero exact matches
against the pinned three-type project factory and nine doubling-cube types.

Header and structural exclusions also separate the candidates from the
named Vinberg `P9`, Tumarkin `n+2`, Bredon `P*`, and Allcock `F`, `Q`,
minus-two, and `D7D9` stocks recorded in the audit. These are finite,
reproducible comparisons. They do not turn the CoxIter bundle or the cited
named families into a complete classification of all noncompact nonsimple
`H9` Coxeter polytopes.

## Scope boundary

This hostile audit certifies hyperbolic geometry, bounded literal-diagram
lineage, and arithmetic overclosure. It makes no surface-automorphism claim,
does not identify these cones with full nef cones, and makes no worldwide
novelty claim for the four unmatched diagrams.

## Replay

```bash
PYTHONPATH=.:src python3 -u scripts/audit_alexeev_parked_rank10_h9_examples_lineage.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_alexeev_parked_rank10_h9_examples_hostile_lineage.py
```
