# Two finite-volume `H10` sections from Roulleau's 59-wall chamber

## Theorem

Roulleau's complete 59-wall nef chamber for

`NS(X) = U + D4(-1)^2 + A1(-1)^2`

has weighted wall group `C2^2 x S6`, of order 2,880.  All **914** subgroup
conjugacy classes have been enumerated.  Exactly two proper classes have
fixed rank at least 11, and both give finite-volume noncompact hyperbolic
Coxeter polytopes in `H10`:

| semantic ID | pinned GAP class | exploratory alias | facets | finite vertices | cusps | elliptic edges |
|---|---:|---:|---:|---:|---:|---:|
| `p67-H10-38` | `C4` | `C3` | 38 | 526 | 141 | 4,518 |
| `S6transposition-H10-45` | `C12` | `C9` | 45 | 1,048 | 93 | 7,332 |

GAP class ordinals depend on the generating presentation.  The semantic IDs
are primary; `C4/C12` and `C3/C9` denote the same two conjugacy classes in the
pinned small-generator and exploratory all-element presentations.

## Exact finite-volume gate

For both sections:

- every orbit normal is spacelike and every multiwall orbit is internally
  disjoint;
- every pair has an allowed exact Coxeter relation;
- projection of the subgroup-averaged ample class `D40` gives a strict
  timelike witness for every displayed facet;
- every elliptic edge has exactly two finite/ideal endpoints; and
- every maximal affine component has a valid full zero-neighborhood.

The 38-facet endpoint histogram is `{(0,2):920, (1,1):1936,
(2,0):1662}`; the 45-facet histogram is `{(0,2):578, (1,1):3028,
(2,0):3726}`.  Each has zero bad edges and zero bad parabolic
neighborhoods.  Thus “finite volume” is a literal Vinberg conclusion, not an
inference from a finite list of vertices.

## Lattice and maximal-group lineage

Both saturated fixed lattices are even of rank 11 and determinant 64, with
Smith invariants `1^6 2^4 4`.  Exhaustive discriminant-module closure under
all facet reflections gives:

| section | maximal stable overclosures | each overclosure | Coxeter Euler characteristic | index in each `P10` group |
|---|---:|---|---|---:|
| `H10/38` | 15 | index 8, odd unimodular `I_(1,10)` | `-527/2580480` | 2,276,640 |
| `H10/45` | 3 | index 8, odd unimodular `I_(1,10)` | `-527/2903040` | 2,023,680 |

Vinberg's classical 12-wall `P10` chamber, equivalently
Scharlau--Kirschmer `List10` entry 1, has characteristic
`-1/11147673600`.  The two stated indices are the exact Solomon ratios.
Consequently these are new literal subgroup chambers only in the bounded
sense; their arithmetic commensurability class is classical.

Neither facet count occurs in the pinned ordinary CoxIter `H10` sizes
`{12,13,14,15,16,19,20,28,44}` or doubling sizes
`{15,16,18,22,30,46}`, and neither appeared in the pre-catch project ledger.
Worldwide literal priority remains open.

## Surface-action qualification

Neither selected action is an automorphism action on the general K3 surface.
The independent source/action audit proves:

- the first double-plane deck action is the central fixed-rank-4 action that
  swaps `B1_ij <-> C1_ij` and `B2_ij <-> C2_ij`;
- the second `F4` deck action is trivial on `NS`; and
- the central fixed-rank-11 `p67` action, which swaps `A6 <-> A7` and all
  `B1/B2`, `C1/C2` pairs, is neither deck action.

Accordingly both results are exact integral divisor-lattice/fixed-nef
constructions.  For the 45-wall transposition, the proposed special-period
Nikulin lift is now ruled out: every primitive integral `E7(-2)` reversal
has the wrong discriminant action.  An order-four or higher lift remains
open; see the
[`H10/45` primitive-glue boundary](r12_h10_45_e7_embedding_search_boundary.md).

## Reproducibility

- full certificate:
  [`roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json`](../data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json)
- certifier:
  [`certify_roulleau_r12_ud4x2a1x2_59_subgroup_folds.py`](../scripts/certify_roulleau_r12_ud4x2a1x2_59_subgroup_folds.py)
- independent surface-action certificate:
  [`roulleau_r12_ud4x2a1x2_surface_actions.json`](../data/enumerations/roulleau_r12_ud4x2a1x2_surface_actions.json)
- regression test:
  [`test_roulleau_r12_ud4x2a1x2_59_subgroup_folds.py`](../tests/test_roulleau_r12_ud4x2a1x2_59_subgroup_folds.py)

Replay with:

```text
PYTHONPATH=src:. python3 -u scripts/certify_roulleau_r12_ud4x2a1x2_59_subgroup_folds.py --check-archive
```
