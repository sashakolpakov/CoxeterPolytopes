# TOOT: eight ADEG polytopes return through complete K3 nef stocks

## Result

The enlarged exact denominator gives **eight polytope-first complete-nef
divisor-lattice returns** among the 30 remaining stable ADEG triality
packages.  Five are new in this pass; the last three were already certified
from the rank-9/10 stocks.

| ADEG package | orientation | complete K3 nef stock | Gram isomorphisms | fixed rank | invariant ample square |
|---|---:|---|---:|---:|---:|
| `N2-01-path-66` | 1 | rank 5, `[4]+D4(-1)`, 5 walls | 12 | 3 | 198 |
| `N3-01-path-633` | 0 | rank 6, `U+D4(-1)`, 6 walls | 6 | 4 | 630 |
| `N3-04-star-6` | 0 | rank 6, `U(2)+D4(-1)`, 6 walls | 120 | 4 | 198 |
| `N3-05-triangle-leaf-6` | 1 | rank 8, `U+A2(-1)^3`, 10 walls | 48 | 4 | 162 |
| `N5-01-path-633336` | 0 | rank 10, `U+D4(-1)^2`, 11 walls | 72 | 6 | 504 |
| `N5-02-path-triangle` | 0 | rank 8, `U+A2(-1)+D4(-1)`, 9 walls | 12 | 6 | 504 |
| `N6-01-diamond-tail` | 1 | rank 9, `U+D4(-1)+A3(-1)`, 10 walls | 12 | 7 | 504 |
| `N7-02-path-star` | 0 | rank 10, `U+D4(-1)^2`, 11 walls | 72 | 8 | 504 |

For each row an exact weighted-wall permutation extends to a unimodular
integral Neron--Severi isometry, transports clone-cycling to an order-three
nef-wall symmetry, fixes a strict integral ample class, and returns the
original folded ADEG Gram by orbit sums.  The complete published wall list
then proves that the fixed section has no hidden nef wall.

Thus the certified transition is

```text
ADEG arithmetic polytope
  -> complete divisor/root system of a published K3 nef cone
  -> the same ADEG polytope as the triality-fixed nef section.
```

These are **positive inverse calibrations** at the marked divisor-lattice
level: the search starts with each ADEG polytope and identifies a complete
surface carrier afterward.  They are not yet the requested flagship
polytope-to-surface theorem.  This finite 30-package comparison is relative
to the declared 55-stock local denominator; it proves neither class-wide
essential surjectivity onto surface carriers nor an algebro-geometric
pushback theorem capable of sharpening dimension bounds.  Both endpoints
are published objects, so no new abstract polytope, K3 lattice, or reflection
group is claimed.

## 1. The new lower-rank source compiler

Roulleau's primary atlas ([arXiv:2003.08985](https://arxiv.org/abs/2003.08985))
contains exactly 17 complete rank-5--8 wall systems whose `(Picard rank,
wall count)` header occurs among the previously uncovered ADEG packages:

| header | `5/5` | `6/6` | `7/7` | `8/8` | `8/9` | `8/10` | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| stocks | 2 | 3 | 2 | 3 | 5 | 2 | **17** |

The source compiler freezes every ordinary and parallel intersection, an
integral wall basis, Smith form, exact inertia, a published strict ample
witness, source line range, and complete-curve lineage.  It also runs an
independent exhaustive Vinberg audit.  All 17 polytopes are finite-volume
noncompact: every elliptic codimension-one edge has exactly two finite/ideal
endpoints, every affine zero-neighborhood is full parabolic, and each has at
least one cusp.  In particular, none is being promoted merely because its
Gram has Lorentzian signature.

The compiler is
[`certify_roulleau_rank5_8_matching_full_nef_stocks.py`](../scripts/certify_roulleau_rank5_8_matching_full_nef_stocks.py),
with [machine archive](../data/enumerations/roulleau_rank5_8_matching_full_nef_stocks.json)
and [source audit](roulleau_rank5_8_matching_full_nef_stocks.md).

## 2. Exact local denominator and gates

The comparison denominator now has 55 complete K3 nef-root stocks:

| Picard rank | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| stocks | 2 | 3 | 2 | 10 | 10 | 9 | 4 | 4 | 3 | 3 | 1 | 1 | 1 | 1 | 1 | **55** |

The 17 lower-rank rows are all published Roulleau stocks meeting the target
headers.  The 38 rank-9--19 rows are all complete stocks frozen in the
project's existing machine denominator.  This remains a declared local
denominator, not a worldwide or all-rank K3 classification.

The exact gate sequence is:

1. unfolded rank and complete wall count;
2. full nonzero Smith diagonal;
3. exact weighted-Gram isomorphism;
4. unimodular wall-lattice basis change;
5. integral order-three action and invariant strict ample class;
6. complete-nef/no-hidden-wall provenance; and
7. exact orbit-sum return to the original folded Gram.

Across the 30 packages there are 61 target--stock pairs after the header
gate, nine after Smith form, and eight after weighted-Gram isomorphism.  The
single Smith survivor rejected by the graph gate is
`N5-02-path-triangle -> U+A1(-1)+A5(-1)`; the other determinant-12 stock
`U+A2(-1)+D4(-1)` is the exact match.

## 3. Why these are full-nef same-polytope returns

Write `G` for an unfolded positive-normal ADEG Gram and `C` for the negative
intersection Gram of the complete K3 curve stock.  Each promotion freezes a
permutation `pi` satisfying

```text
G_ij = -C_pi(i),pi(j)
```

entrywise, together with a determinant-`±1` change between integral wall
bases.  Clone-cycling transported through `pi` is therefore an integral
order-three isometry of the entire divisor lattice, not just a rational
symmetry of a graph.

Averaging the published ample class over triality gives an integral fixed
class of positive square; its minimum pairing with every wall is `3` in all
eight rows.  Hence it lies strictly inside the full nef chamber.  It also
forces the coinvariant lattice to be rootless: a square-minus-two class
orthogonal to an ample class would contradict K3 Riemann--Roch.

On the fixed subspace, inequalities in one triality orbit coincide.  Their
orbit-sum normals reproduce the original folded Gram exactly.  Because the
Roulleau list is the complete irreducible `(-2)`-curve set, no omitted wall
can further cut this section.  Therefore

```text
complete K3 nef cone intersect fixed divisor subspace
    = original ADEG Coxeter chamber
```

at the marked divisor-lattice level.

## 4. Honest remainder

The 30-package ledger is now

```text
P = 8,    R = 15,    O = 7.
```

`P` is an exact complete-nef return.  `R` means every stock in the declared
denominator with that header fails Smith form or exact weighted Gram.  `O`
means no denominator stock has the header; it is not a nonexistence result.

Only three headers remain open for lack of a local stock:

| unfolded rank / roots | packages |
|---:|---:|
| `10/12` | 1 |
| `14/16` | 5 |
| `16/17` | 1 |
| **total** | **7** |

The lower-rank compilation therefore closes all 15 formerly header-empty
packages in ranks 5--8: five promote and ten become exact relative
rejections.  Further positive fishing should now target complete `10/12`,
`14/16`, and `16/17` carriers, or enlarge the carrier category beyond K3
nef-root stocks.

## 5. Geometric-action boundary

All eight transported order-three divisor-lattice actions are nontrivial on
the discriminant group.  The archive consequently certifies **zero**
order-three automorphisms of the named K3 surfaces.  A geometric lift would
need a compatible transcendental Hodge action, glue, and Torelli argument.
Until then these are exact complete-nef divisor-lattice fixed sections, not
quotient surfaces or automorphism fixed loci.

This distinction is essential: it preserves the polytope-to-surface result
without silently turning a lattice symmetry into a surface automorphism.

## 6. Reproducibility

- certifier:
  [`certify_bredon_adeg24_remaining_k3_stock_promotions.py`](../scripts/certify_bredon_adeg24_remaining_k3_stock_promotions.py)
- archive:
  [`bredon_adeg24_remaining_k3_stock_promotions.json`](../data/enumerations/bredon_adeg24_remaining_k3_stock_promotions.json)
- 30-package ledger:
  [`bredon_adeg24_remaining_k3_stock_promotion_ledger.md`](bredon_adeg24_remaining_k3_stock_promotion_ledger.md)

The frozen promotion payload SHA-256 is
`91a4867d958bb05af27ef08fbbc752779c65657ee7e94350fb1cb87db40302fa`.
