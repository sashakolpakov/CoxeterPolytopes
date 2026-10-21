# Cubic-seven H22 `B21 x A1` adaptive Vinberg checkpoint

This checkpoint extends the exact search for

\[
L_{22,t}=\operatorname{diag}(-\theta,1^{21},t),\qquad
\theta=2\cos(2\pi/7),\qquad t\in\{2,2-\theta\},
\]

through Vinberg heights 20, 50, and 100.  It does **not** find a
finite-volume chamber and makes no compactness, termination, novelty, or
nonreflectivity claim.

## Complete arithmetic input

For a primitive root `r`, let

`I(r)=(-theta*r0,r1,...,r21,t*r22)`.

Reflection integrality gives `(q(r)) | 2I(r)`.  Localizing this divisibility
and using primitivity leaves exactly these possible norm ideals:

- `t=2`: `(1)` and `(2)`;
- `t=2-theta`: `(1)`, `(2)`, `p=(2-theta)`, and `(2)p`.

Here `(2)` is inert.  For `t=2`, the only nonunit value of `I(r)` can occur at
`(2)`; if all non-final coordinates are even and the final coordinate is a
unit, the diagonal equation gives `v_2(q)=1`.  For `t=2-theta`, the same
argument gives `v_2(q)<=1`; at `p`, the exceptional primitive case gives
`q=(2-theta) r_final^2 mod p^2`, hence `v_p(q)=1`.  No other prime can divide
the norm.

This is also complete at the element level modulo totally-positive unit
squares.  Dirichlet's theorem gives `|O*/O*^2|=8`, and the units `-1`,
`theta`, and `1+theta` have independent real sign vectors `---`, `+--`, and
`++-`.  The signature map is therefore an isomorphism, so every totally
positive unit is a square.  The canonical norm elements are consequently
`1,2` and `1,2,2-theta,2(2-theta)` respectively.

## Exhaustive Weyl quotient

Height zero is exactly `B21 x A1`, not a sampled sector.  Signed permutations
of the 21 ordinary coordinates give the `B21` action, and the independent
sign of the distinguished final coordinate gives `A1`.  After orienting a
positive-height root by its control coordinate, taking absolute values,
sorting the ordinary coordinates, and making the final coordinate
nonnegative gives the unique orbit representative.  The ordered
sum-of-squares recursion enumerates every such representative, including all
multiplicities and zero padding.

The three height-zero orbit patterns expand to `2*21^2+2=884` roots.  At each
positive height, all candidates are tested against only lower-height walls;
the complete surviving packet is proved pairwise acute before being added.
Thus an arbitrary tie-break cannot change a same-height wall set.

## Exact H22 checkpoints

The ray columns list `(finite, cusp, hyperideal)` projective extreme rays.
Every reported H20/H50 cone passed exact full-rank, irredundancy, and
time-sheet checks.

| seed | height | pool | accepted positive walls | total walls | exact rays | finite volume |
|:--|--:|--:|--:|--:|:--:|:--:|
| inert-2 (`t=2`) | 20 | 11 | 4 | 26 | `(9,0,132)` | no |
| inert-2 (`t=2`) | 50 | 30 | 10 | 32 | `(19,0,554)` | no |
| inert-2 (`t=2`) | 100 | 94 | 24 | 46 | capped / not archived | untested |
| ramified-7 (`t=2-theta`) | 20 | 15 | 4 | 26 | `(9,0,151)` | no |
| ramified-7 (`t=2-theta`) | 50 | 37 | 8 | 30 | `(16,0,712)` | no |
| ramified-7 (`t=2-theta`) | 100 | 109 | 22 | 44 | not run by policy | untested |

The H100 norm splits are:

- inert-2 pool: 87 norm-2 and 7 norm-1 roots; all 24 accepted roots have norm
  2;
- ramified-7 pool: 90 norm-2, 11 norm-`2(2-theta)`, 6 norm-1, and 2
  norm-`2-theta` roots; the accepted packet has 20 norm-2 and 2
  norm-`2(2-theta)` roots.

An exact incremental inert-2 H100 ray run was stopped after exceeding the
requested ten-minute checkpoint budget without completing.  No partial ray
set is treated as a certificate.  The implementation now checks time/ray
limits within facet-pair passes and reuses projective keys for retained rays;
the frozen H100 claim remains arithmetic walls only.

Both trajectories move away from prefix closure between H20 and H50: the
hyperideal counts grow from 132 to 554 and from 151 to 712.  Inert-2 is only
the marginally better continuation candidate (fewer H50 hyperideal rays and
more finite vertices), not a plausible closure claim.

## Exact zero-padding audit, dimensions 22--29

At H50, endpoint enumeration proves for **both** seeds that the complete H29
root pool is exactly the H22 pool with seven zeros inserted before the final
coordinate.  The accepted H29 positive walls are exactly the similarly
padded H22 walls.  All complete-pool roots use at most five nonzero ordinary
coordinates.  The counts are:

| seed | complete pool patterns | accepted patterns | maximum ordinary support |
|:--|--:|--:|--:|
| inert-2 | 30 | 10 | 5 |
| ramified-7 | 37 | 8 | 5 |

Zero-padding injects the ordered-partition solutions from dimension `d` into
dimension `d+1`.  Since the complete endpoint sets at dimensions 22 and 29
are equal after padding, no new pattern can first occur in an intermediate
dimension.  This proves the same H50 pool and accepted-pattern stabilization
for every dimension 22 through 29.

This bounded stabilization does **not** yield a persistent hyperideal
witness.  A root above H50 can have a new control, coordinate support, and
wall inequality, and may cut every current hyperideal ray.  No inequality
against the complete unbounded root family has been proved, so the correct
scope remains a finite exhaustive prefix rather than a nonreflectivity
theorem.

## Reproduction

- Engine: `src/coxeter_polytopes/search/cubic7_bna1.py`
- Generator: `scripts/search_cubic7_bna1_prefixes.py`
- Replay tests: `tests/test_cubic7_bna1_search.py` (9 passed)
- Archive: `data/enumerations/cubic7_bna1_h22_adaptive_h20_h50_h100.json`
- Archive SHA-256:
  `af631899484c0b823d163fe1bbd191e271d5b486e3a75ae4d31f98747e9871b5`

```console
PYTHONPATH=src python3 scripts/search_cubic7_bna1_prefixes.py \
  --dimensions 22 --heights 20,50,100 --seeds both \
  --geometry-through-height 50 --geometry-engine incremental \
  --geometry-ray-limit 5000 --geometry-time-limit-seconds 600 \
  --output data/enumerations/cubic7_bna1_h22_adaptive_h20_h50_h100.json
PYTHONPATH=src pytest -q tests/test_cubic7_bna1_search.py
```
