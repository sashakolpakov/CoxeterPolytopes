# A changed conductor-nine H9 host

The new host `L = <-2T> + B² + A4 + A1`, where `T³−3T−1=0`, `T=2cos(π/9)`, and `B=[[2,−T],[−T,2]]`, gives a certified 34-wall height-100 Coxeter prefix with infinite volume. It is not a new finite-volume polytope. Its full reflection chamber remains unresolved.

This is a different quadratic space up to similarity from the four conductor-nine `B⁴+<1>` hosts already searched. Writing `β=4−T²`, its determinant is `−20Tβ²`. Every determinant ratio to those four hosts has odd valuation at the prime above 5, whereas their determinant support is at 2 and 3. Equivalently, the independent audit records a nonsquare rational field norm for each ratio. In dimension ten, multiplying a form by any field scalar changes its determinant by a square. This distinction is relative to those four explicitly compared hosts; it is not a worldwide novelty claim.

## Root classes and the zero-height chamber

The field arithmetic and class-number/unit facts are reviewed in [the conductor-nine source audit](fresh_cubic9_source_review_20261013.md). In particular `O_K=Z[T]` is principal, 2 is inert, and every totally positive unit is a square. The new Gram matrix is even and has dual exponent `10β`. For a primitive integral reflective vector, its norm divides twice that exponent. A complete safe normalized norm pool is therefore

`{2,4,10,20,2β,4β,10β,20β}`.

These are possible classes, not a claim that every class occurs. A vector divisible by 2 is discarded individually. This suffices for primitivity: the prime factors 5 and β have exponent one in every candidate norm, so divisibility of the whole vector by either is impossible. No unproved exclusion of the norm-four classes is used.

The full dominant zero-height census over all eight classes has six vectors: two norm-two binary orbits, two norm-`2β` binary orbits, one norm-two A4 orbit and one norm-two A1 orbit. The resulting finite chamber is exactly `I₂(18)² × A4 × A1`, with 94 oriented primitive roots and nine simple walls. The census is complete even for roots that might mix blocks: every zero-height root can be moved into this product chamber by its known finite reflection subgroup, and the exhaustive dominant census leaves no extra orbit.

For a fixed norm `q` and time coordinate `x`, spatial energy is `E=q+2Tx²`. A Cartan covector `w` satisfies `w²≤2E` in every embedding and reflection integrality gives `w=(q/2)k`, `k∈O_K`. The producer enumerates this exact field alphabet. A4 covector prefixes are pruned by `w_prefixᵀ A_prefix⁻¹ w_prefix≤E`, the minimum possible full energy with that prefix fixed. The independent verifier instead uses a Cartesian join with the exact A4 integrality condition `w₁+2w₂+3w₃+4w₄=0 mod5`. Binary blocks and the A1 coordinate are joined by exact energy equality. Time coordinates use all three embedding bounds and are processed in exact height order.

Every pair passes the finite-angle gate: for intersecting walls, `4 cos²(angle)` must be one of `0,1,2,3,T²,2+T`, corresponding to labels `2,3,4,6,9,18`. Arithmetic reflection integrality and an acute angle alone are not substituted for this check.

## Bounded results and independent replay

| Height bound | Time/norm controls | Walls | Exact pair checks | Producer time | Independent census time |
|---|---:|---:|---:|---:|---:|
| 30 | 12 | 20 | 190 | 1.54 s | 3.7 s |
| 100 | 23 | 34 | 561 | 4.99 s | 7.18 s |

The independent verifier imports only a separate standard-library Fraction polynomial implementation. It repeats every time/norm control, every primitive candidate census, chronological acceptance, all wall angles, the strict timelike interior, and the exact escape inequalities. Its root alphabet uses an unpruned certified coefficient box, rather than the producer's embedding-interval sieve. Thus this is a second bounded enumeration, not just a check of the accepted roots. Three focused regression tests also pass, including an actual primitive norm-four root, the full zero census, and rejection of the acute non-Coxeter angle 2π/9.

The height-100 escape has time coordinate 1, zero binary and A1 coordinates, and A4 coordinates

`[(2+4T−2T²)/5] · (1,2,3,4)`.

Its norm is `(-48−122T+80T²)/5 > 0`. The audit also records a strict spacelike perturbation inside every wall, so the infinite-volume conclusion does not depend on a floating-point optimizer or a boundary-only witness. Further roots might cut this escape; this is not a full-host nonreflectivity theorem.

Artifacts:

- [Height-30 roots and geometry](fresh_cubic9_mixed_h9_h30_20261013.json), [independent census](fresh_cubic9_mixed_h30_independent_audit_20261013.json).
- [Height-100 roots and geometry](fresh_cubic9_mixed_h9_h100_20261013.json), [independent census](fresh_cubic9_mixed_h100_independent_audit_20261013.json).
- [Producer](../scripts/fresh_cubic9_mixed_search_20261013.py), [independent verifier](../scripts/fresh_cubic9_mixed_independent_20261013.py).

```sh
PYTHONPATH=src:. python scripts/fresh_cubic9_mixed_search_20261013.py --height 100 --output results/fresh_cubic9_mixed_h9_h100_20261013.json
PYTHONPATH=scripts python scripts/fresh_cubic9_mixed_independent_20261013.py --source results/fresh_cubic9_mixed_h9_h100_20261013.json --output results/fresh_cubic9_mixed_h100_independent_audit_20261013.json
```

## Why the initial B³+A3 proposal was set aside

For `α=2T` and `2Tβ`, the A3 factor has a unique integral index-`(2)` overclosure I3. The exact A3 embedding and the dyadic residue check are saved in [the triage certificate](fresh_cubic9_mixed_a3_triage_20261013.json). The full integral isometry group of `<-α>+B³+A3` equals that of this overclosure: uniqueness gives one inclusion, and the original even lattice is intrinsically `{v:q(v)∈2O_K}`, giving the other.

This does **not** identify either integral group with a previous `B⁴+<1>` group. Only the weaker field-space comparison is proved. For `α=2T`, the overclosed form is field-isometric to the old `α=2Tβ` form. For `α=2Tβ`, it is field-isometric to β times the old `α=2T` form. Exact PARI Hilbert symbols at the unique primes above 2 and 3 agree after those choices; the determinant ratios are squares and real signatures agree. At every other finite prime the displayed diagonal entries are units and the local Hilbert symbols are trivial. Hasse–Minkowski then gives the asserted field isometry or similarity. Integral old-host identification remains unresolved, so no positive enumeration of this ambiguous lane was launched.

[Python triage wrapper](../scripts/fresh_cubic9_mixed_a3_triage_20261013.py) and [PARI source](../scripts/fresh_cubic9_mixed_a3_triage_20261013.gp) reproduce the calculation:

```sh
PYTHONPATH=scripts python scripts/fresh_cubic9_mixed_a3_triage_20261013.py
```
