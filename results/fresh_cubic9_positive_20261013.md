# Changed-coefficient prime-five hosts and targeted A4 completion

The new H9 route uses `L=<-T>+B²+A4+A1`, where `T³−3T−1=0`, `T=2cos(π/9)`, and `B=[[2,−T],[−T,2]]`. A height-30 full search gives 13 walls. A targeted search through height 1000, restricted to two free A4 Cartan covectors and the A1 scalar, finds nine additional mutually Coxeter-compatible walls. Several of these cut the earlier exact escapes. The resulting **22-wall cone still has infinite volume**, certified by an exact strict spacelike witness. No finite-volume polytope is claimed.

## Why this is a different arithmetic route

This host has determinant `−10Tβ²`, where `β=4−T²`. Its determinant ratio to the previous prime-five host `<-2T>+B²+A4+A1` is `1/2`, with nonsquare field norm `1/8`. Since the ambient vector dimension is ten, determinant squareclass is unchanged by scaling the quadratic form. Thus the new coefficient is not a disguised field similarity of that previous host. Its odd determinant valuation at 5 also distinguishes it from all four earlier conductor-nine `B⁴+<1>` spaces. These are comparisons to the explicitly searched hosts, not a worldwide novelty assertion.

Removing the two spherical `I₂(18)` components from a hypothetical compact H9 chamber forces a compact H5 Coxeter face with lattice `<-T>+A4+A1`. This face is cheaper to inspect and its missing walls identify useful A4 directions. I also inspected the alternative coefficient `Tβ=T−1` in this H5 face. Both lower searches have seven-wall height-100 prefixes with exact escapes; these prefixes do not prove full-host nonreflectivity. The previously proved nonreflectivity of `<-(T−1)>+B+<1>` does not automatically apply to these changed positive blocks.

## Complete primitive classes and zero-height chambers

The field/class-number/unit facts are those independently checked in the [conductor-nine source review](fresh_cubic9_source_review_20261013.md). The H9 dual exponent is `10β`, so primitiveness and reflection integrality initially allow norms dividing `20β`. Norms divisible by 4 are imprimitive here. Indeed, α is a unit at 2, as are the determinants of B and A4. Reflection integrality at 2 then forces the time, B and A4 coordinates to be even. Their contributions to the norm are divisible by 4. The remaining A1 term is `2z²`, so norm divisibility by 4 forces z even as well, because `O/(2)` is a field. This excludes `4,20,4β,20β` without discarding valid norm-four roots in the previous **even-time** host.

A complete safe primitive norm pool is therefore

`{1,2,5,10,β,2β,5β,10β}`.

These are possible classes, not an assertion that each occurs. Every remaining ideal is squarefree, so an integral vector of one of these norms is primitive. For the H5 coefficient T, β drops from the dual exponent and the pool is `{1,2,5,10}`. For coefficient `Tβ`, it remains in the pool.

The exhaustive dominant zero-height census gives `I₂(18)²×A4×A1` in H9, with six dominant orbits and 94 oriented roots, and `A4×A1` in H5, with two dominant orbits and 22 oriented roots. Mixed-block zero roots are included in that enumeration, rather than excluded by assumption.

## Positive targeted completion

The initial H9 escape lies entirely in A4, on the end fundamental-weight ray `Tω4`; the H5 computation also exposes a middle-weight escape. This suggests changing spatial freedom instead of increasing the full H9 cutoff.

The first five slices have time, one of `5ω1,5ω2,5ω3,5ω4,ω2+ω3`, and the A1 scalar free. At height 1000 their complete root census finds one new compatible wall. It cuts a later mixed A4/scalar gap but does not cut the original pure-A4 escape.

Next, all six pairs of A4 Cartan covectors are allowed independently, with scalar zero. The exact integrality condition is

`w1+2w2+3w3+4w4 = 0 mod5`.

This searches the saturated intersection with the lattice; it does not require integral multiples of the earlier weight rays. Four new mutually compatible norm-two walls appear, including a reversal pair that cuts the original end-weight escapes.

Finally, the scalar is also freed in those same six slices, at the same height bound. This includes both preceding families and yields nine compatible new walls in total. Every one is retained: there are no cross-angle conflicts among them. The first 13 walls and all nine additions pass the full exact Coxeter-angle gate, including rejection of acute non-Coxeter angles. The complete targeted pool is `{1,2,5,10}`: if a norm had a factor β, reflection integrality in the unimodular-at-β time/A4/A1 slice would force every coordinate divisible by β, contradicting its norm valuation one.

| Search | Controls | Walls | Result |
|---|---:|---:|---|
| H5, coefficient T, full height 100 | 15 | 7 | Exact escape |
| H5, coefficient Tβ, full height 100 | 23 | 7 | Exact escape |
| H9, coefficient T, full height 30 | 12 | 13 | Exact escape |
| Five time/weight/scalar slices, height 1000 | 53 | 14 | One new wall; exact escape |
| Six saturated time/two-weight slices, height 1000 | 53 | 17 | Four new walls; exact escape |
| Six saturated time/two-weight/scalar slices, height 1000 | 53 | 22 | Nine new walls; exact escape |

The final escape has norm

`2446/1369 − (11857/20535)T − (596/20535)T² > 0`.

It satisfies the entire compatible bounded dictionary, since all nine compatible roots occur in the final cone. This excludes completion by this saved dictionary, not by arbitrary later roots in the same slices or by unrestricted roots. The full host remains open.

The scalar-free symmetric-middle plane has an additional local obstruction. Its root norm equation would be `4t²−Tx²=q`. Reflection integrality excludes primitive q with factors 5 or β, leaving q=1 or2. Reducing modulo 4 excludes both: for q=2, x is even and the left side is divisible by 4; for q=1, reduction modulo 2 forces `x=T+1`, whose square lift gives `−Tx²=−1−2T² mod4`, never 1. This explains why freeing the scalar or a second independent A4 direction matters. It is only a statement about this single binary plane.

## Exact verification and artifacts

The independent Fraction verifier repeats the full bounded census for all three initial hosts using unpruned coefficient boxes and a residue-indexed Cartesian A4 join. A separate slice verifier independently repeats the saturated two-weight and two-weight-plus-scalar censuses, their exact chronological selections, every wall angle, strict timelike interior and strict spacelike perturbation. No numerical failure is used as a volume certificate. All three initial-host reenumerations and all three targeted-family audits pass. The final cone has 231 exact pair checks. Three focused regressions also pass, covering the time-parity norm-pool change, saturation beyond integral five-weight rays, and a new wall cutting the original end-weight escape.

- [H5 coefficient T](fresh_cubic9_positive_20261013_h5_theta_h100.json), [H5 coefficient Tβ](fresh_cubic9_positive_20261013_h5_theta_beta_h100.json), [H9 initial chamber](fresh_cubic9_positive_20261013_h9_theta_h30.json).
- [Initial-host independent census](fresh_cubic9_positive_20261013_independent.json).
- [Five ray/scalar slices](fresh_cubic9_positive_20261013_target.json), [saturated two-weight slices](fresh_cubic9_positive_20261013_weights.json), [final two-weight/scalar slices](fresh_cubic9_positive_20261013_rank4.json).
- [Independent saturated-slice census and geometry](fresh_cubic9_positive_20261013_slices_independent.json).
- [Generalized producer](../scripts/fresh_cubic9_positive_20261013.py), [targeted producer](../scripts/fresh_cubic9_positive_20261013_rank4.py), [independent full census](../scripts/fresh_cubic9_positive_20261013_independent.py), [independent slice census](../scripts/fresh_cubic9_positive_20261013_slices_independent.py).

```sh
PYTHONPATH=scripts python scripts/fresh_cubic9_positive_20261013_independent.py results/fresh_cubic9_positive_20261013_h5_theta_h100.json results/fresh_cubic9_positive_20261013_h5_theta_beta_h100.json results/fresh_cubic9_positive_20261013_h9_theta_h30.json --output results/fresh_cubic9_positive_20261013_independent.json
PYTHONPATH=scripts python scripts/fresh_cubic9_positive_20261013_slices_independent.py results/fresh_cubic9_positive_20261013_weights.json results/fresh_cubic9_positive_20261013_rank4.json --output results/fresh_cubic9_positive_20261013_slices_independent.json
```
