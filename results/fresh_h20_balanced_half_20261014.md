# A radial gap beyond the complete balanced endpoint families

For the unit-inradius affine `A6 × A6 × A7` cusp seed, let `z` be any alcove vertex, so `|z|²=196`. The point `x=2z/5` has ball power at least `164/25` against **every cap giving a finite or ideal endpoint on an original cusp edge**, provided it retains all 22 seed walls at Coxeter angles or dotted/parallel relations. This includes every real radius and every omitted-wall angle label. It is not a cutoff on a sampled cap bank.

The same point lies strictly outside every ball contained in the Euclidean ball of radius five about the incenter. Consequently, endpoint caps together with any independently established family of such central-contact bridges still leave a strict spacelike escape. Arbitrary interior caps and changed cusp seeds are outside this statement.

The exact checker is [fresh_h20_balanced_half_20261014.py](../scripts/fresh_h20_balanced_half_20261014.py), and the certificate is [fresh_h20_balanced_half_20261014.json](fresh_h20_balanced_half_20261014.json). It uses the frozen [complete endpoint census](fresh_h20_balanced_endpoint_20261014.md) and [first all-radius bounds](fresh_h20_balanced_endpoint_20261014_bounds.json). The endpoint census has 41 profiles and 392 cyclic transports, or 16,072 literal profile occurrences. It includes all four positive-semidefinite endpoint profiles; an ideal endpoint alone is not asserted to have a bounded link.

## Sharper bounds with no finite-label cutoff

An endpoint cap has center `w=z−Ry`, with `|y|²≤1`. The three omitted seed ratios are `A_f/R−B_f`, where `A=(7,7,8)` and `B_f` is the sum of prescribed cosines in that factor. A Coxeter-compatible ratio is zero, a finite cosine `cos(π/m)`, or at least one. In particular, it belongs to `{0,1/2} ∪ [1/√2,∞)`.

The earlier certificate therefore bounded the continuous part by `min A_f/(B_f+1/√2)`, with only the finitely many exceptional radii `A_f/B_f` and `A_f/(B_f+1/2)` left to test. The present checker eliminates twelve retained exceptional occurrences using a necessary arithmetic property of the angle itself: **`2cos(π/m)` is an algebraic integer**. This property holds without any arithmetic assumption on the polytope. The checker constructs the monic minimal polynomial from the distinct Galois conjugates in `Q(√2,√3,√5)`; a nonintegral coefficient excludes every finite label at once.

The resulting necessary upper bounds are:

| Endpoint profile | Radius upper bound |
|---|---:|
| All right | `8` (in fact strictly less) |
| One `π/3` contact in either `A6` | `7` |
| One `π/4` contact in either `A6` | `14(√2−1)` |
| One `π/3` or `π/4` contact in `A7` | `16(√2−1)` |
| One `π/3` contact in each `A6` | `7` |
| One `π/3` contact in `A6` and one in `A7` | `7(√3−1)` |
| Both endpoint nodes of one `A6`, each `π/3` | `7` |
| Both endpoint nodes of `A7`, each `π/3` | `16/3` |

Two rows use additional short all-label arguments. For the all-right profile, `R>8` would give an `A6` cosine less than `7/8<cos(π/7)`, hence label 3, 4, 5, or 6. The corresponding `A7` cosine is `8/7` times this value. In each case its doubled value is not an algebraic integer. At `R=8`, the `A6` doubled cosine is the noninteger rational `7/4`.

For the mixed `A6/A7` profile, `R>7(√3−1)` would put the touched `A6` cosine below `cos(π/6)`. Its only possibilities are zero or labels 3, 4, 5. The touched `A7` cosine would be `8c/7+1/14`, and again all four doubled values fail the algebraic-integer test. The exact minimal polynomials are saved in the certificate. The remaining table rows follow the continuous bound and the exhaustive exceptional-radius check.

## Exact checks of all transported directions

Write `I` for the inverse norm-two Gram of the retained `A6+A6+A7` chains. If `d_i=1−n_i·x` and `c_i` is the profile cosine, then

- `|z−x|² = 2dᵀId`;
- `(x−z)·y = −2dᵀIc`;
- `|y|² = 2cᵀIc`.

All entries of `I` are nonnegative. For `x=λz₀` with the three tested values of `λ` below, every `d_i` is positive. Hence the ball power

`|x−w|²−R² = 2dᵀId −4R dᵀIc + (|y|²−1)R²`

is nonincreasing for positive `R`. Substituting the rigorous upper bound from the table gives a lower bound for every actual radius. The checker evaluates all 16,072 expressions exactly in the multiquadratic field. It archives their digest, their minimum, and every nonpositive case.

| Probe | Result |
|---|---|
| `2z₀/5` | Every lower bound is positive; minimum `164/25` |
| `49z₀/100` | Only the own-vertex all-right profile and five inner `A7` profiles can cut; no zero lower bounds |
| `z₀/2` | The same six possible strict-cutter classes; 44 further zero lower bounds are retained explicitly |

The positive gap at `2z₀/5` is uniform, and cyclic symmetry gives the corresponding statement at every old vertex.

For the Lorentz form `|x|²−2DU`, use

`X=(2z₀/5, 3111/200, 1)`.

Its norm is `1/4`, its spatial square is `784/25`, and every seed gap is at least `3/5`. An exterior sphere inequality at `X` subtracts `1/4` from its Euclidean ball power, so every endpoint inequality is strict. Also `X` pairs negatively with the future timelike vector `(0,16,1)`.

If a further ball is contained in the radius-five ball, then `|x−w|−R≥3/5`, and its power is at least `9/25>1/4`. Thus the same `X` survives such additions strictly. A strict future spacelike feasible vector gives an open set of ideal directions and therefore infinite hyperbolic volume. This is a theorem about the stated families; it does not say all central-contact caps have the radius-five property without its separate hypotheses.

## Which endpoint caps can cut half-vertex probes?

The exact screen also sharply restricts genuine strict cutters of `z₀/2`. Every one is centered from the same old cusp vertex `z₀`, and is one of:

1. The all-right profile, with `R=7/cos(π/m)` for an integer `m≥7`.
2. A single `π/3` contact at one of the five inner `A7` nodes, with the unique radius `r=16(√2−1)` and omitted `A7` angle `π/4`.

To see uniqueness in the second case, put `k=i(8−i)∈{12,15,16}`. Half-probe power is `49−kR/2+(k/16−1)R²`. A strict cut forces `R>49/8`; thereafter the omitted `A7` cosine is at least `1/√2` and below `79/98<cos(π/5)`. Its only finite value is `1/√2`.

There is also a compatibility obstruction between two adjacent half-probes whose old vertices differ by a nearest-neighbor `A6` step, of squared distance 84. Two all-right strict cutters have positive normal pairing because both radii exceed seven. An all-right cutter and an inner-`A7` cutter have normalized negative pairing strictly between zero and `1/2`, which is not a Coxeter cosine. For two inner-`A7` cutters, equal node indices give a positive pairing; unequal indices give a cosine strictly between zero and `1/2`. The latter follows from `42<r²<44` and the exactly checked direction-difference squares `{7/16,3/4,15/16,1}`. Thus no mutually Coxeter-compatible pair of endpoint caps can strictly cut both probes.

This last statement is deliberately about **strict cuts**. The 44 zero-bound occurrences at the exact half-probe were not discarded. The universal strict escape at `2z₀/5` is the direct volume obstruction and does not depend on this secondary two-probe argument.

Run:

```sh
python scripts/fresh_h20_balanced_half_20261014.py --verify
```

The producer uses the previous exact multiquadratic arithmetic. It introduces no numerical approximation in a sign or angle decision, and does not claim an independent second arithmetic implementation.
