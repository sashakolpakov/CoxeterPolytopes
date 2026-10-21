# A continuous limit on the paired-support extension

**Changing the parameters of the four paired-cap curves used by the frozen 92-wall lift cannot complete the chamber.** This holds for every real parameter, both sheets of each normalized hyperbola, all four height reflections, and arbitrary compatible selections or replacements of the paired caps. A single exact strict spacelike point survives the entire admissible family.

This is a restriction on these four curves, not on the H18-to-H20 construction or on arbitrary normals perpendicular to their supports. Each selected parameter must contribute **both** orthogonal caps and retain all sixty predecessor facets with their acute Coxeter pair gates. A single unpaired cap, a rotation of the pair within the three-dimensional support complement, a changed support, or replacement of predecessor walls is outside the theorem.

The [exact checker](../scripts/fresh_h20_structured_extension_20261014.py) and [certificate](fresh_h20_structured_extension_20261014.json) use Fraction arithmetic and outward rational bounds for positive square roots, with no producer imports or numerical acceptance. They bind the [92-wall construction](fresh_h20_coupled_lift_20261014.md) and its [60-wall predecessor](fresh_dimension_construction_routes_20261014.md). The separately archived [numerical scout](../scripts/fresh_h20_structured_extension_20261014_scout.py) and [proposal output](fresh_h20_structured_extension_20261014_scout.json) suggested the inequalities but supply no premise of the proof.

The [separate review](fresh_h20_structured_extension_20261014_review.md) and fresh primary saved replay passed. The review independently checks the four rational frames, support ranks, predecessor center/escape signs and symbolic normalization/derivative identities; it does not claim a second complete interval census.

## Complete parameterization of the stated curves

Let `v` be the predecessor's rational boundary escape, with square 25. For each of the four saved PD18 supports, let `h` and `w` be the frozen rational complement vectors. Define

```text
E = v/5,
Z = h − (h·E)E,
T = Z/sqrt(−Z²),
W = w/sqrt(w²).
```

The checker verifies `E²=W²=1`, `T²=−1`, all cross pairings zero, and that T is future-directed. All three vectors are perpendicular to the support. Every unit spacelike vector in the plane spanned by v,h has one of the two forms

```text
M_eta(t) = eta*(t+t^−1)E/2 + (t−t^−1)T/2,
eta in {+1,−1}, t>0.
```

The paired normals are `r_plus=M+W` and `r_minus=M−W`; both have square two and their mutual pairing is zero. Thus this parameterization includes all real values of the original dyadic-height recipes for this fixed W, not merely an interval around the saved choices. The final two coordinate sign changes preserve the form and permute all sixty old walls; these isometries are checked explicitly.

## Necessary bounds from old-wall acuteness and the angle gap

The negative-E sheet is impossible. For each support, a saved old wall has positive coefficients in both terms of `(-E+T)t+(-E−T)/t`. Its pairing with the pair's mean is therefore positive for every t>0, contradicting simultaneous acuteness of both caps.

On the positive-E sheet, simultaneous acuteness with an old wall r implies

```text
((E+T)·r)t² + 2|W·r|t + (E−T)·r <= 0.
```

Two explicit concave quadratics per support are strictly positive on consecutive closed intervals. Endpoint positivity and concavity certify the whole intervals, including singleton boundaries. They imply the following relaxed necessary bounds:

| Support | First necessary bound |
|---|---:|
| 0, 2 | `t > 3689/2000 = 1.8445` |
| 1, 3 | `t > 11157/2500 = 4.4628` |

For supports 0 and 2 an additional old norm-two wall gives a stronger all-label bound. Its pairing with `r_minus` is strictly decreasing for t>1. At the first lower bound it is negative; at t=5 it is still greater than −1. Consequently the whole intervening interval has pairing in `(−1,0)`, hence normalized cosine in `(0,1/2)`. No finite Coxeter label has such a cosine, and the pair is neither parallel nor dotted. Every admissible pair therefore satisfies **t>5** for supports 0 and 2. This argument has no label cutoff and does not use the previous source-section angle theorem.

The certificate also retains an instructive counterexample to a stronger acuteness-only assertion. On the height-19 reflected support 0, t=37/20 supplies a cap acute to all sixty old walls that cuts the saved escape. It has a forbidden old-wall cosine in `(0,1/2)`. Thus acuteness alone does not prove the all-reflected-family obstruction; the Coxeter angle gap is essential.

## One witness for the complete continuous family

Set the relaxed lower parameters to `(5,11157/2500,5,11157/2500)`. At each lower parameter write the mean as M and its forward derivative as U. For every larger parameter, with `a=log(t/t_lower)>=0`, a cap is

```text
r(a) = cosh(a)M + sinh(a)U ± W.
```

Use the exact rational strict escape X already saved for the 92-wall chamber. For all sixteen support/reflection families the checker proves

```text
(M,X) + |(W,X)| < 0,
(U,X) < 0.
```

Because `cosh(a)>=1` and `sinh(a)>=0`, both cap pairings stay strictly negative for every a>=0. The weakest certified upper bound is approximately −0.04735, so this is not a zero-pairing or limiting argument. All sixty old inequalities are independently strict at X, its time orientation is future, and

```text
X² = 1411253767129583933/200000000000000000 > 0.
```

The same old strict timelike point H is strict for every admissible pair: `E·H<0`, `T·H<0`, `W·H=0`, and t>1. The checker verifies these identities and that all four reflections fix H. Thus every finite selection of admissible paired caps from these curves defines a nonempty chamber, retains a strict spacelike escape and has infinite volume. Allowing deletion or replacement of the previous thirty-two caps does not change that conclusion.

## No finite contacts between different curve families

The same necessary parameter tails also have a stronger local separation property. At their lower endpoints, every cross pair from two different named support/reflection families has pairing below −2; since both normal squares are two, it is strictly dotted. Eight signed monotonicity inequalities for each of the 120 family pairs show that this pairing can only decrease as either parameter increases. All 480 boundary pairings and the full monotonicity conditions are certified.

Within a single family, plus caps at two distinct parameters have pairing `cosh(a−b)+1>0`, so two complete pairs at different parameters cannot coexist acutely. This explains why simply replacing the four dyadic recipes or accumulating more parameters cannot join the finite components by new finite cross contacts.

A finite anchor event inside one curve can still improve its local edge geometry; the theorem says it cannot remove this global escape. The next constructive change must leave the certified family, for example by changing the pair orientation within a PD18 complement, changing the support, or retaining only a differently coupled individual normal. No global nonexistence statement for the sixty-wall core follows.

```sh
python scripts/fresh_h20_structured_extension_20261014.py --verify
```
