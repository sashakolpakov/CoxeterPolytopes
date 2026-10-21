# A complete own-contact gate for the saved T8 block

For the saved case-8 incidence, the fixed mandatory contacts below admit **no compact hyperbolic Coxeter realization in dimension nine**, for any finite internal labels `m >= 7`, `n >= 4`, and for **any four finite own-contact labels** of walls 11 and 17 with walls 12 and 13. This closes all nine nonempty own-support patterns, including the previously saved asymmetric state. It does not exclude the whole case-8 incidence.

The [exact checker](../scripts/fresh_h9_t8_asymmetric_20261014.py) and [certificate](fresh_h9_t8_asymmetric_20261014.json) use six continuous polynomial blocks, 408 cases and 11,736 strictly positive Bernstein coefficients over `Q(sqrt(2))`. Their smallest coefficient is `7/50`. No dotted distance is fixed, and there is no upper cutoff on an allowed finite angle label.

The [independent review](fresh_h9_t8_asymmetric_20261014_review.md) reconstructs all 408 polynomial cases with separate sparse Fraction arithmetic, checks the literal domains of all nine patterns and the eighteen fresh actual-vertex failures, and matches the saved certificate. Its generation and saved replays pass; it does not import the producer's determinant functions or SymPy.

## Fixed data and the nine domains

The source is the realizable [case-8 incidence](fresh_h9_tree_cube_wedges_20261014.json) and the [saved mandatory assignment](fresh_h9_t8_continuation_20261014_next_case_08.json). All its mandatory contacts stay fixed except

```text
(1,2), (12,13), (11,12), (11,13), (12,17), (13,17).
```

The first two labels are `m` and `n`. The positive eight-wall support is

```text
S = [1,2,4,5,7,8,12,13].
```

Twenty-four literal actual-vertex witnesses force all cross contacts between its four two-node blocks right. For a contact involving the first block, the three-node path determinant is at most `3/4 - cos(pi/m)^2 < 0`. For the other blocks, the relaxed `4-3-4` path has positive null vector `(1,sqrt(2),sqrt(2),1)`. Larger contacts only decrease its quadratic value. Consequently `S` is the orthogonal sum `I2(m) + B2 + B2 + I2(n)` and has determinant `(1-t)(1-c^2)/4 > 0`, where `t=cos(pi/m)^2` and `c=cos(pi/n)`.

The missing triples `{11,12,13}` and `{12,13,17}` must be Lannér in a compact realization. Neither anchor can therefore be right to both 12 and 13. Each has precisely three possible own supports, yielding nine ordered patterns. Every pattern has its own actual-vertex witnesses in the certificate; the argument does not assume that swapping roots or anchors is an incidence automorphism.

Four eight-wall budgets retain either own root and delete its partner. A nonright own contact costs at least `1/4` in the Schur projection. Each nonright contact with a foreign `B2` costs at least `1/2`, while a nonright contact with `I2(m)` costs more than one. Thus each anchor has at most one foreign contact, at one of `4,5,7,8`, and that contact must have label three. If it exists, every own nonright contact of that anchor also has label three: an own label at least four would give cost at least `1/2 + 1/2` in an actual positive-definite vertex.

Two further budgets fix wall 14's support projection: it has contacts of label three with `2,4,5`, and right contacts with the other members of `S`. For each of the nine own-support patterns, two fresh actual-vertex minors exclude label four between 14 and either anchor. Hence both these meeting labels belong to `{2,3}`. If the anchors share an own root, a literal three-node affine triangle forces their mutual contact right. If their own supports are disjoint, their mutual contact remains any finite label; the polynomial check allows its cosine to vary over the larger closed interval `[0,1]`.

## One residual identity

Write the two own cosine vectors on roots `(12,13)` as `(e,a)` and `(b,d)`. Absent contacts have value zero. Put

```text
B = 1-c^2,
P = e^2 + 2cea + a^2,
C = b^2 + 2cbd + d^2,
R = eb + ad + c(ed+ab),
h = delta+p.
```

Here `delta` is the original mutual-anchor cosine. The foreign projection costs `u,v` are either zero or `1/2`. Their mutual overlap `p` is zero, `1/2`, or `sqrt(2)/4`, according to the two foreign choices. The third-wall residual pairings `ac,bc` equal minus the original cosine (`0` or `1/2`) minus `1/2+sqrt(2)/4` if the relevant foreign contact is at 4 or 5.

The residual Gram of walls `(11,17,14)` is

```text
[ 1-u-P/B       -h-R/B       ac ]
[ -h-R/B        1-v-C/B      bc ]
[ ac             bc         -Q ],
Q = sqrt(2)/2 + 1/[4(1-t)].
```

The checker expands its determinant directly as a polynomial identity. The two-anchor determinant is `A/B`, and the three-wall determinant is `N/B`, where

```text
A = ((1-u)(1-v)-h^2)B - (1-u)C - (1-v)P - 2hR + (ed-ab)^2,
N = -QA - B[(1-u)bc^2+(1-v)ac^2+2h ac bc]
    + C ac^2 + P bc^2 - 2R ac bc.
```

An exact Chebyshev calculation gives `cos(pi/7)^2 > 81/100`: the polynomial `64x^3-80x^2+24x-1` is increasing above `3/4` and has value `-559/15625` at `81/100`. Therefore `Q > 25/19+sqrt(2)/2`. The certificate proves `N>0` and `A<0` at that lower value. Increasing `t` only increases `N/B`, since its derivative is `-A/[4B(1-t)^2] > 0`.

## Complete continuous boxes

Local renaming of the two anchor vectors and the two coordinates of `I2(n)` reduces the nine already-verified literal domains to the following polynomial blocks. These renamings concern the residual algebra only.

| Block | Own supports | Internal range | Continuous own bounds | Cases |
|---|---|---|---|---:|
| asymmetric | single / double | `n >= 4` | each active cosine in `[1/2,1]` | 100 |
| double | double / double | `n >= 4` | each active cosine in `[1/2,1]` | 100 |
| same_high | same single root | `n >= 7` | both in `[1/2,1]` | 100 |
| opposite_high | opposite single roots | `n >= 7` | both in `[1/2,1]` | 100 |
| same_low | same single root | `n = 4,5,6` | both in `[sqrt(2)/2,1]` | 4 |
| opposite_low | opposite single roots | `n = 4,5,6` | both in `[sqrt(2)/2,1]` | 4 |

For the first two blocks, `c` ranges over `[sqrt(2)/2,1]`; for the high single-root blocks it ranges over `[9/10,1]`. In the low blocks it again ranges over `[sqrt(2)/2,1]`. Indeed, a single own label three with `n <= 6` gives a nonnegative triple determinant `3/4-c^2`, so it cannot be Lannér. Both own labels must then be at least four, and the vertex budgets exclude foreign contacts. This explains why the low blocks have four cases instead of one hundred.

In the other blocks, a foreign contact fixes all active own cosines of that anchor to `1/2`. The two foreign choices and the two third-contact labels give `5*5*2*2=100` cases. Shared-root blocks have `delta=0`; opposite-root blocks use its entire closed interval `[0,1]`.

Every box is mapped affinely to a unit cube. All tensor Bernstein coefficients of both `N` and `-A` are strictly positive. The certificate saves the variable order, exact interval endpoints, power coefficients, degrees and Bernstein coefficients for every case, including constant-polynomial cases. Closed endpoints are retained; the actual finite internal labels have `c<1`, so division by `B` is legitimate.

Finally, a Gram matrix of hyperplanes in `H9` has rank at most ten. The positive rank-eight support leaves residual rank at most two. The strictly positive three-wall determinant contradicts that requirement. This proves the stated family exclusion, without assumptions on the remaining dotted distances.

## Reproduction and remaining scope

```sh
python scripts/fresh_h9_t8_asymmetric_20261014.py --verify
python scripts/fresh_h9_t8_asymmetric_20261014_search.py --case 8 --limit 100000 --verify
```

Generation and saved replay pass for both commands. The new theorem uses the frozen incidence and arithmetic helper as explicit hashed dependencies. It is a complete necessary-metric obstruction for the specified family, not a classification of all labels on this incidence or of all compact nine-dimensional polytopes.

## The preserved next state

The separate [continuation](../scripts/fresh_h9_t8_asymmetric_20261014_search.py) applies the complete family gate, without substituting any one angle for the finite-label category `m >= 7`. It stops at the next necessary assignment after 1,029 search nodes: 11,969 spherical-vertex failures, 595 forced-orthogonal negative-block failures, and nine completed mandatory assignments removed by the new family theorem. This is a first-survivor search, not an exhausted incidence search; its 100,000-node limit was not reached.

The [saved next assignment](fresh_h9_t8_asymmetric_20261014_search_case_08.json) changes the middle internal label `(7,8)` from four to five, takes both anchors' own contacts to be `(3,3)`, and has last internal label `(12,13)=4`. The first internal label still denotes every finite `m >= 7`. All other mandatory contacts equal the previous state. It passes only the necessary triangle, actual-vertex and bridge conditions; no Lorentzian metric is claimed.

A useful next reduction is to vary that middle internal label coherently. If it is at least six, a nonright anchor contact into that block alone has projection cost at least one, so the positive own-contact cost excludes it. The previously proved label-four case and a separate exact label-five block calculation would therefore cover every label at least four, provided the other literal domain guards are rechecked. That extension has not been proved in this package. The current state is preserved as an open metric frontier.
