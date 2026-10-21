# Complete endpoint locus on the H4+H4 edge [7,30]

The retained 55-wall chamber admits a complete, explicitly parametrized
family of Coxeter endpoint walls at its open S6 edge `[7,30]`. The exact
answer is a countable sequence with label `m>=8` to old wall 51, followed
by a continuous strictly dotted tail. Thus this edge is not a local
nonexistence obstruction. No new wall is appended by this certificate,
and compactness of the entire H3 face remains unresolved.

The [source](../scripts/fresh_nonarithmetic_h3_edge730_domain_20261017.py)
and [saved certificate](fresh_nonarithmetic_h3_edge730_domain_20261017.json)
use the frozen [complete endpoint domains](fresh_nonarithmetic_h3_endpoint_domains_20261017.md).
Every old native contact is included. Finite labels are treated without
a cutoff or a coefficient-field restriction.

## The conic and its exact admissible parameters

The eight old endpoint roots are
`[0,1,2,3,5,6,7,30]`, with diagram `H4+H4`. A finite terminal wall must be
right to all eight: the positive inverse of either H4 component has every
diagonal greater than 2. Their perpendicular space has signature `(1,1)`.

Let `k0²=-1`, `k1²=1`, `k0.k1=0` be the frozen exact kernel construction,
and set

```
A=(k0-k1)/sqrt(2), B=(-k0-k1)/sqrt(2),
r(t)=t*A+B/t, t>0.
```

Then `A²=B²=0`, `A.B=1`, and `r(t)²=2`. This is one entire component
of the norm-two hyperbola. On its other component the product with old
wall 47 has two strictly positive coefficients, so it is always obtuse.

On the viable component, every nonframe old product is exactly
`g_j(t)=alpha_j*t+beta_j/t`. All 47 coefficients `alpha_j` are strictly
negative. For wall 51, `alpha<0<beta`, so its product decreases strictly
from positive infinity to negative infinity. Prescribing product `-c`
there determines the unique positive parameter

```
t(c)=(-c-sqrt(c²-4*alpha*beta))/(2*alpha).
```

Write `t_m=t(2 cos(pi/m))` and `t_infinity=t(2)`. The complete coherent
Coxeter locus is

```
{t_m : integer m>=8} union [t_infinity,infinity).
```

For orientation, `t_8` is approximately `5.2454339474` and
`t_infinity` is approximately `5.3022486405`. Their authoritative values
are the algebraic recipes and rational isolating intervals in the JSON.

At every `t_m`, only wall 51 has a nonright finite contact, of label `m`.
All other 46 nonframe products are strictly less than −2. At
`t_infinity` wall 51 is parallel; above it every nonframe product is
strictly dotted. The eight defining old products remain exactly zero.

## Completeness of the finite-label decision

At `t_8`, the numerator `alpha_j*t²+2t+beta_j` of `g_j(t)+2` is
strictly negative for every nonframe wall other than 51, and its derivative
is strictly negative. These inequalities remain true for all `t>=t_8`.
The monotonicity of `2 cos(pi/m)` therefore proves acceptance for every
integer `m>=8`, including arbitrarily large labels.

Below `t_8`, wall 51 can have a Coxeter product only at labels 2 through 7.
The complete rejection table is:

| Label to 51 | Old blocker | Exact reason |
|---:|---:|---|
| 2 | 49 | Positive product |
| 3 | 49 | Strict gap between labels 2 and 3 |
| 4 | 49 | Strict gap between labels 4 and 5 |
| 5 | 50 | Strict gap between labels 6 and 7 |
| 6 | 50 | Strict gap between labels 10 and 11 |
| 7 | 50 | Strict gap between labels 24 and 25 |

The last row is particularly close to label 25, but its rational interval
is strictly separated. The new cosine isolators use Machin's formula for
pi with alternating rational arctangent bounds, followed by an alternating
cosine series with an explicit next-term remainder. They do not assume a
finite algebraic angle alphabet or use floating-point equality decisions.

## Actual finite endpoints and coherent interiors

The new endpoint has the explicit formula

```
x(t)=t*A-B/t.
```

It has square −2, is future directed, and is right to the eight old roots
and `r(t)`. Every other old product is strictly negative for `t>=t_8`:
the certificate checks `alpha_j*t_8²-beta_j<0`, and this numerator
decreases thereafter. Thus every admitted wall gives an actual finite
endpoint of the selected geodesic edge, with positive definite diagram
`H4+H4+A1`.

All admitted normals also pair strictly negatively with the same actual
old finite vertex supported by `S6+[7,8,11]`. Perturbing this vertex a
sufficiently small distance toward the old strict center supplies a common
strict future interior. For any finite collection of such normals, one
perturbation works for the entire collection. This statement concerns
interior feasibility; mutual acuteness is a separate condition.

Indeed, two normals from this component have product
`r(t).r(s)=t/s+s/t>=2`. No two can occur together in an acute completion.
A coordinated repair must choose one endpoint wall and use other supports
for its remaining facet neighbors.

## Parallel boundary and the incidence distinction

The parallel parameter is coherent as a Coxeter wall, but cannot belong to
a finite-volume acute completion. At `t_infinity`, the positive sum
`nu=r+old51` is a nonzero future null vector. Both summands are right to the
retained pi/5 pair `[0,1]`. The existing
[persistent forbidden-angle cusp theorem](fresh_h9_fixed_core_cusp_obstruction_20261013.md)
therefore applies. The finite-volume-compatible candidate locus is the
displayed locus with `t_infinity` removed. This is a necessary filter, not
a finite-volume sufficiency claim.

Old wall 51 is **not** an intrinsic facet of the S6 face. It has labels 6
and 4 to A2 nodes 5 and 6, and its projected square is
`2-(2/3)(5+sqrt(6))<0`. Consequently, the finite contact `r.51` does not
create an edge on the S6 H3 face. It does yield a different positive
support `H4+I2(m)`. On that different H3 face, old walls 7 and 30 are
perpendicular to the support and retain their pi/5 wedge. Closing that
face would be a separate task.

## Replay and scope

Generation and fresh saved replay pass. Source hashes, syntax and report
links are checked. The implementation shares the frozen physical interval
engine; no independent arithmetic reconstruction is claimed here.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_edge730_domain_20261017.py --verify
```

This replaces arbitrary sampled endpoint parameters by a complete exact
domain. It resolves all 55 old Coxeter gates for this one conic and proves
actual finite endpoint geometry. Mutual compatibility with walls on other
supports, compact closure of the S6 face, and ambient finite volume remain
open.
