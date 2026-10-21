# An exact cap cutting the dominant residual axis

One exact cap extends the
[153-wall chamber](fresh_h20_pi8_triangle_face_20261016.md) to **154 actual
facets**. It has19 right contacts, one pi/5 contact, and133 strictly
dotted contacts. It creates exactly two new finite vertices, bringing
the certified total to **at least176**, while preserving all five proper
cusps and the compact triangle face.

This cap cuts a large residual exterior region rather than the current
escape E6. That strict escape survives, so the whole chamber still has
infinite volume. A separate fixed-chart numerical comparison estimates
about76% reduction of the exterior Euclidean volume. This percentage
is a numerical region score, not a hyperbolic volume or exact cone claim.

The [source](../scripts/fresh_h20_dominant142_cap_20261016.py) and
[certificate](fresh_h20_dominant142_cap_20261016.json) pass generation and
a fresh saved replay. Independent review of this new package is pending.
All predecessors remain frozen; the optional all-dotted152 branch is
not retained.

## Exact construction

Use the original mutually orthogonal residual roots n141,...,n144 and
their finite corner C. Put

```
q = q142 = 309/200,       d = -C^2 > 0,
e142 = n142/sqrt(q),     e0 = C/sqrt(d),
gamma = cos(pi/5) = (1+sqrt5)/4,
c = [q*gamma + 5*sqrt(25-q+q*gamma^2)]/(25-q),
s = sqrt(q)*(c+gamma)/5,
r153 = sqrt2*(-c*e142+s*e0).
```

The positive radical solves

```
(25-q)c^2 - 2q*gamma*c - (25+q*gamma^2) = 0.
```

Consequently c^2-s^2=1 and r153^2=2. Numerically c=1.107721126697,
s=.476493540912. The old affine completion gives

```
n149 = -n142 + 5*(e0-n144/sqrt(q144)).
```

Thus `r153.n149=-sqrt(2q)*gamma` exactly. The root is right to
`S16+{141,143,144}`. Reconstructing every old physical normal and checking
all153 contacts certifies pi/5 to149 and strict dotted separation from
all remaining133 roots. No finite label cutoff is used in this
certificate: it verifies this explicit root.

The inherited strict future center has negative pairing with every
wall. Acute facet projections certify all154 facets are actual. The
five old proper cusp rays have strictly negative new pairing, so their
bounded links remain unchanged.

## Complete new finite-vertex census

There are exactly20 finite old neighbors,

```
S16 + {141,143,144,149}.
```

Among their20 subsets of size19, exactly two give a finite vertex
with r153. Omitting149 gives `A10+A4+A2+4A1`; omitting144 gives
`A10+A4+A2+2A1+I2(5)`. The literal20-node Gram matrices are positive
definite, and the projected common-center points are future timelike
with all154 required pairings checked. Every other subset contains
the dotted pair144/149 and cannot be spherical.

Every old finite vertex persists under an acute extension sharing the
strict future center: the inverse of its positive-definite Coxeter
Gram is entrywise nonnegative, and the center decomposition forces the
new wall strictly negative there. In particular, all three vertices of
the old compact triangle persist, so convexity preserves that entire
face. The two new vertices are the complete new finite census; the
inherited ambient finite census remains only a certified lower bound.

## Remaining escape and monotone region score

The old strict E6 has exact norm `4750061/5000000000>0` and remains
strictly negative at all154 walls. Its new pairing is approximately
`-.812822587371`. The fixed future-chart sign is also checked. Thus the
ambient volume is still infinite.

In the fixed C chart

```
x = C - sum_i z_i*sqrt(d/q_i)*n_i,
```

the new wall imposes the particularly simple exact inequality

```
z142 <= s/c = .4301565885384587...
```

The previous largest exterior region reached z142 about4.02. The
[numerical source](../scripts/fresh_h20_dominant142_cap_20261016_residual_scout.py)
and [saved region score](fresh_h20_dominant142_cap_20261016_residual_scout.json)
retain every old inequality and use the same chart and sampled directions
for both chambers:

| Quantity | Old153 | New154 |
|---|---:|---:|
| Maximum vertex squared radius | 16.3345753 | 6.98324997 |
| Numerical residual vertices | 163 | 168 |
| Clearly finite | 16 | 18 |
| Near unit sphere, tolerance1e-6 | 19 | 19 |
| Clearly exterior | 128 | 131 |
| Estimated exterior angular fraction | .961934 | .492572 |
| Estimated exterior Euclidean4-volume | 1.28986 | .305608 |

The fixed-seed2^20 Sobol estimate gives76.3068% reduction; the2^18 prefix
gives76.3036%. These are numerical integration diagnostics without
rigorous error bounds. The new cap cuts seven old exterior vertices,
then introduces two finite and ten exterior vertices. The larger
exterior vertex count does not contradict the substantial monotone
region reduction.

The remaining largest directions are near the141 and143 axes. A future
coordinated construction should address those regions while retaining
all old Coxeter gates. This package makes no claim of their closure.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_dominant142_cap_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_dominant142_cap_20261016_residual_scout.py --verify
```
