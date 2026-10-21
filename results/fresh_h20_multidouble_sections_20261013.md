# H20 sections of successive ambient doubles

**This construction leaves the previously excluded finite spherical stars.** Successive reflections produce exactly certified convex, finite-volume ambient unions of up to 32 Borcherds chambers. A first batch gives six new actual right-angle repairs; a second changes the cusp quotient direction and gives eight more finite-volume sections, six with actual target-angle repairs. Every tested section still has an actual forbidden ridge. No Coxeter H20 survivor is claimed.

The [producer](../scripts/fresh_h20_multidouble_sections_20261013.py), [exact archive](../data/enumerations/fresh_h20_multidouble_sections_20261013.json), [independent Fraction verifier](../scripts/fresh_h20_multidouble_sections_verify_20261013.py), and [PASS record](../data/enumerations/fresh_h20_multidouble_sections_verify_20261013.json) record the literal construction and its limits.

## Leaving finite reflection stars

Start with the original 210-facet Borcherds H21 chamber. Double successively across original walls 169, 194, and 71. Their Gram matrix is

```
[[ 4,-4,-4],
 [-4, 4,-4],
 [-4,-4, 4]].
```

The first two mirrors are parallel. On their normal span, the product of their reflections is `I+N`, where

```
N=[[2,-2],[2,-2]],    N!=0,    N^2=0.
```

Its powers are distinct. Thus the generating subgroup is infinite, and these finite chamber unions are not among the finite spherical stars excluded in the earlier report. The third seam also removes the remaining original seam implicated in that report's simplest three-ridge obstruction.

At each step, the selected outward seam normal has nonpositive pairing with every other current boundary normal. This sufficient convexity condition is checked exactly; the whole ambient union need not be acute. Every reflected boundary normal is a nonnegative combination of its predecessor and the seam normal, so its plane supports both halves. Conversely, the sign of the seam pairing selects a half, and reflecting the other half recovers the preceding inequalities. The new intersection therefore equals the doubled union. Each listed boundary plane contains a preceding full facet, and exact deduplication gives true facets.

## Literal sections and adaptive continuation

The first tests use the irrational cusp family from the preceding round,

```
u=r34+r46+r204+r210,
e0=r46/2-u/36 + sqrt(2)(r34-r46-r204+r210)/4,
e(tau)=e0+tau*u.
```

One old control plane is tested in the four- and eight-chamber unions. Six new values of `tau` impose orthogonality to newly reflected boundary mirrors and give actual right ridges in the eight-chamber section. They differ from the previous twelve-plane and single-double samples. Each nevertheless has an exact forbidden ridge.

To repair these displayed defects, the eight-chamber union is doubled across its boundary normal 396, then across original wall 22. These two norm-two seams are orthogonal, are eligible for convex doubling, and are not incident to `u`. The enlarged unions remove the old displayed seam facets. The same section planes are then checked again.

| Ambient chambers | True ambient facets | Literal section checks | Result |
|---:|---:|---:|---|
| 4 | 744 | 1 | Actual forbidden ridge |
| 8 | 1,456 | 7 | Six actual right repairs; all seven still forbidden elsewhere |
| 16 | 2,625 | 1 | Actual obtuse local repair |
| 32 | 3,353 | 7 | Seven actual obtuse local repairs |

Every seam is strictly nonincident to `u`, and the exact set of cusp-incident boundary normals remains the original set. All sixteen section checks retain a complete rank-19 Coxeter cusp, including the positive affine-kernel identities and the original redundant-wall identity. A strict future timelike interior certifies full dimension. Intersecting the finite-vertex Klein model of each finite-volume ambient union gives a finite-volume H20 section; the certified cusp makes it noncompact.

The last local failure has a simple explanation. In this family `(e,r1)=0`. Each repaired defect uses a norm-two seam `z` with `(r1,z)=-1`. Its reflection sends `r1` to `r1+z`, so

```
(r1,s_z(r1))=1.
```

Projection into the section preserves that positive pairing because `e` is orthogonal to `r1`. The new angle is therefore obtuse if the pair is a ridge. The archive supplies actual timelike ridge certificates for all eight local repairs in the last two rows; this is not merely a positive entry in a raw Gram matrix.

## Changing the cusp quotient direction

The next batch leaves that orthogonality condition. Four directions from the earlier 26-direction cusp-bridge census had not been materialized as global sections. They satisfy `(e0,r1)=2`, with `e0^2=4`, and remain orthogonal to their respective cusp vectors. The three seams 169, 194 and 71 are strictly nonincident to each cusp, so the eight-chamber union preserves its whole projected cusp link.

For each direction, a bounded scout selects twelve distinct angle functions involving reflected boundary mirrors and proposes roots for labels 2 through 6. Two selected algebraic parameters per direction are checked exactly. The literal inputs are preserved in the [proposal archive](../data/enumerations/fresh_h20_multidouble_sections_bridge_proposals_20261013.json); the numerical proposal counts are not exhaustive geometric claims.

| Original cusp bridge | Exact section parameters checked | Actual target ridges | Remaining actual bad ridge |
|---|---:|---:|---|
| 1/22 | 2 | 2 | Original 1/23 |
| 1/24 | 2 | 1 | Original 1/23 |
| 1/49 | 2 | 1 | Original 1/2, obtuse |
| 1/50 | 2 | 2 | Original 1/3 |

All eight are finite-volume, noncompact H20 sections with a complete rank-19 Coxeter cusp. Six target equations correspond to certified actual ridges, at angle pi/3 or pi/5. The other two are only exact angle identities between projected planes, and are not claimed as dihedral repairs. All eight rejection ridges are actual, with exact timelike witnesses.

An additional constructive improvement is an exact interior-point method. The projection of the original center did not lie inside these sections, so the cusp equations were used directly. If a cusp component has positive marks `m_i` and

```
sum_i m_i (q*r_i - (r_i,e(t))*e(t)) = mu(t)*u,    mu(t)>0,
```

solve the linear equations

```
(x0,r_i) = -mu(t)/(q*sum_i m_i),
(x0,u)=-1,    (x0,h)=0.
```

For these four pencils the equations have an exact solution in the parameter field, and `(x0,e(t))=0` is verified. Then `x=x0+L*u` for a sufficiently large positive integer `L` leaves all cusp inequalities strict, makes every nonincident inequality strict, and has square `x0^2-2L<0`. The archived shift gives a future timelike section interior without depending on an interior-point numerical scout.

The [second producer](../scripts/fresh_h20_multidouble_sections_bridges_20261013.py), [eight exact sections](../data/enumerations/fresh_h20_multidouble_sections_bridges_20261013.json), [independent verifier](../scripts/fresh_h20_multidouble_sections_bridges_verify_20261013.py), and [PASS record](../data/enumerations/fresh_h20_multidouble_sections_bridges_verify_20261013.json) preserve all details. The independent checker reconstructs each original bridge and full affine link, the exact cusp center, target-angle branch, and actual forbidden ridge; a forged center is rejected.

## Replay and scope

```sh
python scripts/fresh_h20_multidouble_sections_20261013.py
python scripts/fresh_h20_multidouble_sections_20261013.py --continuation
python scripts/fresh_h20_multidouble_sections_verify_20261013.py
python scripts/fresh_h20_multidouble_sections_bridges_20261013.py
python scripts/fresh_h20_multidouble_sections_bridges_verify_20261013.py
```

The independent checker uses only standard-library arithmetic and earlier standalone verification helpers. It reconstructs every doubling and boundary list, checks the infinite subgroup, strict interiors, full cusp links, six new right ridges, sixteen forbidden ridges and eight actual obtuse local repairs, and rejects a forged ridge certificate.

The first batch concerns the declared unions and seven literal section planes; the second checks eight further literal sections in four changed cusp directions. The local angle identity explains failure of its particular doubling repair wherever the actual ridge persists. These results do not exclude arbitrary successive doublings, all parameters in the four new pencils, other cusp quotient directions, nonreflection gluings, or a finite-volume Coxeter polytope in H20.
