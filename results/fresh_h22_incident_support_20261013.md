# Walls perpendicular to the whole exposed H22 cusp diagram miss the escape

An exact conic identity rules out another local way of cutting the fixed
global escape. This is separate from the coupled affine-A20 bridge exclusion.
It does not rule out other incident pairings or walls away from this cusp.

Let `delta4` be the previously certified null ray, let `I` be its 17 incident
old walls, and let `X` be the lower placement's strict spacelike escape. The
[certificate](fresh_h22_incident_support_20261013_orthogonal.json) gives the
exact identity in `Q(sqrt(11))`

\[
X=\tau\delta_4+\sum_{i\in I}\mu_i w_i
 +\lambda_6w_6+\lambda_8w_8+\lambda_{38}w_{38}
 +\lambda_a k_{95611}+\lambda_b k_{95665}.
\]

The five `lambda` coefficients are strictly positive. The `k` normals are
two retained positive-height members of the old null-cap family. The
coefficients `tau` and `mu_i` are unrestricted. No center wall is used.

Consequently every proposed new normal `r` perpendicular to `delta4` and
all 17 old incident walls, and nonpositive against the five displayed
retained walls, satisfies `(X,r)<=0`. This includes every direction in the
four unused cusp tangent dimensions and every real shift along `delta4`,
with no restriction on the normal's norm or coefficient field. Any finite
compatible collection of these walls leaves a spacelike feasible point;
perturbing it toward a common strict center gives a strict escape.

The identity holds uniformly for the old center labels `m>=3`, because its
support does not use those center facets. It only concerns this fixed cusp
and the all-right incident profile. In particular, it does not rule out
transporting new repairs to other cusps.

An LP found the support. Verification uses no optimization: it reconstructs
the 23 exact coordinate identities and checks the coefficient signs. A
second checker independently implements arithmetic on coefficient pairs,
using rational bounds on `sqrt(11)` to prove the five signs. Both pass.

```sh
python3 scripts/fresh_h22_incident_support_20261013.py --orthogonal --verify
python3 scripts/fresh_h22_incident_support_20261013_independent.py --verify
```

The [independent result](fresh_h22_incident_support_20261013_independent.json)
records the positive rational coefficient bounds. A different LP attempting
an identity for **all** acute incident normals was infeasible on its selected
support; that failed numerical scout is not a general exclusion.
