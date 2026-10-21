# The single-cap S61 quadrilateral gate

The [exact certificate](fresh_h20_s61_single_cap_20261014.json) excludes a specific next step from the frozen [107-wall chamber](fresh_h20_coupled_roles_20261014.md): one new wall perpendicular to S61 cannot close its existing three-sided region as a compact quadrilateral while retaining the old walls. The calculation exhausts sixteen finite contact-label pairs and both quadratic branches for each pair. Every branch has a strictly positive pairing with retained wall 0 or retained wall 95, so it fails acuteness.

This is a whole-parameter result for that one-cap quadrilateral family. It does **not** exclude several new closing walls, a cap with nonright S61 contacts, or a changed retained core. The 107-wall construction and its review remain unchanged.

## The actual retained frame

Use the frozen [98-wall edge-closure recipe](fresh_h20_bridge_edge_closure_20261014.json). Its vectors X,W,r95 satisfy

```
X² = −A,  W² = w > 2,  r95² = 2,
X·W = X·r95 = W·r95 = 0,
Δ = A(w − 5/6),  A > 0.
```

The checker reconstructs these identities in the inherited radical tower. The unit Lorentz frame on the S61 plane is

`x = X/√A`, `e = r95/√2`, `f = W/√w`.

After division of each old side normal by the positive A, the S61 projections of walls 96 and 97 are

`L± = √(w−5/6)x − e/√2 ± √w f`.

The checker directly compares the two serialized ambient normals against the frozen formula `A(−ω−r95/2±W)+√Δ X` before deriving these projected coordinates. The full original side-wall proof remains an inherited certificate.

Each has projected square 4/3, whereas each corresponding normalized ambient normal has square 2. The discrepancy is the removed fundamental weight of square 2/3. Thus the existing corners on wall 95 have intrinsic squared cosine 3/8; treating the three old sides as part of a right-angled polygon would impose the wrong metric conditions.

## Why there are exactly sixteen label pairs

Normalize a prospective normal perpendicular to S61 to square 2. If it meets old side 96 at ambient Coxeter angle π/m, the Gram of the twenty incident ambient normals has final Schur pivot

`2 − 3cos²(π/m)`.

An actual finite endpoint requires this to be positive. For m≥6 it is at most −1/4, while m=2,3,4,5 all pass. The same argument applies to side 97. Hence the full finite endpoint problem has precisely the sixteen ordered pairs `(m,n) ∈ {2,3,4,5}²`; this is a proved label bound, not a numerical angle cutoff.

A compact quadrilateral closing the retained three-sided S61 region must meet both of its outgoing sides at finite endpoints. The following calculation covers every normal with those two contacts, including branches that would have the wrong geometric position or orientation. Rejecting this larger family is sufficient.

## Complete conic solution

Write the new normal as `a x + b e + c f`, and put

```
p = cos(π/m),  q = cos(π/n),
s = p + q,  d = q − p,
k = √(w−5/6),  B = 2w−8/3 > 0.
```

The two contact equations are

`ak + b/√2 = s`, `c√w = d`.

The norm equation reduces to

`B a² − 4sk a + 2s² + d²/w − 2 = 0`.

Consequently every branch is explicitly

```
D = 2s² + B(2−d²/w),
a = (2sk ± √D)/B,
b = √2(s−ak),
c = d/√w.
```

All sixteen discriminants are strictly positive. The two real branches for each pair therefore give all 32 possibilities. Substitution into the quadratic and the contact equations proves the identities algebraically; interval overlap checks in the implementation are only sanity checks on their evaluation.

## Exact rejection and next scope

The checker evaluates the needed frame pairings with old wall 0 directly from the physical X,W,r95 vectors. It also uses `r95·N0 = −1` exactly. Positive square roots and all divisions are enclosed by outward rational intervals at 384 bits. The only cosine radicals required are √2 and √5.

Sixteen branches have positive product with wall 95; the other sixteen have positive product with wall 0. The smallest certified lower bound exceeds 0.99999, so these are strict, robust acuteness failures. No numerical proposal is promoted to a certificate, and there is no need to evaluate all 107 pairs once one positive old-wall pairing has been established.

The right/right case is the particularly simple candidate proportional to `X−√Δ r95`. It is perpendicular to both outgoing side walls, but its orientation compatible with wall 95 is obtuse to wall 0. The other fifteen label pairs fail for the same certified two-wall reason.

A next construction preserving these sides must therefore use at least two new walls to close S61, or leave the perpendicular-to-S61 family. This gate says nothing about whether such a larger coherent closure exists and makes no global H20 nonexistence claim.

## Reproduction

```
PYTHONPATH=src:. python scripts/fresh_h20_s61_single_cap_20261014.py --verify
```

The [checker](../scripts/fresh_h20_s61_single_cap_20261014.py) passed generation and saved replay. It binds the frozen geometric inputs, reconstructs the physical frame identities, and proves every rejection by rational interval arithmetic. It reuses the previously frozen exact field helper; it is not a second implementation of the old 107-wall construction.

The [independent review](fresh_h20_s61_single_cap_20261014_review.md) also passed generation and saved replay. Without importing the primary arithmetic, it reconstructs the physical frame and both serialized side normals, solves an alternative quadratic in the e-coordinate, and checks all 32 branches with a positive old-wall pairing greater than 99/100. The source/proof review confirms the complete finite-label bound and the limited one-cap quadrilateral scope.
