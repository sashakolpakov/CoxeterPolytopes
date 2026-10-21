# Effective compact closure of the H2 face with its cutter retained

There is a finite constructive mechanism that retains the exact wall `r95` and closes the entire `S61` face to a compact right-angled polygon. It may require many sides. This note proves existence and termination of a specified doubling algorithm; **no resulting polygon has been computed or certified here**, and the conclusion is not finite volume of the full H20 chamber.

The inputs are the frozen [whole-H2 face certificate](fresh_h20_bridge_faces_20261014_review.md) and [exact cutter](fresh_h20_e7_bridge_face_20261014.md), with its [normal data](fresh_h20_e7_bridge_face_20261014_exact.json). No predecessor arithmetic is replayed in this note. The positive-definite support is

```
S61 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,61,87].
```

Its orthogonal complement `V` has signature `(2,1)`. Every old nonsupport normal has a strictly future timelike projection into `V`. The norm-two vector `r95` lies in `V`, is Coxeter-compatible with all old walls, and cuts the old strict escape. Let `L` be its geodesic and `H−` its retained half-plane.

## Turning old-wall compatibility into finite distance bounds

Write an old normal as `w_i=s_i+U_i`, with `s_i` in the support span and `U_i` in `V`. Put

```
U_i² = −a_i < 0,  w_i² = n_i > 0,
P_i = U_i/sqrt(a_i).
```

Thus `P_i` is a unit future point of H2. For any outward norm-two normal `r=sqrt(2)v` in `V`, where `v²=1`, if its geodesic is at distance `d_i` from `P_i` and `P_i` is in its negative half-plane, then

```
r·w_i = −sqrt(2a_i) sinh(d_i).
```

Consequently the pair is strictly dotted whenever

```
d_i > asinh(sqrt(n_i/a_i)).
```

Only finitely many such constraints exist. It is enough to make every new side other than `L` avoid a sufficiently large compact ball about a chosen point `O` in `H−`, with `O` on its retained side. For example choose

```
R > max_i [d(O,P_i) + asinh(sqrt(n_i/a_i))].
```

The balls involved may cross `L`: the deliberately finite old contacts of `r95` make that unsurprising. The argument below does not require those whole balls to lie in `H−`.

## Explicit initial polygon and doubling rule

In the form `−x0²+x1²+x2²`, a regular compact right-angled pentagon has unit outward normals

```
b = 5^(−1/4),   q = sqrt(1+1/sqrt(5)),
v_j = (b, q cos(2pi j/5), q sin(2pi j/5)),  j=0,...,4.
```

Direct calculation gives `v_j²=1`, adjacent products zero, and nonadjacent products `−(1+sqrt(5))/2`. Its center is `(1,0,0)`. Map its oriented side `v_0` to `r95/sqrt(2)` by a future-preserving Lorentz isometry of `V`. An algebraic such map is obtained from an algebraic unit timelike point on `L` and an orthonormal tangent; the existing exact edge point supplies the former after normalization. Denote the mapped pentagon by `P0`, and its center by `O`.

For a compact right-angled polygon `P` containing `P0`, double across any side `F` other than the distinguished line `L`: replace it by `P ∪ s_F(P)`. This is again compact and convex. At the two seam endpoints the two right angles add to a straight angle; every other angle remains right. The side `L` remains a side, possibly lengthened. The elementary doubling fact also appears in [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*, Lemma 1.3](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf).

For norm-two normals the exact update is

```
s_v(w) = w − (w·v)v.
```

Remove the doubled side, retain the other sides, add their reflected copies, and merge coincident consecutive sides. This uses a fixed algebraic number field after the initial pentagon is embedded. Intermediate added sides need not satisfy the old-wall constraints: they are temporary construction walls, not walls that must be retained in the final chamber.

A **round** means doubling once across each nondistinguished side present at the start of that round. Such a side remains available until it is doubled: doubling another side keeps its original supporting half-plane and original facet segment. Thus this is a finite, unambiguous procedure, with any fixed ordering of the sides.

## Why a finite round suffices, including the clipped-ball issue

The initial pentagon generates a locally finite reflection tiling of H2. Each doubled polygon is a finite union of its chambers. After `k` rounds it contains every chamber reachable from `P0` in at most `k` chamber-adjacency steps without crossing `L`: every exterior neighboring chamber is included when its containing boundary side is doubled. The chamber graph on either side of the mirror `L` is connected, so these polygons exhaust `H−`.

A boundary mirror other than `L` either meets `L` orthogonally or is disjoint from it. In the latter case, since it is an actual polygon side in `H−`, its whole geodesic lies in `H−`. In the former case, its nearest point to any point of the interior of `H−` also lies in that interior. Algebraically, projection to a normal orthogonal to `r95` preserves the sign of the pairing with `r95`. In particular, the closest point of any such boundary mirror to `O` lies in `H−`.

If such a mirror is within distance `R` of `O`, two chambers on opposite sides of its nearest point lie in `H−` and meet a slightly enlarged ball. Once both chambers have entered the doubled polygon, that mirror cannot remain a boundary side. Local finiteness gives only finitely many relevant chambers. Therefore after finitely many rounds **all** boundary mirrors except `L` have full distance greater than `R` from `O`. Their retained half-planes contain `O`, hence contain the whole ball `B(O,R)`, not merely its clipped part. The required dotted inequalities follow.

One can make this termination bound explicit. Let `rho>0` be the initial pentagon's inradius and let `D` bound its diameter. For the regular pentagon one may take

```
rho = asinh(5^(−1/4)),
D = 2 acosh(cot(pi/5)).
```

Any chamber meeting `B(O,R+1)` has its center within `R+1+D`; its disjoint radius-`rho` inball lies in `B(O,R+1+D+rho)`. Thus fewer than

```
K = 1 + ceil[(cosh(R+1+D+rho)−1)/(cosh(rho)−1)]
```

chambers meet that ball. A slightly perturbed geodesic from `O` to either chamber at a relevant nearest point stays inside the ball and `H−`, giving gallery length less than `K`. Hence `K` rounds suffice. This is a deliberately crude bound, not a feasible size estimate. If `f_k` is a facet-count upper bound, `f_0=5` and `f_(k+1)=f_k*2^(f_k−1)` is a valid computable upper bound. In practice one stops as soon as all exact dotted tests pass.

## What the final polygon proves in the ambient chamber

Every added final normal is orthogonal to `S61`, strictly dotted to every old nonsupport normal, and Coxeter-compatible with every other polygon side. The original `r95` is retained with its original finite contacts. The polygon is compact, so its vertices are actual future timelike points with two incident right sides. Together with `S61`, each gives a positive-definite rank-20 incident Gram and an actual finite H20 vertex: all old nonsupport inequalities are strict because their projections are future timelike.

An interior point of the polygon is strict for all nonsupport old and new walls. A sufficiently small perturbation into the old support sector supplies a common strict ambient center. Acute extension then retains the old facets; the new facets are actual by perturbing relative interior points of their polygon sides off the support. The previously cut escape stays cut because `r95` remains.

This proves a genuine local compact-face completion mechanism beyond the bounded same-face cap bank. It does not certify that a small polygon suffices, does not retain temporary construction sides, and does not settle the other faces, cusps, or global volume of the ambient H20 chamber.
