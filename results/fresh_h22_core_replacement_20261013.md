# An exact one-wall replacement in the H22 incidence chamber

Deleting one original wall permits two new Coxeter walls with **spacelike old projections**. They cut an inherited cusp and replace it with an actual finite vertex. The resulting chamber has **813 actual facets in H22**, but a strict spacelike escape proves that its volume is still infinite. This is a changed-core construction, not a finite-volume example or a claim of literature novelty.

The inputs are the 812-facet incidence branch: 210 vertical walls of the finite-volume H21 source and 301 null-projection caps at each of two opposite heights. This construction does **not** retain the 48 additional walls of the later 860-facet branch.

## Exact construction

Use the old rational form `Q` of signature `(21,1)` and the H22 form `diag(Q,2)`. The old strict center is `h`, with `h²=-14`. Cusp orbit 14 in the [complete source census](fresh_h21_source_corner_walk_20261013.md) has a single affine component on 21 walls; its previously classified type is `D20`. Let `δ` be its saved primitive future null vector. Its center pairing is `h·δ=-46`.

Delete old wall 1. Let `I` be the remaining 20 incident walls, and let `p` be the orthogonal projection of `h` to `span(I)⊥`. The standalone checker reconstructs this projection and proves

```
p² = -2484,       p·δ = -46,
u = -p + 350δ,    U = u² = 29716,
r± = (u/sqrt(U), ±1/sqrt(2)).
```

Thus both new normals have square 2 and their mutual pairing is zero. Each is right to the 20 walls in `I`, and strictly dotted to every other retained wall. The complete new-pair census is **41 right and 1,582 dotted pairs**. All 812 original normals and all their pair relations are independently rebuilt in the same standard-library checker before the deletion.

Both new walls cut the old cusp because

```
(r±,(δ,0)) = 46/sqrt(29716) > 0.
```

Their old projections have square 1. Consequently the previous causal-projection and ball-covering arguments that assumed retention of all 210 original facets do not apply to this branch.

The common center `(h,0)` remains strict: its new pairings are `-13616/sqrt(29716)`. Complete nonpositive Coxeter pairings and a strict timelike center prove that all 813 inequalities are true facets: project that center to any one wall; acuteness makes every other inequality strict there.

## An actual finite vertex and an exact remaining escape

The vector

```
v = p + 296δ
```

satisfies `v²=-29716`, `v·h=-16100`, and `v·u=0`. At `(v,0)` exactly the 20 walls in `I` and the two new walls vanish. Every other retained inequality is strict. The incident Gram is the positive-definite 20-wall Gram direct-summed with `2I₂`, so this is an actual finite H22 vertex, not merely a positive-definite subdiagram.

The old integer escape `(100δ₀+h,40)`, where `δ₀` is cusp representative 0, still satisfies all 813 inequalities strictly. Its square is **1586**, its center pairing is **−814**, and its largest retained old pairing is **−1**. For the new caps the old-coordinate numerator is `−692116`; the rational squared margin

```
692116² - 2*40²*29716 = 478929466256 > 0
```

proves both full pairings are strictly negative. Hence an open portion of the ideal boundary remains exposed and the chamber has infinite volume. No global vertex or cusp census of the changed chamber is claimed.

## Why the replacement mechanism works

This construction is not a lucky finite parameter. For any sufficiently large `t`, set

```
u(t) = -p + tδ,
U(t) = p² - 2t(p·δ),
v(t) = p + (t - p²/(p·δ))δ.
```

Then `u(t)·δ=46`, `U(t)>0`, `v(t)²=-U(t)`, and `u(t)·v(t)=0`. Every retained old wall outside `I`, and every one of the 301 parent null projections, has strictly negative pairing with `δ`. Its numerator against `u(t)` therefore tends linearly to negative infinity, while `sqrt(U(t))` grows only as `sqrt(t)`. All those new relations are eventually dotted; the finite-vertex nonincident inequalities are eventually strict as well. The parameter 350 is a convenient exact choice, with no optimality claim.

Deleting the one affine wall is essential to this mechanism. More generally, if an entire affine component at an old cusp is retained, its positive null relation forces every acute-compatible new projection to pair nonpositively with that cusp. To cut that cusp, a deletion set must meet every irreducible affine component. Orbit 14 has only one component, so a single deletion removes this particular obstruction.

The [supplementary plane experiments](fresh_h22_core_replacement_20261013_planes.md) archive an earlier exact five-wall deletion and two all-real, all-label cusp-cutting obstructions. The present construction changes the old projection plane and escapes those scoped obstructions.

## Arithmetic scope

This is deliberately outside the noncompact arithmetic search excluded above dimension 21 by [Esselmann's original theorem](https://www.sciencedirect.com/science/article/pii/S0022314X96901419). That theorem concerns integral Lorentzian lattices of maximal rank 22; it is not the general dimension-29 bound for all arithmetic reflection groups, including compact ones over other totally real fields.

The new norm-two cap and parent wall 210 have length-two cyclic product

```
12645774/7429 - (3554/7429)*sqrt(59432).
```

It is irrational. The 41 retained rational short walls have Gram value **−118** on their all-ones vector. Any finite-volume completion retaining these facets therefore fails the conjugate-semidefiniteness condition for quasi-arithmeticity: a nontrivial embedding on the displayed quadratic cyclic field leaves that rational indefinite block unchanged. This is the usual Vinberg criterion, stated with the quasi-arithmetic variant in [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf). The conclusion is conditional on a cofinite completion retaining these facets; the present infinite-volume chamber is not being described as a lattice.

## Reproduction

- [Standalone exact checker](../scripts/fresh_h22_core_replacement_20261013.py)
- [Frozen certificate](fresh_h22_core_replacement_20261013.json)
- [Independent source and sparse-radical arithmetic review](fresh_h22_core_replacement_20261013_review.md)

```
python scripts/fresh_h22_core_replacement_20261013.py --verify
```

The producer and a fresh deterministic replay both pass. The checker uses only Python integers and `Fraction`, reconstructs the complete parent input, binds all source hashes, and imports no producer or shared project arithmetic. Radical inequalities are verified by sign-controlled rational squaring. The old cusp census supplies the named representative, but all geometric properties needed for this construction are checked directly. A separate sparse-radical checker also passes every new pairing, norm, center, cusp cut, vertex, escape, and field witness; it does not duplicate the full old 812-wall pair census.
