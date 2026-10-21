# Noncentral four-wall sections with three angle constraints

This batch constructs **21 exact finite-volume H20 sections** of the frozen Borcherds H21 chamber. Seventeen have **three certified Coxeter ridge angles** after imposing two independent auxiliary orthogonality conditions and solving a remaining angle equation. All 21 have a different **actual forbidden ridge**, so none is a Coxeter polytope.

The construction moves beyond the earlier three-wall searches: fifteen accepted section normals are irrational and are proved not to lie in the span of **any three original walls**. The result is a bounded construction experiment, not a four-wall impossibility theorem or a new Coxeter example.

## Construction and declared scope

The input is the same [210-wall Borcherds chamber](../data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json) used in the earlier [mixed gluing](fresh_h20_mixed_20261013.md) and [nonelliptic three-wall sections](fresh_h20_bivariate_nonelliptic_20261013.md). Its wall normals, Lorentz form and strict center are denoted `r_i`, `Q` and `h`. All wall labels here are one-based.

For four independent walls `I`, choose two auxiliary walls `r_a,r_b` and solve

```
e ∈ span{r_i : i∈I},    (e,r_a)=(e,r_b)=0.
```

The two independent equations leave a rational normal plane with basis `e₀,e₁`. The sections use `e(t)=e₀+t e₁`; the selected planes are not identically central, and central or lower-support special parameters are excluded from this noncentral batch. Every selected four-wall Gram has signature `(3,1)`. Thus the ambient elliptic-triple face obstruction from the preceding round does not apply directly.

For projected wall normals the scaled Gram entries are

```
g_ij(t) = q(t) (r_i,r_j) − (e(t),r_i)(e(t),r_j),
q(t) = (e(t),e(t)).
```

The two orthogonalities preserve right angles between each auxiliary wall and its originally orthogonal neighbors. A further chosen pair must satisfy `g_ab²/(g_aa g_bb)=cos²(π/m)` for `m=2,3,4,5,6`. These equations are solved exactly, including every real parameter root. The wrong conjugate of the label-five equation is rejected. At the projective parameter infinity, the chosen basis vector uses at most three listed walls, so no full four-coordinate direction is lost by this chart.

The numerical proposal stage examined 167 literal support representatives distinguished by weighted Gram matrices and 2,551 full four-coordinate normal planes. It retained 195 promising pencils and selected twelve. This is **not** a symmetry-orbit classification: equal support Grams were used only to choose proposals. Raw complete-Gram conflict counts ranked proposals; they are not reported as actual ridge counts. Two proposed planes that were contained in other three-wall spans were removed before the final selection.

## Exact outcome

| Gate | Count |
|---|---:|
| Declared pencils | 12 |
| Distinct real roots of their angle equations | 31 |
| Lower-support parameters | 5 |
| Central parameter | 1 |
| Wrong label-five conjugates | 4 |
| Noncentral, finite-volume section candidates | 21 |
| Candidates with all three repaired ridges certified | 17 |
| Candidates with another actual forbidden ridge | 21 |
| Irrational normals proved outside every three-wall span | 15 |

The seventeen full repair certificates consist of fifteen sections with angles `π/2,π/2,π/3`, and two with three right angles. Four other sections satisfy the algebraic target equations, but the target pair was not certified as an actual ridge; they are not counted as three-angle repairs. Their separate actual forbidden-ridge witnesses still reject them. Among the 21 decisive witnesses, eleven are obtuse and ten are acute with a squared cosine below `3/4` that differs from every Coxeter value in that interval.

The parameter polynomials have degrees 1, 2, 3, 4 and 8. These are candidate constructions in algebraic fields, not a rational-only screen. No unverified surviving candidate remains in this literal root list.

For a concrete quartic example, take

```
e(t)=−r₁+r₄₈+t(−3r₁+r₆₁+r₁₁₅),
182t⁴+152t³+57t²−2t−5=0,
−407/1000 < t < −406/1000.
```

The selected root is approximately `−0.40689566551243816`. The normal is orthogonal to walls 76 and 62; actual ridges 76/2 and 62/4 have angle `π/2`, and actual ridge 7/50 has angle `π/3`. But actual ridge 1/4 is obtuse. The archive gives exact field-valued timelike points for all four ridges and a strict timelike interior point. This irrational normal has minimum original-wall support at least four.

## Why the geometry and support claims are exact

Most interiors are the exact projection `C=q h−(h,e)e`. For two candidates this particular point failed an inequality. A linear program proposed a different point, which was rationally approximated and projected exactly into `e⊥`; exact signs then proved its timelike square, future orientation and all 210 strict inequalities. LP success alone is never a certificate.

Finite volume is inherited by these genuine interior sections. In the Klein model, the original finite-volume chamber is the convex hull of finitely many finite or ideal vertices. A hyperplane section again has finitely many vertices. A new point on an edge between distinct original vertices lies inside the ball unless it is an inherited ideal endpoint. Triangulation into finitely many finite/ideal hyperbolic simplices proves finite H20 volume.

Each ridge certificate supplies a timelike point in `e⊥`, two independent projected normals with positive definite Gram, and nonpositive pairings with all 210 walls. Every extra vanishing inequality is checked to be a nonnegative combination of the two ridge normals. All other inequalities are strict. This proves an actual codimension-two face; it does not infer incidence from a pair Gram alone.

The hidden-support check is also complete for the asserted fifteen irrational examples. If a rational three-wall span contains `e₀+t e₁` with irrational `t`, it contains both `e₀,e₁`. Quotient all 210 rational wall vectors by their plane `E`. Any three-wall span containing `E` has quotient dimension at most one, so its nonzero quotient vectors must be proportional. The checker groups those projective vectors, identifying both signs, includes all zero images, and tests every possible triple in those groups. No selected plane is contained in a three-wall span. The rational-root theorem separately proves the fifteen parameters irrational. The six rational specializations are not assigned this stronger minimum-support conclusion.

## Artifacts and replay

- [Producer](../scripts/fresh_h20_fourwall_sections_20261013.py) and [exact section archive](../data/enumerations/fresh_h20_fourwall_sections_20261013.json).
- [Independent Fraction verifier](../scripts/fresh_h20_fourwall_verify_20261013.py) and [PASS record](../data/enumerations/fresh_h20_fourwall_verify_20261013.json).

```sh
python scripts/fresh_h20_fourwall_sections_20261013.py --limit 12
python scripts/fresh_h20_fourwall_sections_20261013.py --replay-pencils
python scripts/fresh_h20_fourwall_verify_20261013.py
```

The second command retains the archived literal pencils and repeats their exact root and geometry calculations. The independent replay uses only the Python standard library and earlier standalone verification helpers; it imports no producer, SymPy, NumPy or LP implementation. It reconstructs the constraints, counts all 31 roots by Sturm sequences, checks all 21 interiors and forbidden ridges, checks the seventeen triple repairs and the fifteen minimum-support claims, and rejects forged ridge points and normal coefficients. Both scripts compile. A second agent reviewed the quotient-span argument and confirmed that grouping opposite projective directions is required; the implementation does so.

The remaining search space includes other four-wall normal planes, independent higher-support coefficients, other target pairs and labels, and different assemblies. The failed local repairs do not exclude those constructions.
