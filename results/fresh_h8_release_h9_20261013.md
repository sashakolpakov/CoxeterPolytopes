# Releasing an H8 wall and constructing algebraic H9 caps

**Outcome:** deleting one old wall removes the previous affine obstruction and gives a genuine rank-ten H9 seed. A targeted search constructs 452 exact Coxeter-compatible individual caps, including sixteen caps with a freely solved dotted distance. Some use quartic fields and cut the earlier common escape. Nevertheless, exact conflicts and sixteen strict spacelike witnesses prove that **no subset of this archived cap bank completes the ten retained seed walls to a finite-volume Coxeter polytope**. This is a finite-bank result, not a general completion obstruction or a new-dimensional example.

The [independent verifier](../scripts/fresh_h8_release_h9_independent_verify_20261013.py) uses only standard-library `Fraction` arithmetic. It imports neither search code nor numerical or symbolic packages. Its [result](../data/enumerations/fresh_h8_release_h9_independent_verify_20261013.json) verifies all accepted normals, rank and affine gates, cap conflicts, and escape witnesses. Replay:

```sh
python scripts/fresh_h8_release_h9_independent_verify_20261013.py --verify
```

## The released seed

Start with the compact eleven-wall Bugaenko H8 example. In the earlier single-coordinate lift, append a positive coordinate of square `2φ`, lift old root 0 in that coordinate, and divide that whole root by φ. This changes its incident label from 5 to 3 and raises the full Gram rank to ten. Keeping all eleven walls forces an affine E7 null ray; adding acute-compatible caps cannot make that system compact.

Now **delete original wall 4**. The remaining ten roots are independent and have signature `(9,1)`. Their diagram is a chain in the original-index order

```
0 — 1 — 2 — 3 — 9 ⋯ 10 — 5 — 6 — 7 — 8
```

All ordinary edges have label 3 except the final 7–8 edge, which has label 5. The dotted 9–10 pairing is `−2φ`; every root has square 2. This changes the retained facet set, so the earlier eleven-wall affine obstruction no longer applies.

The independent checker examines deleting each of the eleven walls **after each of the two previously classified fixed lifts**. Precisely four cases pass the rank-ten and no-affine-subdiagram gates: deletions 4, 5, 6 after the 5→3 lift, and deletion 4 after the 5→4 lift. These necessary gates do not assert compactness. They also do not classify arbitrary lifts made after walls have first been deleted.

The cap experiment uses the connected 5→3, delete-4 seed. If `G` is its Gram matrix, `(G⁻¹)₀₀=φ/2`. Its natural polar cut at the opposite hyperideal vertex would have squared cosine `φ−1` with root 0. This lies strictly between `1/2` and `(1+φ)/4`, hence is not a Coxeter angle. The experiment therefore searches nonorthogonal caps instead of relying on polar truncation.

## Exact cap constructions

For ten independent seed normals, prescribe the new cap pairings `p_i=(s,r_i)≤0`. Its seed coefficients are `G⁻¹p`, and its square is exactly `pᵀG⁻¹p`. Thus each candidate has a direct rank-ten realization; the norm and angle checks are independent of any lattice-reflection enumeration.

| Cap construction | Accepted caps | Exact effect |
|---|---:|---|
| Square 2, at most seven nonzero pairings from `−{1,φ,2,1+φ,2φ}` | 434 | All share a spacelike escape |
| Mixed pairings involving `√2` or `√(2+φ)` | 2 | New labels 4 and 10; old escape survives |
| One freely solved dotted pairing, up to four other nonzero pairings from `−{1,φ}` | 16 | Each cuts the old escape |

The first bank also screened square-1 proposals and found no additional cap. A parallel `√3` test found none. Proposal filters use floating point and are explicitly bounded; exact arithmetic checks every accepted cap. No claim of exhaustive enumeration is made for these broader coefficient or field families.

The initial all-434-cap escape has squared norm

\[
\frac{2908}{8649}-\frac{584}{2883}\phi>0.
\]

The free-distance construction addresses that particular gap. One cap is orthogonal to eight seed walls, has label 5 with original wall 2, and dotted pairing `−x` with original wall 10. Its norm-two equation is

\[
(-39/2+12\phi)x^2+(-15+9\phi)x+(-7/2+3\phi)=0,
\qquad x\approx2.183321849318035.
\]

The discriminant is `−111+69φ>0`, with field norm `−99`, so it is not a square in `Q(φ)`. This realizes a cap over a genuine quartic extension. It is a new normal in this experiment, not a claim of a new finite-volume reflection group or commensurability class.

## Why the finite bank cannot close

All 120 pairs among the sixteen free-distance caps conflict. In 103 pairs the outward pairing is positive. For each of the other seventeen pairs, exact arithmetic in the compositum proves

\[
-\phi<(s_i,s_j)<0,\qquad (s_i,s_j)^2\notin\{1,2\}.
\]

Both norms are 2. The interval permits only Coxeter labels 2, 3, 4; the strict pairing excludes 2, and the square checks exclude 3 and 4. This argument does **not** bound the compositum degree or omit higher-degree angles: every label at least 5 is excluded by the exact interval itself. Hence a Coxeter completion can have at most one of these sixteen normals as an actual facet.

For each free-distance cap, the verifier retains an **overinclusive** set of potentially compatible caps from the 434-wall bank, ranging from 44 to 273. A cap-versus-old-cap pairing belongs to a field of degree dividing four, so the possible finite labels are

```
2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30.
```

The test retains every root of the associated squared-cosine minimal polynomials, including conjugate values that need not be actual Coxeter angles. This can retain extra inequalities; it cannot reject a valid angle. The two additional quadratic-field caps are retained unconditionally in every geometric test, so their compositum degrees never enter a rejection gate.

Each of these sixteen supersystems has an archived strict spacelike witness. Consequently every subsystem also has infinite volume. If no free-distance cap is selected, the original common escape survives the 434 base caps and both quadratic caps. Together with the pair-conflict proof, this covers every subset of the 452 archived cap normals.

Possible redundant caps do not invalidate the argument: discard them and apply it to the actual facet normals. The ten seed walls themselves remain true facets. A common strict timelike point exists, and projecting it orthogonally to any seed wall makes every other inequality stricter because every new cap pairs nonpositively with that seed wall.

## Exact volume gate and artifacts

Use coordinates `y_i=−(x,r_i)≥0`. The normalization `Σy_i=1` is a bounded simplex, cut by the cap inequalities `(G⁻¹p)·y≥0`; the Lorentz square is `yᵀG⁻¹y`. This supplies an explicit bounded affine section if a future compact candidate requires a complete vertex census. No compactness claim here is based only on sampled vertices.

For every negative result, the stored witness has all coordinates and cap constraint values strictly positive, positive Lorentz square, and negative pairing with the strict timelike center. The segment from that center crosses the future light cone inside every inequality, giving an open ideal-boundary patch and therefore infinite volume. LPs only propose points; rational reconstruction and exact sign checks certify them. The independent replay needs no LP solver.

Producer and certificate files:

- [Pairing-bank search](../scripts/fresh_h8_release_h9_caps_20261013.py), [434-cap archive](../data/enumerations/fresh_h8_release_h9_caps_expanded_20261013.json).
- [Common escape gate](../scripts/fresh_h8_release_h9_volume_20261013.py), [exact common escape](../data/enumerations/fresh_h8_release_h9_volume_expanded_20261013.json).
- [Quadratic-angle search](../scripts/fresh_h8_release_h9_quadratic_caps_20261013.py), archives with suffixes `_2_`, `_3_`, `_2phi_` in `data/enumerations/`.
- [Free-distance construction](../scripts/fresh_h8_release_h9_free_caps_20261013.py), [sixteen caps, 120 conflicts and sixteen escape certificates](../data/enumerations/fresh_h8_release_h9_free_caps_20261013.json).

Other released wall sets, altered old projected directions, multiple free dotted parameters, other labels and other cap normals remain outside this finite-bank result.
