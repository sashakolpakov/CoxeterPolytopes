# A parameterized explanation of the two D9 obstructions

The two verified transports have the same explicit form. This gives a **conditional structural test**, not a theorem that every coefficient fails and not a reason to continue scanning small coefficients blindly. The [exact replay](../scripts/fresh_h9_tripled_d9_20261013_mechanism.py) reproduces both archived matrices from the formulas below; its [certificate](fresh_h9_tripled_d9_20261013_mechanism.json) passes for `a=T−1` and `a=3(T−1)`.

Suppose

```text
Z²−aX²=1,       W²−3aY²=1,
Z,W>1, X,Y>0,
Z,X,W,Y in O_K, (W+1)/3 in O_K.
```

Let `B_j` be the negative Lorentz boost on time and spatial coordinate j, with block `[[Z,−X],[−aX,Z]]`, for j=1,2,3. Let P reverse those three spatial coordinates. Define C on the time-plus-three-coordinate subspace by

```text
C = [ W       Y*(1,1,1)       ]
    [ aY*1    −I+(W+1)J/3    ],
```

where J is the three-by-three all-ones matrix. Extend every map by the identity on the other six coordinates. Then

```text
F = C P B_3 B_2 B_1.
```

The two Pell identities prove that C and the boosts preserve `diag(−a,1,1,1)`. C acts as minus the identity on the spatial sum-zero plane and as a Lorentz boost on its orthogonal complement. Thus `det(C)=1`, `det(F)=−1`, and F preserves the future cone. The displayed integrality assumptions make F and its ambient-coordinate inverse integral.

This transport maps the structured three-step axis corner to the corner whose new simple root is

```text
(2Y, (2W+2)/3, (2W+2)/3, (2W−4)/3, 0,...,0).
```

For an actual full-chamber symmetry, these must be the **complete first corners**, with the oriented sectors connected through the initial B9 corner. A merely compatible Pell cap does not suffice. The two completed hosts establish that hypothesis through their independent all-norm, all-profile and all-unit-orbit censuses.

## Four exact lattice congruences

Since F fixes the last six spatial coordinates, preservation of the D9 sum condition reduces to four congruences. Define

```text
c0 = a*(3Y Z³ − W X*(1+Z+Z²)),
c1 = W*(Z+aX²*(1+Z)) − 3aXY Z²,
c2 = W*(Z+aX²) − 3aXY Z,
c3 = W Z − 3aXY.
```

These are the coefficients of the sum of the first three output spatial coordinates. Therefore F preserves D9 exactly when

```text
c0 in 2O_K,   c1−1,c2−1,c3−1 in 2O_K.
```

The last six coefficients force the parity multiplier to be one. Conversely these congruences give an integral matrix in the D9 basis; determinant minus one then gives an integral inverse. Both concrete certificates satisfy all four conditions. This explains why checking only ambient-coordinate integrality would have missed the essential lattice issue.

Finally, checking `abs(trace(F))>10` proves infinite order. Combined with the complete-corner hypotheses, this gives nonreflectivity. The matrix factorization, four congruences and trace formula all replay exactly in the two completed examples. A new coefficient must satisfy these hypotheses afresh; in particular new prime norm classes could introduce earlier facets and destroy the corner pattern.

## A more discriminating next host

A concrete alternative is

```text
<−3(T−1)> + E8 + <1>.
```

The repository search found no completed conductor-16 E8 spatial-block lane. This proposal changes the zero-height chamber from B9 to **E8+A1**, so the three-coordinate B9 staircase above does not start. Total definiteness makes the zero group finite; the rank-nine finite reflection classification cannot enlarge a subgroup E8+A1 to an irreducible rank-nine group. This is a structural reason to test the host, rather than another scalar in the same D9 family.

Its unimodular spatial block and time coefficient `3(T−1)` give safe norms dividing six, hence ten dyadic/three classes. The edge obtained by dropping the A1 wall is fixed by E8. A first finite extension must be E8+A1, so the initial problem is an orthogonal binary norm equation in the already certified field `K(sqrt(3(T−1)))`. An E7+A1 edge then offers a small finite-type profile list for a genuine second step.

The determinant ideal differs from the tripled D9 lattice, but its rational determinant squareclass is the same. **Integral overclosure and full-group duplication checks remain necessary before calling this new reflection-group coverage.** The suggestion is a tractable structural pivot, not an executed search or a claimed new example.

The subsequent [exact structural triage](fresh_h9_tripled_d9_20261013_e8_triage.md) proves that the displayed maximal lattices and reflection sets differ and that the saved D9 transport fails E8 preservation. It does not settle possible conjugacy by another indefinite isometry. No expensive edge search has begun.

The next round completed [two exact first edges](fresh_h9_e8_first_20261013.md) while preserving that conjugacy caveat. Its additional mirror also proves that the standard maximal overclosure has extra mirrors inside the closed D9 chamber, so saturation alone does not transfer nonreflectivity.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_mechanism.py --verify
```
