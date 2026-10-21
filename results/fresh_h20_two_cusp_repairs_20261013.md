# Two simultaneous angle repairs in H20 sections

Nine explicit sections of the Borcherds H21 chamber have finite volume and two prescribed **actual** Coxeter ridges each. All nine also have an exact actual forbidden ridge, so none is a Coxeter polytope. Seven of their normals require at least five of the original wall normals to span them. This is a new construction family within this repository, not a new Coxeter example or an exclusion of all H20 sections.

The [independent replay](../scripts/fresh_h20_two_cusp_repairs_verify_20261013.py) uses only the Python standard library and earlier Fraction verification helpers; it imports neither the numerical scout nor the exact producer. Its [result](../data/enumerations/fresh_h20_two_cusp_repairs_verify_20261013.json) is PASS: nine strict section interiors, 18 actual repaired ridges, nine actual forbidden ridges, seven support exclusions, and rejection of forged ridge-point and target-angle data.

## A genuinely two-parameter section family

Take the literal cusp rays numbered 0, 3 and 4 in the [earlier exact cusp archive](../data/enumerations/h20_gluing_independent_cusps_20261013.json), with indices starting at zero. Their Gram matrix and pairings with the original strict center h are

```
       u0   u3   u4
u0      0  -24   -2
u3    -24    0  -14
u4     -2  -14    0

(h,u0) = -30, (h,u3) = -30, (h,u4) = -18.
```

Use E(x,y) = 7u0 − x u3 − 12y u4, with independent positive x and y. Its square is 336(x+y−xy). Seventeen original walls are orthogonal to the whole three-dimensional normal span: 1, 28, 29, 30, 31, 32, 35, 39, 40, 41, 46, 47, 49, 50, 51, 64 and 210, using one-based wall indices. Their Gram rank is 17.

For an orthogonal wall i and another wall j, write v_j=(E,r_j) and q=E². The projected squared cosine is

```
G_ij² / (G_ii (G_jj − v_j²/q)).
```

Prescribing a label m gives v_j² = delta_j q, where delta_j = G_jj − G_ij²/(G_ii cos²(pi/m)). Two angle equations can be solved together: v_j = ±sqrt(delta_j/delta_k) v_k reduces the system to a line and a quadratic; a zero delta gives a linear equation directly. This uses two independent freedoms, rather than increasing a cutoff in an old one-parameter family.

The numerical [scout](../scripts/fresh_h20_two_cusp_repairs_20261013.py) used labels 3, 4, 5 and 6 and five specified interior probes. Its [archive](../data/enumerations/fresh_h20_two_cusp_repairs_20261013_scout.json) records 43 profile types, 11,427 two-equation templates, 5,485 distinct numerical roots with positive parameters, 5,057 numerical strict-section tests, and 61 proposals passing both target-ridge scouts. These are numerical proposal counts, not a completeness theorem. The exact batch consists of the first rational control and all proposals involving a target label 5 or 6.

## Exact selected sections

All wall indices below are one-based. “Acute forbidden” is certified by a positive definite projected two-wall Gram matrix, a strict timelike ridge point, and exact exclusion of every Coxeter angle; “obtuse” is also certified at an actual ridge.

| Scout ID | x | y | Two actual target ridges | Actual forbidden ridge |
|---|---|---|---|---|
| 0 | 35/11 | 28/33 | 1/185: pi/4; 29/180: pi/4 | 2/4, obtuse |
| 31 | (791+112sqrt(19))/93 | (91−14sqrt(19))/186 | 210/194: pi/3; 47/202: pi/6 | 1/185, acute forbidden |
| 41 | 63−28sqrt(5) | 7/6 | 210/202: pi/3; 210/34: pi/5 | 1/199, acute forbidden |
| 42 | (133+56sqrt(7))/29 | (119−14sqrt(7))/174 | 47/202: pi/6; 32/9: pi/3 | 1/185, acute forbidden |
| 43 | 7/4+21sqrt(5)/8 | 7/9 | 47/202: pi/5; 30/27: pi/3 | 2/202, obtuse |
| 44 | 651/61+1008sqrt(5)/305 | 105/122−84sqrt(5)/305 | 47/202: pi/5; 39/6: pi/3 | 2/202, obtuse |
| 45 | 7 | 77/18−14sqrt(5)/9 | 47/202: pi/5; 64/18: pi/3 | 2/202, obtuse |
| 46 | 7 | 35/54 | 47/202: pi/6; 64/18: pi/3 | 2/4, obtuse |
| 60 | quartic root described below | 7/12−x/4 | 49/5: pi/3; 46/204: pi/5 | 3/8, acute forbidden |

For case 60, x is the unique root in the archived rational isolating interval

```
4237530754/1971844737 < x < 1336721595/622014941
```

of 396x⁴−84x³−490x²−1372x−2401. Numerically x is about 2.14901846706605; this approximation is not used in verification.

The [exact producer](../scripts/fresh_h20_two_cusp_repairs_exact_20261013.py) saves [algebraic certificates](../data/enumerations/fresh_h20_two_cusp_repairs_exact_20261013.json) for every strict interior and every claimed ridge. The replay checks the defining polynomial, real-root interval, exact section orthogonality, future timelike norm, all wall inequalities, precise zero sets or explicit redundant-zero cone identities, projected Gram positivity, and target or forbidden angle identities. Thus these are actual dihedral-ridge statements, not a screen of arbitrary wall pairs.

Finite volume follows by taking the section in the Klein model. The ambient chamber is the convex hull of finitely many finite and ideal vertices. A full-dimensional hyperplane section has vertices inherited from those vertices or from intersections with edges; an intersection in the interior of an edge is strictly inside the ball. Its hyperbolic convex hull therefore has finite volume. The strict interior certificates establish full dimension for all nine chosen sections. This batch makes no compactness claim: it does not classify their complete ideal-vertex sets.

## Difference from earlier section families

The two rational controls have no minimum-support claim. For each of the seven irrational sections, the replay extracts independent rational normal vectors e0,e1 and an irrational parameter alpha with E proportional to e0+alpha e1. A rational wall span containing this normal must contain both e0 and e1.

The complete four-wall support filter then works modulo 1009. The two normal vectors retain rank two. Rational containment in four walls would force the combined rational rank to be at most four; reduction cannot increase that rank, so the four quotient wall images would have rank at most two. The filter includes quotient ranks zero, one and two; there are no surviving four-wall sets for any of these seven planes. This proves minimum original-wall support at least five, without assuming the four wall vectors preserve rank modulo 1009.

All seven rational normal planes have nonzero Gram determinant (respectively 7007, −784, 3920, −12544/9, −3136, −1016064 and 15680). They therefore also differ from the earlier cusp-translation planes spanned by a positive vector e and an orthogonal null cusp ray u, whose Gram determinant is zero. This distinguishes literal construction mechanisms, not isometry classes of polytopes.

## Initial six-chamber gluing attempt

The same exact archive records an initial attempt to remove the old persistent ridge 1/22 by taking the full six-chamber star generated by their A2 reflections. The rebuilt union has 416 distinct supporting boundary planes and is convex with finite volume. This is verified by reflection-orbit enumeration and nonnegative combinations of the two seam normals, together with folding into their fundamental chamber.

For every section normal E=r60+s r99−t r125 with s,t>0, however, that union has an actual ridge with boundary normals r24 and r24+r1. Their Gram matrix is [[2,1],[1,2]], giving angle 2pi/3, and both are orthogonal to the whole section-normal span. Two archived future timelike face points have opposite section-normal signs for all positive s,t; positive interpolation places a point of this ridge in every such section, strictly below every other boundary plane. The independent replay checks all these claims. This rules out this particular six-star throughout that positive A2 family, not arbitrary ambient unions.

The cusp incidence matters here: wall 1 is incident to the chosen A11×D7×A2 cusp, while wall 22 is not. Treating both seams as nonincident would give an incorrect cusp-preservation argument.

## Reproduction and scope

From the repository root:

```sh
python scripts/fresh_h20_two_cusp_repairs_verify_20261013.py
```

To regenerate the proposal and certificate archives:

```sh
python scripts/fresh_h20_two_cusp_repairs_20261013.py
python scripts/fresh_h20_two_cusp_repairs_exact_20261013.py
python scripts/fresh_h20_two_cusp_repairs_verify_20261013.py
```

The verified result is a bounded rejection of nine new finite-volume sections after two successful simultaneous angle repairs, plus one all-positive-parameter obstruction for a specified six-chamber union. It leaves other parameters, other triples of cusp rays, other ambient unions and the general H20 construction problem open.
