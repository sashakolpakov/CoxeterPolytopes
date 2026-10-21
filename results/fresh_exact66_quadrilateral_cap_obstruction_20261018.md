# Exact obstruction to one quadrilateral cap on the 66-wall seed

18 October 2026. **Source and saved certificate complete; independent source/report review pending.**

A single new wall cannot satisfy the four finite-corner duties below, for any of the seven retained white profiles and any finite Coxeter labels. This is a fixed-parent necessary-norm obstruction. It does not exclude a split cap, changed incidence, moving a parent wall, or another white framework. It does not assert that the infinite-volume seed has no finite-volume completion.

The unchanged seed is [the exact 66-wall intersection](fresh_exact65_escape33_cap_20261018.md), with 205 finite vertices and a certified strict spacelike escape. There are no wall deletions, no lost finite vertices, and no added wall in this packet. The inherited compact old60 triangle is a face of the fixed-six-white hyperbolic three-dimensional section; no ambient facet is asserted compact.

## Geometric question

The remaining exterior regions at 30/52/76 and 47/52/76 suggest closing old52 and wall76 simultaneously. The proposed new wall K must make all four genuine finite corners

- 30/52/K;
- 52/47/K;
- 47/76/K;
- 76/30/K.

The fixed cyclic frame is `(30,52,47,76)`, with edge labels `(4,2,2,2)`. Both diagonals, 30/47 and 52/76, are dotted. No finite duty involving either diagonal is imposed. The four proposed corners would be the fan on this cyclic boundary; this is only a proposed local incidence, not a verified chamber boundary or cone.

The original white roots are `W=(0,1,2,3,5,6)`, of type H4+A2. Each frame wall is orthogonal to W. The seven inherited profiles for K are profile0 and a single bond of label 3, 4, or 5 to either endpoint of the A2. This seven-profile restriction is retained as an explicit hypothesis, not reproved or enlarged here.

## Physical Gram matrix and necessary norm equation

The producer reconstructs the frame products from the exact original physical root intervals and the frozen exact wall76 products. Original native products are isolated between adjacent native label values before intersecting with the corresponding exact identity. This uses the previously certified Coxeter nature of those original walls. It does not infer exact labels from floating approximations of a new wall.

All frame white covectors vanish. Thus the physical and projected frame Gram matrices agree. Writing `a=−(30·47)>2` and `b=−(52·76)>2`, the exact matrix is

\[
G=\begin{pmatrix}
2&-\sqrt2&-a&0\\
-\sqrt2&2&0&-b\\
-a&0&2&0\\
0&-b&0&2
\end{pmatrix}.
\]

Directed rational bounds give `−5.281 < det G < −5.280`. In the inherited projected space of signature (3,1), nonsingularity makes these four walls a basis. Consequently their four prescribed products determine a projected normal uniquely; there is no omitted second norm branch or free perpendicular parameter.

Put `theta_m=2 cos(pi/m)` and let `c=(-theta_m0,...,-theta_m3)`. The projected square must be `c^T G^{-1} c`. For profile0 its required value is 2. For a nonzero profile with label `l`, it is

\[
q_l=2-\frac23\theta_l^2.
\]

The coefficient 2/3 is the endpoint diagonal of the inverse A2 Gram matrix. Thus every candidate must satisfy the exact scalar equation

\[
c^TG^{-1}c-q_l=0.\tag{1}
\]

The cap fails this necessary equation before any common-center, retained-wall, or full-cone tests are needed.

## Complete finite-label domain

For labels 2 through 5, the complete seven-profile four-label cube has 1,792 recipes. The finite-corner classifier is a combinatorial finite-diagram test on A2 plus the three indicated black roots; the independent checker builds the catalogue explicitly. Its connected components are among A1–A5, B2–B5, D4/D5, F4, H3/H4, and I2(5). The fixed disconnected H4 factor is already finite.

Exactly 1,708 recipes fail at least one finite duty. For each of the remaining 84, rational intervals strictly separate the left side of (1) from zero. The nearest residual occurs at profile0 and labels `(2,3,2,2)`, and lies strictly between `−0.019744` and `−0.019743`. The full certificate retains every rejected recipe and all four duty classifications, not just survivors.

The infinite high-label domain also reduces completely. In a finite Coxeter diagram an edge labelled at least 6 is an isolated I2 factor. The existing 30/52 edge labelled4 therefore forces K/30 and K/52 to have labels 2 or 3, with at least one of these labels2. A high label can occur only at K/47 or K/76. Isolation in the two adjacent finite corners forces K to have profile0 and forces both neighboring K contacts right. The only four symbolic families are

\[
(2,2,m,2),\quad(3,2,m,2),\quad
(2,2,2,m),\quad(2,3,2,m),\qquad m\ge6.
\]

There cannot be two high contacts. The standalone checker independently confirms these four forms using an isolated-high-edge category, rather than sampling large labels.

In each family the left side of (1) is `A theta^2+B theta+C`, where `theta=2 cos(pi/m)` lies in the complete closed overinterval `[sqrt(3),2]`. Directed intervals prove `A>0`, `B>=0`, and a positive value and derivative at `sqrt(3)`. The lower bounds for the four starting values are respectively `0.863`, `2.097`, `6.213`, and `15.485`. Hence every polynomial is strictly positive throughout the closed interval. This covers every finite integer `m>=6`, without a label cutoff, and also covers the unused endpoint limit2 as an overdomain.

## Reproduction and independence

The producer [fresh_exact66_quadrilateral_cap_obstruction_20261018.py](../scripts/fresh_exact66_quadrilateral_cap_obstruction_20261018.py) binds the exact64, exact65, and exact66 source/data/report triples. It saves the reconstructed Gram matrix, inverse, all low recipes, all strict norm residuals, and all four high-family polynomials in the matching JSON.

The standalone [check_exact66_quadrilateral_cap_saved_20261018.py](../scripts/check_exact66_quadrilateral_cap_saved_20261018.py) imports no project modules. It reconstructs the physical frame and native-isolation intersections, rebuilds the finite catalogue and the complete low/high domain, and checks all saved arithmetic. Its directed rational arithmetic follows the earlier standalone-checker design. In addition to the adjugate reconstruction, a separate Gaussian inverse gives the same strict norm-exclusion signs and signature. This is independent verification of these saved identities and domains, not a fresh proof of the original scaffold, its seven-profile theorem, or the original finite-vertex census.

Run from the repository root:

```sh
python scripts/fresh_exact66_quadrilateral_cap_obstruction_20261018.py --verify
python scripts/check_exact66_quadrilateral_cap_saved_20261018.py
```

The final generation `52172/43d5d5` and standalone saved check `d4f5d2/8c3af8` exited zero and passed. The separate source-bound producer replay `28810/dd924d` also exited zero and passed. The initial successful generation `45564/999c3b` preceded only a scope-key correction: ideal or dotted corner duties are untested, while right finite labels are included.

Root completed a full mathematical producer/report/standalone-helper read. The fixed projected basis and uniquely necessary norm equation, seven-profile domain, all1,792 low-label finite-corner tests, four symbolic high-label families with strict monotonicity on the closed interval, and the directed adjugate/Gaussian verification paths are coherent. The reviewer did not execute a program or reprove the original scaffold or seven-profile theorem. This packet is final/frozen after those author-owned generation, replay and saved-checker outcomes and the completed independent source/report review.

No numerical cone is generated because no norm-two cap survives. The result changes the next construction choice: the requested four finite duties cannot be met by one wall in this framework, so a future attempt must split that cap or alter a parent/incidence. No new finite-volume example, compact ambient polytope, or global completion impossibility follows.
