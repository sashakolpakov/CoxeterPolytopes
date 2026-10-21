# An exact noninducing H20 cap and a finite edge

The retained chamber has **43 actual facets**. A new norm-two cap cuts the saved escape of the direct A1²/H18 coupling and closes an actual finite edge. The new cap contains exactly one finite vertex, as certified by a complete census on that cap. A renewed exact strict spacelike escape proves that the resulting chamber still has **infinite volume**.

This is a geometric construction outside the preceding direct-inducing-profile obstruction. It retains the complete 40-wall H18 source and both ordinary A1 support walls. It does not assert a finite-volume H20 example or an arithmetic lattice.

The [exact producer](../scripts/fresh_h20_joint_route_20261015_exact.py) and [certificate](fresh_h20_joint_route_20261015_exact.json) passed generation and saved `--verify`. They depend on the frozen [A1²/H18 source certificate](fresh_h20_a1a1_h18_coupling_20261015.json) and the separate [Q(√2) arithmetic helper](../scripts/fresh_h20_b2_h18_coupling_20261014_review.py), and have no numerical-scout dependency. The [independent review](fresh_h20_joint_route_20261015_exact_review.md) also passed generation and saved replay: separate four-field Fraction arithmetic reproduced all 903 pairs, all 43 facet witnesses, all 1,540 supports using full matrix LDL, both edge endpoints and the strict escape, without producer or helper imports.

## Exact recipe and retained face

Use the source's physical 19-dimensional Lorentz form $Q_0$, and append two positive coordinates with form $2I_2$. Number the 40 source normals $0,\ldots,39$, the two support normals 40 and 41, and the new normal 42.

The vector $Y\in\mathbb Q(\sqrt2)^{19}$ is uniquely specified by zero products with

\[
1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,21,36
\]

and products $-1$ with source normals 12 and 18. Exact linear algebra gives $Y^2=-3/2$, with $Y$ future directed. Set

\[
R=(Y,-\sqrt3/2,-1),\qquad R^2=-3/2+3/2+2=2.
\]

The 42 new pair relations are 18 right, three π/3, one π/6, one parallel, and 19 dotted. All 903 chamber pairs are checked: 649 right, 30 π/3, six π/4, one π/6, 32 parallel, and 185 dotted. All scalar data use the explicit basis $1,\sqrt2,\sqrt3,\sqrt6$, with four rational strings per entry. No floating-point acceptance enters this certificate.

The frozen strict center remains timelike and strict at all 43 walls. For each wall $N_i$, its orthogonal projection

\[
C_i=C-\frac{(C,N_i)}{N_i^2}N_i
\]

is timelike, lies on that wall and is strict at the other 42 walls; the producer checks all 1,849 facet signs. Thus every proposed wall is an actual facet.

Since $Y$ is strictly future timelike, $(Y,x)<0$ for every nonzero future causal source vector $x$. Hence the new inequality is redundant on the entire source face defined by walls 40 and 41, including its ideal points. This **noninducing** behavior is precisely outside the direct inducing-wall families previously excluded.

## A finite edge and a complete census on the new cap

The new cap has exactly 22 finite-angle neighbors. All \(\binom{22}{19}=1540\) possible 20-root supports containing it are checked. Pairwise parallel/dotted supports are rejected; a nonright cycle cannot be positive definite; the remaining forests are tested by exact rational leaf-Schur pivots of their normalized Grams. There is exactly one positive support:

\[
S=\{1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,20,21,36,40,42\}.
\]

Replacing 42 by 41 gives another positive support $S_0$, belonging to the retained 42-wall base. Both physical vertex vectors are stored explicitly, normalized by $(C,v)=-1$. Their squares are negative; each has exactly its 20 stated zero products and is strict at every other wall. In particular, $R$ is strict at the old endpoint and wall 41 is strict at the new endpoint, proving distinctness.

The shared 19-root support $S\cap S_0$ is positive definite. Its feasible geodesic therefore has these two finite endpoints, giving an actual bounded edge. The census is complete for vertices incident to $R$; it is not a complete census of all ambient finite vertices, edges or cusps.

## Escape and arithmetic scope

The old direct-coupling escape has new-wall product

\[
-\frac{1252829}{1150000}+\frac{12431}{57500}\sqrt2+\frac{\sqrt3}{2}>0,
\]

so $R$ genuinely cuts it. Keep its first 19 source coordinates and replace its two transverse coordinates by $-3/4,-1/1000$. The resulting exact vector $E'$ satisfies all 43 wall inequalities strictly, has $(E',C)<0$, and has square

\[
E'^2=\frac{20028484177299}{661250000000}
-\frac{87506613029}{4132812500}\sqrt2>0.
\]

The strict spacelike escape gives an open ideal-boundary region, so the new chamber remains infinite volume. The certificate rechecks every sign, rather than inheriting the old escape.

The retained source pairs also preserve a **conditional nonquasi-arithmeticity** witness: walls 0 and 36 have Cartan two-cycle $12+8\sqrt2$, while walls 38 and 39 have a rational Gram minor of determinant $-32$. Any finite-volume completion retaining these four wall reflections fails the nonidentity-embedding positive-semidefiniteness condition. The current infinite-volume chamber itself is not asserted to be a lattice.

## Bounded diagnostics and stopping scope

All numerical files below passed their own saved replays. They are proposal diagnostics, not exact exclusions or all-label searches.

| Diagnostic | Explicit bounded work | Result |
|---|---|---|
| [Corner scout](fresh_h20_joint_route_20261015_scout.json) / [producer](../scripts/fresh_h20_joint_route_20261015_scout.py) | 594 fixed LP queries; 361 distinct finite source corners; 89,528 tests with support labels 3–24 and source angle fits 2–60 | 43 numerical compatible proposals; none cut the old escape |
| [One-release test](fresh_h20_joint_route_20261015_release.json) / [producer](../scripts/fresh_h20_joint_route_20261015_release.py) | 90 shared 17-wall supports; 900 parameter cases; five anchor labels and two transverse label pairs | One compatible proposal; no cutter. Eighty near-null equations were treated linearly, including 16 degenerate cases skipped; no exhaustive branch claim |
| [Two-anchor test](fresh_h20_joint_route_20261015_two_anchor.json) / [producer](../scripts/fresh_h20_joint_route_20261015_two_anchor.py) | Same 90 supports; 810 source-anchor cases; 2,224 transverse tests; anchor labels 3/4/6 | 1,856 accepted numerical recipes, 1,844 after numerical deduplication; two cutters sharing $Y$. The π/6–parallel cutter is the exact $R$ above |
| [Selection diagnostic](fresh_h20_joint_route_20261015_selection.json) / [producer](../scripts/fresh_h20_joint_route_20261015_selection.py) | Only the saved bank; finite angle fits 2–60 at tolerance $10^{-7}$ | The two cutters have positive mutual product about 1.96081. The renewed $E'$ is numerically strict for all 1,844 bank vectors, with largest pairing about −0.48264, even before coherence restrictions |

The exact construction is therefore a finite-edge milestone, not a larger coherent bank or a volume completion. Further progress must use contacts outside this saved bounded bank or change the geometry in another way; the numerical common escape is not a theorem about all such future choices.

Reproduce the exact result with `python3 scripts/fresh_h20_joint_route_20261015_exact.py --verify`. For the numerical saved replays, use `OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_joint_route_20261015_<scout|release|two_anchor|selection>.py --verify` separately for each named file.
