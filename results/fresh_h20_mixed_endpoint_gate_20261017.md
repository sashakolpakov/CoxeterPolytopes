# Finite and ideal endpoint caps do not close this dotted-other family

**No finite coherent Coxeter batch from the complete finite-or-parabolic endpoint family below gives finite volume.** This closes the ideal-endpoint loophole in the earlier finite-endpoint theorem. It is a restriction on the three selected old edges and their specified old contact pattern; it is not a nonexistence result for H20 polytopes.

The batch retains all168 walls. Each new square-two normal r extends one of the spherical rank-19 supports W = T+{159,160}, T+{159,161}, or T+{160,161}; its W+r Gram is either positive definite or positive semidefinite of corank one. Every unprescribed old168 contact other than158 is strictly dotted; all Coxeter contacts with158 are allowed. The new walls are mutually Coxeter-compatible and have a common strict future timelike interior. There is no bound on angle labels, heights, or the finite number of added walls. Additional finite old contacts and intermediate facets outside these endpoint profiles are not covered.

The [source](../scripts/fresh_h20_mixed_endpoint_gate_20261017.py) and [certificate](fresh_h20_mixed_endpoint_gate_20261017.json) pass exact generation and a fresh saved replay and are final/frozen. They bind the complete [finite-endpoint gate](fresh_h20_all_endpoint_separation_20261017.md), [parabolic deletion census](fresh_h20_parabolic_endpoint_profiles_20261017.md), and [continuous parabolic height gate](fresh_h20_parabolic_endpoint_gate_20261017.md). Arithmetic uses the existing outward-rational interval engine; this is a primary certificate, not an independent arithmetic review.

## Complete mixed contact reduction

The finite census has159 fixed profiles and three singleton151 families. The parabolic census has20 profiles, five of which have positive projected normal square perpendicular to T. The other15 have projected square zero and cannot meet the hyperbolic T-space. Their positive affine relation makes the nonzero projected null vector future-pointing; hence their inequalities are strictly satisfied throughout that hyperbolic space. They cannot add sides to a target-face facet.

For each positive parabolic profile the preceding gate supplies r(x)=b+x nu, x=-r.158, with nu the unique old-feasible future null ray and nu.158=-1. Its continuous necessary lower bounds combine with the exact Coxeter magnitude list

\[
0,\ 1,\ \sqrt2,\ \varphi,\ \sqrt3,\ldots,\ 2,\quad\text{or a value greater than }2.
\]

Consequently the edge159/160 family has x>=varphi; the two middle-edge nonzero-bound families have x>=1; the remaining middle-edge family has x>=0. The edge160/161 family splits into the single value x=1 and the continuous tail x>=sqrt2. These are safe relaxations of every finite, parallel, and dotted contact with158, not sampled angle domains.

For the fixed finite profiles and the singleton151 labels3,4,5, the other wall's158 contact uses the inherited discrete labels2..6 plus its complete continuous conic tail. Multiplying by the positive conic denominator reduces each product r.s+2 to a polynomial. The certificate checks coefficientwise strict negativity for every mixed continuous case. Of508 discrete pairs,506 are strictly dotted. The other two have the following exact interval certifications:

| Parabolic profile | Finite profile | Contact magnitudes with158 | Product, approximately | Forbidden label gap |
| --- | --- | --- | --- | --- |
| edge160/161, profile0 | edge159/160, profile38 | 1,1 | -1.59092250348 | 4–5 |
| edge160/161, profile0 | edge159/160, profile50 | 1,0 | -1.97077696724 | 18–19 |

The consecutive-angle brackets use rigorous rational Machin/cosine bounds. The decimals in the table are explanatory only.

For a finite singleton151 label m>=6, an actual neighbor in the T-face would make T+r+s spherical. The I2(m) component is isolated, so the parabolic wall must be right to151 and to the finite cap. The one parabolic profile with a nonright151 contact is immediately excluded. For the remaining profiles the finite cap divided by d=2cos(pi/m) has its largest pairing at d=sqrt3, holding x158/d fixed: the coefficient of the future old-corner vector increases with d, and every coherent wall is strictly negative on that old finite vertex. The inherited m6 comparison therefore applies. All ten mixed high-label continuous/discrete comparisons are strictly negative and exclude the necessary orthogonality.

The case counts are546 fixed tail/tail,2542 fixed tail/discrete,30 low tail/tail,146 low tail/discrete,8 high tail/tail,102 fixed discrete/tail,6 low discrete/tail,2 high discrete/tail, and508 discrete/discrete. All continuous coefficient upper bounds are negative; the least margin is greater than0.4203 in absolute value. The certificate saves the exact bounds and stream digest.

## Same-edge mixtures are impossible

This part does not depend on angle gaps. Write a parabolic cap as

\[
s=p_s+h\nu,\qquad p_s=-WJd_s,\quad h>0.
\]

Write a spherical endpoint cap on the same edge as r=p_r-a e+b V, with p_r=-WJd_r. Here e is the normalized projection of158 perpendicular to W, V is the old future unit finite vertex, and nu=(V-e)/sqrt(sigma). Coherence selects b>0. Since r.158<=0 and p_r.158>=0, a>0. Its positive spherical Schur residual gives a²-b²>0, hence a>b. Therefore

\[
r\cdot\nu=(a-b)/\sqrt\sigma>0,\qquad
p_r\cdot p_s=d_r^T Jd_s\ge0,
\]

and r.s>0. Such two normals cannot belong to an acute chamber. Equivalently, the positive affine relation of s would force r.nu<=0, contrary to the finite cap cutting that old open edge.

Two different parabolic profiles that meet T on the middle edge have the strictly positive constant products sqrt2,1,sqrt2. Two different heights of the same profile have product2. Thus at most one parabolic facet meeting T can occur on each edge. At most one spherical endpoint cap per edge was already proved by persistence of the old finite endpoint, acute finite-vertex simplicity, and convexity of the one-dimensional edge. These conclusions allow zero-projected-norm cusp walls; those walls supply no target-face sides.

## Why a finite batch still has an open ideal region

For any added cap meeting T, the Gram of T+r is positive definite, including for the parabolic endpoint profiles. Projection of a common strict center onto (T+r)-perp gives a hyperbolic point strictly inside every other inequality: the spherical inverse Gram is entrywise nonnegative, while all cross pairings are nonpositive. Thus r defines an actual facet of the target H3 face. In particular no wall outside T union {r} can vanish identically on this facet.

The old compact158 triangle remains strictly inside every cap, because its three old finite vertices persist strictly under coherent acute additions. All other unprescribed old162 walls are intrinsically disjoint from r by the dotted assumption and the nonnegative spherical projection correction. The six later cusp-completion walls miss T by the frozen height gate. The mixed arithmetic above, together with the same-kind gates, excludes every added neighbor from another old edge; the same-edge arguments exclude extra neighbors there as well.

Consequently an added facet in T has at most the two prescribed old partner sides. Its nonempty hyperbolic region is a wedge with an open ideal arc; the parabolic endpoint case may make the two sides asymptotic, which still leaves such an arc. The158 inequality cannot introduce a hidden side inside this wedge: such a boundary point would lie in the actual old158 triangle, already strict for the cap. Other old planes are disjoint, and other new planes cannot meet the facet inside T.

Choose an ideal point on the open arc avoiding the finitely many boundary zeros of walls outside T union {r}. Perturb it slightly to a spacelike vector within (T+r)-perp, preserving all these strict inequalities and its future affine chart. Then perturb toward the common strict future center. The resulting vector is still spacelike and is strict for every wall, including T and r. A finite-sided chamber with such a strict spacelike cone direction has an open ideal-boundary patch and infinite hyperbolic volume.

If the batch contains no cap meeting T, the same conclusion follows without assuming a numerical target-face census. Each certified old open-edge ray nu is strict for all retained old walls outside W. The two partner normals are spherical after projection off T, so nearby future null directions in T-perp can be moved strictly inside both partner inequalities. All remaining old inequalities stay strict. The finitely many new projected null normals vanish at at most one ideal point each; avoid those points. This again supplies a strict null patch, and the spacelike/interior perturbation above gives an ambient strict escape.

The obstruction therefore includes finite and ideal direct endpoints and arbitrary finite cusp-completion batches within the stated contact pattern. A useful next construction must leave that pattern: introduce an additional finite old contact, an intermediate facet not completing one of these three edges, or change the retained old walls. Local bounded cusp links alone do not close the remaining target-face sides.
