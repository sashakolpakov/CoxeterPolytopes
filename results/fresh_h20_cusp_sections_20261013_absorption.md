# A failed A1-absorption construction at a Borcherds cusp

The proposed local construction does not produce a Coxeter H20 section. The exact result is stronger than the initial width test: **no interior section through cusp 0 can be Coxeter if both inequalities of its affine A1 factor are redundant in the section’s Euclidean cusp link.** This leaves sections retaining one or both A1 link facets open. It does not exclude all sections through this cusp, other cusps, or H20 Coxeter polytopes.

The source is the complete [H21 corner-orbit census](fresh_h21_source_corner_walk_20261013.json). This cusp has affine factors D5, A7, A7, and A1, of total Euclidean rank 20. Its A1 walls are zero-based indices 162 and 164, both of norm four, with

\[
r_{162}+r_{164}=2\delta.
\]

The [standard-library checker](../scripts/fresh_h20_cusp_sections_20261013_absorption.py) reconstructs these data from the original 210 roots, verifies all positive affine relations and principal minors, enumerates the sign-support cases below, and saves the [exact certificate](fresh_h20_cusp_sections_20261013_absorption.json). Both generation and `--verify` pass. No numerical optimizer, finite angle cutoff, or earlier search implementation is used.

An [independent proof review and separate Fraction replay](fresh_h20_cusp_sections_20261013_absorption_review.md) also pass, including independently reconstructed D5 norm forms, every A7 width, and the actual-ridge interpolation. The source and certificate are frozen; this review link does not change either.

## Why the apparent mixed-factor opening becomes small

Work in the Euclidean horosphere quotient at the cusp and normalize the section normal to have square one. A shift by a multiple of δ changes the affine position without changing the projected Gram matrix. Write the normal as a sum of its orthogonal factor components.

If its A1 component vanished, both bounds of that independent interval would remain necessary. Thus absorption of both bounds requires a nonzero A1 component. Projection from the section onto the other three factors is then an affine isomorphism. When the two interval bounds are redundant, its image is the entire original D5 × A7 × A7 product. Every one of the other 22 facet inequalities is consequently an actual facet of the section link; none can be discarded to remove a bad pairing.

In any active remaining factor, the positive affine relation forces both signs among the values `(e,r_i)`. If two remaining factors were active, opposite signs in different factors would give a positive projected pairing, since their original pairing is zero. Hence at most one of D5, A7, A7 can be active. Its normal component has square strictly less than one because the A1 component is nonzero.

Let `v_i=(e,r_i)` in that factor. Its projected Gram is `G_ij−v_i v_j`. Every positive-support node must be joined by an original edge to every negative-support node. Both the D5 tree and the A7 eight-cycle contain no four-cycle, so one sign support is a singleton and the other lies among its neighbors. The checker enumerates every such signed support, rather than assuming a sparse normal.

One elementary all-label bound is sufficient. If a changed node has an unchanged neighbor, their old angle is π/3 and the new squared cosine is

\[
\frac{1}{2(2-v_i^2)}>\frac14.
\]

A Coxeter angle must therefore have label at least four, which forces `v_i²≥1`. On the other hand, the component square being less than one gives `v_i²<2`. No upper label cutoff or crystallographic-angle assumption enters this argument.

## The A7 attempt: the interval is too narrow

For an opposite-sign adjacent pair, the affine marks are all one, so the two values are `a,−a`. Each has an unchanged neighbor, giving `a²≥1`. Acuteness on their common edge gives `a²≤1`. Thus `a²=1`.

A three-node star is impossible: its two outer values `p,q` have unchanged neighbors, so `p,q≥1`; the central magnitude is `p+q≥2`, contradicting its square being less than two.

The unique surviving Gram update changes affine A7 into affine C7: the chosen edge disappears, its two endpoints have norm one, the new path has two π/4 edges and five π/3 edges. The checker reconstructs all 16 choices across the two A7 factors. Each A7 normal component has square `7/8`; the A1 component must have square `1/8`.

Normalize the horosphere by `(δ,x)=−1`. The A7 wall values are nonpositive and sum to −1. The reconstructed coefficients of its normal have range `7/8`, so the A7 functional ranges over an interval of width `7/8`. Meanwhile `r162−δ` has norm four and ranges from −1 to 1 on the A1 interval. A component of square `1/8` therefore contributes width `√2/4`.

An affine translation cannot fit the first interval inside the second:

\[
\frac{(7/8)^2}{(\sqrt2/4)^2}=\frac{49}{8}>1.
\]

Thus the Coxeter C7 metric exists abstractly, but both A1 inequalities cannot be redundant in this actual cusp. This is a width obstruction, not a failed bounded parameter search.

## The D5 attempt: every star violates the norm budget

The D5 factor has four norm-two leaves of affine mark one and two adjacent norm-two centers, 29 and 48, of mark two. Up to overall sign, the following cases exhaust the 26 signed supports in the certificate.

- A leaf–center pair has leaf magnitude twice the center magnitude. The center has an unchanged neighbor, so its magnitude is at least one; the leaf magnitude then violates `v_i²<2`.
- A center–center pair has equal magnitudes. The unchanged-neighbor bound and their mutual acuteness force both magnitudes to be one. The exact inverse-Gram computation gives component square one, leaving no nonzero A1 component.
- A center opposite its two outer leaves has `p+q=2c` and `c≥1`. Orthogonality of the leaves gives component square at least `(p²+q²)/2≥c²≥1`.
- A center opposite one leaf and the other center has `p+2b=2c`. Both centers have unchanged neighbors, so `b,c≥1`; the relation gives `c>b`, while their mutual acuteness requires `bc≤1`.
- A center opposite all three neighbors has `p+q+2b=2c`. The other center has an unchanged neighbor, so `b≥1` and again `c>b`, contradicting `bc≤1`.

This rejects every D5 rank-one update with a nonzero A1 component. Together with the product argument, the only remaining candidate normal plane is the pure A1 plane.

## The remaining plane has an actual obtuse ridge for every position

Set

\[
e(t)=\frac{r_{162}-\delta}{2}+t\delta.
\]

Its square is one. Exactly `−1/2<t<1/2` meets the strict interior: the two endpoint normals are `−r164/2` and `r162/2`. Exact projections of the chamber center onto those two facets give interior intersections for every intermediate parameter. The endpoints and exterior parameters cannot meet the strict interior; the pure δ direction is null and is not a hyperbolic section normal.

The original short walls 45 and 63 are orthogonal. Their pairings with the base normal are `−1/2,+1/2`, and both pair with δ as −1. After projection into the section their mutual pairing and two-by-two Gram determinant are

\[
\frac14-t^2>0,\qquad 3-4t^2>0.
\]

This is an **actual ridge**, not merely an intersecting hyperplane pair. The checker projects the strict center onto the two elliptic ambient faces `{164,45,63}` and `{162,45,63}`. The resulting future timelike points `zL,zR` have precisely those zero inequalities. Their positive interpolation

\[
(1/2-t)Bz_L+(t+1/2)Az_R,
\quad A=-(\delta,z_L)>0,\quad B=-(\delta,z_R)>0
\]

lies in `e(t)⊥` and has exactly walls 45 and 63 active. Both projected walls are therefore facets meeting at an obtuse angle, for every interior parameter. This closes the pure A1 plane and completes the stated absorption obstruction.

```sh
python3 scripts/fresh_h20_cusp_sections_20261013_absorption.py --verify
```

The conclusion depends on the two A1 inequalities being redundant in the **local cusp link**. It does not impose an unproved rule that all Coxeter sections use only one affine factor: products of simplexes can have Coxeter sections mixing factors when appropriate original facets disappear. The width and actual-ridge checks above are what rule out this particular proposed construction.
