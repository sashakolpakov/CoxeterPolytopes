# New H20 cusp planes and angle rotations — 2026-10-13

This search produced four new rational supporting planes and eight new angle-rotation planes. It did **not** produce a finite-volume Coxeter polytope. All 24 tested compatible selections involving one rotation have exact strict spacelike escapes. The exclusions apply to those recorded selections; they do not exclude further rotations, mixed assemblies, or inward cuts.

Indices below are zero-based in the 580-wall order of `data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json`. The Lorentz form, section normal, and common strict timelike center are the same as in the independent six-tile and deletion certificates.

## Exact local cusp replacements

The seven saved null vectors were checked exactly against all 580 walls and decomposed by their incident Gram graphs. Their component sizes/ranks are respectively:

| Cusp | Component size/rank |
| --- | --- |
| 0 | 9/8, 5/4, 8/7 |
| 1 | 11/10, 8/7, 3/2 |
| 2 | 13/11, 9/8 |
| 3 | 8/7, 9/8, 5/4 |
| 4 | 8/7, 10/9, 4/3 |
| 5 | 13/11, 9/8 |
| 6 | 5/4, 9/8, 8/7 |

Thus the links do not all split into three simplex factors. Each total Euclidean rank is 19.

At cusp 0, the non-Coxeter rank-seven factor on walls `26,29,30,31,32,38,100,101` is a section of an affine E8 alcove. An exact reflection closure gives 240 finite E8 roots, of which 126 are orthogonal to the section wall. Their closest affine translates select an affine E7 alcove. It keeps walls `26,29,30,31,32,38` and replaces the last two by

\[
e_1=(r_{100}+r_{101})/2,\qquad
e_2=2r_{26}+2r_{29}+r_{31}+r_{38}+r_{101}.
\]

Both new roots have square 2. The resulting Gram graph is the affine E7 tree with arms of lengths 3, 3, and 1. The archive contains a nonnegative old-wall expansion for every new root, and a nonnegative new-root expansion for all 126 translated E7 roots. These verify the proposed Euclidean alcove directly.

At cusp 6, rescaling the factor on `0,1,2,3,41` gives affine F4. Its 24 short roots select an affine D4 alcove. The replacement keeps `2,3,41` and introduces

\[
d_1=r_1+r_2,\qquad d_2=r_0+r_1+r_2.
\]

The scaled Gram is the affine D4 star, with diagonal 18 and edge entries −9. All four added inequalities are nonnegative combinations of old inequalities, so their halfspaces contain the original six-tile polytope. They describe genuine local cusp coarsenings, but do not by themselves close the global construction.

In fact both exact mandatory-wall escape vectors from the earlier deletion obstruction satisfy all four new inequalities strictly:

| Escape | Pairings with `e1,e2,d1,d2` |
| --- | --- |
| Delete wall 0 | −46, −98, −38, −14 |
| Delete wall 74 | −178461556, −390384662, −217500034, −518653914 |

Consequently the earlier two-wall obstruction still applies to the 584-plane dictionary consisting of the old walls and these four planes.

## Rotating the mandatory bad ridge

Let `a=r0`, `b=r74`. Their exact Gram is

\[
(a,a)=4,\quad (b,b)=6,\quad (a,b)=-4.
\]

For `m=2,3,4,5`, put `t_m=√2 cot(π/m)`; these values are `0`, `√(2/3)`, `√2`, and `√(2+4√5/5)`. Two orientations were tested:

\[
w_m=3a+(2-t_m)b,\qquad
z_m=2b+(2-t_m)a.
\]

The first replaces `a` and retains `b`; the second replaces `b` and retains `a`. Each pair has angle `π/m`. Their squared norms are `12+6t_m²` and `8+4t_m²`. Since `0≤t_m<2`, each plane is a positive combination of old normals and its halfspace contains the known polytope.

This genuinely bypasses the old obstruction: for example, `w4=3a+(2−√2)b` cuts the simple old escape `x0=c+10a`, since `(w4,x0)=−28+50√2>0`.

For each orientation and target angle, the search kept all individually compatible old walls available, forced the anchor wall, then generated up to six greedy maximal compatible selections. Every pair was checked exactly. The declared angle alphabet was `{2,3,4,5,6,8}`; parallel and disjoint pairs were also allowed. This was a bounded selection procedure, not an exhaustive search over all subsets or all labels.

| Replaced wall | Target angle | Compatible old walls available | Distinct selections | Total inequalities |
| --- | --- | ---: | ---: | ---: |
| 0 | π/5 | 443 | 1 | 142 |
| 0 | π/4 | 504 | 1 | 160 |
| 0 | π/3 | 548 | 1 | 160 |
| 0 | π/2 | 562 | 1 | 174 |
| 74 | π/5 | 174 | 5 | 150–152 |
| 74 | π/4 | 425 | 6 | 275–281 |
| 74 | π/3 | 523 | 6 | 334 |
| 74 | π/2 | 569 | 3 | 362 |

All 24 selections have saved integer vectors `x` with `(x,x)>0`, `(x,c)<0`, `(x,f)=0`, and strictly negative pairing with every selected wall, including the algebraic new wall. The line from the strict center to such a vector reaches an open set of ideal boundary inside all inequalities, certifying infinite volume. Numerical LPs only proposed these vectors; every sign and squared norm is checked exactly. The π/5 field is handled with the polynomial `5t⁴−20t²+4` and rational isolating intervals.

An independent direction now being explored is **inward** replacement. For example `v=3a−√2b` makes angle π/8 with `b` and still contains the strict center. The inequalities `v≤0,b≤0` imply `a=(v+√2b)/3≤0`, so replacing `a` in this way preserves finite volume by containment, provided the other old inequalities are retained or exactly proved redundant. This report makes no claim that the remaining boundary angles of that inward-cut construction are Coxeter.

## Artifacts and validation

- `scripts/fresh_h20_cusp_repair_20261013.py`
- `data/enumerations/fresh_h20_cusp_repair_20261013.json`
- `scripts/fresh_h20_cusp_repair_rotations_20261013.py`
- `data/enumerations/fresh_h20_cusp_repair_rotations_20261013.json`

The cusp replay reconstructs the finite root systems, all cone coefficients, and the saved decompositions. The rotation replay recomputes every selected Gram pair and every escape without invoking LP optimization. Both commands passed:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_cusp_repair_20261013.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_cusp_repair_rotations_20261013.py --verify
```

A corrupted rotation certificate substituting the timelike center for the first escape was rejected at the exact positive-norm assertion.
