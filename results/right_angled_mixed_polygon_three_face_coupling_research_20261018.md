# Polygon attachment and a three-face exception budget

Date: 2026-10-18. New written universal incidence inequalities, with exact
small coefficient checks. No new optimizer run, claimed new LP optimum,
Lean verification or dimension reduction.

## Results

For a finite-volume right-angled polytope, use f,c,a2,a3,B,T as in the
reviewed weighted packet. The new universal rows are

    n=12: 101B-5130a2+69080f+22049280c >=0,
    n=11: 172B-7866a2+88275f+15367680c >=0.

Combining them with the existing weighted and canonical P3 sum gaps gives

| Dimension | Lower bound for 2B/f | Lower bound for C/f | Fraction of finite corners containing K4 |
| --- | --- | --- | --- |
| 12 | 184460/1619 | 40422/40475 | 6737/21047 |
| 11 | 117825/1379 | 50254/62055 | 25127/79982 |

All three bounds in each row are strict. The twelve-dimensional triangle
lower bound is approximately 113.93453, still below the universal
reciprocal-square upper cap 165.43744. These are improved necessary
densities, not a dimension exclusion.

There is also an independent polygon/corner-graph coupling. If J2 counts
finite corners of nonbad two-faces whose two neighboring vertices are
ideal, then in dimension twelve

    649T >21846f+1430J2.

If s(G_q) counts nonisolated directions in a finite corner graph,
J2>=sum_q[binomial(s(G_q),2)-e(G_q)]. This rules out the specific nonempty
multipartite mixture displayed in the earlier exact LP certificate as a
literal polygon/corner realization. It does not by itself exclude every
point of that optimal face.

## A polygon corner carries more information than its size

Let a nonbad right-angled polygon have u finite and v ideal vertices, and
let j count its finite vertices with two ideal neighbors. If v<=1, j=0.
Otherwise j<=min(u,v): each such vertex uses two finite/ideal polygon
edges, and each ideal vertex supplies at most two. Put

    lambda=2123/1280, alpha=331/3840, beta=143/384.

Then

    u+lambda v-22/5 >=alpha u+beta j.

The necessary minimal cases are (u,v,j)=(0,3,0),(1,3,1),(5,0,0),
(3,1,0),(2,2,2). All pass, with equality for the last two. On the region
u,v>=2, adding a finite or ideal vertex increases the left minus right
side by at least 1-alpha-beta or lambda-beta, both positive. The separate
u=0,u=1,v=0,v=1 boundary regions have the displayed minimal cases. This
proves the inequality for every permissible nonbad polygon; the sole
excluded positive-area triangle with finite vertices is the counted
(u,v)=(1,2) case.

A counted triangle has deficit 318/3840. Summing, using the strict weighted
two-face inequality and the finite corner sum 66f, gives

    331(66f-T)+1430J2 <318T,

which is the stated twelve-dimensional row.

At q, a nonisolated direction i belongs to a bad triangle, so its actual
edge ends at the ideal point q-n_i. Each nonedge between two such
directions belongs to a nonbad two-face through q whose two neighbors of
q are ideal. Thus it contributes to J2. In the earlier LP candidate the
nonempty graph has all twelve directions nonisolated and exactly three
nonedges. Its weight 662/1121 gives J2>=1986/1121; the new row then has
negative slack -262152/1121. This is a genuine attachment obstruction.
The statistic s is not a clique count, so the old 77-vector reduction may
not simply be reused for it.

## The zero-excess eight-facet three-face

Call a right-angled three-polytope with F=8 facets and s=2 cusps a Q8
three-face; uniqueness up to isometry is unnecessary. Its Euler data give
p=8 finite vertices. The total polygon excess u+2v-5 is zero, so each
facet is a bad triangle, a (3,1) quadrilateral or a finite pentagon. A bad
triangle would have a bad triangle across its ideal/ideal edge; the
reviewed attachment lemma forces the whole polytope to be P3, a
contradiction. There are therefore no triangles. Its thirty-two polygon
corners on eight facets force every facet to be a (3,1) quadrilateral.

At a finite Q8 vertex there is exactly one incident edge with ideal other
endpoint. There cannot be two, because their quadrilateral would contain
two ideal vertices. No edge joins the two cusps, for the same reason.
The two degree-four cusps therefore supply eight finite/ideal edges.
The eight finite vertices have at most one apiece, so each has exactly one.

## Every finite corner requires exceptional three-faces

Call a three-face exceptional if it is neither P3 nor Q8. At a finite
ambient vertex q let k be the number of its n incident edges ending at
ideal vertices, and t=t(G_q). Exactly t of its incident three-faces are P3.
A Q8 can use only a triple with exactly one of the k ideal directions, by
the preceding lemma. Thus the exceptional three-face count X(q) satisfies

    X(q) >=binomial(n,3)-t-k binomial(n-k,2).

Every graph triangle uses three ideal directions, so t<=binomial(k,3).
The following finite integer inequalities hold for every k=0,...,n:

    (1-gamma)[binomial(n,3)-binomial(k,3)] >=k binomial(n-k,2),
    gamma=13/27 for n=12; gamma=11/23 for n=11.

Equality occurs at k=4. They are checked exactly for all 25 values across
the two dimensions. Consequently

    X(q) >=gamma[binomial(n,3)-t(G_q)].

If X counts exceptional finite-vertex/three-face incidences, summing gives
X>=gamma[binomial(n,3)f-2B]. This couples actual ideal edge endpoints and
the types of incident three-faces; it uses no graph symmetrization.

## An exact budget for exceptional incidences

For any right-angled three-face use F,s,p and nu=F+s. The established
three-dimensional facts are F>=6 with equality only for P3, nu>=9 with
equality only for P3, F+2s>=12, and p=4F-2nu-4. For every exceptional face,

    p <=22(nu-10)+4(8-F).

Here is the complete case proof. If nu=10, F+2s>=12 implies F<=8, and the
six-facet classification excludes F=6. The remaining F=8 is Q8, while
F=7,s=3 has p=4, exactly the displayed budget. If nu=11, then F<=10 and
the right minus left side is 24nu-184-8F>=0. For nu>=12, F<=nu makes this
difference at least 16nu-184>=8. The remaining nu=9 case is P3, excluded.
The budget is zero for Q8. No additional lower-cusp theorem is used.

Define the nonnegative nu-excess E and signed facet defect R by

    E=sum_non-P3(nu-10)
     =B-10a3+(n-2)a2+binomial(n-1,2)2^(n-3)c,
    R=sum_non-P3(8-F)=8a3-(n-2)a2-2B.

R need not be nonnegative; it is not a face count. Summing the budget gives
X<=22E+4R. Combined with the preceding local lower bound,

    22E+4R >=gamma[binomial(n,3)f-2B].

For n=12, substitute a3=5a2-55f-14080c and simplify to the first row in
the results. For n=11, substitute a3=(9/2)a2-(165/4)f-5760c to get the
second row. The checker verifies every coefficient.

## Eliminating a2 and c

The twelve-dimensional weighted row a2>15f+4246c makes the new inequality
strictly imply

    101B>7870f-267300c.

Combine it with B>50f+880c (the P3 sum gap plus the weighted row). Multiply
the latter by 267300 and the former by 880; the cusp terms cancel. This
gives B>92230f/1619, the first bound in the table.

Similarly in dimension eleven, a2>(25/2)f+(21230/11)c gives

    172B>10050f-186300c.

Combine with B>(75/2)f+360c, with respective weights 360 and 186300, to
obtain B>117825f/2758. Finally use t<=64+10k4 on twelve vertices and
t<=49+9k4 on eleven vertices, together with sum k4=5C, for the C bounds.
The finite-corner fractions follow from the respective total triangle
maxima 220 and 165. All eliminations use positive multipliers and preserve
strictness.

## Scope and reproduction

Run `python3 scripts/check_right_angled_mixed_polygon_three_face_coupling_research_20261018.py`.
The compact JSON records every finite case and rational coefficient.
The supporting geometry consists of reviewed P3 attachment facts and the
written arguments above. No additional optimizer or Lean process is run.
The old exact LP minimum remains correct for its earlier, weaker model;
the new rows are not asserted to attain their new lower bound.
