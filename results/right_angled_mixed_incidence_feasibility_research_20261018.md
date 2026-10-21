# Exact feasibility of the current mixed-dimension aggregate inequalities

18 October 2026. **Exact diagnostic obstruction to this aggregate route.**
The new canonical-face gaps, their coface/vertex incidences, the used
Nikulin–Khovanskii inequalities, and all simple-polytope g-vector
conditions listed below have a simultaneous integral numerical assignment
in dimension12. They therefore do not, by themselves, exclude dimension12.
This is a statement about the explicitly listed aggregate conditions.
No boundary complex, actual hyperplane arrangement or polytope is constructed.

## Exact assignment

Normalize the cusp count c to1 and prescribe

\[
 (a_6,a_7,a_8,a_9,a_{10},a_{11})
 =(990,100,5,1/20,1/100000,1/10^{11}).
\]

Determine f and a_1,...,a_5 from the exact dimension12 Dehn–Sommerville
rows in the [weighted packet](right_angled_mixed_weighted_triangles_research_20261018.md).
In particular f=4830020999979/50000000000, approximately96.60042.
Set the numbers of one-finite/two-ideal triangles and canonical three-,
four- and five-faces to

\[
 T=5590,\qquad B=9000,\qquad C=8000,\qquad D=3000.
\]

Multiply every count by N=396 times10^22. All face, polygon-class and
finite-corner graph counts below become integers. The resulting abstract
facet and cusp counts also exceed the absolute published lower bounds
used in the preceding packets. Scaling does not claim geometric gluing.

The [compact certificate](right_angled_mixed_incidence_feasibility_research_20261018.json)
records all rational counts and exact minimum slacks. The
[standalone checker](../scripts/check_right_angled_mixed_incidence_feasibility_research_20261018.py)
recomputes them with Python's standard-library rational/integer arithmetic.
No optimizer was needed.

## Precise conditions satisfied

The assignment satisfies the dimension12 Dehn–Sommerville identities,
positive face counts, the vertex-degree identity, and the three-face
Euler identity. In addition it satisfies all the following.

1. The Nikulin–Khovanskii inequalities with2<=k<=6 and0<=l<k, strictly.
   For positive l their exact form is
   a_l < binom(6,l)a_k/binom(6,k); for l=0 the incidence numerator is
   binom(12,k)f+beta_k c, with beta_k=binom(11,k-1)2^(12-k).
2. The weighted two-face bound with ideal weight2123/1280, both strict
   vertex-ratio bounds, and the strict triangle/B/C density consequences
   proved in the weighted packet.
3. The facet minima for dimensions3,...,11:
   6,10,16,27,46,82,218,1923,183968. Their summed inequalities are
   (13-k)a_(k-1)>=minimum_k a_k.
4. The facet-plus-cusp minima in dimensions3,...,11:
   9,15,26,45,81,217,1922,183967,1674688396. Their summed inequalities are
   (13-k)a_(k-1)+beta_k c>=nu_k a_k.
5. Both facet-count and minimum-sum canonical-face gaps:

\[
\begin{aligned}
 B&\ge7a_3-10a_2,&B&\ge10a_3-10a_2-28160c,\\
 C&\ge11a_4-9a_3,&C&\ge16a_4-9a_3-42240c,\\
 D&\ge17a_5-8a_4.
\end{aligned}
\]

The minimum-sum B gap follows from the already proved uniqueness of
facet-plus-cusp sum9 in dimension3. The C and D gaps use the
[reviewed low-facet rigidity packet](right_angled_mixed_low_facet_rigidity_research_20261018.md).

6. All canonical face upper bounds0<=B<=a_3,0<=C<=a_4,0<=D<=a_5;
   their finite/cusp incidence capacities, using respective vertex counts
   (2,3),(5,5),(16,10); and all the coface inequalities

\[
 6B\le10T,\quad10C\le9B,\quad16D\le8C,\quad
 30C\le45T,\quad160D\le120T,\quad80D\le36B.
\]

7. The full g-vector conditions for the truncated simple twelve-polytope.
   After the integral scaling, its computed h-vector is symmetric with
   h_0=1, and g_1,...,g_6 are positive integers satisfying all five
   Macaulay growth bounds. These are checked exactly, not by floating
   point or merely by nonnegative h-coefficients.

The published scalar minima in items3–4 are the statements already
used in the [weighted packet](right_angled_mixed_weighted_triangles_research_20261018.md)
and recorded in [Alexandrov2023, Section6](https://arxiv.org/html/2303.09533v1#S6).
The simple-polytope g-vector test is an additional numerical consistency
check; it supplies no prescribed simplex/cube links for the untruncated
object.

## Polygon and finite-corner counts are also explicit

The polygon types, listed by finite/ideal vertex counts, are
(1,2),(0,3),(3,1),(5,0). Put

\[
 E=66f+22528-5a_2,
 \quad Q=11264-2T-3E,
 \quad R=a_2-T-E-Q.
\]

They occur with counts T,E,Q,R respectively. All four counts are
positive, sum to a_2, and give exactly66f finite and11264 ideal corners.
Thus feasibility did not come from forgetting the two different corner
counts or the polygon area condition.

For the graph on twelve edge directions at a finite vertex, the exact
weighted graph multiset is

|Graph|Number before scaling|
|---|---:|
|K12|8000/99|
|Complete tripartite K4,4,4|125/36|
|Complete bipartite K6,6|5/2|
|Empty graph on twelve vertices|f-8000/99-125/36-5/2|

Every weight is positive and their sum is f. Their edge, triangle and
four-clique totals are respectively T,2B,5C, exactly as required by the
actual-face correspondences proved in the weighted packet. There are
also enough five-cliques to supply the necessary16D finite-corner
incidences of canonical five-faces:792 times8000/99=64000 exceeds16D=48000.
This is only the necessary implication from a canonical five-face to
its corner five-cliques. A five-clique is never identified with such a
face in this diagnostic.

## Replay and the missing geometric condition

```sh
python3 -B scripts/check_right_angled_mixed_incidence_feasibility_research_20261018.py
```

Generation and separate-process replay passed in terminal chunk0b80fe.
The checker runs in a fraction of a second; no LP, Lean or heavy job ran.

This witness includes many K12 finite-corner graphs. A graph records
which particular triangles occur; it does not provide the outward
normals, consistent opposite facets at shared cusps, actual face gluing,
or widths of the cubical links. Those conditions may exclude the
assignment. Even a complete K5 corner does not force the standard
five-polytope without an additional minimum-facet argument. The next
useful step is therefore a geometric attachment obstruction for the
K12/K5-rich corner systems, rather than another elimination of these
same aggregate inequalities. The published ceilings12 and6 remain
unchanged.
