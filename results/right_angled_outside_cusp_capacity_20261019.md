# Outside-cusp capacity and stronger entirely ideal right-angled counts

19 October 2026. New written proof, awaiting independent review. No calculation,
optimizer, Lean process or TeX build was run for this report. These are necessary
conditions for existence, not dimension exclusions or claims of published
priority. Previously frozen reports and certificates are unchanged.

The new geometric input is a capacity of two: an ideal vertex outside one
facet of an entirely ideal right-angled polytope is adjacent to at most two
vertices in that facet. With the new Figure7 exclusion this gives:

| Entirely ideal class | Earlier internal necessary counts | New proposed necessary counts |
|---|---:|---:|
| Right-angled dimension5 | Facets at least56; vertices at least170 | Facets at least59; vertices at least180 |
| Right-angled dimension6 | Facets at least227; vertices at least4462 | Facets at least240; vertices at least4980 |

All polytopes here are full-dimensional, convex and finite-volume. Write V(X)
for the number of vertices of X and F(X) for its number of facets. Entirely ideal
means that every vertex is ideal. Positive dimensional faces inherit finite
volume and right angles.

## 1. Acute projection, including walls outside the chosen cusp

Use the hyperboloid model with outward unit spacelike normals n_i and
P={x in H^n: (x,n_i)<=0}. The standard acute-chamber normal fact gives
(n_i,n_j)<=0 for distinct actual facet walls of a right-angled polytope,
including walls which do not meet inside P. This is an inherited foundational
acute-polytope fact, not inferred from a finite list of examples.

Here is the projection argument used below. Let selected walls have positive
definite Gram matrix M whose inverse is entrywise nonnegative, and let q be
a strict interior future timelike point of P. Put b_i=(q,n_i)<0,
u=M^{-1}b and p=q-sum_i u_i n_i. Then u_i<0, (p,n_i)=0, and

    (p,p)=(q,q)-b^T M^{-1}b<0.

For every other actual wall n_j,

    (p,n_j)=(q,n_j)-sum_i u_i(n_i,n_j)<0,

because every summand subtracted is nonnegative. Also (q,p)<0, so p is
future timelike. After normalization p is an actual point of the selected
face, with every other wall strict. Thus other noncusp walls cannot remove
the face produced by this projection.

Two consequences will be used. First, n mutually orthogonal actual walls
give an actual finite vertex. Second, if two actual unit normals had product
in (-1,0), their Gram matrix would be positive definite with nonnegative
inverse. The projection would give an actual codimension-two face with
nonright angle, impossible here. Consequently distinct actual unit normals
have product either0 or at most-1. This conclusion concerns the actual global
facets, even when their hyperplanes might intersect outside the initially
visible chamber.

## 2. An outside ideal vertex has at most two neighbors in a fixed facet

Let P be entirely ideal and right-angled of dimension n>=3, let F be a facet,
and let w be an ideal vertex not in F. Send w to infinity in the upper-half
space model. A horospherical link at w is a bounded rectangular (n-1)-box;
write its incident vertical walls as

    x_i=0 and x_i=L_i,     1<=i<=n-1,     L_i>0.

The side of F containing P also contains w. Therefore F's supporting
hyperplane is a hemisphere of Euclidean center c=(c_1,...,c_{n-1}) and
radius R>0, with P on the exterior-of-ball side. Its unit-normal products
with the two i-th vertical walls are

    -c_i/R and (c_i-L_i)/R.

Section1 makes them nonpositive, so 0<=c_i<=L_i. Its zero-or-at-most-minus-one
alternative gives the stronger coordinate conditions

    c_i in {0} union [R,infinity),
    L_i-c_i in {0} union [R,infinity).                 (2.1)

Every edge from w is vertical over a vertex xi of this box. If its other
endpoint belongs to F, that endpoint is ideal, has Euclidean coordinate xi,
and lies on the hemisphere's boundary sphere. Hence

    sum_i (xi_i-c_i)^2=R^2.                           (2.2)

Each absolute coordinate difference in (2.2) is either0 or at leastR by
(2.1). Exactly one difference is nonzero, and it equalsR. Thus c lies on at
least n-2 boundary coordinate hyperplanes of the box.

The center cannot be a box vertex. If it were, F and those n-1 incident
vertical facets would be n mutually orthogonal actual facet walls. The
projection in Section1 would give an actual finite vertex of P, with
every other facet inequality strict. This contradicts entire ideality; it
does not assume that an untested hemisphere apex is feasible.

Consequently c lies in the relative interior of one box edge. Only the two
endpoints of that edge can satisfy (2.2). There are therefore at most two
edges from w to vertices of F. Distinct edges from w give distinct box
vertices and distinct other endpoints. We have proved

    every w outside F has at most two graph-neighbors in F.              (2.3)

If there are no such neighbors the conclusion holds without placing c on
a box edge; that placement was deduced only after assuming a solution of
(2.2). No claim is made that every permissible sphere or box direction
produces an actual neighboring vertex.

## 3. Vertex propagation across an ideal facet

At a vertex of P the link is an (n-1)-box. Its number of vertices, and hence
the graph degree in P, is 2^(n-1). At a vertex of F the degree inside F is
2^(n-2). Thus exactly 2^(n-2) edges leave F at each vertex of F. An edge
whose two endpoints lie in F lies entirely in F, since F is a face, so each
of these edges has an outside endpoint. Counting with (2.3) yields

    2^(n-2)V(F) <= 2(V(P)-V(F)),
    V(P) >= (1+2^(n-3))V(F).                          (3.1)

The link degree is the number of box vertices, not its number of facets.
This distinction matters from dimension4 onward. In particular, for a
three-facet R of an ideal four-polytope Q,

    V(Q)>=3V(R).                                     (3.2)

## 4. The new large three-facet carrier, with two separate counts

The [new Figure7 exclusion](ideal_ra4_ten_vertex_facet_exclusion_research_20261019.md)
proves that every ideal right-angled four-polytope all of whose facets have
at most ten vertices is an O/D chain. Its numeric cutoff uses the stated
complete published classification through ten vertices, while its metric
exclusion and chain steps are separate internal written proofs. Each such
chain has kappa=2-2m<=0, where m is its number of24-cell pieces.

The [five-dimensional carrier theorem](right_angled_five_face_carrier_20261018.md)
gives an actual four-facet Q with positive
kappa(Q)=3V(Q)-4F(Q)+24 in every ideal right-angled five-polytope. It follows
that Q has a three-facet R with m=V(R)>=11. Equation(3.2) proves

    V(Q)>=3m>=33.                                    (4.1)

The facet count needs a distinct argument; (4.1) alone would not prove it.
For an ideal right-angled three-polyhedron, degree4 and Euler give
F(R)=V(R)+2=m+2. The already reviewed actual-facet injection gives

    F(Q)>=1+F(R)+V(R)=2m+3>=25.                       (4.2)

For completeness, the injection consists of R itself, one distinct adjacent
facet for every ridge of R, and the unique facet opposite R in the cusp-box
link for every vertex of R. An opposite has hyperplane asymptotic to R's
hyperplane at that one cusp, so cannot serve a second cusp. Nor can it be a
ridge neighbor, whose hyperplane is orthogonal to R. All entries are actual
facets; no census or extension conjecture is used. This is exactly the
injection already written in the [five-face density report](right_angled_mixed_five_face_density_20261018.md)
and [ideal-six propagation report](right_angled_ideal6_carrier_propagation_20261018.md),
now applied one dimension lower.

Applying the same injection to Q in an ideal five-polytope P gives

    F(P)>=1+F(Q)+V(Q)>=1+25+33=59.                   (4.3)

## 5. Strict ideal-five Euler count

Write v=V(P), f=F(P), and g=f_3(P) for an ideal right-angled five-polytope.
Here is the box-link derivation of the needed Euler identity. A four-box
link has16 vertices and24 two-faces. Thus f_1(P)=8v and the sum of the
vertex counts of all three-faces is24v. Every ideal right-angled three-face
has V+2 facets, and every two-face has exactly three three-cofaces.
Consequently 3f_2(P)=24v+2g. Substitute this and f_1(P)=8v into
v-f_1(P)+f_2(P)-g+f=2 to obtain

    g=3v+3f-6.                                       (5.1)

The published facet-minimum theorem for ideal right-angled four-polytopes
gives at least24 facets for every four-facet. At least one actual four-facet
has at least25 by (4.2). Thus the coface incidence is strictly larger:

    2g=sum_Q F(Q)>24f.

Together with (5.1) this gives v>3f+2. Integrality and (4.3) imply

    V(P)>=3F(P)+3>=180.                              (5.2)

This strictness uses only the published minimum24 and the actual large
carrier (4.2). It does not use the uniqueness of the minimum, the second
volume theorem, or a claim that every non24-cell four-polytope is a chain.
Equivalently, sum_Q(F(Q)-24)=6(v-3f-2) is a positive multiple of6.

## 6. Ideal-six facets and the sharp reciprocal transport

Every five-facet S of an ideal right-angled six-polytope has at least59
facets and180 vertices by (4.3) and (5.2). The actual-facet injection gives

    F(P)>=1+59+180=240.                              (6.1)

Let V=V(P), F=F(P), and a_j=f_j(P). The already written and independently
reviewed [sharp reciprocal transport](ideal_ra6_sharp_reciprocal_20261018.md)
proves the universal relation

    a_4-F <= (136/99)V.                              (6.2)

Its unchanged inputs are the ideal-six face identities, triangle counts,
good-edge/cubic-cusp transport, and the exact formal reciprocal maximum
103/132. No reoptimization or changed counting convention is needed here.
Each five-facet has at least59 four-facets and every four-face has two
five-cofaces, so 2a_4>=59F. Consequently

    (57/2)F <= a_4-F <= (136/99)V,
    V >= (5643/272)F.                                (6.3)

At F>=240 this gives

    V>=ceil(1354320/272)=4980.                        (6.4)

Indeed 272*4979=1354288<1354320 and 272*4980=1354560.
The older4462 count used the same relation(6.2) with only56 facets per
five-facet and F>=227. The new conclusions depend on the new large-facet
theorem and Sections1-5, not a numerical improvement to transport.

## 7. Consequences for mixed vertices in dimensions eleven and twelve

Retain the exact published Nikulin--Khovanskii and Dufour inputs and their
source-access qualifications from the [five-face density report](right_angled_mixed_five_face_density_20261018.md).
In dimensions11 and12 the strict average number of facets of a five-face
is less than20. Every finite-volume right-angled five-face has at least16
facets and every entirely ideal one now has at least59. If I_5 counts ideal
five-faces and N_5 counts all five-faces, then

    16N_5+43I_5<20N_5,
    I_5/N_5<4/43.                                    (7.1)

For dimensions9 and10 the same substitution into the strict average30
gives I_5/N_5<14/43. All these fractions are strict.

For six-faces in dimensions11 and12 the strict average number of facets is
less than42. The published general finite-volume minimum is27 and (6.1)
gives the entirely ideal minimum240. Therefore

    27N_6+213I_6<42N_6,
    I_6/N_6<5/71.                                    (7.2)

These improve the earlier necessary fractions1/10 and3/40. They do not
by themselves exclude any ambient dimension.

One can also update finite/ideal vertex ratios without changing any
incidence convention. Let v_f,v_i count finite and ideal vertices. Every
nonideal five-face contains a finite vertex, each finite vertex belongs
to binom(n,5) five-faces, and (7.1) gives

    binom(n,5)v_f>(39/43)N_5.                         (7.3)

The unchanged edge/NK inequalities in the old density report are

    N_5>(7/2)v_f+(3584/11)v_i,   n=11,
    N_5>6v_f+1024v_i,           n=12.

Substitution in (7.3), with binom(11,5)=462 and binom(12,5)=792, yields

    20669v_f>13312v_i,          n=11,
    5637v_f>6656v_i,            n=12.                 (7.4)

These are necessary joint inequalities, also meaningful when v_i=0;
ratio notation would require v_i>0. They are not contradictions.

## 8. Verification and remaining gap

The new capacity lemma and all arithmetic above are written here in full.
They have not been executed, formalized in Lean, externally reviewed or
assigned a publication-priority claim. Independent mathematical review is
pending. The inherited Figure7 theorem, carrier theorem, acute-normal fact,
published small-three-polyhedron classification, published four-dimensional
facet minimum, and sharp reciprocal transport remain explicit dependencies.
The mixed consequences additionally use the published NK/Dufour inputs
with exactly their previously documented access limits.

In particular, the universal ideal-six upper dimension bound remains6 and
the mixed finite-volume right-angled bound remains12. These arguments
increase necessary facet and vertex counts but supply no opposing upper
count or sufficient density contradiction. A dimension exclusion still
needs a quantitative restriction not satisfied by arbitrary large boundary
data; no such final restriction is claimed here.

