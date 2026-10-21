# Forced mixed triangles and a canonical four-face in dimensions eleven/twelve

18 October 2026. **New written necessary conditions; independent review pending.**
The ceiling12 is unchanged. This packet derives a substantially narrower
finite/ideal vertex ratio and forces a specific mixed four-face in either
dimension11 or12. It uses no restricted four-facet spectrum, arithmetic
field, facet census, numerical optimizer or Lean compilation.

For dimension12, write f for the number of finite vertices, c for the
number of ideal vertices, and a_j for the number of j-faces. Then

\[
 \boxed{\frac{1474}{21}c<f<\frac{1298}{9}c.}                 \tag{1}
\]

If T counts triangular two-faces having exactly one finite vertex, then

\[
 \boxed{T>\frac{331}{437}a_2,
        \qquad T>\frac{1986}{59}f.}                    \tag{2}
\]

Every hypothetical polytope in dimension11 or12 has a four-face isometric
to the ten-facet right-angled polytope with five finite and five ideal
vertices, described explicitly below. This is a forced-face theorem,
not a nonexistence argument. In particular a clique of four triangles at
a finite vertex is geometrically possible.

## External geometric inputs and exact scope

We use the standard Lorentz normal description of an acute-angled
polytope, simplex links at finite vertices, rectangular box links at
cusps, and inheritance of right angles and finite volume by faces.
Truncating ideal vertices gives a simple Euclidean polytope, to which
Dehn–Sommerville applies. These facts and this truncation method are in
[Alexandrov2023, Sections2–5](https://arxiv.org/html/2303.09533v1).
We also use the published lower bound82 for the number of facets of a
finite-volume right-angled eight-polytope, recorded there in Section6.
The absence of compact right-angled polytopes above dimension4 ensures
that this noncompact lower bound applies to every eight-face.
The remaining numerical identities and arguments are derived here.

The old [five-face density bound](right_angled_mixed_five_face_density_20261018.md)
gave only f/c>512/437 in dimension12. Its hypotheses are compatible with
(1); the new lower bound is1474/21. The ideal5/ideal6 carrier size bounds
alone are not used in the new proof.

## The dimension12 weighted identity

Set lambda=2123/1280. The exact identity is

\[
 66f+11264\lambda c-\frac{22}{5}a_2
 =\frac{-22a_7+88a_8-198a_9+308a_{10}-308a_{11}}5<0.    \tag{3}
\]

Here is a compact coefficient certificate for its Dehn–Sommerville step:

\[
\begin{aligned}
f&=a_6-5a_7+14a_8-28a_9+42a_{10}-42a_{11}-462c,\\
a_2&=15a_6-74a_7+206a_8-411a_9+616a_{10}-616a_{11}-2684c.
\end{aligned}                                                   \tag{4}
\]

To reproduce (4), put b_0=f+2048c, b_11=a_11+c and
b_k=a_k+2^(11-k) binom(11,k)c for1<=k<=10. Form
H(t)=(t-1)^12+sum_(k=0)^11 b_k(t-1)^k and equate the coefficients
of H(t) and t^12 H(1/t). Solving the six resulting independent linear
equations gives (4). Thus the displayed certificate is a finite rational
identity, not a claim based on a feasible optimizer assignment.

For the strict sign in (3), count incidences between seven/eight-faces:
5a_7>=82a_8. Every ten-polytope has at least11 facets, so
3a_9>=11a_10. Consequently

\[
 -22a_7+88a_8\le-\frac{1364}{5}a_8<0,
 \qquad -198a_9+308a_{10}\le-418a_{10}<0.
\]

The remaining term is negative. Hence (3) is strict.

A polygon with a finite and b ideal vertices satisfies

\[
 a+b\ge3,\qquad a+2b\ge5.                               \tag{5}
\]

The second inequality is the strict positivity of its hyperbolic area,
with finite angles pi/2. Summing (5) gives, with N_f=66f and N_i=11264c,

\[
 a_2\le\frac{N_f+N_i}{3},\qquad
 a_2\le\frac{N_f+2N_i}{5}.
\]

Combining each with N_f+lambda N_i<(22/5)a_2 gives (1). In particular,
these inequalities require both vertex types; neither a large finite
vertex count nor a large cusp count on its own is a contradiction.

## Triangle density

Give each finite polygon vertex weight1 and each ideal one weight lambda.
The triangles counted by T have weight1+2lambda. Every other polygon
has weight at least3+lambda. Thus

\[
 (3+\lambda)a_2-(2-\lambda)T
 \le N_f+\lambda N_i<\frac{22}{5}a_2,
\]

which gives the first inequality in (2).

A slightly stronger counting device retains the finite corners. For any
polygon other than the counted triangle,

\[
 a+\lambda b-\frac{22}{5}\ge\frac{331}{3840}a.           \tag{6}
\]

This follows by checking the four minimal nonnegative integer pairs
(0,3),(2,2),(3,1),(5,0); all other nonexceptional pairs satisfying (5)
dominate one of them. A counted triangle has deficit106/1280. Summing
(6) and using the strict weighted bound yields

\[
 331(N_f-T)<318T,
\]

which gives T>(1986/59)f. The average number of counted triangles at a finite vertex is therefore
greater than1986/59; the three-face argument below gives the stronger
clique conclusion.

## Three counted triangles force an actual standard three-face

At a finite vertex q, normalize q squared=-1 and let n_1,...,n_n be its
orthonormal outward facet normals. Define a graph G_q on the n incident
edge directions: ij is an edge when their two-face is a triangle with
one finite vertex. Its ideal endpoints are q-n_i and q-n_j, and its
opposite supporting wall has the uniquely determined unit normal

\[
 m_{ij}=q-n_i-n_j.                                      \tag{7}
\]

Suppose ij,ik,jk are all present. Their three normals are mutually
orthogonal, and

\[
 q'=2q-n_i-n_j-n_k,
 \qquad (q')^2=-1.                                     \tag{8}
\]

This is an actual finite vertex. To check every other supporting wall,
write its normal as r=a q-sum_l c_l n_l. Since q is strict for every wall
not through q, a>0; acute products with the old n_l give c_l>=0.
For r distinct from the three m in (7), acute products with them imply
c_i+c_j<=a, c_i+c_k<=a, c_j+c_k<=a. Therefore

\[
 r\mathbin\cdot q'=-2a+c_i+c_j+c_k\le-a/2<0.
\]

The three selected m normals and the n_l outside {i,j,k} vanish at q'
and are orthonormal. The original n_i,n_j,n_k are strictly negative there.
This checks both feasibility and the full finite vertex star.

The three-face cut out by the n_l outside the triple is exactly the
triangular bipyramid with finite vertices q,q' and ideal vertices
q-n_i,q-n_j,q-n_k. Indeed the six known inequalities in its Lorentz
four-space have this finite-volume convex hull; all five displayed
vertices belong to the ambient polytope, so no additional wall cuts it.
Its six polygon faces are precisely triangles with one finite vertex.

Let B count such three-faces. Every one has two finite corners, and every
triangle in G_q supplies exactly one such face. Thus

\[
 \sum_{q\ {\rm finite}} e(G_q)=T,
 \qquad \sum_{q\ {\rm finite}} t(G_q)=2B.              \tag{9}
\]

## A universal three-face gap makes four triangles meet

Every finite-volume right-angled three-polytope has at least6 facets.
Equality holds only for the preceding triangular bipyramid. Here is an
elementary proof, avoiding a classification input. If it has p finite
vertices, s cusps and F facets, Euler and vertex degrees3/4 give
p+2s=2F-4. Summed polygon areas give p+6s>=20; with p>=0 this excludes
F<=5. If F=6, the same inequality gives s>=3. The total number of polygon
corners is24-2s, so at least2s of the six facets must be triangles; hence
s<=3. Thus s=3,p=2 and all six facets are triangles. Each has at most one
finite vertex, while the total finite-corner incidence is3p=6, so each
has exactly one. At either finite vertex (7) applies to all three pairs,
forcing precisely the previously derived geometric bipyramid.

All other three-faces consequently have at least7 polygon facets. In
ambient dimension12, their Euler identities and actual incidences give

\[
 10a_2=2a_3+110f+28160c,
 \qquad B\ge7a_3-10a_2.
\]

Equation (3) gives a_2>15f+4246c, and hence

\[
 B>7590c-10f>\frac{2515}{59}f.                       \tag{10}
\]

The second step uses f<(1298/9)c. Consequently (9) forces
sum_q t(G_q)>(5030/59)f>64f.

A graph on12 vertices without a four-clique has at most64 triangles.
For a short proof, maximize the weighted sum of triangle products for
nonnegative vertex weights summing to1, choosing a maximizer with
smallest support. Two nonadjacent supported vertices occur in no common
triangle, so concentrating their combined weight on the better one
preserves or improves the sum and contradicts minimal support. The
support is therefore a clique of size at most3. Its triangle product
is at most1/27. Uniform weights give t<=12^3/27=64.

Some G_q therefore contains a four-clique. Quantitatively more than
209/1534 of the finite vertices have one: otherwise bounding their
triangle counts by220 and all remaining counts by64 contradicts (10).

## The forced four-face, with all its vertices

Let S={i,j,k,l} be the four-clique. In the Lorentz five-space spanned by
q and these four n normals, the four old walls and the six m_ab from
(7) define ten walls. With x=tq-sum_(a in S)y_a n_a, their inequalities
are y_a>=0 and y_a+y_b<=t.

The vertices of this projective polytope are precisely

* q and the four finite vertices 2q-sum_(a in A)n_a for the three-element
  subsets A of S;
* the four cusps q-n_a and the cusp 2q-sum_(a in S)n_a.

The former are actual vertices by (8). The first four cusps already
belong to the counted triangles. The last cusp is m_ij+m_kl; its product
with any other outward normal is nonpositive by acuteness, including the
selected normals themselves. It therefore also belongs to the ambient
polytope. The ten known halfspaces have the convex hull of these ten
points as their entire finite-volume region. Convexity shows that no
additional ambient wall cuts the four-face. This proves the claimed
canonical ten-facet face with five finite and five ideal vertices.

The conclusion is not that the entire ambient polytope has this local
facet spectrum. Additional facets outside this four-face remain
unrestricted. In particular this proof makes no inference that a larger
clique is itself impossible.

## The four-face minimum and a quantitative occurrence bound

The same elementary three-dimensional inequalities give F+2s>=12 and
F>=6, hence F+s>=9. In a four-polytope, any three-facet with F polygon
facets and s cusps forces at least1+F+s ambient facets: its neighbors and
one distinct cusp-opposite per cusp. Thus every right-angled four-polytope
has at least10 facets. This recovers the published minimum in
[Potyagailo–Vinberg2005, Section4, equation(7), p71](https://ems.press/content/serial-article-files/43032).

Equality is rigid here. If there are10 facets, then every three-facet has
F+s=9. Its inequalities F>=6 and F+2s>=12 force s=3,F=6. Every such
facet is therefore the geometric bipyramid proved above. All polygon
faces have one finite vertex; choosing any finite vertex gives a complete
four-clique in its triangle graph. The preceding explicit ten-wall
argument identifies the whole four-polytope. Consequently every other
four-polytope has at least11 facets; uniqueness of the ten-facet equality
case has been proved here, rather than inferred from attainment of the
published minimum.

Let C count the canonical four-faces. Double counting yields the new
necessary inequality

\[
 C\ge11a_4-9a_3.                                      \tag{11}
\]

There is also an unconditional positive frequency. For a graph on12
vertices, delete an edge of a four-clique until none remains. At most
k_4(G) deletions are required, and each deletes at most10 triangles.
The remaining graph has at most64 triangles, so

\[
 t(G)\le64+10k_4(G).
\]

Each canonical four-face has exactly five finite vertices, with one
four-clique at each, and every four-clique supplies one such face. Thus
sum_q k_4(G_q)=5C. Equations (9) and (10) give

\[
 \boxed{C>\frac{627}{1475}f.}                           \tag{12}
\]

Neither (11) nor (12) is a contradiction. They quantify the attachment
patterns which any proposed dimension12 boundary must support.

The minimum-sum equality in dimension4 and the sixteen-facet equality in
dimension5 are classified in the separate
[low-facet rigidity packet](right_angled_mixed_low_facet_rigidity_research_20261018.md).
They are not premises of the weighted-density proof above.

## Dimension11 coefficient certificate and consequences

The same calculation in dimension11 gives lambda=2123/1280 and

\[
\begin{aligned}
f&=2a_6-12a_7+42a_8-112a_9+252a_{10}-252c-504,\\
a_2&=25a_6-148a_7+515a_8-1370a_9+3080a_{10}-1220c-6160,\\
55f+5120\lambda c-\frac{22}{5}a_2
 &= -\frac{44}{5}a_7+44a_8-132a_9+308a_{10}-616<0.
\end{aligned}
\]

For the sign use4a_7>=82a_8 and2a_9>=11a_10. Repeating the preceding
polygon steps yields

\[
 \frac{268}{7}c<f<\frac{236}{3}c,
 \quad T>\frac{331}{437}a_2,
 \quad T>\frac{1655}{59}f.
\]

Here the three-face identities give

\[
 a_3=\frac92a_2-\frac{165}{4}f-5760c,
 \quad B\ge7a_3-9a_2>
 3105c-\frac{15}{2}f>
 \frac{7545}{236}f.
\]

Thus sum_q t(G_q)>(7545/118)f, whereas a four-clique-free graph on11
vertices has at most11^3/27 triangles. The same four-face is forced.
The edge-deletion argument, using the integer bound49 for the remaining
triangle count, also gives C>(1763/5310)f. The four-face gap reads
C>=11a_4-8a_3 in this dimension.

## Verification boundary and remaining obstruction

The [compact certificate](right_angled_mixed_weighted_triangles_research_20261018.json)
and [standalone exact checker](../scripts/check_right_angled_mixed_weighted_triangles_research_20261018.py)
recompute the Dehn–Sommerville equations by rational Gaussian elimination,
derive the weights directly from their finite/cusp coefficients, and check
every displayed ratio and density. They enumerate every independent
active-wall set for the canonical three- and four-polytopes,
checking all20 and210 wall choices, with respectively2/3 and5/5
finite/ideal vertices.
The checker uses only the Python standard library and performs no proof
compilation. Replay it without modifying the certificate with

```sh
python3 -B scripts/check_right_angled_mixed_weighted_triangles_research_20261018.py
```

The displayed rows and rational eliminations were additionally checked
with a separate small symbolic calculation. The independent standard
library calculation caught and corrected a draft finite/cusp convention
error before this certificate was produced; the final weights are
2123/1280 in both dimensions. The certificate does not encode the universal
Lorentz geometry or constitute an independent mathematical review.
No heavy numerical job, proof compilation, frozen-result modification or
new chamber construction was performed. The geometry, equality case and
local-to-global counting are written proofs awaiting independent review.

These conditions leave the intervals in (1) and its dimension11 analogue
nonempty. The forced four-face is a known realizable four-polytope; the
remaining question is whether it can occur in a complete right-angled
boundary of dimension11 or12 with the required adjacent finite/cusp stars.
The packet therefore supplies a precise universal attachment obstruction
for further work, while retaining the current published ceilings12 and6.
