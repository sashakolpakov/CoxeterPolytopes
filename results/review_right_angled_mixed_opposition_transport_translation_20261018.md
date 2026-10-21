# Independent review of the opposition-transport translation

18 October 2026. **PASS for the complete mathematical/source translation.**
I read all1,217 lines of
[the manuscript](../manuscripts/right_angled_mixed_opposition_transport.tex),
[its proof map](../docs/right_angled_mixed_opposition_transport_proof_map.md),
the mapped actual-support, symmetric-quad, nu10-transition and
reciprocal-square source notes, and the separate small-carrier,
carrier-fraction, actual-support and quad-charge written reviews.
I checked all five figure constructions and captions. No correction to
the TeX is required by this review. I ran no checker, optimizer, Lean
process or TeX build and left the TeX unchanged.

The verdict retains the stated boundaries: the178 interval inequalities
and nine quadratic-field rows use their saved exact arithmetic evidence;
this review did not execute them again. The geometric premises listed
in the manuscript are explicit. The translation proves two conditional
dimension-eight theorems and necessary dimension12 incidence conditions;
it makes no new unconditional dimension bound.

## Actual-wall feasibility and the complete small-three-face list

The projection lemma proves actual feasibility with the outward-normal
signs used in the paper. Starting from a strict future interior point,
projection off orthogonal normals adds positive multiples of outward
normals, decreases every other wall product, preserves negative square
and future orientation, and leaves the selected products zero. The
two-normal positive Z-Gram argument then excludes products in(-1,0).
Thus the stated distinct-normal gap is derived from actual walls; an
intersection outside the chamber is not silently declared a ridge.

The bad-wall equivalence checks its null endpoints against the selected
walls as well as all other walls. The P3 proof establishes strictness
against every unselected wall using the three pair-sum inequalities,
then uses the actual five feasible model vertices and convexity. It
does not infer an actual opposite from a formal Gram matrix alone.
The same distinction is maintained for the symmetric quad and its D3
carrier. Their extra ideal endpoint is a future null sum of actual
normals, with the selected inequalities checked. The general quad
parametrization has positive factors a-c,a-d, and its discarded minus
branch really cannot give both c,d>0.

The attachment proof correctly normalizes both shared null rays, retains
the positive parameter s, and uses the two products-s,-1/s to force
s=1. Its final null point and six-wall convex hull prove the whole P3
conclusion. This supplies the geometric attachment statement needed in
the later combinatorial counts.

The nu<=10 argument is complete, rather than an assumed census:

* Euler, finite/cusp degrees and polygon area give the displayed p,c,s
  identities and nonnegative integer polygon excess.
* The minimum s=6 is proved. At s=6 the corner count and excess force
  six bad triangles and P3.
* Outside P3, s+c<=10 and s+2c>=12 leave exactly(8,2) or(7,3).
  Attachment excludes a triangle in the zero-excess(8,2) case; the total
  corner count then forces eight(3,1) quads and ideal degree one at every
  finite vertex.
* At(7,3), at least three bad triangles have distinct ideal--ideal edges.
  Only the single excess-one all-ideal triangle can contain all three.
  Corner counts force the remaining three quads. The six finite--ideal
  edges and three bad triangles force ideal degrees0,2,2,2 and the
  finite-edge star. No metric uniqueness of Q8 is used.

In R7 the three bad facets are opposite pairs in their rectangular cusp
links. Their products-1 give the three positive sine products equal to
one, hence cosh length sqrt2 on every finite edge. The ambient small
carrier argument correctly finds an R7 leaf or center around each
finite edge. It does not require an extra hypothesis on its endpoints.

## The first conditional theorem and the fraction estimate

Finite edges and P3 oppositions are actual undirected relations. The
P3 coordinate formula gives distinct opposite vertices for distinct
triples; the reversed move belongs to the same three-face. Their
weighted displacement sums therefore cancel with the stated
multiplicities. The Delta expansion uses C(k,3)=k*C(k-1,2)/3 correctly.
The nine base rows and their increasing first differences cover k<=8;
the nonnegative higher-k tail covers all remaining n>=9. Its null
k9,m0 row has positive time coefficient21/5 and is not zero. Evaluation
against one common future vector yields the contradiction.

For the fraction theorem, g=1/H gives reciprocal energy ratios, so the
global row-minus-degree sum has the stated nonnegative sign. Replacing
signed coordinates by absolute values and raising their norm to one
enlarges every reciprocal term; the finite and triple denominators
remain positive. The selected row is complete at a small vertex, while
the unrestricted row excess bound uses only the selected bounded move
lengths. The resulting b>=f/385 and p(F)<=14(nu-10) incidence conversion
give E>=f/5390 with no leftover energy weighting.

The physical-variance sign is correct. At an actual finite-group point,
y^2<=x/m<=h/m, and at an actual ideal-group point ky^2<=1-x<=1-ell.
Only there does the paper lower positive denominator parameters and
inflate variance to obtain the rational upper function. It enlarges the
mean interval afterward. A negative numerator at an added nonphysical
point cannot invalidate maximization over a superset. I independently
differentiated both displayed functions: their derivative numerators
are exactly the stated decreasing affine expressions, so clipping the
critical point gives the global interval maximum. The grid, sqrt2
lower bound, all176 ordinary boxes and two degenerate boxes are specified.
This matches the reviewed178-box interface and its strict1/8 deficit.

## Actual opposition support and the inherited cap

The cross-wall identity is an equality of actual outward normal vectors,
so it preserves bad pairs in both directions. Internal triples persist;
outside pairs explicitly require a new wall. The manuscript does not
infer that wall merely from a five-clique. Incoming small triples
therefore give true-twin triples at the recipient, with the correct
binomial capacity. The invariant h=3q-sum(n_i), its square n-9 and
its componentwise constancy follow from the displayed new normals.

The K9 exit expansion is correct. In a graph on at most12 vertices,
c_i<=27,d_j<=min(M,55),O>=M give the displayed upper-2M. Every
nontrivial actual component must have a noncausal displacement row, so
it reaches a K9-free corner. No quantitative exit density is inferred.

The reciprocal-square section retains exactly the reviewed
516992/3125 cap. In dimension12, h is spacelike, permitting a future
null z perpendicular to h separately on each component. Its coordinates
have sum3 and squared norm1; centering at1/4 gives sum zero and squared
norm1/4. The zero-sum bound |d|<=3/4, the repeated-index moments and
the quartic reciprocal identity have the stated signs and constants.
The remainder multiplier is bounded by13 at r>=1/2. Restricting the
positive formal row to actual neighbors gives an upper bound. The
componentwise choices of z are compatible with summing boundary
inequalities; the paper does not try to choose one z for unrelated
components. Isolated vertices are treated separately.

Internal reciprocal pairs contribute at least2 and boundary ratios at
least1/4. Degree190 and incoming-complete capacity20 consequently give
102344/3125 and39844/3125 with the correct signs. The clean-edge status
switch, exceptional empty-corner count, and disjoint incidence sets give
the stated42E budget. The cap is not replaced by the later sharpening,
and no infeasibility claim about an older scalar relaxation is made.

## Quad reversal, unrestricted errors and the second conditional theorem

The unrestricted quad rule uses one actual quad per missing pair,
weighted by d_ij/2. The two exact cross-wall identities preserve its
common bad-neighbor set on reversal, so this weight is symmetric.
Different D3 carriers deliberately contribute through d_ij; they do not
create different opposite vertices. The resulting time contribution is
one per backed path and its spatial contribution is one at each missing
endpoint. Thus D_actual, its cancellation, and the sqrt2*M error bound
have the correct unrestricted multiplicities.

The second dimension-eight theorem uses a different, explicitly stated
weight: each actual quad once with1/2, independent of carrier count.
Its hypotheses include a finite vertex and **global absence of every
finite--finite edge**, as well as nu<=11 on every three-face. The
nu11 classification uses3p<=4c globally within that carrier, obtains
(s,c,p)=(7,4,2), and puts both finite vertices opposite in its quad.
All carriers are then P3 or D3, and the corner-graph complement is a
matching. The triangle and quad contributions give precisely the
displayed a,b_un,b_pair, with nonnegative margin at n9 and positive
minima975,3135,7384 at n10,11,12. Every row is nonzero future causal.
The reduction of n>12 to a twelve-face through a finite vertex preserves
both global hypotheses. The caution following the proof correctly
forbids replacing the global edge condition by a condition at one corner
or one opposition component.

For the unrestricted error charge, distinct bad triangles use distinct
ideal--ideal edges outside P3, giving4Z<=nu-2. The special integer
nu12 step is retained before the larger-nu linear estimate. At nu11,
the unbacked(8,3,6) case really has all three cusps neighboring the
chosen corner. Simplicity and convexity put an edge between two of
those cusps in their common polygon, forcing exactly their triangle.
This proves E_II=2,E_FI=8,E_FF=5 and the one-corner restriction. The
factor binom(n-1,2)/5 counts actual ambient finite-edge incidences;
it assumes neither globally ideal directions nor bounded larger carriers.
The companion K3,3 and larger-carrier estimates preserve the same scope.

The maximum-height restrictions remain local. Their norm costs, the
strict289>288 comparison and the5/9 common-K5 cost are correct. The
K4-joined-independent8 example satisfies all triple inequalities with
norm8/9, yet each alleged quad pair leaves two disjoint actual P3
triples and is impossible at that maximum. Its112 paths are therefore
unbacked. Neither a general global triangle average nor a universal
triangle-count56 claim is inferred.

## Figure audit and final scope

1. The P3 drawing has the correct five vertices and nine incidence
   edges; the finite-apex opposition is dashed and explicitly not an
   edge. The three ideal labels and opposite formula are correct.
2. The R7 panel draws only its finite-edge star, with three equal
   cosh lengths sqrt2. The companion causal cone is explicitly schematic
   and the caption states the common-future-vector test.
3. The cross-edge panel marks directions rather than hyperbolic
   vertices, preserves the displayed blue cross pairs, and warns that
   outside pairs require new walls.
4. The symmetric-quad panel has the correct alternating finite/ideal
   incidence cycle, formula and cosh diagonal. Its second panel is a
   direction graph and retains the actual-quadrilateral hypothesis.
5. The maximum-height example is correctly labeled K4 joined to eight
   independent directions. Its4/9 and1/9 coordinates, norm8/9, two
   disjoint complementary triples and112 unbacked paths agree with the
   proof. It is explicitly not asserted to be a realizable polytope.

The translation preserves every reviewed hypothesis and does not promote
an abstract graph, formal Lorentz candidate or selected corner into an
actual global configuration. Its two conditional theorems and necessary
dimension12 estimates coexist with the unchanged unconditional mixed
ceiling12. No mathematical translation issue was found, and no source
or figure edit was made. Typesetting and any later exact replays remain
separate coordinated, resource-gated work.
