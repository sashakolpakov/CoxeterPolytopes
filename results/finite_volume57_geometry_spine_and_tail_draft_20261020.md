# Finite-volume 57: geometric spine, higher tail, and completeness audit

Staging text for the existing paper, 20 October 2026. No primary TeX, frozen
source, certificate, or build was changed. This document supplies prose-ready
geometric arguments and identifies the precise remaining exposition burden.
It is not a new dimension theorem, an execution report, or a fresh independent
verification of every historical atlas case. The complementary finite/cusp
algorithm fragment is being prepared by `batch_guard` in
`results/finite_volume57_algorithm_specifications_draft_20261020.md`.

The intended order is: common geometric counting framework; the actual-face
payment theorem; separate finite and ideal budgets; their exact computational
proofs; the ambient58 contradiction; an explicit higher-tail lemma. The
higher tail cannot be supplied by an unstated inheritance of Coxeter angles
to arbitrary metric faces.

## 1. Actual walls and actual corners

Let P be a finite-sided finite-volume Coxeter polytope in H^n. Work in a
Lorentz space of signature (n,1), on one sheet of the unit hyperboloid, and
normalize each outward wall root to square2. The chamber inequalities are
`(x,e)<=0`; distinct wall roots have nonpositive products. Fix a strict
interior point p. Thus `(p,e)<0` for every wall root e.

We use the standard finite-volume local-structure theorem: every finite
vertex is simple, every positive-dimensional k-face is contained in exactly
n-k facets, and the link of an ideal vertex is a product of Euclidean
Coxeter simplices. The last statement allows several affine factors and
does not say that an ideal vertex is simple. These hypotheses are stated
in the primary sources cited by
`fresh_finite_volume_simplex_truncation_20261017.md` and
`fresh_simple_cusp_distance_bound_20261017.md`.

In particular, if F is an actual two-face, its n-2 supporting roots form
a positive definite set W. Each side of F is cut out by one further wall,
and W together with its side root is positive definite: its common edge
has an interior hyperbolic point. At a finite corner W together with its
two adjacent side roots is positive definite. At an ideal corner it is
positive semidefinite. Every principal Gram has negative index at most1.
These facts concern the original chamber, not a truncation at its cusps.

The full star of F is indefinite. Indeed, projection perpendicular to the
positive span W takes its side walls to the side walls of a nondegenerate
finite-area polygon in a Lorentz three-space. Their normals span that
three-space. This is the dimension-independent geometric input for both
triangle and quadrilateral payments.

## 2. Two elementary feasibility lemmas

**Elliptic projection.** A positive definite collection of original wall
roots has an actual finite common point in the chamber closure.

To prove this, write its Gram as B and put `b_i=-(p,e_i)>0`. A positive
definite symmetric matrix with nonpositive off-diagonal entries has an
entrywise nonnegative inverse. Set

`q=p+sum_i (B^(-1)b)_i e_i`.

Then `(q,e_i)=0`, `q^2=p^2-b^t B^(-1)b<0`, and q is future pointing.
Every other wall e has `(q,e)<0`, since its product with p is negative
and the added coefficients are nonnegative. Normalize q to the hyperboloid.
This proves feasibility against **all** chamber walls, not just a displayed
principal subdiagram.

**Affine feasibility.** A connected positive semidefinite singular Coxeter
block supplies a nonzero feasible future null vector.

Its kernel has a strictly positive vector u. The root combination
`nu=sum_i u_i e_i` has square0. It is nonzero because `(p,nu)<0`, and it
pairs nonpositively with every other chamber wall. It therefore defines
an ideal point of P. A point of the ideal sphere is extreme in the closed
Klein ball, hence also in P. It is an ideal vertex. In its full cusp
diagram the connected affine block is an entire irreducible affine factor:
a proper principal block of an irreducible affine diagram is positive
definite. This also explains why a connected displayed affine block cannot
quietly grow into a larger actual affine factor.

**Closed-opposite-side lemma.** If U is either opposite pair of side roots
of an actual quadrilateral, `W+U` is indefinite. Positive definiteness
would, by elliptic projection, give a common finite point of its opposite
closed sides. A singular positive semidefinite block would, by affine
feasibility, give a common ideal point. Both are impossible: the four
vertices of the closed quadrilateral are distinct and opposite sides are
disjoint, including when some vertices are ideal. Compactness is not used.

## 3. Why a finite native-core classification applies

Choose an inclusion-minimal indefinite principal block in each of the two
opposite-pair blocks. It is connected. It contains both marked black roots,
since deleting either of them leaves a subdiagram of an actual positive
definite edge diagram. Thus both marked deletions are strictly elliptic.

If a minimal indefinite block has a proper affine block A, connectivity
gives an outside node joined to A. Its nonpositive contact row pairs
strictly negatively with the positive kernel of A. The block obtained
by adjoining that one node is already indefinite. Minimality makes it the
whole block. Otherwise all proper principal blocks are positive definite.
The published Lanner/quasi-Lanner classification consequently applies:
the order is at most10, and the branch with all proper blocks elliptic
has order at most5. An order-two dotted pair is included in this statement.

The finite reconstruction of orders6--10 has a particularly simple complete
domain. Each such diagram is one extra node attached to an irreducible
affine source. An edge of label at least5 could be extended to a connected
five-node proper principal block, impossible by finite/affine classification.
A parallel pair with a joined third node would already be indefinite;
a dotted pair is itself indefinite. Hence every attachment is in the
alphabet `{2,3,4}`. The 24 affine sources and **all** their nonzero
attachment vectors give146262 tests; exact classification and isomorphism
give26 diagrams. The order-ten ceiling is a published input, not a result
of stopping a search at ten. These domain and correctness arguments are
written in `fresh_high_rank_quasilanner_census_20261017.md`.

For a quadrilateral the two minimal cores cannot be disjoint and
orthogonal, which would give two negative directions. They therefore
overlap or are joined. A shortest path between their opposite marked sets
has only white internal roots and ends at an adjacent pair of sides.
This is an actual corner path. The complete higher-order catalogue has
marked radius at most7; this first gives distance15, before the sharper
local reductions are used.

For a triangle every indefinite core contains all three side roots,
since deleting one lies in a positive semidefinite actual corner. A
minimal subtree joining the three marks gives two distinct corner paths
avoiding the third mark. For low orders their length is at most5. The
1025 marked triples of orders7--10 have two retained payments except for
four explicitly classified triples, each with one elementary payment and
a distinct protected E8 alternative. The current triangle theorem uses
finite paths at most7, or the protected whole exceptional affine factor
of order at most9. These alternatives are proved in
`fresh_finite7_ideal5_triangle_payments_20261017.md`; the later finite6
change retains its exceptional length6 alternatives through a protected
nearby whole leaf. This is a theorem about distinct actual corners, not
the sum of three marked-pair weights inside one diagram.

## 4. Restoring every supporting white

A calculation on a displayed core has no whole-component conclusion
until arbitrary additional supporting whites have been controlled. The
following tools account for direct contacts and indirect tails.

**Lorentz projection.** Let G be a nonsingular core Gram with one negative
direction. Its span contains a timelike vector, so its perpendicular space
is positive definite. A prospective square-two root with core contact
column q has projection square `q^t G^(-1)q` and residual square
`2-q^t G^(-1)q>=0`. A projection above2 excludes the root. At equality
the root lies in the old span, so its products with every other old-span
root are fixed exactly. It cannot be replaced by a freely chosen extension.

**Radical seal.** If an indefinite displayed core has a Gram radical z,
then `sum_i z_i e_i=0` as an actual root relation. This combination is
null and perpendicular to a timelike vector in the core span, and the
perpendicular space is positive definite. Consequently every additional
white has `sum_i z_i (e,e_i)=0`. If all its permitted contact coordinates
have nonzero z-coefficients, a nonzero contact must meet both signs.
The appropriate short-path test either pays that branch independently
or forces all contacts to vanish. A component orthogonal to every added
root is whole: an indirect extension would have a first root contacting it.

This zero root relation must be distinguished from an affine positive
kernel, whose actual root combination is a **nonzero** null vector.
An extra root lying in a positive semidefinite corner with that affine
block is orthogonal to all its roots. Indeed a zero quadratic-form vector
of a positive semidefinite matrix is in its kernel, and all the individual
products being summed have the same nonpositive sign.

**Connected inverse bound.** If B is a connected positive definite
square-two Coxeter Gram of rank k, every inverse diagonal is at least
`k/(k+1)`. Delete the chosen root and let the remaining connected
components have ranks k_j. By induction, the unit-or-larger contact into
each component has energy at least `k_j/(k_j+1)`. Inverse positivity and
the inequality

`sum_j k_j/(k_j+1) >= (k-1)/k`

give a root Schur denominator at most `(k+1)/k`. The rank-one base is1/2.
Applying the same argument separately to a forest of total rank M, with
a unit-or-larger contact in every component, gives energy at least
`M/(M+1)`.

**Attached-forest bound.** Suppose additional whites Y contact a
nonsingular Lorentz core only through one position p, with inverse square
`r=(G^(-1))pp>1`. Take **all** new roots of the growing component, including
ones orthogonal to the core. Its induced Gram B is positive definite;
each component of Y meets the old small block. The full enlarged Gram
has at most one negative direction. Eliminating the old core therefore
gives

`B-r u u^t >=0`, hence `u^t B^(-1)u<=1/r`.

Together with the inverse bound this yields `M/(M+1)<=1/r`, or
`M<=1/(r-1)`. A bound on immediate neighbors alone would not prove this.
For the final58 cases `(old small rank,r)=(1,4/3),(2,3/2),(3,2)`, the
whole small component has rank at most4.

These are the common proofs behind the extension certificates. A matrix
table need only specify the permitted contacts, inverse/radical values,
and selected paying role; it must still establish that every possible
additional white falls under one of these alternatives.

## 5. The universal local payment theorem that must precede the counts

The paper needs an intrinsic definition of its one weight function w,
not the instruction “run the preceding version”. Its arguments are the
actual finite forest with its marked pair, or the actual varying affine
factor with its marked pair. Whole-component presence conditions belong
to the same corner; several witnesses do not multiply a presence coupon.
The definition can consist of classical all-rank formulas and finite
marked-orbit tables, with explicit priority for any contextual alteration.
An arbitrary-label rank-two component requires an analytic or dominating
type clause, not an unstated finite label cutoff.

For that defined rule the required theorem is:

> Every actual triangle has total weight at least2 and every actual
> quadrilateral has total weight at least1, for arbitrary supporting
> roots and without a simple-cusp hypothesis.

Here is the complete logical decomposition of its quadrilateral proof.
First use an independently valid short-path or actual-ideal-corner route
when available. Otherwise choose the two minimal opposite cores. A core
of order at most5 belongs to the complete low/high reduction, including
the unbounded black-black label in order3. Two disjoint higher cores are
excluded by the all-label index obstruction. Overlapping higher cores
are enumerated by every injective partial identification preserving all
edges and nonedges. Finally restore all supporting whites using Section4.
The outcomes are finite whole components, protected leaves/features,
whole affine factors, or named independent replacements.

For clarity, in the distance8 predecessor these exhaustive domains are
written, rather than inferred from the saved survivors:

* low/high joins include all possible native contacts after derived
  inverse bounds; the order-three family is excluded symbolically for
  every allowed real black-black product;
* higher overlaps include all overlap sizes and all consistent induced
  identifications (82197 before the strict edge tests); every omitted
  cross-copy edge is tested at its smallest nonzero magnitude;
* nonnegative nonpositive-square witnesses remain valid when a bond is
  increased or more acute contacts are added;
* every surviving marked core has an arbitrary-support extension proof,
  not merely a search for one extra node.

The resulting 2614 overlap survivors and the low/disjoint branches are
proved in `fresh_distance8_low_core_20261017.md`,
`fresh_distance8_high_disjoint_20261017.md`, and
`fresh_distance8_high_overlap_20261017.md`, then composed in
`fresh_finite8_ideal5_combined_payments_20261017.md`.

Later transfers form a directed proof of retention. For each changed role,
one must either show it was never an output of a particular earlier
universal route, or give an independently paid replacement. A census of
selected finite maps proves only its stated selected branch. In particular
the B13 endpoint has two genuine older affine-B consumers outside that
census, and the final D13 reduction has the genuine affine-D12 consumer.
Their repairs are necessary. They are not optional checks on sample faces.
The final independent short/ideal priorities are used before the selected
finite route, and no replacement may refer back to the payment it removes.

For a standalone presentation the full definition and exhaustive routing
statement are unavoidable. The current manuscript's historical names and
row counts identify the evidence, but do not replace those mathematical
objects. This audit has read the stated native-domain and support proofs;
it has not newly rederived every later transfer in the recursive atlas.

## 6. The final58 transfer in geometric terms

Only the following three role changes are new at ambient58: halve each
whole D13 long-end/short-leaf unit when three other whole A15/D15 factors
are present; add a half at D28(0,12) when two such factors are present;
add a half between an A12 endpoint and a D16 long endpoint under the same
two-factor guard. The separate B13 theorem deletes its endpoint unit.
All ideal weights stay unchanged.

The negative guard uses all `13+3*15=58` ranks. If the selected native
support has a sealed additional small component, it is impossible under
the guard. The complete affected selected contexts divide into radical
seals, impossible proper completions, inverse-isolated small factors, and
the attached-forest cases of Section4. In the last class the small whole
factor has rank at most4, again incompatible with the rank guard. This
shows how rank capacity converts a dimension-specific statement into
an arbitrary-support proof; it does not truncate the attached forest.

The older affine-D12 fork-leaf theorem requires a separate argument.
With its elliptic spectator forest Z and disjoint nonnegative contact
vectors x,y, its exact four-black Schur determinant gives, in the extreme
one-sided boundary,

`(x-y)^t G_Z^(-1)(x-y)=1`.

The guard makes Z exactly three whole A15/D15 factors. Strict edge
positivity leaves A15 simple contacts only at0,1,13,14 and doubled
contacts only at0,14; D15 permits only a simple long-end contact.
Every nonzero allowed component energy is at least15/16. The complete
disjoint two-sided difference table has no A15 energy1 and no additional
energy1 profile. Two active factors would already cost at least15/8.
Thus exactly one D15 has a single unit long-end contact, and both other
rank15 factors are inactive.

Those inactive factors are orthogonal to the full star and stay whole.
At the distinct opposite paying corner the active extension gives either
whole D28 at(0,12), or whole A12 endpoint against whole D16 long endpoint.
Its new half and the quiet D13 half sum to1. Nonextreme equality profiles
use an unchanged A1/D12 unit, the independent affine-B13 route, or two
distinct guarded D13 halves. The repaired B13 route uses its retained
D12 cross unit, so the argument is not circular. The full source is
`fresh_d13_three_rank15_context_transfer_20261018.md`.

The finite arithmetic is now exceptionally short. The predecessor's
only eight formal forests above377 are B13/D13 plus three A15/D15 factors,
each of old cost378. They all drop to377. New positive guards force exactly
the six forests D28+two rank15 factors or A12+D16+two rank15 factors.
Their costs are370.5 and375.5. All other forests are pointwise no larger
than before. This is the place for the finite algorithm theorem from the
companion fragment; no geometric realization of a formal maximizer is
needed.

## 7. The finite/ideal split and the shared global lower bound

At an ideal vertex write the affine factor ranks as `r_i>=1`, with
`sum_i r_i=d=n-1`. Its link is `L=product_i Delta^(r_i)`. Let

`V=product_i(r_i+1)`,
`H(t)=product_i(1+t+...+t^(r_i))=sum_j h_j t^j`,
`G=sum_j h_j min(j,d-j)`.

An actual two-face corner corresponds to an edge of L. Such an edge
varies in exactly one simplex factor. A fixed pair of its facet labels
determines one edge in that factor and a vertex in each other factor.
It therefore occurs exactly `product_(j!=i)(r_j+1)=V/(r_i+1)` times.
The marked pair restores the full varying affine factor; the other
factors appear with one node deleted. Consequently a factor pair budget
`(r+1)c(r)` gives the **actual** cusp-charge bound

`C_actual<=V sum_i c(r_i)`.

This is why the finite and ideal cases must be stated separately. A
nonsimple cusp does not have n incident roots, and summing one unweighted
copy of each affine pair table gives the wrong charge.

Here is a complete shape-curvature proof usable by every tail row.
Choose a generic linear height on the bounded Klein closure of P. At
each vertex v let E_v count the edges of its simple vertex figure and
E_v^0 the edges not cut by the induced level section. The published
incidence identities give

`I=sum_v E_v`, `2f_2(P)=sum_v E_v^0`,

where I counts all polygon vertices with multiplicity. Hence

`D5=sum_F(5-s(F))=sum_v[(5/2)E_v^0-E_v]`.

Orient each simple vertex figure by its induced generic linear function.
At a link vertex of index j, assign each uncut edge below the section to
its upper end and each uncut edge above it to its lower end. Its assigned
count is either j or d-j, and is at least their minimum. The index
distribution is the h-vector, so `E_v^0>=G_v`.

At the distinct global maximum and minimum of P the link is wholly
uncut. The index formula and symmetry give

`E_v-G_v=sum_(j>d/2) h_j(2j-d)`.

For a simplex or product of simplices every coefficient h_j is at least1.
Thus each extreme supplies at least
`s_d=sum_(j>d/2)(2j-d)=floor((d+1)^2/4)` extra uncut edges. The two
extrema together contribute `5s_d` to D5. A finite vertex has
`H=1+t+...+t^d`, `V=n`, and `G=floor(d^2/4)`; its curvature is

`beta_n=[5 floor((n-1)^2/4)-n(n-1)]/2`.

At an ideal vertex put `kappa_v=(5G_v-dV_v)/2`. We have proved

`D5>=beta_n v_f+sum_(ideal v) kappa_v+5s_(n-1)`.

This supplies a quantitative constant and remains strict when there are
no ideal vertices. It counts the two-extrema contribution once. The exact
incidence/index inputs and this derivation are in
`fresh_cusp_link_shape_curvature_20261017.md`.

At n58, `G_simplex=812`, `beta_58=377`, and `s_57=841`.
The independent ideal algorithm proves `kappa-C_actual>=36` for every
rank57 link. Nonnegative payments, with triangle total at least2 and
quadrilateral total at least1, bound D5 above by the total vertex charge;
larger polygons have nonpositive deficit. Thus

`377 v_f+sum C_actual >= D5 >=377 v_f+sum kappa+4205`,

which implies `0>=4205+36v_i>0`. This excludes ambient58.

## 8. A genuine higher-tail lemma

The conclusion above excludes one dimension. It does **not** imply that
dimensions59 and above are impossible. A section of a Coxeter chamber
by a metric face need not itself have Coxeter angles; neither that claim
nor an unproved dimension-descent is available here.

The minimal logical tail assertion is:

> For every integer n>=59, an explicit valid payment comparison, or the
> all-rank narrow-facet comparison below, contradicts finite volume.

The existing reports cover exactly the following intervals. Names in
the right column are report stems under `results/`, with suffix `.md`.
An entry is an internal theorem requiring its local-rule and complete
budget proofs; the table by itself is only the coverage argument.

|Ambient dimensions excluded|Existing internal endpoint/report|
|---|---|
|59|`fresh_a15_d15_dimension58_20261018`|
|60|`fresh_depth6_dimension59_20261018`|
|61|`fresh_order11_order16_dimension60_20261018`|
|62|`fresh_rank15_rank16_dimension61_20261018`|
|63|`fresh_sparse_coupon_order11_13_dimension62_20261018`|
|64|`fresh_order11_order17_dimension63_20261018`|
|65|`fresh_deep_finite8_dimension64_20261018`|
|66|`fresh_ideal5_weighted_dimension65_20261017`|
|67--68|`fresh_ideal5_weighted_dimension66_20261017`|
|69|`fresh_ideal5_weighted_dimension68_20261017`|
|70|`fresh_ideal5_weighted_dimension69_20261017`|
|71|`fresh_ideal5_weighted_dimension70_20261017`|
|72--74|`fresh_ideal5_weighted_dimension71_20261017`|
|75|`fresh_ideal5_piecewise_dimension74_20261017`|
|76--77|`fresh_ideal5_piecewise_dimension75_20261017`|
|78--80|`fresh_ideal5_piecewise_dimension77_20261017`|
|81--86|`fresh_ideal6_piecewise_dimension80_20261017`|
|87--91|`fresh_order16_piecewise_dimension86_20261017`|
|92--97|`fresh_ideal7_cross_dimension91_20261017`|
|98--108|`fresh_distance8_ideal_dimension97_20261017`|
|109--116|`fresh_selective_distance9_dimension108_20261017`|
|117--126|`fresh_selective_distance10_dimension116_20261017`|
|127--136|`fresh_selective_distance11_dimension126_20261017`|
|137--146|`fresh_selective_distance12_dimension136_20261017`|
|147--156|`fresh_selective_distance13_dimension146_20261017`|
|157--163|`fresh_distance14_unrestricted_dimension156_20261017`|
|all n>=164|`fresh_nikulin_two_budget_frontier_20261015_dimension163`|

The first rows illustrate the common scalar proof. At n59 the finite
budget is783/2, the ideal surplus is at least91/2, and the extrema
constant is4350. At n60 the finite budget404 is below curvature405,
the ideal surplus is at least54, and the constant is4500. At n61 the
finite comparison is again equality and the ideal surplus is at least64,
with constant4650. Thus the contradictions are respectively

`0>=4350+(91/2)v_i`,
`0>=4500+v_f+54v_i`,
`0>=4650+64v_i`.

All are strictly impossible, including the all-finite case. Later interval
rows use the same Section7 comparison with their stated simultaneous
weight rule and exact complete finite and ideal domains. The endpoint
name is not an additional published assumption.

For the unbounded last row, one needs a distinct geometric argument,
not extrapolation of a finite rank search. The published narrow-facet
theorem supplies a facet whose framing walls have pairwise product
greater than-14 in the norm-two normalization. If it is compact, the
independent distance17 proof bounds n by143. If it contains a cusp,
write the factor ranks with total N=n-1 and use two outgoing actual
edges from that cusp. Their terminating framing walls have positive
cusp heights with common ratio M>=1 across all factors. Eliminating the
cusp span gives the exact scale-invariant energy identity

`-(x,y)=(1/2) sum_a sigma_a tau_a E_a - q - 1/q`.

The retained elliptic edge supports give separate endpoint Schur budgets
at most2. The all-rank classical inverse formulas and the finite
exceptional/continuous-support certificates bound the total target loss
by4 in the near range `M<=21/20`. In the far range the older ordered
budget bound applies. Their precise resulting lower bounds are

`-(x,y)>=M(N-37)/9` for `M<=21/20`,
`-(x,y)>=M(N-43)/9` for `M>=21/20`.

At `N=163+t`, the first is at least `14+t/9` and the second at least
`14+7t/60`. Both contradict the strict framing-wall inequality for
every integer t>=0. Actual terminating walls, ideal second endpoints,
the shared cusp-height ratio, and the compact-facet alternative are
proved in `fresh_dimension225_hostile_20261015.md`; the sharper ordered
budget is the dimension163 source in the table. The finite enumeration
and all-rank tails of that budget are a separate computational lemma,
not established by the inequality at N163 alone.

Once every interval theorem is supplied, ordinary integer case splitting
proves the tail n>=59. Combining it with the ambient58 contradiction
then proves n<=57. No face-inheritance theorem is used.

## 9. What remains to make the six-paper project genuinely self-contained

The following are exposition/proof-boundary obligations, not newly found
counterexamples to an internal theorem.

1. **Define the full current rule.** The final D13/B13 changes are explicit,
   but their predecessor is still introduced by artifact name. A reader
   needs the complete finite/affine marked-orbit tables, all-rank formulas,
   cross-component roles and same-corner guards. The same object must feed
   the local geometric proof and the counting algorithm.
2. **State and prove the universal routing theorem.** Native selected maps
   do not cover short, low-core, disjoint, protected/adaptive, actual-ideal,
   and later replacement routes by themselves. The written domain proofs
   and arbitrary-white arguments must be included or re-presented as
   explicit mathematical lemmas with a well-founded retention order.
   Counts801/862/33/76 are outputs of finite domains, not that theorem.
3. **Give each finite computation its mathematical certificate statement.**
   The companion fragment handles definition, finite termination,
   correctness and complete leaf coverage for finite forests and ideal
   products. In particular, the final binder does not reconstruct the
   historical3298/748 support cover or1549 residual tables. Their own
   statements and proofs remain necessary.
4. **Supply the tail as an actual internal lemma.** The table removes the
   circularity risk and exposes every interval. To call the57 paper
   standalone, the row rules/local theorems and complete budget results
   must appear in an appendix or a mathematically specified shared atlas.
   Merely citing “completed58” still leaves a separate theorem dependence.
   The unbounded n>=164 row also requires the ordered continuous-support
   budget proof; no finite partition checker replaces it.
5. **Keep the external boundary narrow and exact.** Finite/affine
   classification, the Lanner/quasi-Lanner order ceilings, finite-volume
   local structure, the index/incidence identities, and the narrow-facet
   theorem are published inputs. The factor multiplicity, arbitrary-white
   seals, chosen-output retention, contextual finite budgets, positive
   concentration and every new dimension conclusion are internal.

The generated atlas is74616 TeX lines. It preserves evidence, but that
length and its source hashes are not a substitute for this mathematical
spine and the missing complete definitions. Conversely, no new computation
or proof-assistant claim is needed merely to write the already supplied
geometry coherently. This staging text deliberately leaves the unresolved
exposition blocks visible rather than declaring a full standalone audit
finished on the basis of endpoint replay.

## Source identity of the main geometric reads

|Source under `results/`|SHA256|
|---|---|
|`fresh_cusp_factor_payment_framework_20261017.md`|`c65cfc5a4409a8cc83fcf0a083a6c1a7a29bb7bc7d649c11435773877e70ccab`|
|`fresh_cusp_link_shape_curvature_20261017.md`|`0946759e18de93d13ae4729b95026b6cb92c59e3c48c0f2ef33c61032b8ebd3e`|
|`fresh_high_rank_quasilanner_census_20261017.md`|`deee1eee4ab3e01e508f6d23c995d3d66557d0724b80d533fabdcf6fca065b56`|
|`fresh_distance8_low_core_20261017.md`|`634e921bd21f0d5c39c1da343f21389368a522a9e49489b493d6226e11bcf8b0`|
|`fresh_distance8_high_overlap_20261017.md`|`84e0c3cd3dbd16545038ccdf0161e4561eaeecd241532e877d38c294643c41d1`|
|`fresh_finite8_ideal5_combined_payments_20261017.md`|`bb75329fd65933fe5a784fe21cecc64b5464523db6ffa8693b94c88392fd19d1`|
|`fresh_finite7_ideal5_triangle_payments_20261017.md`|`00d3f334ffefc3b926edfdbe7b319a7f6792d83e23fe94511ede7dc60b71bf75`|
|`fresh_d13_three_rank15_context_transfer_20261018.md`|`91e93c3f7ff5eb4e64b435eea2892bf20af08711aed131110f40e898a9322014`|
|`fresh_nikulin_two_budget_frontier_20261015_dimension163.md`|`a0b39b448d543ade48327fd5ab78e8c698834f008ced600aaa5f5f80071643f2`|
|`fresh_dimension225_hostile_20261015.md`|`45b5115b8157a342015abb3338dcc609fa48f85cfd1d2b4b71161bfd01426d7d`|
