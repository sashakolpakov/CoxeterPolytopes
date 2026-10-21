# Independent review of the unrestricted dimension156 proof

**Mathematical and source review PASS: no gap identified.** The
[primary proof](fresh_distance14_unrestricted_dimension156_20261017.md)
excludes dimensions157 through163 for every finite-sided finite-volume
real hyperbolic Coxeter polytope. Together with the inherited unrestricted163
theorem, this proves dimension at most156. No arithmetic, compactness,
simple-cusp, minimum-rank, factor-count or facet-count hypothesis is used.

This review independently checks the new physical Gram obstructions,
the graph-radius cases, the convex-order reduction and all924 reduced
cusp profiles using a [separate Fraction-only checker](../scripts/fresh_distance14_unrestricted_dimension156_20261017_review.py).
Its [certificate](fresh_distance14_unrestricted_dimension156_20261017_review.json)
binds the final producer source, data and report. It imports no producer
arithmetic and does not use SymPy. The complete high-rank minimal-diagram
classification, actual cusp-fiber interpretation and unrestricted shape
theorem remain explicitly inherited proof inputs; their large censuses
are not repeated here.

## Distance14 at actual polygon corners

The actual-face incidence argument extends to ideal corners without
assuming that the ambient cusp is simple. An edge interior is hyperbolic,
so W plus one side is positive definite. Opposite closed sides of a
genuine quadrilateral have no common finite or ideal point. The acute
projection and positive affine-kernel arguments therefore exclude PSD
for either opposite block. Each inclusion-minimal indefinite core in
that block contains both black roots, and both black deletions are PD.
This is exactly the marked domain of the inherited catalogue.

The reviewer separately recomputes384 marked graph-distance arrays from
the catalogue edges, finding the same four radius-seven cases and their
unique farthest nodes. If no white path of length at most14 joins the
two opposite black sets, the minimal cores are disjoint. A shared white
node would give a path of length at most14. Two disjoint orthogonal
indefinite cores would already contradict ambient negative index one.
Every cross bond must join radius-seven nodes, since either smaller
distance yields a path of length at most14. Thus the only possible cross
bond is between the two unique farthest nodes; all other cross products
are zero. Choosing a shortest marked path removes any internal black
node, so a surviving path lies in an actual adjacent corner diagram.

For L23 the reviewer uses a rational congruent Gram: scale node8 by
sqrt2, changing its squared norm to4 and its last bond to−2. The
inverse diagonal in the original norm2 normalization is twice the
rational inverse diagonal. Deleting black0 or1 gives respectively9/2
and4. The exterior farthest node is white and would extend a finite
edge-support subdiagram with a nonright product of magnitude at least1.
Its Schur residual is negative. This excludes L23 without a field or
label restriction and without relying on a common-center choice for
the hypothetical new node.

For two L25 cores, their A3 tails7--8--9 consist entirely of white
nodes. The union of the tails is therefore PD. Eliminating the two
inner nodes on each side leaves diagonal4/3 at the far endpoints. A
bridge of magnitude at least sqrt2 is impossible; any nonright Coxeter
bridge must consequently have label3. Dotted and parallel bridges are
covered by the same inequality.

The literal joined20-root Gram is reconstructed separately. Removing
the two fork centers leaves the positive18-root block
A12+2A2+2A1. Rational elimination gives the Schur complement

    [[−7/78, −1/13], [−1/13, −7/78]],

with two negative eigenvalues−1/6 and−1/78. Thus this last possible
distance15 configuration also cannot occur. The triangle argument
continues to give two paths of length at most8; no new claim about
maximal-simplex incidence at ideal vertices is involved.

## Costs at finite and arbitrary ideal vertices

The inherited all-rank finite and affine graph formulas bound the
distance-at-most14 pair count by14 times component order. The small
exceptional diagrams have at most4 times component order pairs in
total, which is enough here. Finite forests therefore cost at most14n.

At an arbitrary cusp, a paid marked pair lies in one affine factor.
For a factor of rank r and order r+1, its number of paid pairs is at
most min(14(r+1),r(r+1)/2). The actual factor-fiber multiplicity divides
this by r+1 when normalized by the link vertex count V. This gives
exactly sum min(14,r/2). It does not substitute total incident wall
count n for a nonsimple cusp, or charge pairs from different affine
factors as though they were adjacent.

## Concentration and complete partition coverage

The transfer(a,b) to(a−1,b+1), 1<=a<=b, preserves the total rank and
mean. On the common integer support0,...,a+b, the new-minus-old mass is
positive for j<a or j>b and negative for a<=j<=b. On the central
interval the difference is−1/((b+1)(b+2)); on either tail it is the
old rising or falling numerator times
(b−a+1)/[(a+1)(b+1)a(b+2)]. These formulas are checked independently
for all6,561 positive pairs with a+b<=162, which includes every
transfer needed in the target dimensions.

The same formulas prove the all-rank convex-order statement. Center
the support and take the chord of an arbitrary convex test function
through the points ±(b−a+1)/2. The function lies below that chord on
the negative-mass central interval and above its extension on the
positive-mass tails. Mass and first moment cancel the chord, so the
expectation change is nonnegative. Convolution with the other
independent uniforms preserves this comparison. Hence normalized
curvature cannot increase under the transfer.

The cost-preserving endpoints matter. For low factors the receiving
rank must start below28: use1<=a<=b<=27 and stop at0 or28. For high
factors both start above28 and the smaller stops at28. Thus the cost
sum min(14,r/2) is unchanged at every step. Strict increase of the
integer sum of squared ranks proves termination. Only repeated28s,
at most one rank1..27, and at most one rank>=29 can remain. Zero
ranks contribute a constant uniform and zero cost and may be dropped.
This proves domination for every initial rank partition, regardless
of its number of factors.

The independent enumeration chooses the optional large rank first;
Euclidean division of the remaining total by28 uniquely supplies
the repeated blocks and small remainder. Counts for link dimensions
156 through162 are129,130,131,132,133,134,135, totaling924. Integer
sliding-window convolution is checked coefficient by coefficient
against an independent bounded-composition inclusion-exclusion
formula. All curvatures, costs and margins agree with the final
producer certificate. The smallest margin is

    1476616/20511149 > 0,

at ranks(16,28,28,28,28,28). The finite budgets beta_n>14n are separately
checked throughout n=157,...,163. At n=156 the displayed local profile
(15,28,28,28,28,28) instead has margin−4172075/82044596; this correctly
limits this payment comparison without asserting an ambient example.

## Composition and validation status

All actual two-faces receive the required5−k payment, and the complete
partition reduction bounds every ideal-vertex contribution. The
unrestricted shape inequality has finite coefficient beta_n, ideal
curvature kappa and the additional positive5(n−1) term. It therefore
contradicts the distance14 upper payment throughout dimensions157–163.
The already complete unrestricted163 proof covers every dimension
above163. There is no finite upper-rank sampling gap or conditional
cusp hypothesis in this composition.

Independent generation52415 and fresh saved replay87174 both exited0
with the identical saved certificate. Source syntax, direct manifests,
the producer's direct source bindings and local report links pass.
This review's source, data and report are final and frozen.
Source-integrity checks are separate from the mathematical proof reads.
The completed result is an internally checked computer-assisted proof;
external refereeing, publication priority and proof-assistant
formalization are not claimed.
