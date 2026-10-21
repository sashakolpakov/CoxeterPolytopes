# Ideal cost seven and a half with finite cross-component half-payments

This local argument replaces whole order17 ideal payments selectively and
preserves the total finite-vertex cost at most9n for n<=108. Every triangle
receives at least2 and every quadrilateral at least1. Every irreducible
ideal factor of order m costs at most min(15m/2,m(m-1)/2).
An ambient dimension conclusion requires the separate curvature comparison.

The [source](../scripts/fresh_ideal7_finite_cross_payments_20261017.py) and
[certificate](fresh_ideal7_finite_cross_payments_20261017.json) inherit the
complete [distance9/order17 classification](fresh_selective_distance9_dimension108_20261017.md)
and its [independent review](fresh_selective_distance9_dimension108_20261017_review.md).
They do not presume that the known H17 chamber is the only ambient chamber
containing one of its face stars.

## The new weights

At finite vertices retain weight1 for connected marked pairs at distance
at most9 or in a component of order at most17. Add weight1/2 for a pair
in different irreducible components in either of these cases:

* one component is A1 and the other is D16;
* one component is A2 and the other is D14.

In either case the mark in D must be its unique degree3 fork center or
one of the two short fork leaves. Either A2 node is allowed. All these
definitions refer to whole irreducible components of the actual finite
vertex diagram. Other disconnected pairs retain weight0.

At ideal vertices use the actual varying affine factor, as justified by
the [factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md).
For order at most16 pay every pair. For order at least18 pay distance<=7,
and distances8 or9 if either mark lies within7 steps of a feature. A
feature has degree>=3 or meets a bond of label>=4.

At order17 the rules are:

| Affine type | Paid pairs | Count |
|---|---|---:|
|A16 cycle|distance<=7|119|
|B16|distance<=9, plus a fork leaf paired with the spine node nine steps from the fork|111|
|C16|distance<=9|108|
|D16|all pairs except distance>=13 and the pair of fork centers|127|

The B exception consists of two pairs. The D deletion consists of nine
pairs: four at distance14, four at distance13, and the two fork centers
at distance12. These rules are intrinsic under diagram automorphisms.
Every ideal weight is0 or1.

## Why all faces remain paid

The [mixed finite9/ideal7 theorem](fresh_ideal_distance7_mixed_replacement_20261017.md)
gives two length<=7 triangle paths, so their two unit payments survive.
It also handles an ideal quadrilateral corner at distance8 or9 in any
factor of order>=17 outside the feature neighborhood. Its replacement
has finite distance<=9 or ideal distance<=7 and hence a unit payment.

Consequently a quadrilateral of total weight<1 can have no adjacent
white-interior path of length<=9. At a finite corner such a path would
already pay1. At an ideal corner it either pays directly or invokes the
mixed replacement. Thus the complete inherited no-short-path
classification applies unchanged, including all its cross-bond exclusions.

It produces a preferred connected marked component of order at most17
after every supporting-white extension. A finite preferred corner or an
ideal factor of order<=16 still pays1. The sole issue is an ideal factor
of order17. The preferred core component must embed as an induced labelled
subdiagram of that actual factor, with its two marks retained.

Only80 of the949 inherited high-order contexts have upper bound17;
all remaining contexts and the complete low-core supplement have smaller
bounds. The new checker reconstructs each of these80 preferred components
and enumerates every injective induced labelled embedding into all four
classical affine order17 diagrams. It uses explicit backtracking: every
unused target vertex consistent with every previously mapped edge and
nonedge is tried. Necessary degree pruning does not exclude an embedding.
The domain is complete by induction on the number of mapped vertices.

Taking all preferred marked images together with the basic short/feature
rule gives costs119,111,108,136 on A,B,C,D. In D, the nine deleted pairs
can arise only in context indices383,393,400,647,662,675,813,828,841.
Every one is two L25 diagrams sharing node9, with the exact19-root Gram
already fixed by the inherited cross-bond tests. No other core, label,
embedding or indirect extension is discarded.

## Sealing the nine exceptional stars

In each of those stars the actual corner obtained by omitting black
roots3 and13 contains the entire17-node affine D16 factor F. An additional
supporting white e is orthogonal to every root of F: its products are
nonpositive and the positive affine kernel pairs to zero with e.

The19-root core is indefinite, so a Gram-radical vector is an actual
linear relation of its roots in Lorentz space. Its coefficients at3 and13
are nonzero with opposite signs. Since e is orthogonal to all other
core roots, the relation forces e to contact both3 and13 or neither.
Contacting both gives the adjacent white-interior path3-e-13 of length2,
a unit payment at that actual corner. Thus in an otherwise unpaid
quadrilateral every extra supporting white is orthogonal to the entire
19-root core. This seals all its finite components; no outside tail can
enlarge them or turn a finite corner into an ideal one.

The checker inspects all four corner decompositions in each of the nine
literal stars. Five contexts already have a finite unit payment in the
retained connected-component rule. The other four, indices383,400,813,841,
each have two finite half-payments of the stated A1/D16 or A2/D14 kind.
The pair of fork-center marks uses two A2/D14 payments; two opposite
fork-leaf marks use two A1/D16 payments; the mixed case uses one of each.
Every exceptional quadrilateral therefore receives at least1.

The separate [half-payment certificate and review](fresh_h17_half_coupon_transfer_20261017_review.md)
reconstructs the nine literal stars, their Gram signatures and kernels,
and every finite corner over the rational numbers. It also gives a
separate written proof of the support-sealing argument. It does not
replace the complete80-pattern embedding enumeration above.

This is a universal local argument about actual face stars. Its support
sealing is essential: merely observing the same finite component types
in the known H17 chamber would not prove the transfer.

## Every-rank ideal-factor cost

Whole factors of order<=16 have m(m-1)/2<=15m/2 pairs. The explicit order17
costs above are also<=15m/2. In order m>=18, direct spine counts give:

| Type | Distance<=7 pairs | Extra distance8/9 pairs |
|---|---:|---:|
|A|7m|0|
|C|7m-28|min(m-8,18)+min(m-9,18)|
|B|7m-27|4+min(m-10,17)+min(m-11,17)|
|D|7m-26|8+min(m-12,16)+min(m-13,16)|

For C, the simple path has m vertices. The feature neighborhoods at its
two ends have nine positions each. At distance k=8 or9, the starting
positions paid from the two ends form intervals whose union has size
min(m-k,18). For B the spine has m-2 vertices and end neighborhoods of
sizes8 and9. For D it has m-4 vertices and both neighborhoods have size8.
The union count is always min(number of positions,sum of end lengths).
Each fork leaf contributes one extra spine pair at each distance; this
gives4 in B and8 in D. Opposite fork leaves are farther than9 in these
orders. Their same-fork pairs were already included in the short count.

These are all-rank formulas. The checker also corroborates them by444
literal BFS counts through order128; that finite range is not the proof
of the infinite tail. Each sum of minima has three linear branches.
For C they are9m-45,8m-19,7m+8, breaking at26 and27. For B they are
9m-44,8m-17,7m+11, breaking at27 and28. For D they are
9m-43,8m-15,7m+14, breaking at28 and29. Comparison with15m/2 at the
finite interval endpoints, and the increasing slack on each final branch,
proves the desired bound for every m. The smallest D slack is1/2 at m=29.
The finite exceptional affine types all have order<=9.

## Total finite cost remains9n in the relevant dimensions

Let a,b,c,d be the numbers of whole components A1,A2,D14,D16 at a finite
vertex. Their retained costs are0,1,91,120 respectively, and the inherited
cost on every other component of rank r is at most9r. The new coupon
cost is at most(3/2)ad+3bc. Thus

    cost <=9n-9a-17b-35c-24d+(3/2)ad+3bc,
    a+2b+14c+16d <= n.

For n<=108, d<=6 makes a((3/2)d-9)<=0. If c<=5 then
b(3c-17)-35c<=0. For c=6, b<=12 gives at most12-210<0; for c=7,
b<=5 gives at most20-245<0. The remaining term-24d is nonpositive.
Hence total finite cost is at most9n. This analytic argument permits
arbitrary other finite components and all their Coxeter labels.
The source additionally enumerates all admissible a,b,c,d in dimensions
92 through108 as an exact finite corroboration.

## Verification status

Primary generation passes the complete80-pattern embeddings, all nine
sealed replacements,444 large-factor graph counts, the small affine
catalogue, nine unbounded linear branches and the finite-budget checks.
Fresh saved replay also terminates exit0 with the identical certificate.
A separate full mathematical/source read checks the embedding domain,
the reduction to the inherited no-short-path classification, sealing,
finite9n budget and every-rank ideal counts. That review independently
compares the saved order17 paid sets with the intrinsic rules above;
it does not duplicate the80-pattern embedding enumeration or444 BFS runs.
The separate rational half-payment certificate independently reconstructs
the exceptional stars and their finite payments.
The new package does not rerun the inherited complete high-core census;
it binds that completed proof and uses its explicit order and extension
bounds. External refereeing and formalization are not asserted.
