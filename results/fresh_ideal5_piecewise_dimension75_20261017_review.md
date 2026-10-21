# Independent review of the dimension75 comparison

The independent exact comparison and complete written composition pass. The
[primary argument](fresh_ideal5_piecewise_dimension75_20261017.md) proves that
every finite-sided, finite-volume real hyperbolic Coxeter polytope has dimension
at most75, without an arithmetic, compactness or cusp-shape restriction.
This review uses the completed dimension77 theorem for the higher-dimensional
tail. It finds no open mathematical dependency in the new composition.

The [checker](../scripts/fresh_ideal5_piecewise_dimension75_20261017_review.py)
imports only Python's standard library. Its
[certificate](fresh_ideal5_piecewise_dimension75_20261017_review.json)
independently reconstructs all2,164 reduced scalar profiles, compares every
primary row and checks the new local interfaces. It does not rerun the
inherited2,614-context overlap or1,153-pattern embedding searches.

## Independent arithmetic and complete reduction

The normalized ideal-factor budget at affine rank r is

    r/2 for1<=r<=12;
    87/14 at13;20/3 at14;111/16 at15;126/17 at16;
    163/23 forr>=17.

The checker chooses the optional low rank1,...,11, optional13 and15 and optional
rank at least18 first, then recursively enumerates all multiplicities of
12,14,16,17 at the remaining total rank. This differs from the primary's
anchor-first loops and imports no producer list.

Sliding-window convolution forms H(t)=product_i(1+t+...+t^r_i). For
d=sum r_i, V=H(1) and G=sum_j h_j min(j,d-j), a separate inclusion-exclusion
calculation expands product_i(1-t^(r_i+1))=sum_l a_l t^l and obtains

    G=floor(d/2)*V
      -2*sum_l a_l*binom(floor(d/2)-l+s,s+1).

Both exact values agree for every profile. The checker also verifies polynomial
symmetry, H(1)=V and the absolute-deviation identity

    (5G-dV)/(2V) = 3d/4-(5/2)*E|sum U_i-d/2|.

Subtracting the new ideal budget gives the following exact minima.

|Ambient dimension|Profiles|Minimum ideal margin|Finite slack beta_n-8n|
|---:|---:|---:|---:|
|76|1053|11615/448188|57|
|77|1111|41558/485537|68|

The first minimum occurs at(11,12,12,12,12,16), the second at(12,12,12,12,12,16).
The ambient75 boundary profile(12,12,12,12,12,14) has margin-16873/171366.
It is a formal obstruction to this comparison, not a realized polytope or
a dimension74 theorem.

The inherited convex-order theorem permits transfer from(a,b) to(a-1,b+1)
for1<=a<=b without increasing normalized curvature. Concentrating two ranks
at most11 preserves their linear cost until one vanishes or the other reaches12.
Replacing duplicate13 by12 and14 raises cost by5/21. Replacing duplicate15
by14 and16 raises it by83/408. Concentrating two ranks at least18 until one
reaches17 preserves their constant cost. Each move strictly increases the
bounded integer sum of squared ranks, so even moves that recreate a middle
rank cannot cycle. The terminal domain is exactly the independently enumerated
one. Positivity there therefore covers every partition of the same total rank.

The checker additionally reconstructs the15 all-rank B/C/D affine cost branches
and their exact ratio maxima160/23,75/11,163/23 from spine-window formulas.
Each branch has cost/m=a+b/m, so its endpoints and limiting value5 exhaust its
extrema. This is an all-rank calculation, not extrapolation from sampled ranks.
All16 small-factor cost and elementary-triangle inclusion interfaces are checked.

## The new local transfer

The full source and written proof of the
[D17 singleton transfer](fresh_d17_singleton_pair_coupon_20261017.md) were read.
Here D17 counts nodes: the standard affine type is D16. Removing only the
automorphism-fixed pair(1,11) reduces its cost127 to126. Independent breadth-first
search gives its distance10. The checker validates all eight saved automorphisms,
the invariance of the new paid set and retention of all elementary triangle
pairs; it does not claim an independent automorphism enumeration.

The complete saved induced-map data contain precisely one affected context,
2518, the two L25 copies marked{3,4} and identified at node9. The reviewer
reconstructs its19-root Gram from the literal edge list using rational arithmetic.
Independent congruence elimination gives signature(17,1,1), a positive-definite
15-root white diagram, four positive-definite16-root edge diagrams and the
four stated corner signatures. Direct multiplication verifies the saved affine
positive kernel and full-core radical. The two finite coupon corners have
whole A4 and D13 components with the correct endpoint and long-arm-neighbor
roles. These small checks are independent of the local producer's helpers.

The geometric handling of additional whites is essential. At the actual ideal
corner containing the17-root affine diagram, its positive kernel and nonpositive
contacts force every extra supporting white to be orthogonal to every affine
root. No simplicity assumption is needed. The full-core Gram radical represents
a vector orthogonal to a timelike core vector; its zero square in that positive
orthogonal complement makes it an actual zero relation. Pairing it with an
extra white equates that white's contacts with the two omitted roots.

If those equal contacts are nonzero, the opposite corner has a paid path of
length2. Otherwise all additional whites are orthogonal to the entire core.
Their Gram is positive definite as part of the actual face support, so the
two displayed finite corners stay finite and their A4/D13 components stay
whole. The two new half-payments then supply the required unit. This proves
the transfer for arbitrary additional whites and nonsimple ideal vertices.

The literal whole-core hypothesis is supplied by the completed overlap theorem,
which excludes extra cross-core bonds. An arbitrary affine embedding alone
would not justify sealing the finite components. The proof and its saved-map
interface retain this distinction. The distance10 omission changes no short8
reduction, elementary triangle payment or immediate small-factor payment.

## Finite budget and global composition

The full source and proof of the
[new finite budget](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.md)
also pass the narrow review. A whole A4/D13 component pair adds one unit:
two A4 endpoints, one specified D13 node and half-weight each. If a is total
small-component rank and b is the number of D13 components, this is at most ab/4.
Each D13 has unit-pair cost78 and slack26 against8r, supplying18 more than
the previously reserved8 and consuming two additional ranks beyond the
generic large-component bound11. Thus

    E_new = E_old+b(a/4-18).

For b>0, the rank bound in dimensions at most80 gives a<=67, so the added
term is at most-5b/4. For b=0 it vanishes. The old finite8 budget therefore
continues to cover every old coupon and the new one. This argument was read
in full; the3,465 endpoint enumeration was not independently rerun here.

The completed [finite8/ideal5 theorem](fresh_finite8_ideal5_combined_payments_20261017.md)
supplies all remaining triangle and quadrilateral payments, using the full
2,614-context classification rather than the old no-short9 domain. Only the
one saved D17 pair changes. Its new transfer, all16 elementary triangle
inclusions and existing removed-pair transfer interfaces preserve the required
face payments. Actual nonsimple-cusp multiplicities multiply a factor cost by
the link-vertex count divided by its order; no simple-cusp restriction enters.

The finite coefficient exceeds8n by57 and68 in dimensions76 and77. Every ideal
profile has the positive margin above, and the inherited shape inequality
retains positive global terms375 and380. This yields the contradiction in
both dimensions. The completed dimension77 proof excludes the higher tail.
The unrestricted conclusion also gives an effective compact-facet ambient
ceiling75; the separate compact arithmetic21 and arbitrary compact29 results
are unchanged.

## Validation scope

Independent generation and a separate fresh saved replay both terminated
exit0/PASS. All28 direct source/data/report bindings, syntax and all6 local
report links pass. The checker independently recomputes every scalar row,
the boundary, all-rank cost branches and the stated small local/interface
checks. The full primary and local mathematical composition was read without
an identified gap. A separate complete source and mathematical peer read
passed the alternate enumeration, sliding convolution, cumulative identity,
absolute-deviation check, concentration argument, independent inertia and
literal core, saved-map and automorphism interfaces, all-rank branches and
final bindings. It also read this full report. That peer read did not
duplicate the2,164-row arithmetic or inherited censuses.

The inherited diagram classifications, full overlap and map searches retain
their frozen validation records. They are not independently reconstructed by
this package. All reviews are internal; external refereeing, publication
priority and proof-assistant formalization are not claimed.

Source, data and report are final/frozen. No frozen predecessor was modified.
