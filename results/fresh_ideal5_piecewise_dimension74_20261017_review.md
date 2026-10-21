# Independent review of the dimension74 comparison

The independent exact reconstruction and full written composition pass.
The [primary proof](fresh_ideal5_piecewise_dimension74_20261017.md) establishes
dimension at most74 for every finite-sided, finite-volume real hyperbolic
Coxeter polytope, without arithmetic, compactness or cusp-shape restrictions.
The completed dimension75 theorem supplies all higher dimensions; the new
comparison excludes dimension75.

The [reviewer](../scripts/fresh_ideal5_piecewise_dimension74_20261017_review.py)
imports only standard-library arithmetic. Its
[certificate](fresh_ideal5_piecewise_dimension74_20261017_review.json)
independently enumerates999 reduced profiles, compares every primary row and
reconstructs the seven new literal rational cores. The complete inherited
overlap and all-label supporting-white censuses remain explicit frozen inputs.

## Independent scalar reconstruction

Only the affine-rank14 cost changes, from100/15 to97/15. The envelope is

    r/2 for1<=r<=12;
    87/14 at13;97/15 at14;111/16 at15;126/17 at16;
    163/23 forr>=17.

The reviewer chooses optional low, middle and high ranks first, then
recursively enumerates all multiplicities of12,14,16,17 at the remaining
total rank. It imports no producer profile list. Sliding-window convolution
forms H(t)=product_i(1+t+...+t^r_i), with V=H(1) and
G=sum_j h_j min(j,d-j). A separate inclusion-exclusion calculation expands
product_i(1-t^(r_i+1))=sum_l a_l t^l and obtains

    G=floor(d/2)*V
      -2*sum_l a_l*binom(floor(d/2)-l+s,s+1).

The two values agree for every profile. Polynomial symmetry, coefficient sum
and the absolute-deviation identity

    (5G-dV)/(2V)=3d/4-(5/2)*E|sum U_i-d/2|

are also checked. All999 rows at total rank74 match the primary certificate,
including every coefficient, curvature, cost and margin. The minimum ideal
margin is63/74698 at(10,12,12,12,12,16). The finite slack is95/2 and the
positive global term is370. At ambient74, the boundary profile(1,12,12,12,12,12,12)
has margin-168143/742586; it is a limitation of this comparison, not an
asserted polytope.

The inherited convex-order move(a,b) to(a-1,b+1) cannot increase normalized
curvature. Low ranks concentrate to12 or zero at unchanged linear cost.
Duplicate13 can become12 and14 with cost increase4/105; duplicate15 can
become14 and16 with increase7/2040. Both increments remain strictly positive
under the changed envelope. Two ranks at least18 concentrate until one
reaches17 at unchanged constant cost. Every move increases the bounded
integer sum of squared ranks, so the process terminates even if middle
ranks are recreated. The terminal domain is exactly the one independently
enumerated. The unchanged15 all-rank affine cost branches and their maxima
are independently recomputed and matched to the frozen branch certificate.

## Three removed pairs and seven complete finite outcomes

The full source and mathematical proof of the
[local affine-B reduction](fresh_order15_affine_b_pair_reduction_20261017.md)
were read without an identified gap. It removes only(0,11),(11,13),(11,14)
from the order15 affine B14 target. Independent breadth-first search gives
distances11,12,12. The leaf swap preserves their union. Their complete saved
marked-image support consists exactly of1284,1623,1641,1646,1666,1933,2326.
All16 elementary-triangle inclusions remain valid, including the prior
dimension75 singleton deletion in order17 affine D.

The reviewer reconstructs every local Gram from the literal L25 tree and
the saved identifications, rather than importing a producer matrix helper.
Rational congruence elimination and Gauss-Jordan inversion independently
check the signatures, white and displayed-corner positive definiteness,
marked components, contact projections and endpoint extensions. These
small checks verify the seven listed outcomes; they do not reconstruct the
complete all-label contact classification that supplies their candidate lists.

In six contexts the core has signature(16,1,0), and the preferred component
is D14. Every allowed nonzero direct contact is at node3 or13, outside that
component. No extra-white path can reach it without a first direct contact.
The complement in the actual corner consists of supporting whites and is
positive definite. Thus the entire corner is finite and D14 is whole,
so the already paid whole-component rule supplies a unit.

In context1641 the core has signature(15,1,0), and the preferred component
is A12. Only its two endpoints can receive relevant direct contacts. Each
endpoint-contact root has projection square2. Its residual lies in the
positive-definite perpendicular to the Lorentz core and has square zero,
so the residual is exactly zero. It is therefore orthogonal to every
core-orthogonal white. Against another allowed core-contact white, its
actual product is a verified projection product: positive, or at most-2.
The former violates acuteness for distinct actual walls; the latter violates
strict positive definiteness of two supporting whites. Identical contact
and zero residual fix the same normal, rather than another occurrence.

If no endpoint-contact root exists, A12 remains whole. If one exists,
it gives a whole A13. Both exact extension matrices are checked. The
remaining corner components are again positive-definite white blocks.
This establishes actual finite corners, not merely finite subdiagrams
inside an ideal corner. No additional finite coupon is required.

The argument uses the complete frozen no-extra-cross and all-label support
domain for each literal core. The removal is not justified by an arbitrary
abstract affine embedding. Basic, triangle and short8-first payments remain
unchanged because the three omitted pairs have distance greater than9.

## Global composition and validation

The complete finite8/ideal5 theorem and the dimension75 refinement retain
all other local payments and the whole A4/D13 finite coupons. Since the new
seven-case replacements use existing finite units, the finite budget remains
8n through80. At order15 the other type costs75,95,96 are below97, so the new
rank14 envelope is valid. Actual nonsimple-cusp factor costs retain the
link-vertex multiplicity divided by factor order; no simple-cusp hypothesis
is introduced.

Every ideal profile at dimension75 has positive margin, the finite margin
is95/2 and the inherited global shape inequality retains370. This contradicts
the summed local face payments in dimension75. The completed dimension75
theorem excludes dimensions at least76, proving74. The compact arithmetic21
and arbitrary compact29 ceilings are unchanged. No new construction or
dimension73 result is asserted.

Independent bound generation and separate fresh saved replay both terminated
exit0/PASS. Every999-row scalar field and both comparison interfaces match
the frozen primary. All25 direct bindings, syntax and4 local links pass.
Root completed a full source and mathematical read of the independent
scalar algorithms, all-rank branches, seven literal-core reconstructions,
local interfaces, bindings and this entire report, finding no gap. That
peer read did not duplicate the999-row arithmetic, seven rational checks
or inherited censuses. Source, data and report are final/frozen; all
frozen predecessors are unchanged.

This is an internal computer-assisted proof review. External refereeing,
publication priority and Lean or other proof-assistant formalization are
not claimed. The inherited2,614-context overlap and1,153-pattern embedding
censuses are not rerun by this package.
