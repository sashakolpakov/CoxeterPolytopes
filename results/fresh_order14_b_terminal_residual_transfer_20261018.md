# Removing the last two extra order14 affine B payments

In a whole affine B13 factor of order14, delete the pairs(9,12) and(9,13).
The total falls from84 to82, exactly its basic ideal5/feature5 pair count.
No finite weight, finite coupon or other ideal payment changes. The target
has spine0 through11, doubled bond10–11, and leaves12/13 at0. The two
deleted pairs form one automorphism orbit of distance10. Every short8
and elementary triangle pair remains a unit.

Together with the separately proved
[C14 singleton deletion](fresh_order14_c_inner_endpoint_projection_transfer_20261018.md),
the four order14 totals are A70/B82/C83/D82, with maximum83. This is a
local payment improvement. The ambient62 finite obstruction is separate;
no new ambient dimension bound follows here.

The [source](../scripts/fresh_order14_b_terminal_residual_transfer_20261018.py)
and [certificate](fresh_order14_b_terminal_residual_transfer_20261018.json)
start from the [previous B84 theorem](fresh_order14_affine_b_fork_projection_transfer_20261017.md)
and the [completed62 rule](fresh_sparse_coupon_order11_13_dimension62_20261018.md).
They reconstruct all600 eligible high-overlap selected patterns, including
failed searches, compare207 successful patterns and506 induced maps with
the frozen complete certificate, and check all108 changed maps in53
contexts. The affected maps are not inferred to be every quadrilateral;
the separate nonselected interfaces are treated below.

## Complete old-core contact control

Every affected old core H is nonsingular Lorentzian, simply laced and
physically normalized to square2. The inherited full contact checker
recomputes its inverse J, signature, actual selected corner, and entire
necessary contact domain for any additional supporting white. It also
checks that old corner components outside the selected component are
orthogonal to it and that the omitted blacks are absent from this corner.

The contact domain splits into two end sets with entrywise-positive
inverse blocks. Any extra white meeting both end sets gives an actual
opposite-black path of length at most8. The present simultaneous payment
route for that case is supplied below. In the remaining branch, every
extra white's contacts with H lie in one end. If a is its nonnegative
physical contact vector, then

    projection_H(e)^2 =a^T J a <= e^2=2.

The positive complement of the Lorentz span proves this inequality.
Positivity of the end block gives a^T J a>=q²J_ii for every specified
contact q at i, regardless of arbitrary other same-end contacts. In
particular a label4 contact has q²=2 in this normalization.

In92 affected maps, the actual target forces a missing root adjacent to
the displayed component with q²J_ii>2. This root must be outside H:
omitted blacks are not in the corner, and other old corner blocks are
orthogonal to the displayed component. These are immediate projection
obstructions even when the old classification calls the context a
relevant-component case. They do not require a whole completion or a
finite replacement coupon.

## The sixteen two-root obstructions

The remaining16 maps occur in the eight contexts

    630,638,1648,1668,1849,1935,2242,2328.

Each displayed image contains every target node except10 and11. The
actual whole B13 factor therefore needs a root e at target10, simply
joined to the old target9 at core node i, followed by a root f at target11,
joined to e by label4. Both actual roots have square2. The source checks
these exact incidences for every map.

The first root e is outside H by the immediate-neighbor argument above.
For these contexts, the full end blocks also satisfy

    J_aa+J_bb+2J_ab>2

for every distinct pair a,b in one end. Native nonzero contact magnitudes
are at least1; block positivity therefore excludes two contacts on the
same end. Contacts on both ends already give the short8 route. Thus e
has exactly the one old-core contact at i. Moreover J_ii is3/2 or5/4,
so a label at least4 there would have projection square at least2J_ii>2.
Its contact is exactly the simple magnitude1 fixed by the target. Its
residual in the positive complement has square

    rho=2−J_ii =1/2 or3/4.

The terminal f is genuinely new as well. It cannot be an omitted black.
If it were an old root in another displayed corner component, e would
have an additional old-core contact to f through the required doubled
bond, contradicting the just-proved uniqueness. Old roots of the selected
component already have their fixed distinct images and cannot supply a
missing target node. Thus the complete extra-white contact domain applies
to f; it is not restricted to a guessed subset of spectator roots.

If f contacts both old end sets, its own contacts give the earlier short8
alternative. Otherwise its contact vector b lies on one positive inverse
block, giving b^T J b>=0 and therefore

    ||f_perp||²=2−b^T J b <=2.

This bound includes equality when f is orthogonal to H. We do not need
to assume that f has one old contact, or that it has no old contacts.
For every possible contact node j, the source checks the sign of J_ij.
If J_ij<0, the actual path from an opposite black through i,e,f,j to the
other black has length at most8; in all eight contexts its saved bound
is3. This path uses the actual doubled e–f edge and the actual nonzero
f–j contact, independently of their magnitudes. That case pays by the
short-path route. After excluding it, every actual term in

    <projection_H(e),projection_H(f)>=sum_j b_j J_ij

is nonnegative. Zero coefficients and an empty support are permitted.
Since the physical e–f bond is−sqrt(2), their residual product satisfies

    <e_perp,f_perp> =−sqrt(2)−sum_j b_j J_ij <=−sqrt(2).

Its square is at least2. Cauchy in the positive complement instead gives

    <e_perp,f_perp>² <=rho*||f_perp||² <=2rho <2.

This is a strict contradiction, with gap1 when rho=1/2 and gap1/2 when
rho=3/4. Arbitrary additional contacts of f only retain this necessary
contradiction or supply a short path. No unspecified Gram entry has been
set to zero. No residual-rank bound, optional exterior-root enumeration,
or assumption excluding indirect spectator bridges is required.

Thus all108 changed selected occurrences are impossible unless an
independently retained short-path payment already resolves the face.
There is no surviving affine completion needing a new typed coupon.

## Complete domain and simultaneous present weights

The low-core theorem has whole marked-component bound at most12, so
cannot require either changed order14 pair. The all-label high-disjoint
branch is excluded. Protected/basic ideal choices have distance at most9
and stay paid; every finite whole-component, adaptive or protected-leaf
payment is unchanged. The remaining high-overlap selected cases are
precisely the complete600-pattern domain reconstructed here. Both removed
pairs are outside the basic set; all82 basic pairs remain units.

The universal tree and older D13/D14 immediate affine outputs have orders
at most13. The cycle12, D11 and order12 B/C successors use finite payments
or smaller/shorter ideal factors. Later order15/16/17 D transfers use their
own D-factor roles and typed finite coupons. The newer C11/B12/D13 and
coupled order11 theorems do not output these order14 B distance10 pairs.
Neither the separate C14 singleton proof nor B16 first-entry proof uses
them. Hence the selected domain has no unaccounted later long-pair use.

To invoke the short8 alternative under the present reduced rule, use
actual-corner priority. First apply independent C11-zero, B12-zero and
D13-half or D13-zero theorems whenever those changed actual corners occur.
Their outputs are short units, retained order11/12 center units, whole
D12 leaf9, or the stated whole finite coupons. Next apply the independent
cycle11 and coupled B/C/D11 half theorems to any affected actual corners.
They give distinct retained halves or finite protected units without
calling for a long order14 B pair. These routes do not recursively spend
an ideal weight that is being removed here.

Otherwise the explicit current short8 route from the
[completed63 composition](fresh_order11_order17_dimension63_20261018.md)
remains: unaffected small-tree units; the established cycle6/7/8 and
radius-five tree replacements; and, whenever a finite distance8 output
is deep, the independent depth11 corridor supplies its protected short7
replacement. The sparse finite deletions are at distance at least10.
The order17 zero refinement only removes its own distance10 orbit and
uses this same independent short8 route or a distinct whole E6/D11 unit.
Thus neither proof requires the other's changed orbit. The short7 route
uses retained finite paths or protected roles, rather than another deep8
zero, so it is not circular.

After these priorities, an old selected payment at either new B14 pair
is excluded by the complete local obstructions above. Every other old
selected payment is retained. On a face with no actual changed B14 corner,
weights agree pointwise with the prior rule. The separate C14 deletion
has exactly the same independent short-path scope and removes a different
factor type. Its affected occurrences are likewise excluded by their
own first-entry proof. The two deletions therefore compose without new
coupons or any assumption about which old witness would be selected.

Triangle payment2 is unchanged: all elementary retained pairs are checked,
and the removed distance10 orbit is absent from the finite7/ideal5,
protected6/7 and whole small-factor triangle routes. Quadrilateral payment1
follows from the complete cases above. The current finite forest rule
and its maximum are unchanged. Cusp factor-fiber multiplicity remains
V/14, and deleting two unit pairs changes the factor total by exactly2.
The separate C14 theorem gives the stated all-family maximum83; a global
scalar/budget successor must bind both local theorems explicitly.

## Verification scope

Generation486d85 and separate fresh replay3c04d0 both terminated exit0/PASS.
Audit080f4c passes35 direct bindings, six local links and Python syntax.
Root and a separate peer completed full source/report mathematical reads,
including the inherited complete contact checker at their recorded scopes.
They passed the92 direct all-contact bounds and16 residual cases: e's
unique old contact, f's genuine newness, its complete one-end projection
bound, every negative cross coefficient's actual shortcut, and strict
Cauchy including the core-orthogonal terminal case. Both also checked
the exhaustive nonselected/current short8/triangle interfaces and absence
of new coupons. Neither repeated maps, inverses, Cauchy arithmetic or the
upstream census. Source, data and report are final/frozen.

The source re-enumerates all600 selected patterns and506
successful maps, checks53 full necessary contact domains, and verifies all
92 first-entry inequalities and16 two-root Cauchy contradictions. It binds
the inherited all-label overlap reduction, basic/triangle rules and current
simultaneous short-path composition. It does not repeat the entire original
core census or prove those inherited classifications anew. No finite
optimization, external refereeing, formal verification or global dimension
bound is asserted.

    python3 scripts/fresh_order14_b_terminal_residual_transfer_20261018.py
    python3 scripts/fresh_order14_b_terminal_residual_transfer_20261018.py --verify
