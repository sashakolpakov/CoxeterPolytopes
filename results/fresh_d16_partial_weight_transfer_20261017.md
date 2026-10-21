# Ten half-weights on the order16 affine D factor

Starting from the completed finite8/ideal5 payment theorem and its successors,
ten further pair weights in affine D15 fall from1 to1/2. Its diagram has16
nodes, and its total pair weight falls from111 to106. Every affected
quadrilateral still receives at least1. No new finite coupon is needed.
The order16 maximum consequently becomes107, attained by the unchanged B
rule; the unchanged A and C costs are80 and103. This is a local weighted
transfer, with no ambient dimension assertion.

The [source](../scripts/fresh_d16_partial_weight_transfer_20261017.py) and
[certificate](fresh_d16_partial_weight_transfer_20261017.json) cover every
affected map from the frozen complete
[small-factor census](fresh_distance8_ideal5_small_factor_20261017.md).
They do not assume that a displayed finite subdiagram is the whole affine
factor. The proof below controls every possible completion to the actual
order16 factor, including contacts with roots outside the displayed component.

## The intrinsic weighted rule and its complete support

The canonical affine-D15 diagram has spine0–1–…–11, leaves12,13 at0, and
leaves14,15 at11. Its previously deleted nine pairs remain at weight0.
All still-paid pairs at distance greater than9 receive weight1/2:

|Representative|Distance|Pairs in its orbit|Selected contexts|
|---|---:|---|---:|
|0,10|10|(0,10),(1,11)|8|
|1,14|11|(1,14),(1,15),(10,12),(10,13)|12|
|2,14|10|(2,14),(2,15),(9,12),(9,13)|36|

The first and third context lists overlap in two cases, leaving54 distinct
contexts. The complete328 affected induced maps produce66 distinct literal
completions after identifying maps that give the same Gram with the old
core labels fixed. The source checks every map, its marked image and all
induced bonds, then records its completion index through the grouped output.
All three pair orbits are invariant under every automorphism of the target.

The new weights sum to `111-10/2=106`. Every marked path of length at most8,
and every elementary ideal triangle payment, retains weight1. The old nine
zero weights retain their existing replacement proofs. In particular no
half-weight is being treated as a deleted pair.

## Residual control before completing a factor

The frozen [complete overlap theorem](fresh_distance8_high_overlap_20261017.md)
applies when no retained short8 payment has already resolved the face. Its
cross-bond exclusions determine each literal old core H. Of the54 affected
contexts,46 use its ordinary Lorentz-residual argument and8 use its
relevant-component argument. Every H here is simply laced, has norm-two
roots and signature `(order(H)-1,1,0)`.

Write J for the inverse Gram of H. The actual preferred corner determines
the necessary contact domain for an additional supporting white root. The
source reconstructs that corner, including all its finite blocks and inverse
diagonals or affine kernels, and checks equality with the saved domain.
The two end sets cover this domain. A root touching both ends supplies a
marked path of length at most8. Within either end, J is entrywise positive
and every two distinct native contacts have projection square greater than2.
Therefore every remaining extra root has at most one contact with H.

For a contact at i, let p=J_ii and rho=2-p. The projection and positive-space
residual squares of a simple-contact root are p and rho. All relevant contact
nodes have p greater than1, so the next possible Coxeter magnitude, sqrt(2),
already has projection square `2p>2`. Thus these contacts must be simple;
there is no label cutoff. The possible p=1 contacts outside the relevant
component, including their doubled-contact alternatives, receive the
separate exact tests inherited from the residual proof.

In each ordinary context all candidate nodes have p>1. For two candidate
roots the exact test is

`J_ij^2 > (2-J_ii)(2-J_jj)`.

If J_ij is positive, the two roots have positive actual product and cannot
coexist. If it is negative, their actual product is negative and their new
bond supplies a marked path of length at most8. The source recomputes every
test and shortcut length. Hence at most one additional supporting white
root can contact H in an otherwise unresolved face. Any such root has a
single simple contact. All later whites are orthogonal to H. The inherited
inverse bound controls entire attached components, so this statement is
not restricted to the immediately adjoining nodes.

In each of the eight relevant-component contexts, the marked component
has exactly one permitted entry node i. When its attaching root G is
present, there can be at most one further core-contact root C. C contacts
a node outside the marked component and is orthogonal to G. The only two
residual possibilities needed here are:

|p_G|p_C|Projection cross product|Residual determinant|Maximum attached tail rank if C exists|
|---:|---:|---:|---:|---:|
|4/3|4/3|2/3|0|0|
|5/4|1|1/2|1/2|1|

These entries, and every other excluded contact type and doubled alternative,
are recomputed exactly. The complete
[residual/no-bridge proof](fresh_distance9_supporting_extensions_review_20261017.md)
applies with the now-checked shortcut threshold8. For determinant0 the two
residuals have a positive zero relation, so no core-orthogonal white can
attach to either. For determinant1/2 the inverse GG entry is2; a white
touching G must be the unique equality root in their residual span, and
acuteness seals it against every further core-orthogonal white. Thus the
tail bounds0 and1 include indirect chains and prohibit bridges to C.

The marked base component has order13 in the first row and12 in the second.
With C present, its whole component can therefore have order at most14 in
either case. An actual order16 affine factor excludes C. Without C, G is
the sole root contacting H, as in the ordinary case. This also shows that
no other old component of H can enter the actual affine factor through
an indirect bridge.

## Every completion is either impossible or a sealed whole factor

Fix one of the328 maps of the marked base component into an actual affine
D15 factor F. The no-bridge argument just given ensures that the remaining
nodes of F are new roots outside H: they cannot be old roots from another
component of the preferred core corner. Consequently each missing target
node is a new supporting white, with its bonds to the base fixed by the map.

If two or more missing nodes must touch H, the one-contact-root conclusion
excludes the map. If its forced entry node has J_ii>2, projection already
excludes the norm-two new root. In every other proper completion exactly
one new root G touches H. Its single simple contact is fixed by the target.
Every other new node is orthogonal to all of H; the entire Gram on H and
the missing F nodes is therefore determined, including its zeros against
old components outside the marked base. This is where the arbitrary-support
argument is needed: merely drawing the missing target edges would not
justify these other zeros.

Exact rational inertia of every remaining proper completion is `(17,2,0)`.
Two negative directions are impossible in the ambient Lorentz space.
Thus no proper completion actually occurs. The exhaustive outcome is:

|Outcome|Distinct completions|Affected maps|
|---|---:|---:|
|Requires multiple core-contact roots|20|64|
|Forced single contact has projection square greater than2|14|112|
|Forced full Gram has two negative directions|22|72|
|Displayed factor already whole; sealed|10|80|
|Total|66|328|

The22 inertia exclusions arise in18 contexts; four contexts have two
different literal completions. All completed matrices and their exact
signatures are saved. The source also checks the positive affine kernel
on the target F before using the full-matrix obstruction.

For the ten surviving cases, H consists of two L25 diagrams sharing a
reversed edge, has18 roots and signature `(17,1,0)`, and its preferred
component is already the entire16-node affine F. The necessary contact
domain has no feasible projection candidate. Hence every additional
supporting white is orthogonal to all of H. The remaining white Gram is
positive definite, since it supports an actual hyperbolic face. Each
other corner is an orthogonal sum of its displayed finite core corner
and that positive-definite remaining diagram. Every displayed finite
component used below is therefore a whole actual component. This does
not assume a simple cusp or any upper bound on the number of spectators.

## Actual finite halves in the ten surviving contexts

The checker reconstructs all40 literal corners, the white positive-definite
principal and all four strict edge principals, the affine positive kernel,
and every whole finite component and marked role. The needed finite halves
are:

|Contexts|Whole small component|Whole large component and marked role|
|---|---|---|
|1312,2495|A4, endpoint|D12, fork center|
|1975,2368,2505,2512|A4, endpoint|D12, short fork leaf|
|1979,2372,2536,2540|D5, marked node|D11, short fork leaf|

These are all already permitted by the
[whole-small-component coupon budget](fresh_finite8_a17_d11_coupon_budget_20261017.md).
The D5 mark is allowed to be any node of that whole small component.
Each finite half is at a different actual corner from the selected ideal
half. Their sum is1. Four cases additionally have a whole finite D16 unit;
the proof does not need to spend that extra payment.

## Composition and cost

The complete finite8/ideal5 theorem first handles the short-path cases,
then its low-core, disjoint-core and complete overlap domains. The present
changes concern only distance10 or11 pairs in order16 D factors, so every
short8 unit and elementary triangle unit is unchanged. The immediate
whole-affine alternatives in the inherited local transfers have order at
most13; their rules are unchanged here. The other domains are retained for
specific reasons: the low-core theorem gives whole component order at most12;
the all-label high-disjoint theorem excludes its entire domain; finite
whole-component, protected-leaf and adaptive finite outcomes retain their
finite rules; and protected/basic ideal outcomes have marked distance at
most9, so retain their full weight. The complete1,153-pattern overlap map
certificate contains all remaining possible order14–17 factors. Its54
contexts listed here are exactly those with an image among the ten changed
pairs. Their rejected maps cannot represent actual factors, and each surviving
one now receives an ideal half and an existing finite half. Thus a list of54
contexts is used only after the complete domain split, not as a census of
all quadrilaterals.

Thus triangles still receive at least2 and quadrilaterals at least1. The
finite budget is unchanged. At an ideal vertex with link-vertex count V,
the contribution of this factor is now `106*V/16`, using the existing
factor-fiber multiplicity, which is linear in the pair weights. The
unchanged B cost107 is the largest order16 cost. A global scalar comparison
must use107/16 and must check its own concentration moves; this package
makes no dimension assertion.

## Verification scope

Primary generation62629 and separate fresh saved replay45499 both terminated
exit0/PASS. Syntax, all22 direct bindings and all6 local report links passed.
An additional calculation using SymPy characteristic polynomials and exact
Sturm root counts independently checked the nine distinct full Grams underlying
all22 Lorentz-index exclusions and all10 sealed cases; their signatures agree.
It did not import the producer's inertia routine. The checker revalidates all direct input bindings,
the54 affected support cases and328 maps, every residual inequality used
here, all66 completed Grams or projection/contact obstructions, and the40
surviving corners. The frozen complete overlap and induced-map enumerations
remain inputs; they are not independently enumerated again. No external
refereeing or formal proof-assistant verification is claimed.

Two separate internal full source/report mathematical reads passed. The
root review also reread the inherited residual/no-bridge proof and checked
the new actual-corner reconstruction. The second review independently
checked why the completion zeros are forced before the index-two test,
why an exterior contact would bound the relevant whole component by14,
and why the surviving finite components are whole. Both checked the full
local-domain composition. Neither review independently repeated the frozen
overlap or induced-map census. Source, certificate and report are final/frozen.
