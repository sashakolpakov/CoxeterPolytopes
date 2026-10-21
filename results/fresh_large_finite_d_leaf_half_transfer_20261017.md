# Half payments at short fork leaves in large finite D components

In the completed dimension66 payment rule, replace the two unit payments
at distance9 from the short fork leaves of a **whole finite D_n, n>=20**
by halves. Add a half-payment between any marked root of a whole finite
component of order at most5 and a short fork leaf of a separate whole
finite D10. All other payments, including ideal weights, remain unchanged.
Then triangles still receive at least2 and quadrilaterals at least1.
The same-component cost of D_n becomes8n−33 for n>=20.

The [source](../scripts/fresh_large_finite_d_leaf_half_transfer_20261017.py)
and [certificate](fresh_large_finite_d_leaf_half_transfer_20261017.json)
check the affected exact local interfaces. A finite-forest maximum and
any dimension consequence are separate obligations. The review status
is recorded at the end; no new ambient dimension bound is claimed here.

## The complete affected domain

All finite marked paths of length at most8 remain paid, so triangles and
the initial short-path branch retain their finite payments. In the
[universal cycle transfer](fresh_finite8_cycle_replacement_20261017.md),
a second-node attachment gives a whole D_m and marked distance L−1 on
an arc of length L. To need a distance9 payment, L=10. The shorter
cycle arc P is6,7,8 or9, hence m=P+10<=19. These D payments are unchanged.
An endpoint4 gives whole B, also unchanged. An endpoint3 gives whole A
if its new marked leaf has no spectator contact; if it has a contact,
that mark is no longer a leaf and the existing spectator correction
applies. The other mark has two fixed path neighbours. Thus this branch
cannot need a changed short-fork-leaf payment in D_n with n>=20.
The radius-five tree branch and the cycle12/13/14 successors use shorter
finite paths, small whole components or unchanged separate coupons.

After the short-path branch, the complete low-core closure has whole
component bound12 and the disjoint high-core branch is excluded. In the
[complete high-overlap census](fresh_distance8_high_overlap_20261017.md),
2,531 of2,614 selected contexts have whole marked-component bound at
most17. Their payments are unchanged. The nine adaptive contexts protect
a simple marked leaf with its fixed degree-two neighbour, or give a
whole component of order at most16. A short fork leaf of D_n has a
degree-three neighbour, so none of their chosen payments changes.

The remaining74 contexts originally use a protected marked leaf at
distance9. Eighteen protect a simple leaf and its degree-two neighbour;
again they cannot use a changed short fork leaf. Reconstructing every
corner of each of the other56 contexts gives:

|Outcome|Contexts|Retained payment|
|---|---:|---|
|Two distinct protected leaf-nine corners|40|At least one half at each|
|One protected corner plus the new coupon|16|One half plus one half|

Every protected ideal alternative retains its old unit payment. A finite
alternative either keeps its unit or is precisely one of the new halves.
The first40 contexts therefore need no further geometric inference.

All later finite alternatives in the
[rank17 interface audit](fresh_finite_d17_long_pair_support_20261017.md)
have whole order at most17, shorter finite paths or separate coupons.
The change introduces no dependence on a subsequently removed ideal pair.

## The sixteen remaining contexts

Use the norm-two full old core H, with exact inverse J. Each of these
cores is nonsingular Lorentzian. Write a,b for the selected marks, with
b its short fork leaf, and a',b' for the respective opposite black roots.
The selected base component K is already the complete labelled finite
D_s, with13<=s<=17. It contains every old core root except a',b'.
The old root a' is a simple leaf attached only to a in H.

If the actual selected component is not D_n with n>=20, its original
unit is unchanged and the face is covered. Suppose it is such a D_n.
An induced copy of D_s with both short leaves and its fork present can
extend inside D_n only along its existing long arm. Thus all new nodes
of that component lie beyond its old long endpoint j. The checker
reconstructs both canonical D20 maps, while this graph argument covers
every n>=20 without extrapolating a finite rank sample.

Delete the mark a and insert a'. The fork-side part left by deleting a
is exactly D10, because a was at distance9 from the short leaf b. All
ten of these roots already lie in K, and a' is orthogonal to them.
This D10 is a **whole** component of the new actual corner: all other
roots of the selected D_n are separated from it by deleting a, and
every supporting root outside that whole D_n was orthogonal to it.
No unexamined new node or indirect bridge can enlarge this D10.

It remains to bound the whole component containing a'. This requires
the residual argument below; being isolated in the displayed corner
alone would not suffice.

## Complete contacts and residual bounds

Let j be the old long endpoint of K, and p its preceding node. The
complete necessary contact domain of any additional supporting white
is exactly {a',p,j,b'}. Its two end sets are {a',p,j} and {b'}.
The inverse on each end is entrywise positive. A root meeting both
ends already gives a retained marked path of length at most8.

The exact inverse entries, reconstructed separately for all16 contexts,
are

    J[a',a'] = r in {5/4,3/2,7/4,9/4},
    J[p,p] = 2,  J[j,j] = 1,
    J[a',p] = 1,  J[a',j] = 1/2.

Every multi-node support on one end has unit quadratic value greater
than2. Entrywise positivity excludes it for arbitrary native contact
magnitudes, each at least1. Every singleton has inverse diagonal at
least1; labels at least5 have squared contact magnitude greater than2
and are excluded. The complete remaining singleton profiles have labels
3 or4 and are explicitly enumerated. These are all-label necessary
tests, not a finite-label search assumption.

If r=9/4, no extra root can meet a', so its whole component is A1.
Otherwise the only possible root meeting a' is a unit-contact root G,
with residual square rho=2−r<=3/4 in the positive-definite orthogonal
complement of H. Two such roots have a necessarily positive actual
product, so G is unique.

For each other surviving core-contact profile the exact residual Cauchy
test gives one of three outcomes. A root at p, a double-contact root at
j, or a second G would have positive actual product with G and is
impossible. A root meeting b' either is already excluded by its inverse
diagonal or has forced negative product with G; its native bond gives
the marked path a'–G–root–b' of length3 and the face is already paid.
The only profile that can coexist with G in an otherwise unpaid face
is a unit-contact root B at j.

The native G·B must be zero: a product at most−1 would require residual
magnitude at least3/2, larger than the available Cauchy bound. At most
one B can coexist with G. Two B residuals would have to be opposites,
while each has product−1/2 with the G residual, a contradiction.

If B is absent, every additional root reachable from G is orthogonal
to H. The established connected inverse lower bound m/(m+1) controls
the entire attached finite forest. Two nonempty components cost at
least1>rho, and a component of rank4 costs at least4/5>rho. Therefore
there are at most three such roots and the whole a' component has order
at most1+1+3=5. This controls indirect chains, not merely neighbours.

If B is present, their residual Gram is

    S = [[2-r,-1/2],[-1/2,1]].

For r=5/4 its inverse is [[2,1],[1,3/2]]. A core-orthogonal white meeting
G costs at least2 in the inverse form. Equality forces its sole contact
to be the unit G contact and the root to have residual E=−2g−b, with
no additional orthogonal component. A further root meeting B but right
to G would then have positive product with E, impossible; roots right
to G and B are also right to E. This prevents an indirect bridge into
an unbounded B-tail. There is at most one such E.

For r=3/2 the inverse G diagonal is4, so no core-orthogonal root can
meet G. For r=7/4, S has positive kernel(2,1); every core-orthogonal
white must be orthogonal to both G and B. These are the same exact
no-bridge mechanisms as the
[inherited residual theorem](fresh_distance9_supporting_extensions_review_20261017.md),
with the literal matrices checked here. In every case the whole a'
component has order at most5 and stays disjoint from the whole D10.

The alternate actual corner is finite. If it were ideal, the actual
factor-fiber theorem would require both varying marks in the same
irreducible affine factor, whereas b lies in the whole finite D10.
Its new half-coupon and the original corner's retained half therefore
sum to one. This completes all16 remaining configurations.

## Counting and verification scope

There are two short fork leaves in D_n, each with one root at distance9
along the long arm. The two changed pairs are distinct and exchanged
by the diagram automorphism. Replacing their units by halves saves1,
giving8n−33. Between a component of order s<=5 and one D10 the new
coupons have total cost at most s: there are two short leaves and s
possible marks, each with half weight. Existing coupons are unchanged.

Generation85418 and fresh separate-process replay74699 both terminated
exit0/PASS. The checker reconstructs all74
protected contexts, every alternate protected corner, the16 literal
Lorentz inverses, complete necessary contact supports and native label
domains, all residual inequalities and both D20 interface maps in each
coupon case. It records the nine adaptive and universal cycle interfaces.
The geometric completeness of the original census and the inverse-rank
lemma are inherited inputs.

A separate full source/report mathematical review passed. It reread the
complete contact/end-domain helper, protected-leaf helper, connected
inverse-rank and no-bridge proof, universal finite8 cycle/local proofs,
adaptive nine-case proof and later rank17 interface audit. The whole
attached-forest bound, exclusion of other core-contact profiles, residual
equality and positive-kernel sealing, all-rank D_n long-arm restriction,
whole exposed D10, and finite alternate corner all passed. That review
did not duplicate the16 inverse calculations,74-corner enumeration or
original geometric census.

A second reviewer independently inverted all16 saved Grams using rational
arithmetic and checked their Lorentz signatures without producer imports.
All160 end-block entries,128 nonempty end subsets,64 multisupport
exclusions,72 singleton profiles,60 G-interaction rows and12 residual
matrices agree. The check reconstructs all16 exposed D10 components,
their short fork marks and orthogonality to a', and every edge and nonedge
of all32 D20 embeddings. Its complete source/report read also passed the
all-rank extension and no-bridge arguments. That check did not reproduce
the original census, protected-leaf classification, inverse-rank lemma
or a global finite budget. Source, certificate and report are final/frozen.
No finite-forest maximum or ambient bound follows from this local package
alone.

    python3 scripts/fresh_large_finite_d_leaf_half_transfer_20261017.py
    python3 scripts/fresh_large_finite_d_leaf_half_transfer_20261017.py --verify
