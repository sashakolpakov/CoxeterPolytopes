# Surviving-end refunds remove 1,650 compact rank29 heavy types

18 October2026. **New paper proof, awaiting independent review.** No
computation, optimizer, finite census, Lean process or TeX build was run.
All previously settled sources remain unchanged.

A globally valid modification of the unconditional Vinberg angle rule
lowers every original-heavy rank29 vertex containing B4,F4 or H4 to
charge at most84. The possible heavy vertex types consequently reduce
from2,126 to476. The earlier E7³+E8 facewise refund composes with this
rule, retaining its bound171/2 and its active-state bound81. Other
remaining heavy types can still have charge87, so the unconditional
compact dimension ceiling remains29.

## 1. Original rule features used here

Let w0 be the unconditional Vinberg planar-angle weight. Its values are
0,1,2. Marked roots in different finite components have weight0. Every
pair in a connected component of rank at most4 has weight2. Adjacent
marked roots always have positive weight. Pairs contained in a classical
end segment of order at most6 have positive weight. An end segment is a
connected subdiagram containing a terminal node, with the stated whisker
invariance in type D.

These are literal features of the published rule, confirmed in
[Vinberg §6.1–2, printed pp80–81](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
The original face theorem is `sum_F w0>=5-k` on a k-gon. We use the
finite Coxeter classification, compact simplicity and actual finite
corners as the same explicit inherited inputs as the existing compact29
packets. No E8 omitted-pair positions or arithmetic field assumption
enter the new reduction.

## 2. The surviving B3/H3 end lemma

**Lemma.** Let C be a whole rank4 component at a finite vertex of a
compact hyperbolic Coxeter polytope. Suppose b is a leaf of C such that
K=C-minus-b is B3 or H3. Then, for every other root a of C, an actual
triangle at the marked pair(a,b) has another corner of positive original
Vinberg weight.

**Proof.** Let c be its actual third wall. At the corner replacing b,
the full K remains, together with c and all untouched components.
If c contacts K, the finite connected component containing K and c is:

* B_n or F4 when K=B3, since it contains a4-bond;
* H4 when K=H3, since it contains a5-bond and has at least four roots.

In F4 or H4 the marked pair(a,c) has weight2. In B_n, the old B3 is
the short end containing its terminal4-bond. Its old nodes have no
contacts outside C, and the only new root is c. A finite tree therefore
allows contact only at the long end of B3; the four roots K+c form an
actual end segment of B_n. Any outside continuation can attach only
through c. Thus(a,c) lies in an end segment of order4 and has positive
weight. If the finite component has rank at most7, positivity also
follows directly from the full-component rule. These arguments retain
all untouched components; they do not assume c is isolated from them.

Consequently zero weight at this corner forces c orthogonal to all of K.
If the other new corner also had weight0, c would be orthogonal to b:
a nonright product c.b occurs in a finite corner and makes the two marks
adjacent, hence positively weighted. We would then have c orthogonal to
all of C. The finite corner at which b is omitted splits as
`K + (outside+c)`, so outside+c is positive definite. All old roots plus
c would split as `C + (outside+c)`, a positive-definite root set of
order n+1 in a space of positive index n. This is impossible. At least
one other triangle corner therefore has positive weight. QED.

The proof is independent of the dimension, the other whole components,
the label on any new permitted finite bond, and the arithmeticity of the
chamber. It requires the three corners of the actual triangle to be finite,
as holds here by compactness.

## 3. Explicit angle selections and a simultaneous face proof

Use these path numberings and selections:

|Whole component|Non-simple bond|Leaf deletion retaining B3/H3|Selected pairs|Saving|
|---|---|---|---|---:|
|B4|0--1 has label4|b=3|{0,3},{1,3},{2,3}|3|
|H4|0--1 has label5|b=3|{0,3},{1,3},{2,3}|3|
|F4|1--2 has label4|b=0 or b=3|Every pair meeting0 or3|5|

Unspecified path edges have label3. There are five selected F4 pairs;
the unselected pair is{1,2}. Its path reversal preserves the selection.
Each selected angle has old weight2. Lower it to1. The surviving-end
lemma applies at either chosen leaf and proves the required other positive
corner on every triangle containing a changed angle.

This modification may be applied on any chosen subset of vertices, and
at every displayed selected pair in all their whole B4/F4/H4 components,
simultaneously. Here is the complete consumer proof:

* An affected triangle has one changed weight1 corner and another
  originally positive corner. Every originally positive weight remains
  at least1 under this rule, including if that second corner is changed.
  Its sum is therefore at least2.
* An affected quad has a changed corner of weight1, which by itself
  supplies its required payment. No restrictions on its other corners,
  opposite products or spectator contacts are needed.
* All longer faces retain nonnegative weights.
* An unaffected face retains the old inequality.

This deals with all simultaneous consumers; it does not refund several
corners on the assumption that each separately has an unmodified partner.
It also does not need a new quadrilateral classification.

## 4. Apply only at original-heavy vertices

In dimension29 call a vertex original-heavy if its original charge
sigma0 is greater than84. Apply the selected reductions of Section3
only at these vertices. Let the resulting weights be w1. The preceding
proof permits this full-forest selector, so w1 is a valid global face
rule. All positive weights remain positive; in particular its zero set
is identical to that of w0.

If a heavy forest contains B4 or H4 it loses at least3; if it contains
F4 it loses at least5. Original charge is at most87. Thus every such
vertex has new charge at most84. Multiple occurrences add their exact
savings3,3,5, since their marked pairs are distinct.

Every original vertex of charge at most84 is unchanged. This restriction
is useful for composing the earlier E7-triple refund: its protected
weight2 source corners all have original total at most84.

## 5. Exact count of the remaining heavy types

The existing complete heavy-type reduction has five allowed rank4
types A4,B4,D4,F4,H4. Removing every forest containing B4,F4,H4 leaves
only A4,D4 at rank4. A multiset of m such rank4 components has exactly
m+1 possibilities. Other type counts remain four at ranks7,8 and three
at ranks9,10, with repeated factors unordered.

The following counts are elementary products of these multiplicities,
not a new computer census:

|Rank partition|Original allowance|Original types|Remaining types|
|---|---:|---:|---:|
|4+4+4+4+4+9|87|378|6*3=18|
|4+4+4+7+10|86|420|4*4*3=48|
|4+4+4+8+9|87|420|4*4*3=48|
|4+4+7+7+7|87|300|3*binom(6,3)=60|
|4+7+8+10|86|240|2*4*4*3=96|
|4+7+9+9|87|120|2*4*binom(4,2)=48|
|4+8+8+9|87|150|2*binom(5,2)*3=60|
|7+7+7+8|87|80|binom(6,3)*4=80|
|9+10+10|85|18|3*binom(4,2)=18|
|Total||2126|476|

Hence1,650 of the2,126 original heavy types are controlled at84 by an
actual valid face rule. The remaining476 have allowances314 at87,
144 at86 and18 at85 before adding the separate E7-triple refund.
The list still includes B7/B8/B9/B10; only B4 has been removed.

This is a necessary-type refinement, not exclusion of a dimension.
The original350-heavy-vertex lower count remains valid for the modified
rule because its global strict h-vector balance is unchanged and its
maximum vertex excess is still3.

## 6. Composition with the E7³+E8 refund

Now use w1 as the background in the prior facewise refund. At T=E7³+E8,
the eligible E7 angles still have weight1. For each actual face let m
count such corners and B1 sum all other w1 weights. Set their weight to
`max(0,(5-k-B1)/m)`, leaving all other weights fixed. The same algebra
proves every face inequality simultaneously. No weight is increased.

For an active T vertex, the six forced triangles' weight2 source vertices
have these changing forests, in addition to the unchanged2E7+E8 of
weight66:

|Changing source forest|Original full charge|
|---|---:|
|A2+H4+A1|80|
|A4+H3|84|
|D5+I2(5)|78|

None is original-heavy, so all protected weight2 source corners remain
unchanged under Section4. All six E7 refunds still apply, and active
T vertices retain the upper81. This is why the heavy selector matters:
applying the B/F/H4 changes at every vertex would reduce one H4 source
angle and would not retain the old six-refund proof unchanged.

For an inactive T vertex the prior essential-quad theorem still applies:
the zero set of w1 is unchanged. A chosen quad with zero other weight
is still an original essential quad and forces an actual(4,5) triangle
with an unmodified I2(7) weight2 source; that source has original full
charge78. Every other chosen quad has either an unmodified positive
payment or at least two eligible T payments and hence refunds at least
1/2. Long faces refund the whole eligible unit. The three factors give
distinct refunds totaling at least3/2. Thus every T vertex still has
charge at most171/2, with81 in the active case.

The simultaneous resulting unresolved allowances are therefore:

* 313 other types at87;
* 144 types at86;
* 18 types at85;
* the T type at171/2, with only inactive T vertices capable of excess>0.

No density of active T vertices is assumed. Even a density result inside
that one subclass would not control the other475 types.

## 7. Precise next blockers: A4/D4 and E8 zero positions

The surviving-end proof almost applies when deleting a leaf of A4 or D4
leaves A3. Classical and rank-at-most7 finite continuations are handled
by their small terminal segment, but an A3 continuation can now sit in
E8. The exact weight of its marked pair is indispensable.

The current original certificate imports only the number of four E8
omitted pairs, not their positions. No general positivity for all E8
pairs at distance at most3 is assumed here. In fact the primary paper's
§6 Lemma2 expressly permits an exceptional E8 pair involving a terminal
node at distance3 or4. Thus a guessed short-distance extension to A4/D4
would be unsafe. These exceptional continuations are a concrete remaining
consumer family, rather than an omitted label cutoff.

Resolving the exact A3-to-E8 exceptional routes, restoring their lost
triangle units at source vertices with sufficient actual capacity, or
showing incompatible simultaneous exception walls could extend the rule.
If all rank4 types were covered, only the80 types of7+7+7+8 and18 of
9+10+10 would remain. This last98 count is a conditional target, not a
proved further reduction. The present unconditional new bound is the
476-type necessary frontier with the allowances stated above.

This report is stable for independent mathematical review. Its new
claims are the surviving-end triangle lemma, the explicit3/5/3 reductions,
their all-consumer proof, the1,650-type count and the compatible E7-triple
composition. No executable verification or Lean status is claimed for
this new paper argument.
