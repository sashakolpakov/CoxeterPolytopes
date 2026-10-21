# Complete arithmetic Lannér input over Q(sqrt5)

The [exact checker](../scripts/fresh_sqrt5_lanner_catalogue_20261015.py) and
[certificate](fresh_sqrt5_lanner_catalogue_20261015.json) give the complete local
Lannér catalogue compatible with compact arithmetic ground field
\(K=\mathbb Q(\sqrt5)\): an all-parameter dotted family, **12 triangles, six
order-four types, and four order-five types**. No ambient quadrilateral
classification or dimension improvement is claimed.

Set \(\phi=(1+\sqrt5)/2\), so \(O_K=\mathbb Z[\phi]\) and
\(\phi'=1-\phi\). All roots are normalized to square 2. A finite bond m has
Gram entry \(-2\cos(\pi/m)\). The arithmetic criterion imposes integral cyclic
products in K, a totally real entry field, and positive semidefinite conjugate
Grams at embeddings nontrivial on K. For nonsingular principal blocks this
requires positive definiteness. These are the published conditions recalled in
[Nikulin, section 2](https://arxiv.org/pdf/0708.3991), with the opposite overall
sign convention there. The maximum Lannér order five and the interpretation as
minimal indefinite diagrams are the published inputs from
[Vinberg, section 2.4](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

## All finite labels and dotted pairs

The two-cycle is \(4\cos^2(\pi/m)=2+2\cos(2\pi/m)\in K\). Its real cyclotomic
degree is at most two, giving candidates \(2,3,4,5,6,8,10,12\). Labels 8 and 12
generate Q(sqrt2) and Q(sqrt3), respectively, so the complete finite alphabet is
\[
\boxed{2,3,4,5,6,10}.
\]
There is no finite-label cutoff assumption.

For an order-two Lannér block write its off-diagonal entry as \(-\sqrt t\).
The complete arithmetic conditions are
\[
\boxed{t\in\mathbb Z[\phi],\quad t>4,\quad 0<t'<4.}
\]
These imply total reality of the radical field and conjugate positive
definiteness, and are also necessary. Endpoint equality cannot occur because
an embedding cannot send the nonzero algebraic number t or t minus 4 to zero.
The example \(t=3+\phi\) is verified. Dotted distances are not discretized or
bounded by a search cutoff.

## Exactly 12 triangle types

The relevant contact squares are
\[
0,1,2,\phi+1,3,\phi+2
\]
for labels \(2,3,4,5,6,10\). The classes of 2, 3, and \(\phi+2\) are independent
in \(K^*/K^{*2}\). Indeed 2, 3, and 6 are not squares in K: writing a square as
\((a+b\sqrt5)^2\) forces \(ab=0\), and neither d nor d/5 is a rational square
for \(d=2,3,6\). The norm of \(\phi+2\) is 5, so multiplying it by any product
of 2 and 3 still has nonsquare rational norm. Thus a nonzero triangle cycle
belongs to K exactly when the counts of labels 4, 6, and 10 are each even.
For a right triangle the cycle vanishes and this parity condition is unnecessary.

The checker handles all 56 unordered label triples. Physical hyperbolicity is
\(1/p+1/q+1/r<1\). If the three contact magnitudes are a,b,c, the determinant is
\(8-2(a^2+b^2+c^2)-2abc\). Each conjugate proper two-node block is positive;
the remaining test is the exact positivity of the conjugate determinant.

| Labels | Norm-two determinant |
|---|---|
| 2,3,10 | \(2-2\phi\) |
| 2,4,5 | \(2-2\phi\) |
| 2,4,10 | \(-2\phi\) |
| 2,5,5 | \(4-4\phi\) |
| 2,5,6 | \(-2\phi\) |
| 2,5,10 | \(2-4\phi\) |
| 2,10,10 | \(-4\phi\) |
| 3,3,5 | \(2-4\phi\) |
| 3,5,5 | \(-6\phi\) |
| 4,4,5 | \(-2-6\phi\) |
| 5,5,5 | \(-10\phi\) |
| 5,10,10 | \(-4-12\phi\) |

All listed conjugate determinants are positive. The certificate records every
rejected triple and its exact rejection reason. All cyclic products that pass
field membership are integral, since the individual finite-contact entries
are algebraic integers. The required radical fields are totally real: 2 and 3
are positive, and both conjugates of \(\phi+2\) are positive. Different choices
of radical signs over the same ground-field embedding are vertex sign switches,
because their cycle products agree; testing one resulting Gram suffices.

## Orders four and five

A connected diagram of order at least four with a bond at least 6 contains a
proper three-node path with that bond and another nonright bond. Its relaxed
unit-normal determinant is at most
\(1-\cos^2(\pi/6)-\cos^2(\pi/3)=0\); an extra bond only decreases it.
Parallel and dotted bonds already violate proper two-node positivity.
Consequently every bond here is 2, 3, 4, or 5.

A Lannér graph is connected. If it contains a cycle, the cycle must span all
nodes: an all-3 cycle already has a nonpositive all-ones quadratic form, and
larger labels or chords decrease that value. A proper cycle would contradict
proper-subdiagram positivity. A chord creates a shorter proper cycle. Thus the
only graph topologies to consider are **trees and spanning chordless cycles**.

The checker enumerates all labelled trees using Prüfer words and all labelled
cycles, with every assignment of labels 3,4,5 to their edges. This gives 675
order-four and 13,041 order-five candidates. The finite proper-diagram filter
is independently checked against exact Sylvester minors for all 64 and 4,096
labelled rank-three and rank-four graphs. Full determinants are computed in
the exact basis \((1,\phi,\sqrt2,\phi\sqrt2)\). A tree has no nontrivial cycle;
a spanning cycle has its product in K precisely when its number of 4-bonds is
even. Every accepted diagram is also checked against all proper physical and
conjugate principal minors.

| Order | Diagram | Norm-two determinant |
|---|---|---|
| 4 | Three-arm star with labels 3,3,5 | \(4-4\phi\) |
| 4 | Path 3–5–3 | \(5-4\phi\) |
| 4 | Path 4–3–5 | \(2-2\phi\) |
| 4 | Path 5–3–5 | \(6-5\phi\) |
| 4 | Four-cycle with one 5, otherwise 3 | \(2-5\phi\) |
| 4 | Four-cycle with opposite 5s, otherwise 3 | \(1-7\phi\) |
| 5 | D5 tree with a 5 at the terminal edge of its long arm | \(4-4\phi\) |
| 5 | Path 5–3–3–3 | \(6-4\phi\) |
| 5 | Path 4–3–3–5 | \(2-2\phi\) |
| 5 | Path 5–3–3–5 | \(8-6\phi\) |

There are 78 and 360 labelled arithmetic survivors, respectively. Before the
arithmetic tests the search recovers nine geometric order-four and five
geometric order-five classes. The older repository census restricted to entries
\(0,-1,-\phi\) is a proper subfamily: it does not cover the new triangle labels
4,6,10 or the two displayed paths containing a 4-bond. This package closes that
local-input gap, without claiming novelty for the underlying classical lists.

## Scope and reproduction

Every Lannér subdiagram in a compact arithmetic Q(sqrt5) chamber must belong to
this catalogue. Each listed finite type satisfies the local arithmetic tests.
The list does not classify unions of opposite Lannér pieces, shared white
vertices, cross contacts, or arbitrarily large external white components in a
quadrilateral star. Unlike the Q(sqrt2) order-four/five list, this field permits
several tree types, so a proof using only cycle pieces cannot be transferred.
The existing E7/E8 marked-triangle restrictions likewise do not replace the
missing global quadrilateral analysis.

```sh
python3 scripts/fresh_sqrt5_lanner_catalogue_20261015.py
python3 scripts/fresh_sqrt5_lanner_catalogue_20261015.py --verify
```

Primary generation and saved replay pass, along with source hashes, syntax,
and local links. Independent review has been assigned separately. No dimension
bound is asserted here.
