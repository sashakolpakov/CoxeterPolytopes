# A heavy D17-to-D23 transfer forces a cheaper third corner

21 October 2026. New manual supporting candidate, **pending full root and
outside-author review**. No scientific, symbolic, enumeration, worker-import,
test, Lean or TeX process was run. Frozen papers and rules are unchanged.

## 1. Actual statement and input boundary

Let F be an actual quadrilateral two-face of a hyperbolic Coxeter57-polytope.
Its supporting-white forest is

    W=A13+2A1+D17+D23.

Write the A13 path q1,...,q13 and its two isolated leaves f1,f2. In cyclic
order let the black sides be a,b,c,d. Assume the corner (a,b) has whole
forest 2D17+D23, where a is the long endpoint and b the fork center of
the designated D17. Thus a.q1=-1, while b.q13=b.f1=b.f2=-1, and a.b=0.
All roots have square2. Assume further that c has exactly ONE white
contact: simple3 at the long leaf of the displayed inert D23. Consequently
(b,c) has whole forest D16+D17+D24: finiteness forces c.b=0, since a
nonzero product would join two already branched large D blocks.
No hypothesis is imposed on d except
that it is the actual fourth side and its two incident corners are finite.

**Proposition.** The third corner (c,d) has held finite charge

    C0(c,d)<=2857/8=357.125<364.

This statement exhausts every finite d completion under the displayed
c-contact hypothesis; it is not a selected list of favorable d models.
The three tempting active completions D16, A15+A1 and B15+A1 are all
geometrically impossible here, including every opposite mutual label.
In particular a quadrilateral cannot alternate just between the two
equal-charge372 forests by this simplest rank-transfer mechanism.

Input sources, read as mathematics, are:

| Path | SHA-256 | Used interface |
|---|---|---|
| manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 | Complete held same/cross/context rule, norm-two actual walls and finite-vertex simplicity |
| results/finite57_d16_full_orbit_priority_with_retained_context662_credits_20261021.md | 0b09e32f10721337ddff468dc521aa1a31178d92d7b1e619747ef9113fdd01ae | Exact current C0 convention and distinction between selected financing and global finite charges |
| results/finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md | 3ef4130fc30f557425fdb06144342054e019c6eedda3e54a197c532b158b10ac | Full-white Schur/rank method and actual norm-two native-gap justification; its conclusion alone is insufficient here |

C0 includes the held order11 increments at delta=11/512. They vanish
on every forest appearing below. The author contributed to earlier
Context662 financing and related source audits. This new calculation is
not an independent review of those foundations and spends none of their
reserved source surpluses.

## 2. Complete finite edge alternatives at the initial corner

Deleting a from 2D17+D23 leaves D16+D17+D23. A connected finite
one-node extension of any of these whole large D factors is only a
simple3 extension of its long leaf. Two contacts to one tree form a
cycle; a short-leaf extension would require a nonexistent large E type;
joining two large D factors creates two branch vertices and is not finite.
There are no small exceptional alternatives. Thus, when that edge's
other endpoint is finite, its full forest is one of

| c extension | Full endpoint forest | C0 |
|---|---|---:|
| isolated | D16+A1+D17+D23 |369.5|
| active D16 long leaf |2D17+D23 |372|
| inert D17 long leaf |D16+D18+D23 |365.5|
| inert D23 long leaf |D16+D17+D24 |372|

For example D16+A1 has charge108.5+1+1.5=111, including its
Boolean restoration. This table is a finite-endpoint classification;
an ideal endpoint has not been silently included or ruled out.
The proposition concerns the last row and the actual quadrilateral
continuation, not all four rows.

## 3. All possible d completions, before geometric elimination

At (c,d) the retained white blocks are A13+2A1+D17+D24. If d meets
one of the two large D blocks, it meets its long leaf by simple3. It
cannot meet both large D blocks, or create another fork. The only extra
simultaneous contacts it can have are simple3 to one isolated A1 or to
an endpoint of A13. All other joins would introduce a second fork, a
cycle or a doubled bond in a large branched finite component. Therefore
the complete possibilities involving an inert D block are

| Full forest at (c,d) | Exact C0 |
|---|---:|
| D18+D24+A13+2A1 |344.28125|
| D25+D17+A13+2A1 |350.78125|
| D19+D24+A13+A1 |350.140625|
| D26+D17+A13+A1 |353.640625|
| D31+D24+2A1 |335|
| D38+D17+2A1 |333.5|

These are all below357.125. The formulas include the A13 presence
restoration and its complete A1 cross charge137/64 per A1. No
large/large or large/small D17-or-higher cross term is present.

If d meets neither large D block, these remain D17+D24 with charge
263.5. Finite classification of the whole A13+2A1+d block gives exactly
the following list; T2,T3 denote arbitrary irreducible finite types of
those ranks, including all dihedral labels and H3.

| Active forest | Full active C0, at most |
|---|---:|
| A13+3A1 |84.421875|
| A13+T2+A1 |85.859375|
| A13+T3 |81.5|
| A14+2A1 |92.5|
| B14+2A1 |87|
| D14+2A1 |93.625|
| A15+A1 |102|
| B15+A1 |97|
| D16 |108.5|

A single contact to A13 is endpoint3/4 or next-to-end3; an internal
exceptional completion has rank14>8 and is not finite. A simultaneous
contact to one isolated leaf extends the path to A15 or B15. Contacting
both isolated leaves and the path makes D16. If the path is untouched,
the component on d,f1,f2 has rank at most3. These observations prove
the displayed classification without assuming a finite MAX.

The largest remaining nonexceptional row is D14+2A1: its same charge
is89.5 after the Boolean presence restoration, and each A1 contributes
the ordinary center/two-leaf1.5, mark10 half, and mark9 sixteenth,
total33/16. Thus89.5+2*(33/16)=93.625. The A13+T2+A1 maximum uses
T2=A2 with its full151/32 cross charge. All named guards and order11
increments absent from these small domains stay absent.

It remains to exclude the last three rows, whose d attaches at q13.
The attachment at q1 instead would create, with a, a nonfinite large
exceptional/fork completion at the actual corner (a,d). Thus no reversed
identification of the A13 endpoints is lost in the next calculation.

## 4. Exact native-gap exclusion of all three heavier completions

Put n=14. Schur-eliminate the FULL W, including D23. Its long-leaf
inverse diagonal is1, so the projected c square is1. The projected a
square is15/14, and its product with c is -g, where g=-a.c>=0.
The projected b square is1/14. Set h=-b.d>=0. The native gap for
distinct actual norm-two walls gives

    g^2 in {0,1} union [2,infinity).

This uses actual wall-pair incidence, not an angle assigned to an
arbitrary projected face. All adjacent raw black products in these
three cases are zero. White projections give a.b=a.d=-1/14 and
c.b=c.d=0. There are55 independent positive white roots and four
black roots in the actual58-dimensional Lorentz space. Therefore the
four-by-four Schur matrix must have determinant zero.

For D16, d has the same white contacts as b. Its b-d difference
eigenvalue is2+h>0. Eliminating the positive c coordinate and the
negative b-d sum coordinate forces

    g^2=15/14+2/[196*(h+13/7)],
    15/14<g^2<=14/13<2.

This contradicts the actual native gap, for every h>=0.

For A15+A1, d contacts q13 and one leaf by simple3. The b,d block
after white projection is

    B=[1/14, -h-10/7; -h-10/7, 4/7],
    -det B=h^2+(20/7)h+2>0.

The a couplings are both -1/14. Determinant zero therefore gives

    g^2=15/14+(7/2+2h)/[196*(h^2+(20/7)h+2)].

The ratio in the second term is at most7/4, since subtracting
7/2+2h from (7/4)*(h^2+(20/7)h+2) gives7h^2/4+3h>=0.
Hence15/14<g^2<=121/112<2, again impossible.

For B15+A1, d contacts q13 by simple3 and one leaf by4. Its
projected square is1/14; the b-d off-diagonal is
`-h-1-1/sqrt2+1/14`. The difference eigenvalue is positive, and
the same elimination forces

    g^2=15/14+2/[196*(h+1+1/sqrt2-1/7)].

The denominator's parenthesis is>1. Therefore
15/14<g^2<53/49<2, which is also forbidden. No condition on h's
own finite or dotted label was needed, beyond acuteness h>=0.

All heavier rows are excluded. Adding263.5 to the largest remaining
active charge93.625 proves C0(c,d)<=357.125, as claimed.

## 5. What is still missing for a global discharge

This is a genuine actual-quadrilateral coupling, not an ambient
existence claim or an algebraic realization of 2D17+D23. It shows why
the simplest equal372 rank transfer cannot by itself form a closed
four-corner orbit. It does not imply that every such two-face is a
quadrilateral, that every opposite c takes this role, or that a cheap
third corner accepts credits from arbitrarily many heavy vertices.

A complete global argument still needs the isolated, same-active and
other-D17 rows of Section2, ideal endpoints and higher-sided faces, plus
an actual multiplicity bound for any chosen heavy-to-cheap routing.
The marked pair at the receiving vertex and the entire source face must
be tracked before claiming an injection or assigning any finite credit.
Neither the reserved Context662/D16 source budgets nor an unselected
vertex cap supplies that missing multiplicity control. No finite372
improvement, dimension56 bound, new example or formalization follows.
