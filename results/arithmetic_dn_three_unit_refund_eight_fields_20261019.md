# Three full-D19/D20 units in the eight high-label fields: a standalone consumer proof

19 October 2026. **New written theorem candidate, awaiting independent
review.** No scientific job, proof script, optimizer, Lean process or build
was run. All earlier source-bound evidence is unchanged.

This theorem is relative to the UNMODIFIED A_n baseline. It must not be
combined automatically with the new A_n end-distance-two deletion:
those roles have the same white graph and create additional mixed A/D
quadrilateral consumers. That separate interface remains open.

## 1. Rule, exact scope and inherited geometry

Let P be compact arithmetic simplest-type of dimension n=19 or20,
with ground field one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Start with the distance-at-most-two baseline plus the entire-I2(m>=6)
addition of the eight-field theorem. At each vertex whose ENTIRE
diagram is D_n, write a for the fork, b,e for the two short leaves,
and a-l1-...-lL for the long arm, L=n-3. Delete the three baseline
units at (a,b),(a,e),(b,e), at EVERY whole-D_n vertex.

**Standalone consumer theorem.** The resulting charge still pays
two on triangles and one on quads; larger polygons retain their
required nonpositive payment. In particular the whole D20 cost
would fall from 38 to35. No deletion with an outside spectator
component is asserted, and no global allowance35 follows.

Inputs are compact Coxeter simplicity and finite classification,
acute all-wall projection, absence of actual positive-null principal
diagrams in compact chambers, arithmetic cyclic integrality and
nonphysical positivity, and the explicitly cited local proofs below.
All unit-normal formulas are determinant-ratio or cyclic identities;
they do not assume individual unit entries lie in the ground field.

Triangles remain paid at their unchanged entire high-label I2 corner,
which cannot itself be a whole D_n. All weights stay nonnegative.
Quads containing any selected whisker-pair (b,e) corner are already
covered by the stronger B* proof in
`results/arithmetic_dn_whisker_pair_quad_retention_20261019.md`:
its B* removes all THREE of the present D_n roles at every corner.
Its triangle companion and independent review are
`results/arithmetic_dn_whisker_pair_triangle_retention_20261019.md`
and `results/review_dn_whisker_pair_retention_finite57_20261019.md`.

It remains to cover a newly unpaid quad containing selected fork/leaf
units but NO selected whisker-pair corner.

## 2. Complete original-zero and newly-zero corner inventories

Put N=n-2 in {17,18}, s=(N+2)/(4N), and orient one selected corner
as (a,b). Let the cyclic sides be a,b,c,d and the fixed whites

    W={e,l1,...,lL}=A1+A_L.

An original-zero corner (a,d) has exactly the four finite continuations

    d=I: isolated from W+a;
    d=F: simple contact only with lL;
    d=B: 4-contact only with lL;
    d=T: simple contact only with l_(L-1).

This is the one-root extension list of A_(n-1)=e-a-l1-...-lL:
near endpoint or near next-to-end attachments give marked distance
at most two and are paid. The list is independent of the new c wall.

If (a,d) is newly zero because it is another selected fork/leaf unit,
a already has two white neighbors, so a must remain the fork. Thus

    d is white-isolated,       a.d=-1/2.               (2)

If (b,c) is a selected fork/leaf unit, b has no white neighbor and
must remain a short leaf. Its new fork c must connect simply to e
and to one endpoint of the long A_L. Precisely two possibilities
exist: the NEAR endpoint l1 and the FAR endpoint lL. Then c_bar^2=s.

An opposite selected fork/leaf corner (c,d), with adjacent original
zero corners, similarly forces d=I and c to be the near or far fork.
If d were F,B,T, its one white neighbor and c contact would give
degree two, neither a fork nor a short leaf. These are complete
structural inventories; the far choice is never identified with
the near one.

Because whisker-pair-selected faces were already settled, every
other modified-zero corner is either original-zero or one of these
selected fork/leaf cases. Count the additional selected fork/leaf
corners beyond (a,b). Sections 4-8 cover counts 0,1,2,3 exhaustively.

## 3. Two exact interval tools specific to the eight fields

Suppose A is integral, K=Q(A), 4<A<U, and every other conjugate
satisfies 0<A'<v=8/19. These tools will be used only after their
field-generation and integrality premises are established.

* If U<=92/17, degree at least three is impossible, since
  (92/17)v^2=5888/6137<1. Degree two has integer trace5 and norm1
  or2, so its field is Q(sqrt(21)) or Q(sqrt(17)), outside our eight.
  Degree one gives K=Q, also outside the eight.
* If U<=6, degree at least four is impossible by 6v^3<1. In degree
  three the norm is1; writing the other roots r,t gives A>1/v^2
  and Tr(A)>=A+2/sqrt(A)>361/64+4/5>6. But Tr(A)<6+2v<7,
  impossible. In degree two the trace is5 or6 and norm1 or2,
  giving fields Q(sqrt(21)),Q(sqrt(17)),Q(sqrt(2)),Q(sqrt(7)).
  None is one of the eight. Degree one is again excluded.

These are n19/n20 interval arguments; no extrapolation from an older
n14 small-conjugate bound or assumption that all fields are quadratic
is made.

## 4. No additional selected fork corner: inherited original-zero theorem

In this case all other baseline marks are genuinely zero. The complete
original-zero proof is the composition of

* `arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`;
* `arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md`;
* `arithmetic_dn_far_fork_quad_obstruction_20261019.md`;
* `arithmetic_dn_far_double_quad_obstruction_20261019.md`;
* `arithmetic_dn_far_simple_unpaid_quad_obstruction_20261019.md`.

For d=I, the first source excludes the complete finite c inventory,
including zero, finite and dotted a.c. For F,B,T the finite opposite
product b.d is excluded independently of c; the remaining dotted
cases are covered by the last three sources. The positive-corner
alternative in the far-simple source is c,d as the two short leaves
of D_(n-1)+A1. Its marked distance is two and the unit is RETAINED:
it is not a whole-D_n context. This is exactly the needed original-zero
contradiction, not an assumption that every such quad is geometrically
impossible.

The isolated dotted input may also be replaced in our eight-field
scope by the first interval tool: its integral parameter 4(b.d)^2
has physical upper bound below92/17 and other conjugates below8/n.
The finite c classification is source-bound by the independent review
`review_arithmetic_dn_fork_whisker_classification_finite57_20261019.md`.

## 5. Exactly one additional selected corner: adjacent at the fork

First suppose only (a,d) is selected. Then (2) holds,
c.b=c.d=0, and write b.d=-beta. Both b,d are white-isolated
opposite sides, hence beta>1. Let c_bar^2=t,
a_bar.c_bar=-X, with a_bar^2=s.

The new corank identity is

    (X^2/t-s)(beta-1)=1/2.

The actual simple triangle a-b-d gives 2beta integral in K. Put

    A=2(beta+1)=4+t/(X^2-st).                         (3)

The independent old-D_n-plus-d cost is

    q_d=n*beta^2/2+(n-2)beta+(n-2)/2,
    q_d-1=(beta+1)(n*beta+n-4)/2.

It follows that K=Q(beta)=Q(A) and 0<A'<8/n at every other embedding.

Here one must NOT apply the older isolated-d geometric lemma literally:
that lemma had a.d=0, whereas here a.d=-1/2. The required transport
is instead algebraic. The complete c proof in Sections 3-10 of
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`
uses the following interfaces after the original geometry is reduced:

    W+c finite; c.b=0; the old D_n+c determinant costs/cycles;
    Z integral, K=Q(Z), Z>4, 0<Z'<8/n;
    Z=4+t/(X^2-st).

The new parameter A in (3) satisfies EVERY one of these interfaces.
The W+c inventory, c's contacts with old roots, and all its cycles
are unchanged. Consequently every finite-c rational nonintegrality
or quadratic trace calculation transfers with A in place of Z.

For precision, the sole zero-a.c exception that used a separate
double-whisker lemma also transfers at the formula level. If c only
contacts e with u=cos^2(pi/m), (3) becomes

    A=4(nu-2)/((2n-2)u-n),     0<u'<2/n.

This is exactly the integral scalar interface of
`arithmetic_dn_double_whisker_obstruction_20261018.md`; its label/trace
proof uses only this interface and u's native cosine-square values,
not a.d=0 downstream. Its contradictions therefore remain valid.
The other zero-a.c exclusions are principal affine or positive-index
arguments entirely inside old D_n+c.

For dotted a.c=-x with x>=1, X>=x and t<=1 give

    4<A<=4+1/(1-s)<=264/49<92/17.

The first interval tool excludes this case. Thus the algebraic
transport covers all zero, finite and dotted c patterns, proving
this branch rather than assuming a desired refund theorem.

## 6. Exactly one additional selected corner: adjacent at the leaf

Suppose instead (b,c) alone is selected. Then c is the near or far
fork, c.b=-1/2, while a.d=c.d=0 and d has one of I,F,B,T.
Put a.c=-x. The old D_n+c costs imply

    near: q_c-1=2(x+1)(N(x+1)-1),
    far:  q_c-1=2N(x+1/2)(x+1/2+1/N).

The triangle a-b-c makes 2x integral. In the near case 2(x+1)
has other conjugates in (0,2/N); in the far case 2x+1 has other
conjugates in (-2/N,0). For x zero or finite, physical bounds below
4 or3 give impossible nonzero integer norm in degree at least two.
The independent-block argument gives K=Q(x); degree one is outside
the eight fields. Hence x>1, with the parallel x=1 excluded by compactness.

Write b.d=-beta. In F,B,T the inherited finite-product theorem
excludes finite beta independently of c; beta=1 is parallel and
forbidden. In I, the unchanged white-isolated opposite pair gives
beta>1 directly. Thus beta>1 in every case.

For the NEAR fork c, its white contacts and its contacts with b
equal those of a. Put X=x+1-s. The symmetric a,c Schur equation is

    (X-s)(beta^2-u)=u/2+2rho*beta+2rho^2,             (4)

where u=d_bar^2 and -rho=a_bar.d_bar=c_bar.d_bar.
For I, u=1,rho=0; then 4beta^2=4+2/(X-s)<264/49<92/17 and the
first interval tool excludes it, using the old-D_n-plus-d cost
n*beta^2/2 to establish K=Q(4beta^2),0<(4beta^2)'<8/n.

For F,B,T use u<=9/17,rho<=1/17,X-s>49/34. The left side of
(4) minus its right side increases with beta>=1, and its value
is bounded below by

    (49/34)(1-9/17)-9/34-2/17-2/289=167/578>0.

Equation (4) is impossible.

For the FAR fork, W+c is the path e-c-lL-...-l1. A nonisolated
d=F,B,T would attach at its third or fourth position, producing
a long nonfinite fork (or a branched 4-bond component). Thus joint
finiteness W+c+d forces d=I. Symmetry then gives

    x(beta^2-1)=1/2,       4beta^2=4+2/x<6.

The second interval tool excludes it. Both rooted fork placements
have been covered.

## 7. Exactly one additional selected corner: opposite

Suppose only (c,d) is selected. The inventory gives d=I,
c.d=-1/2, a.d=b.c=0, and c is a near or far new fork. Again
beta=-b.d>1. If c is NEAR, the symmetric corank equation is

    (x+1-2s)(beta-1)=1/4.

For finite x the independent old-D_n+c cost is
1-s+2N(x+1-s)^2. The cycle a-e-c puts 2x in the ground field,
and its other conjugates force
|x'+1-s|<sqrt(s/(2N))<1/10. Apply the cosine-moment lemma in
Section 5 of the inherited isolated-d classification with
kappa=1-s=(3N-2)/(4N), r^2=s/(2N). Its test coefficient is
(11N-26)/(16N)>0, so the degree is at most three and the complete
finite label list is 3,4,5,6,7,9. Label7 has a positive nonphysical
conjugate. Labels 3,5,6,9 fail by their explicit negative conjugates
or the affine a-e-c triangle, as in that source's near row.
Only label4 remains, giving
the ground field Q(sqrt(2)), outside the eight. The value x=0 has
the principal affine cycle a-e-c-l1-a.

For x>1 put Z=4beta^2. It is integral, K=Q(Z), and 0<Z'<8/n.
The displayed corank equation gives

    4<Z<13225/2401<111/20.

Degree at least three is impossible since (111/20)(8/19)^2<1.
Quadratic trace5 and norm1 or2 force Q(sqrt(21)) or Q(sqrt(17));
rational K is also outside the eight. Thus the near branch closes.

For the FAR fork this is exactly the opposed-fork branch completely
classified in the source-bound
`arithmetic_dn_silver_opposed_fork_quad_survivor_20261019.md`, with
rank relation x(beta-1)=1/4. Its only possibilities are the golden
Q(sqrt(5)) and silver Q(sqrt(2)) profiles. Their actual ground-field
equalities, rather than mere subfield containments, exclude them
in the present eight-field scope. This use is limited to that precise
far branch and does not assume it exhausts all quads.

## 8. Two or three additional selected fork corners

With exactly two additional selected corners (three in all), orient
at a selected fork incident to two selected angles; such a fork
exists in the three-edge subgraph of the four-cycle. Then (a,d)
is selected and c is a near/far new fork, with exactly one of c.b,
c.d equal to -1/2 and the other zero. The old-D_n-plus-d cost gives

    -1<beta'<-1+4/n

at every nonidentity embedding. The independent principal block
W+b+c+d has white Schur determinant

    s(1-beta^2)-1/4<0

physically. Hence its conjugate is positive definite and requires
beta'^2<1-1/(4s)=2/n. This contradicts
beta'^2>(1-4/n)^2>2/n, valid for n>8. The eight ground fields
are nontrivial, so such a nonidentity embedding exists.

If all FOUR corners are selected, a,c are forks and b,d short leaves.
The preceding finite-x norm argument again gives x>1, for either
the near or far rooted placement of c. All four adjacent black
products are -1/2. Splitting both opposite pairs into symmetric
and antisymmetric combinations gives

    (X-s)(beta-1)=1,
    X=x+1-s (near),       X=x+s (far).

The integral parameter A=2(beta+1) has K=Q(A) and 0<A'<8/n.
Near it has A=4+2/(x+1-2s)<264/49<92/17; far it has
A=4+2/x<6. The respective interval tools exclude both.

This completes every possible number and placement of the additional
selected fork corners.

## 9. Conclusion and the mixed A/D obligation

Quads with a selected whisker-pair corner are protected by the
stronger B* theorem. Every remaining newly unpaid fork quad was
covered by Sections 4-8. Triangles retain their high-label I2
payment and larger polygons remain nonnegative. Thus the standalone
three-unit deletion has a complete consumer proof, subject to full
independent review of this source and its explicit local dependencies.

At full D20 the cost is 38-3=35. This does NOT cover D19+A1 or
other spectator contexts. It composes with the earlier B reductions:
their white graphs contain a 4-bond, so no face can carry a whole
D_n mark as well. Nonnegative B coupons only help.

It does NOT yet compose with the proposed A_n two-end refund.
The A_n end-distance-two and D_n fork/leaf roles both have white
graph A1+A_(n-3). A mixed corner may become zero through an A_n
deletion even though it had positive original baseline and is not
a selected D_n fork. In particular its adjacent black products need
not satisfy the zero/selected-D_n inventories used here. All newly
unpaid mixed A/D quads must be classified separately before claiming
simultaneous costs35 at A20 and D20.

The standalone local theorem, its field scope, and that remaining
composition gap are separate assertions. There is no global rank20
allowance35, dimension19 claim, new example or formal verification
in this note.
