# Whole A18 silver weak caps: independent pair-test derivation

21 October2026. Manual mathematical specification, before reading the new
pair worker. No computation or whole-A18 exclusion is asserted. This is
support for the existing arithmetic paper, not an additional manuscript.

The precise inputs are the complete weak inventory
`arithmetic_whole_a18_quadratic_weak_cap_inventory_and_remaining_interface_20261021.md`
(`4710025219a3d8a2d5a485de27d3a9f6962fa61289ef9e6120ee880ee4eec618`)
and the ACTUAL-field composition
`arithmetic_whole_a18_actual_silver_field_reduction_20261021.md`
(`fb185ba3e154863520a62729348339ebde51443f5f94c77ff0abaa528c96485a`).
The old actual feasible-point foundation is3732b826; the actual-foot input
8d93c23c is shared earlier work of this author. Neither a field list nor
a whole-D19 theorem is imported.

## Domain and raw normalization

Let s=sqrt2, K=Q(s), O_K=Z[s]. Old roots a_i,1<=i<=18, have square2,
Gram off-diagonal-1 on the chain, and duals Omega_i. Their inverse is

    B_ij=Omega_i.Omega_j=min(i,j)(19-max(i,j))/19.

Every .52-weak central index i5..14 has one of these44 necessary models:

*20 far-helper caps V=(1+s)Omega_i+Omega_e, e=1 or18;
*4 near4 caps V=s Omega_i+Omega_j at(i,j)=(5,6),(6,7),(14,13),(13,12);
*20 isolated caps V=beta Omega_i, beta=3+2s or6+4s, at every i5..14.

For helpers put eta=1; for an isolated cap put eta=beta. Its K-valued
raw normal is c=A p-V, c^2=2eta, A>0, and

    F=A^2=V^2-2eta.

In particular F=beta*(i(19-i)beta-38)/19 for isolated caps. The raw
marked products are-t and helper-1 for a helper cap, and-beta for an
isolated cap. Neither unit-root rescaling nor adjoining sqrt(beta) is
used. The44 profiles are a necessary superset; they are not all asserted
to be geometrically weak. There are946 unordered model pairs;76 pairs
share the same marked index and cannot be two actual outgoing opposites.

## Exact square and mutual-product tests

For two different marked indices, the orthogonal complement of the FULL
old18-space is one-dimensional over K. Thus F_b/F_a must be a K-square.
Conversely a positive physical witness z in K with z^2=F_b/F_a fixes
their required time ratio and their unique physical raw product

    P=V_a.V_b=sum_(i,j) v_ai*v_bj*B_ij,
    R=c_a.c_b=P-z*F_a.

Both F's are positive physically and negative at the other embedding.
The ratio is therefore totally positive. The positive physical choice
of z is mandatory; changing its sign would put one cap on the wrong
time orientation. No separate sign restriction is imposed on z's other
embedding. The conjugate Gram is then a Euclidean Gram with residual
coordinates related by that conjugate z; it needs no invented conjugate
acuteness condition.

Here is a complete elementary square test for q=a+b*s, a,b rational.
If b=0 and q>0, a square exists exactly when a OR a/2 is a rational
square; the second branch produces a pure-radical witness. For b!=0,
its norm must be a positive rational square c^2. Try BOTH signed values
n=c,-c. A root u+v*s has

    u^2=(a+n)/2, v=b/(2u).

Thus test whether(a+n)/2 is a positive rational square; construct u,v,
and verify the square exactly. The equivalent other equation is
v^2=(a-n)/4. The b!=0 case cannot have u=0, while b=0 explicitly
includes that possibility. Finally replace the root by its negative
if necessary to make its PHYSICAL value positive. Rational-square and
quadratic-field sign tests suffice throughout; norm-square alone does
not suffice.

## All cyclic integrality conditions involving the new edge

Put A=R^2/(eta_a*eta_b). The necessary and sufficient conditions for
cyclic integrality of this PAIR GRAM, given the individually calibrated
old-plus-cap Grams, are

    R in O_K,          A in O_K.                         (1)

This includes the two-cycle; R integral by itself is insufficient in
general when beta=6+4s occurs.

To see completeness, use Cartan entries2(c_i,c_j)/c_i^2. All old-chain
entries and old/cap entries are already integral. A helper has outgoing
entries-t,-1; an isolated cap has outgoing entry-1 and incoming-beta.
The new two-cycle has product A. Every longer simple cycle containing
the new edge uses it once and returns along an old-chain path. Its
product is R times outgoing old-contact ratios from each cap, hence
R times elements of O_K. There is also such a cycle whose contact
ratios are both1: use the helper's unit helper contact, or the isolated
cap's sole contact. Its product is exactly plus or minus R, proving
necessity of the first condition even with either raw norm2beta.
Coincident chosen old nodes give a triangle, not a missing case. Cycles
without the new edge were individually integral already. Arbitrary
closed walks decompose into these cycles, proving sufficiency for the
cyclic criterion, not for ambient chamber existence.

## Physical Coxeter filter

Require R<=0. If R=0 the pair is orthogonal. For R<0 the normalized
squared cosine is A/4. When A<4, positivity of the two-wall Gram and
the positive-inverse projection prove an ACTUAL ridge, so its angle is
pi/m. Since4cos^2(pi/m) belongs to Q(sqrt2), the complete nonorthogonal
finite possibilities are

    A in {1,2,3,2+s}       (m=3,4,6,8).

This is a safe complete alphabet; exact field-square constraints may
remove further values. Labels5,10 generate sqrt5 and12 generates sqrt3;
all other finite labels have larger degree. If A=4, the positive null
sum of the actual UNIT normals is feasible against every other wall,
contradicting compactness. Reject this parallel case. If A>4 the pair
is dotted; there is no finite-label cutoff, but BOTH conditions(1)
still apply. Neither a norm test, a graph edge nor all pairwise tests
prove a compatible global chamber or a simultaneous weak set.

## Required survivor and possible energy consumer

Far helpers(5,18),(7,1), t=1+s, must survive. Their times are
F_a=(200+150s)/19 and F_b=(256+192s)/19, with positive ratio witness
z=4s/5. Their projection product is(202+141s)/19, while zF_a is
(240+160s)/19. Hence R=-2-s and A=6+4s>4. Both are integral and
the pair is genuinely dotted; its presence is not an arithmetic error.

A graph can only constrain weak-index sets. A later energy proof must
use ONE reference with zero coefficients at all retained weak indices,
or pay all actual-foot shifts in that same reference. For surviving
dotted helper caps there is an additional legal actual-point constraint:

    x.(V-(2/t)a_i) >= sqrt(F)-sqrt(2-8/t^2),  t>2.

Indeed c+(2/t)a_i is future timelike, has square-(2-8/t^2), is orthogonal
to c, and is feasible against every other actual wall. For near4,
the finite two-wall projection instead gives, with delta=2-t^2/2>0,

    x.(V-(t/2)a_i) >= sqrt(F+delta)-sqrt(delta).

These pressures include the helper coordinate. It cannot be discarded
or separately charged; any subtraction of its old coefficient must be
paid before using strong floors. These are prospective same-vector
inputs, not a proved exclusion of the known pair or wholeA18.
