# Every one-basic-unit D5 fork quadrilateral has the classified Lambert pattern

18 October 2026. **New complete local classification and triangle/quad
incidence restriction.** At a compact arithmetic vertex of type D5+2E8
over Q(sqrt5), suppose its D5 fork-leaf two-face is a quadrilateral.
If the other three marked corners all have zero basic distance-at-most-two
payment, the quadrilateral has exactly the previously classified Lambert
pattern. Consequently none of the eight designated nonadjacent E8 triangles
can occur at that vertex. Four other finite-corner patterns are excluded
for all opposite distances, by arithmetic and exact quadratic identities.

This is a mathematical restriction on **all** quadrilaterals with this
payment property. It does not by itself prove a new globally admissible
weight rule or a smaller vertex allowance. The published finite spherical
classification and the earlier internal Q(sqrt5) field reduction are
inputs; the new finite-pattern reduction and all-distance exclusions below
are explicit proofs. No Lean verification is claimed.

## Setup and scope

Number D5 by edges 01,12,23,24. Its fork leaves are A=3 and B=4; put
W={0,1,2}. The two E8 factors are white spectators throughout. A new side
at either endpoint of an edge of this two-face must be orthogonal to each
whole E8: adjoining a nonorthogonal root to a whole E8 cannot be spherical.
Thus their Gram matrices split off, and the seven relevant roots lie in
a six-dimensional Lorentz space. Eliminating W leaves four black side
normals of rank three. The old pair A,B is one adjacent pair. Write the
cyclic side order A,B,c,d, so the opposite pairs are A,c and B,d.

The basic payment at a marked spherical corner is one when the two marked
nodes are connected by a path of length one or two, and zero otherwise.
No additional D5 or cross-component bonus is included in this definition.
The original marked D5 corner contributes one. The hypothesis is that
the remaining three corners contribute zero under this specific rule.
This wording does not imply that their payments vanish under every
existing global weight rule.

All normalizations below have positive physical scale. The old D5 roots
can be made field-valued and of square two by rescaling the quadratic
form and its roots over k=Q(sqrt5). Put s=sqrt5, phi=(1+s)/2 and
O=Z[phi]. Arithmeticity provides integral cyclic Cartan products and a
positive-semidefinite conjugate Gram at the other real embedding.
These are necessary conditions, not sufficient conditions for a complete
compact chamber. Additional white spectators beyond the two E8 factors
are covered only if they are orthogonal to the four black roots. In
particular this proof does not presume E7 isolation at quadrilateral sides.

An opposite pair together with W is indefinite. If it were positive
definite, the acute M-matrix projection argument in the earlier reviewed
Lambert report would give an actual common point of the two opposite
sides. A singular positive-semidefinite subsystem is excluded by
compactness. We use these strict inequalities below.

## Complete finite-pattern reduction

At the B,c corner, W+B is A4, in path order 0,1,2,B. The finite rank-five
Coxeter classification shows that c either is orthogonal to this A4,
or attaches at exactly one node. All connected rank-five possibilities
are A5, B5 and D5. A new root cannot attach twice to A4, since that
would create a cycle in a spherical diagram. The complete endpoint list is

|Attachment of c to W+B|Bond label|Marked B,c basic payment|
|---|---:|---:|
|None|2 throughout|0|
|White 0|3 or 4|0|
|White 1|3|0|
|White 2|3|1|
|B|3 or 4|1|

There are no omitted high labels: an irreducible spherical diagram
containing A4 and one additional root has rank five, so it is A5, B5
or D5; the rank-two I2(m) family cannot contain the connected A4.
The same list applies to d at its A,d corner. Let N mean no contact,
P03 white-0 contact with label 3, P04 white-0 contact with label 4,
and P13 white-1 contact with label 3. Only these four types can have
zero endpoint payment.

At the c,d corner the new sides must be orthogonal: if their mutual
product were nonzero, their marked graph distance would be one and the
payment would already be positive. The finite diagram W+c+d and the
zero marked payment now leave exactly five unordered patterns:

  N-N, N-P03, N-P04, N-P13, P03-P13.                 (1)

Indeed two contacts at white 0 give marked distance two and positive
payment; two at white 1 give a nonspherical degree-four vertex; a
P04-P13 contact pair would put a label-4 edge in a branched rank-five
diagram, whereas the only spherical rank-five diagram with a 4 is the
linear B5. P03-P13 is D5 and has marked distance three. This proves
the completeness of (1), including all possible nonzero c,d labels
without placing a cutoff on those labels.

The exact script independently checks all seven finite endpoint types
and all sixteen pairs of the four zero-payment endpoint types. Precisely
nine ordered patterns, or the five unordered patterns in (1), survive.

## Field and finite-angle lemmas

When a new side contacts an old root by label 3, its field-valued normal
can be chosen of square two and contact -1. A cycle through old label-3
edges then shows that its opposite product y belongs to O itself,
not just that y squared is in O. A finite opposite product of this
kind is -2cos(pi/m). Membership in k gives real cyclotomic degree
phi_Euler(2m)/2<=2. Thus m is among 2,3,4,5,6. Removing sqrt2 and
sqrt3 leaves only the values 0,-1,-phi.

For an isolated new root c normalize physically to square two and
write its opposite product as -sqrtX. This does not assume that c's
normalization is field-valued. The invariant length-two Cartan product
is X, so X belongs to O. Its white-augmented opposite pair is indefinite
exactly when X>5/2. If it is a finite Coxeter pair, cyclotomic degree
applied to X=2+2cos(2pi/m) gives

  m in {2,3,4,5,6,8,10,12}.

The fields for m=8 and 12 are Q(sqrt2) and Q(sqrt3), so they do not
embed in k. The condition X>5/2 leaves precisely

  X in {phi^2, 3, (5+s)/2}.                          (2)

These are proved finite-angle lists, not sampled high-label cutoffs.
For completeness the elementary totient restriction uses only primes
p with p-1<=4, namely 2,3,5; enumerating their exponents with totient
at most four gives {1,2,3,4,5,6,8,10,12}.

We repeatedly use the following exact nonsquare test. If a+b*s with
b nonzero is a square in k, its rational norm a^2-5b^2 is a square
in Q. A rational positive number r is a square in k only if r or r/5
is a square in Q. All discriminant obstructions below are verified
by these necessary tests.

## Exclusion of N-P03

Choose c of type N, with A,c product -sqrtX, and d of type P03,
with B,d product y. Thus X in O, y in O and X>5/2. The inverse
white Gram is one quarter of

  [3 2 1; 2 4 2; 1 2 3].

After eliminating W and the square-two root c, the residual matrix on
A,B,d is

  [5/4-X/2, -3/4, -1/4;
   -3/4,     5/4,  y-1/4;
   -1/4,     y-1/4,5/4].

Corank one and opposite-pair indefiniteness give

  (5-2X)(3+y-2y^2)+3y-7=0,    y<-1,                 (3)
  X=(4+4y-5y^2)/(3+y-2y^2).

If y<=-2, the denominator is negative and

  X-4=(3y^2-8)/(3+y-2y^2)<0.

Consequently both opposite pairs cannot be dotted or parallel: either
y is finite or X is finite. The only finite y allowed by y<-1 is
-phi. Substitution gives X=(25+12phi)/11, which is not in O.

For finite X from (2), put a=5-2X. The quadratic in y in (3) has
discriminant 25a^2-50a+9. In the order (2) its values are

  134-50s, 84, 134+50s.

The first and third have norm 5456, not a rational square. Neither
84 nor 84/5 is a rational square. Hence no y in k exists in any case.
This excludes N-P03 for all opposite distances.

## Exclusion of N-P13

Use the same X, and let d contact white 1 by -1. Write B,d=y. The
residual matrix after eliminating W,c is

  [5/4-X/2, -3/4, -1/2;
   -3/4,     5/4,  y-1/2;
   -1/2,     y-1/2,1].

Its determinant gives

  (5-2X)(1+y-y^2)+3y-5=0,    y<(1-s)/2,             (4)
  X=y(8-5y)/(2(1+y-y^2)).

For y<=-2 one has X-4=(3y^2-8)/(2(1+y-y^2))<0, again forcing at
least one finite opposite pair. Finite y=-1 or -phi gives respectively
X=13/2 or X=2+5phi/4, both nonintegral. For finite X the discriminant
is 5a^2-14a+9, namely

  26-6s, 28, 34+14s.

Their rational norms are 496,784,176. The first and third are not
rational squares. The rational middle value 28 is neither a rational
square nor five times one. Hence this entire pattern is impossible.

## Exclusion of N-P04, including root scaling

Choose the label-4 new side d of square four and contact -2 with white
0. This is a field-valued normalization: for an original field root
of square a and white contact b, label 4 gives b^2=a; multiplying it
by -2/b gives exactly this norm and contact. Put z=B dot d. The
Cartan cycle through d,0,1,2,B equals z up to sign: its first edge
contributes -2 and its last contributes z/2. Thus z belongs to O,
with no hidden factor-of-two integrality loss.

After white elimination d has square one and old-A product -1/2.
Its residual matrix is exactly that in (4), with z replacing y.
Thus the same equation, the same strict bound z<(1-s)/2, and the same
finite-X discriminant exclusions hold. The normalized opposite
square is z^2/2. If that pair is dotted or parallel, z^2>=8; equation
(4) gives X<4. Therefore at least one of this pair and the X pair
is finite.

For a finite z pair, the finite-angle list for z^2/2 first allows
labels 2,3,4,5,6,10. Label 2 violates the strict opposite bound.
Labels 3 and 5 would require sqrt2 in k; label 6 would require
sqrt6 in k; label 10 would require sqrt(5+s) in k, whose radicand
has nonsquare norm 20. The only remaining finite value is z=-2,
from label 4. Equation (4) then gives X=18/5, which is not in O.
Every branch is excluded.

## Exclusion of P03-P13

Both new roots can now be field-valued of square two. Let A,c=x
and B,d=y; simple label-3 cycles show x,y in O. White elimination
gives the residual matrix

  [5/4, -3/4, x-1/4, -1/2;
   -3/4,5/4, -1/4,   y-1/2;
   x-1/4,-1/4,5/4,   -1/2;
   -1/2,y-1/2,-1/2,  1].

Corank one and opposite indefiniteness give

  (2x-3)((x+1)y^2-(x+2)y)-2(x-1)^2=0,               (5)
  x<-1,    y<(1-s)/2.

Both opposite pairs cannot be dotted or parallel. In fact for
x,y<=-2 the left side F of (5) has negative partial derivative in y.
Hence F(x,y)>=F(x,-2)=2(5x^2+x-13)>=10>0. This contradicts (5).
Thus either x or y is finite. The finite x option is only -phi;
the finite y options are -1 and -phi. Substitution gives quadratics
with the following discriminants:

|Finite substitution|Quadratic in the other variable|Discriminant|
|---|---|---|
|x=-phi|(3phi-1)y^2+(4-phi)y-(6phi+4)|73+89phi|
|y=-1|2x^2+4x-11|104|
|y=-phi|4phi x^2+3x-(9phi+5)|153+224phi|

The two nonrational discriminants have norms 3905 and 7505, neither
a rational square. Neither 104 nor 104/5 is a rational square. Thus
the remaining variable cannot belong to k. This excludes the last
alternative in (1).

## Remaining pattern and coupled incidence budget

Only N-N remains. Its full arithmetic classification, including all
opposite distances and all root scalings, was proved and independently
reviewed in
[the earlier Lambert packet](fresh_compact21_d5_fork_quad_local_witness_20261018.md).
Its unordered opposite squared products are exactly

  {phi^2,(23+9s)/2}.

It forces an actual D5 fork triangle of type L3 or L4. Integral cyclic
products exclude coexistence with either of the two E8 triangle root
types. We have therefore proved the stated classification and exclusion
for **every** D5 fork quadrilateral whose other basic payments vanish,
not just for the one constructed local witness.

At a D5+2E8 vertex let z_1 be the indicator that the fork face is a
quadrilateral with these three zero basic payments. Let t_E be the
total number of designated nonadjacent E8 triangles, among the eight
possible pairs. Then the complete triangle-domain theorem and this
new exhaustive quad theorem give

  t_E+8 z_1<=8.

This is stronger than the previous budget restricted to one prescribed
finite-corner pattern. If z_1=1, all E8 nonadjacent triangles are absent.
If t_E>0 and the fork face is a quadrilateral, at least one of its
other three corners has a positive basic payment.

Turning this incidence restriction into a global allowance reduction
still requires a valid payment rule covering quadrilaterals with more
basic payments and E8 nontriangular faces. The existing bonus families
may already pay some of the zero-basic corners. No existing consumer
has been deleted here, no linear program has been rerun, and no smaller
dimension bound is asserted.

## Exact certificate and reproduction

The [standard-library checker](../scripts/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.py)
checks seven finite A4 endpoint types, all nine ordered zero-extra
quad patterns, exact sparse-polynomial determinant expansions, the
nonintegrality substitutions, and every nonsquare obstruction above.
It does not replace the written finite-classification or cyclotomic
arguments by a bounded numerical sample. Its
[small certificate](fresh_compact21_d5_fork_quad_single_basic_unit_20261018.json)
pins its source, the shared exact-field helper, the prior Lambert
certificate and the field-reduction certificate.

```sh
python3 -B scripts/fresh_compact21_d5_fork_quad_single_basic_unit_20261018.py --verify
```

Generation passed in chunk `04ea1b`, peak RSS 28,295,168 bytes. A
separate-process replay passed in chunk `d1b0a7`, peak RSS 28,327,936
bytes. Both used the shared nonblocking global proof lock, a 59-second
wall alarm, and a checked one-GiB RSS ceiling. All matrices have order
at most five. No opposite-distance search, optimizer, Lean process,
background worker, temporary proof output or bytecode cache was used.
Independent review of this new packet is pending. It is not a claimed
new compact or finite-volume polytope construction.
