# A quadrilateral-only priority payment for D5+2E8

18 October 2026. **New local charging theorem with complete quadrilateral
consumer coverage.** In the compact arithmetic Q(sqrt5) setting, a D5
fork-leaf unit at a D5+2E8 vertex may be removed when its two-face is a
quadrilateral and that vertex has a designated nonadjacent E8 triangle.
Every affected quadrilateral still receives total weight at least one
from the remaining basic units and the standard D5 long-short half-units.
No long-face or triangle weight is changed by this rule.

This is a precise conditional payment improvement, not a uniform vertex
allowance reduction or a new dimension ceiling. The distinction matters:
the fork face may be long, or the necessary E8 triangle may be absent.
Those cases retain their original fork unit. All other existing weights
and nonnegative bonuses are retained.

The proof is a finite-corner reduction to eleven patterns followed by
all-distance arithmetic arguments. The companion script checks the
exact determinant identities, finite-pattern list and cyclic-integrality
obstructions. The unit and finite-angle completeness arguments are
written below, rather than inferred from sampled distances. Independent
review and Lean verification of this packet are pending.

## Rule, inputs and reduction

Use the setup of the
[complete one-basic-unit packet](fresh_compact21_d5_fork_quad_single_basic_unit_20261018.md):
D5 edges 01,12,23,24, white W={0,1,2}, old fork leaves A=3,B=4,
cyclic quadrilateral order A,B,c,d and whole orthogonal E8 spectators.
Write O=Z[phi], phi=(1+sqrt5)/2 and s=sqrt5. The basic corner weight
is one for marked graph distance one or two and zero otherwise.
There is also the existing half-unit for a marked long-short leaf pair
of D5 at distance three. The field-reduction and spherical-classification
inputs, field normalizations and opposite-pair projection argument are
exactly as in that packet.

Call a marked corner a fork corner when its full vertex diagram is
D5+2E8 and its marked nodes are the two short leaves of D5. Define
I_E(v)=1 if at least one of the eight possible nonadjacent E8 triangles
occurs at v, and zero otherwise. The proposed rule subtracts one
precisely at a fork corner (v,F) satisfying

  F is a quadrilateral, and I_E(v)=1.                 (1)

All other basic units and every D5 half-unit are retained. A face with
no subtraction keeps its old payment. For an affected quadrilateral,
any retained basic unit already proves payment at least one. Thus it
is enough to examine quadrilaterals all of whose positive basic
corners are fork corners potentially subject to (1).

The complete A4 endpoint list in the previous packet now leaves five
types, instead of four: N, P03, P04, P13 and P23. Here Pij means
contact with white root i by bond j; N means no contact with its
adjacent A4. The added P23 creates precisely a marked D5 fork at
that endpoint. Any other positive endpoint payment is retained and
was already dealt with. Likewise c and d must be orthogonal, since
a nonzero mutual bond would be a retained distance-one basic unit.
The finite W+c+d diagram must have either zero basic payment or a
marked D5 fork. There are exactly eleven unordered patterns:

|Pattern|Disposition|
|---|---|
|N-N|Earlier classified Lambert pattern; no E8 triangle at its fork corner|
|N-P03, N-P04, N-P13, P03-P13|Excluded by the complete one-basic-unit theorem|
|N-P23|Classified below; no E8 triangle at either fork corner|
|P03-P23, P04-P23, P23-P23|Excluded below|
|P13-P23|Two retained D5 half-units pay one|
|P03-P03|Classified below; no E8 triangle at either fork corner|

This is exhaustive. The only basic-positive new corners admitted beyond
the previous five-pattern list are P03-P03 and P23-P23, each a marked
D5 fork. The other additions involve P23 at an old endpoint. Branches
with a label 4 in a rank-five fork are nonspherical; P13-P13 has a
nonspherical degree-four vertex. The exact script checks all 25 ordered
endpoint pairs, leaving 19 ordered or eleven unordered cases. No
opposite-distance bound or arbitrary high-label cutoff is used.

## N-P23: an inherited complete unit equation

Take c of type N and d of type P23. Normalize c physically to square
two and write A,c=-sqrtX, while the field-valued square-two d has
B,d=y. The Cartan cycles give X,y in O. Opposite indefiniteness gives
X>5/2 and y<-1/2. After eliminating W and c the residual matrix is

  [5/4-X/2, -3/4, -3/4;
   -3/4,     5/4,  y-3/4;
   -3/4,     y-3/4,5/4].

Its determinant factors by 2-y, which is nonzero, leaving

  (5-2X)(2y+1)=9.

Put U=2X-5 and V=-2y-1. Then U,V in O, UV=9 and U,V>0. The
conjugate old-D5 Schur inequality for c gives 0<=X'<=8/5, so
-5<=U'<=-9/5. For d its squared height relative to D5 is

  H(y)=1-3y+5y^2/4.

Conjugate positive semidefiniteness requires H(y')<=0, hence
2/5<=y'<=2 and -5<=V'<=-9/5. Both U and V are one modulo 2.
These are precisely the complete UV=9 hypotheses proved in the
reviewed Lambert packet: inertness of 3 in O, all positive units
phi^r, the parity r=3 modulo 6 and the conjugate interval leave
only U=phi^(-3) or U=9phi^3. Thus the unordered two-branch list is

  (X,y)=(phi^2,-(19+9s)/2),
  (X,y)=((23+9s)/2,(1-s)/2).

The second has -1<y<0, so it is not a finite Coxeter product; it
cannot be dotted or parallel either. It is excluded. In the first,
c has the single old contact -phi with leaf A: this is the previously
classified L_A root, of height (-1+5s)/8. It excludes both E8 triangle
types by nonintegral length-two cyclic products.

There are two fork corners, A,B and A,d. Relative to either D5,
the same c is orthogonal to W and to the other short leaf, and has
contact -phi with A. Thus the same L_A exclusion works separately
at both vertices. No designated E8 triangle can occur at either.
This is a necessary branch classification; no complete compact
polytope realizing the branch is asserted.

## P03-P23: impossible

Both roots are field-valued of square two. Put A,c=x and B,d=y.
The label-three cycles give x,y in O and opposite indefiniteness gives
x<-1, y<-1/2. White elimination and corank one give

  (2y+1)(3+x-2x^2)+3x-7=0,
  y=(x^2-2x+2)/(3+x-2x^2).                           (2)

Both opposite pairs cannot be dotted or parallel: for x<=-2,

  y+2=(8-3x^2)/(3+x-2x^2)>0.

The finite-angle lemma then leaves x=-phi or y=-1,-phi. The first
gives y=-(3+12phi)/11, not in O. For y=-1, equation (2) gives
x^2+x-5=0, whose discriminant 21 is not a square in k. For y=-phi
it gives 2s x^2+(3-s)x-(7+3s)=0, with discriminant 134+50s and
nonsquare norm 5456. Thus every branch is impossible.

## P04-P23: impossible, with the label-four scaling retained

Give c square four and white-0 contact -2; give d square two and
white-2 contact -1. Write x=A,c and y=B,d, both in O by the same
Cartan-cycle cancellation used for P04 in the earlier packet.
Opposite indefiniteness gives x<(1-s)/2 and y<-1/2. Corank gives

  (2y+1)(1+x-x^2)+3x-5=0,
  y=(x-2)^2/(2(1+x-x^2)).                            (3)

If the x pair is dotted or parallel then x^2>=8; equation (3) gives
y+2=(8-3x^2)/(2(1+x-x^2))>0. Again one opposite pair is finite.
For a finite x pair the field restriction from the previous P04
proof leaves only x=-2; it gives y=-8/5, not in O. The finite
y options are -1 and -phi. Their equations in x have discriminants
28 and 34+14s, respectively. The first is neither a rational square
nor five times one; the second has nonsquare rational norm 176.
Thus this pattern is impossible for all opposite distances.

## P23-P23: the four-D5 pattern is impossible

All four black roots have square two and contact -1 with white 2.
Put A,c=x and B,d=y. Then x,y in O and x,y<-1/2. Splitting the
white residual into opposite symmetric and antisymmetric pairs gives
the corank relation

  (2x+1)(2y+1)=9.                                    (4)

If both x,y<=-2, equality in (4) forces x=y=-2. The resulting Gram
is rational and contains an indefinite opposite white-augmented
pair; its conjugate is therefore not positive semidefinite. Thus
one pair is finite. Its value must be -1 or -phi. If x=-1 then
y=-5; the principal old-D5-plus-c Gram is rational of determinant
-21, again contradicting conjugate positivity. If x=-phi then
y=-1/2-9s/10, which is not in O. The exchanged cases are identical.
No branch survives.

## P03-P03: a second complete unit classification

Here c and d contact white 0 by -1; the old A and B contact white
2 by -1. All four black roots have square two. Put A,c=x and
B,d=y, so x,y in O and x,y<-1. The white residual splits under
the simultaneous interchange A with c and B with d into the blocks

  [x+1, -1; -1, y+1],
  [3/2-x, -1/2; -1/2, 3/2-y].

The second is positive definite. Corank one therefore gives

  (x+1)(y+1)=1.

Put U=-x-1,V=-y-1. They are positive algebraic-integer units with
UV=1. Conjugate positive semidefiniteness of the first block makes
U',V' negative. Its second block gives U',V'>=-5/2; their product
one strengthens this to -5/2<=U',V'<=-2/5. Thus U=phi^r with r
odd, using the complete unit lemma proved in the Lambert packet,
and 2/5<=phi^(-r)<=5/2. Only r=1 or -1 is possible. Consequently

  {x,y}={-phi,-phi^2}.                               (5)

The finite opposite wall contacts the long endpoint of the original
D5 by -1 and a short leaf by -phi. Its height relative to that D5 is

  H=1+phi+5phi^2/4-2=(1+9phi)/4.

The two E8 third-wall types have heights 4+2s and 1+s. Their length-two
cyclic products with this wall would be

  H(4+2s)=(19+29phi)/2,
  H(1+s)=(9+10phi)/2.

Neither belongs to O. Thus no designated E8 triangle can occur at
the original fork corner A,B. At the other fork corner c,d the
opposite old wall has exactly the same long-end/short-leaf contacts,
with the white path reversed. Its height is again H. Therefore the
same obstruction excludes E8 triangles at that vertex as well.
As before, (5) is a complete necessary local branch classification,
not an assertion of a full compact chamber realizing it.

## The retained-half pattern and payment preservation

For P13-P23 the B,c corner is D5 with marked long-short leaves at
distance three. The c,d corner is another such marked D5. Their
unchanged half-units total one. These marked pairs are not fork
leaves, so the subtraction rule (1) does not affect them.

We can now prove the charging theorem. Take any affected quad. If
it has a retained basic unit, its remaining payment is at least one.
Otherwise all its positive basic corners would be subtracted fork
corners. The complete eleven-pattern list applies. Seven patterns
are arithmetically impossible. In each of N-N, N-P23 and P03-P03,
every fork corner has I_E=0 by the proved cyclic exclusions, so none
of its fork units is subtracted. The final pattern P13-P23 keeps
the two half-units. Thus in every case the remaining total is at
least one. No argument relies on keeping only one of several
simultaneously deleted fork units.

All other quadrilaterals keep their original payment. By definition
there is no change at triangles or faces with at least five sides.
Hence any prior nonnegative global weight rule containing the basic
and standard half-units retains all its face inequalities under
this explicitly quadrilateral-only subtraction.

The new vertex allowance decreases by one exactly when its fork
face is a quadrilateral and an E8 triangle is present. This is
not a uniform saving on D5+2E8 vertices. It also does not extend to
A1+D5+E7+E8 without a fresh spectator-coupling argument. No new
linear program or dimension claim is contained in this packet.

## Reproduction and verification scope

The [small checker](../scripts/fresh_compact21_d5_fork_quad_triangle_priority_20261018.py)
and [certificate](fresh_compact21_d5_fork_quad_triangle_priority_20261018.json)
pin the earlier one-basic-unit theorem, Lambert certificate, field
certificate and original D5 half-bonus source. Sparse rational
polynomials verify the determinant factors exactly. All matrices
have order at most six; the only finite enumeration is the 25
endpoint pairs. There is no distance/unit-height enumeration.

```sh
python3 -B scripts/fresh_compact21_d5_fork_quad_triangle_priority_20261018.py --verify
```

Generation passed in chunk `48c694`, peak RSS 28,393,472 bytes;
separate-process replay passed in chunk `2a7ab3`, peak RSS 28,131,328
bytes. Both held the shared nonblocking proof lock, used a 59-second
wall alarm and checked RSS below one GiB. No temporary proof files,
bytecode cache, optimizer, numerical library or Lean process was used.
The written completeness and payment argument still require independent
review; these exact checks alone are not advertised as a formal proof.
