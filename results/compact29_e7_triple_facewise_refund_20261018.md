# A universal 3/2 charge saving at E7³+E8 vertices

18 October2026. **New paper proof, awaiting independent mathematical
review. No computation, census, optimizer, Lean process or TeX build was
run for this report.** All earlier proof packets remain unchanged.

There is a globally valid nonnegative modification of the unconditional
Vinberg angle rule with charge at most

`171/2 =85.5`

at **every** actual vertex of type T=E7+E7+E7+E8 in a compact hyperbolic
Coxeter polytope of dimension29. If any designated E7 triangle occurs,
the stronger upper81 holds. All other vertex charges can be retained
unchanged. This is a genuine local allowance improvement from87, not a
compact dimension28 theorem: the target at T is84 and the other2,125
extremal types remain separate obligations.

The new geometric ingredient is a complete necessity statement for an
essential quadrilateral at E7 pair(3,5). It includes contacts with both
untouched E7 factors, arbitrary finite opposite labels and dotted opposite
products. It is stronger than checking the previously selected Lambert
witness. The proof uses only finite Coxeter classification, compact acute
geometry and small exact Schur identities.

## 1. A facewise refund which automatically covers all consumers

Let w0 be the original unconditional Vinberg angle weight. Thus every
E7 pair has weight1; every pair in a connected finite factor of rank2,3
or4 has weight2; pairs in different factors have weight0. The whole
T vertex has charge87. We use the original face inequality
`sum w0 >=5-k` on a k-gon.

Call an angle eligible if its vertex has full type T and its marks lie
in one E7 factor. Its old weight is1. On each actual face F let m be
the number of eligible corners and B the total original weight of all
other corners. If m>0, replace every eligible weight on F by

`q(F)=max(0,(5-k-B)/m)`.

All other weights are retained. The old inequality gives q<=1 and all
new weights are nonnegative. The new face sum is
`B+max(0,5-k-B)>=5-k`. Thus simultaneous modifications are valid on
every triangle, quadrilateral and longer face, regardless of the other
vertex types or their incidences. This proves the full consumer check;
no selected list of quad consumers is silently substituted.

This removes the maximum possible total weight from the eligible corners
of each face while leaving all other weights fixed; it distributes that
refund equally. In particular:

* every eligible angle on a face of at least five sides is fully refunded;
* a triangle with an unmodified weight2 corner fully refunds its eligible
  corners;
* a quad with B>=1 fully refunds its eligible corners;
* a quad with B=0 and m>=2 refunds at least1/2 at each eligible corner;
* the only quad giving no refund is a quad whose sole original positive
  corner is the eligible one (B=0,m=1).

The choice is global but purely facewise. It never raises a vertex charge
and does not treat an abstract incompatible triangle pair as a refund.

## 2. Active vertices already save six units

Number each E7 by `0-1-2-3-4-5`, with edge `2-6`.
The complete designated-triangle theorem proves that a triangle at any
of(2,4),(2,5),(3,5) has third normal f whose only old contact is
`f.r5=-phi`, where `phi=(1+sqrt5)/2` and all roots have square2.
It is orthogonal to the other whole components.

Replacing each j in{2,3,4,5} by f gives finite corners. Their changing
rank7 diagrams, respectively, are

`A2+H4+A1`, `A4+H3`, `D5+I2(5)`, `E6+A1`.

For each pair j<k in that four-set, the replacement-j corner puts k,f
in its rank-at-most4 factor, while the replacement-k corner separates
j and f. The three corner payments are1,2,0. All six triangles are
actual: elliptic projection realizes their vertices, and the resulting
three edges close the polygon. Equivalently these are the six triangles
of the actual simplex four-face in the earlier refund note.

The weight2 corner is not a T corner and remains unmodified. Therefore
the facewise rule refunds all six old E7 units at the original vertex.
Its total is at most81. The cross-factor triangle incompatibility is
consistent with this argument but is not needed to turn those six
actual triangle surpluses into valid refunds.

## 3. Essential quads at pair(3,5): statement

**Lemma.** At an actual T vertex, suppose the face at E7 pair(3,5)
is a quadrilateral and its other three original Vinberg corner weights
are zero. Write a=r3,b=r5 and cyclic side order a,b,c,d. Then c,d are
orthogonal to all old roots except, respectively, a and b. Their products
are

`c.a=-sqrt(X)`, `d.b=-sqrt(Y)`, `c.d=0`,

where

`Y=4 cos²(pi/7)`,
`X=(3Y-4)/(5Y-15)>4`.

Consequently the old pair(4,5) is an actual triangle with payments1,0,2,
and its old unit is fully refunded by Section1.

This lemma is a necessary classification under the stated zero-other-
payment hypothesis. It is not a classification of every quad at this
pair, and it does not claim that the resulting local chamber is globally
compact or finite volume.

## 4. Zero-payment corners and the complete finite profile table

The white diagram inside this E7 is K+S, where

`K=(0-1-2-6)=A4`, `S={4}=A1`.

Root a contacts K at2 and S at4. Root b contacts only S. The white
diagram also contains the two whole untouched E7 factors and the whole
E8 factor.

Both new walls are orthogonal to the old E8 by finite maximality. At
each new corner, two marked roots with zero weight must lie in different
components. Indeed, a component meeting both cannot use a whole
untouched E7: adding one new root already makes E8 and cannot admit the
other. Otherwise its rank is at most7, where every connected marked pair
has positive weight. Thus the exceptional zero pairs inside E8 never
provide a loophole here.

It follows that c.b=d.a=c.d=0. Root c cannot contact S. Root d cannot
contact K or S. A new wall may contact an untouched E7 only by its unique
one-node finite extension to E8; it then has no other white contacts.
The two new walls cannot contact the same untouched E7. A wall cannot
contact two such components, or one such component and K, in a finite
corner. Thus all outside-contact possibilities are explicitly retained.

For c not using an outside E7, its finite K+c diagram is either K+A1,
A5, B5 or D5. This exhausts the domain by finite Coxeter classification;
no label cutoff is guessed. Put eta_j=-c.r_j on K,
`C=2-eta^T K^{-1} eta` and `s=(K^{-1}eta)_2`. In K order0,1,2,6,

```
K^{-1}=1/5 [[4,3,2,1],
            [3,6,4,2],
            [2,4,6,3],
            [1,2,3,4]].
```

The complete table is

|K contact of c|C|s|
|---|---:|---:|
|none|2|0|
|label3 at0|6/5|2/5|
|label3 at6|6/5|3/5|
|label4 at0|2/5|2sqrt2/5|
|label4 at6|2/5|3sqrt2/5|
|label3 at1|4/5|4/5|
|label3 at2|4/5|6/5|

A wall using an outside E7 instead has reduced norm1/2 and no K contact:
the E7-to-E8 extension has projection square3/2. Root d, which has no
K/S contact, has reduced norm2 if it uses no outside E7 and1/2 if it does.

## 5. Schur equation and exclusion of every outside E7 contact

Write `x=-c.a`, `y=-d.b`, `Y=y²`. These are actual native Coxeter
magnitudes: zero, finite `2cos(pi/m)`, parallel2 or strictly dotted>2.
After eliminating all white roots, the residual matrix on a,b,c,d is

```
[ 3/10  -1/2   -(x+s)   0 ]
[ -1/2   3/2      0     -y ]
[ -(x+s)  0       C      0 ]
[   0    -y       0      D ]
```

Here D is2 or1/2 as above. Opposite pairs are strictly indefinite by
the actual-face projection lemma. Complete-star corank one gives

`((x+s)²/C-3/10)*(Y/D-3/2)=1/4`,

with both parenthesized factors positive.

First let d use an outside E7, so D=1/2. Then Y>3/4, hence native
Y>=1. The equation gives

`(x+s)²/C=3/10+1/(8Y-6)<=4/5`.

If c uses the other outside E7, C=1/2,s=0 and opposite indefiniteness
requires x>0, so x>=1 and the left side is at least2: impossible.
Otherwise c has one of the seven K profiles. For any nonzero profile,
x>=1 gives left side at least49/30>4/5, so x must be zero.

At x=0 the two label3 endpoint profiles fail strict opposite
indefiniteness. The label4-at6 and label3-at2 profiles have left side9/5
and fail the upper bound. Only label4-at0 and label3-at1 give4/5,
requiring Y=1. Both contain a forbidden singular PSD principal block:

* label4-at0 uses roots{c,0,1,2,3,6}. Eliminating leaves3,6 reduces
  root2's norm to1; successive elimination along2,1,0 leaves the c norm
  equal to zero, with all preceding pivots positive;
* label3-at1 uses roots{0,c,1,2,3,6}. Eliminating leaves0,c and3,6 leaves
  the matrix `[[1,-1],[-1,1]]` on1,2.

The no-contact profile gives
`X=x²=(12Y-4)/(20Y-15)` and `3/5<X<=8/5`. Thus native X=1, which
forces Y=11/8. Since Y<4, it must be a finite-bond square, an algebraic
integer. A rational algebraic integer is an integer, a contradiction.
This exhausts d's outside contact.

Now let c use an outside E7. The preceding argument excludes an outside
contact of d, so C=1/2,s=0,D=2. Opposite indefiniteness gives X>0 and
Y>3. The equation becomes

`X=3/20+1/(4(Y-3))`.

Every native Y>3 is at least `Y7=4cos²(pi/7)>81/25`. Therefore
`1<=X<143/120<2`, forcing X=1. This gives Y=56/17<4, again a
rational noninteger finite-bond square. Contradiction. Both outside
E7 contacts are excluded, without assuming spectator isolation.

The bound on Y7 is elementary: it is the root in(3,4) of
`p(t)=t³-5t²+6t-1`; p is strictly increasing there,
`p(81/25)=-559/15625<0` and `p(13/4)=1/64>0`. Thus

`81/25<Y7<13/4`.

## 6. The isolated essential quad has only the Lambert7 parameters

We now have D=2 and one of the seven K profiles. The equation is

`((x+s)²/C-3/10)*(Y-3)=1/2`, with Y>3.

If x=0, the no-contact and both label3 endpoint profiles fail strict
opposite indefiniteness. The label4-at0 and label3-at1 profiles force
Y=4, a parallel two-wall affine block forbidden by compactness. The
label4-at6 and label3-at2 profiles force Y=10/3, a rational noninteger
finite-bond square. Thus x>0 and native x>=1.

For either B5 or D5 profile, `(x+s)²/C>=81/20`. Hence
`Y<=3+2/15<81/25<Y7`, impossible. For either A5 endpoint profile,
the minimum49/30 gives `Y<=27/8<2+sqrt2=Y8`. Consequently native
Y=Y7. The displayed Y7 interval then gives

`69/25<(x+s)²<143/50`.

But native x=1 gives at most64/25, while native x>=sqrt2 gives at
least `(sqrt2+2/5)²>81/25`. No native x lies between1 andsqrt2.
Both A5 profiles are excluded.

Only the no-contact profile remains:

`(X-3/5)*(Y-3)=1`, or `(5X-3)*(Y-3)=5`.

Native X=1 gives Y=11/2, but c then has label3 contact only with root3;
the old roots{1,2,6,3,4} with c form the same double-fork PSD block,
so this is impossible in a compact chamber. X=2 and X=3 give the
rational nonintegers26/7 and41/12 for Y<4. If
`X=4cos²(pi/5)=(3+sqrt5)/2`, then

`Y=(21+25sqrt5)/22`.

It lies strictly between Y8 and7/2, while Y9>7/2, so it is nonnative.
These bounds follow from `223/100<sqrt5<224/100`,
`sqrt2<71/50`, and the increasing polynomial
`t³-6t²+9t-1` for Y9 on(3,4), whose value at7/2 is-1/8.

If a finite X has label at least7, then X>=Y7>81/25 implies
`3<Y<3+25/66<Y8`. Hence Y=Y7, which conversely forces
`X=3/5+1/(Y7-3)>23/5>4`, a contradiction to finite X<4.
Parallel X=4 is excluded by compactness. Thus X is strictly dotted;
the same interval forces Y=Y7 and
`X=(3Y-4)/(5Y-15)>4`. All opposite labels and distances have been
covered, including zero, finite, parallel and unbounded dotted branches.

## 7. Its forced actual triangle gives a safe full refund

The necessity result says d's only old contact is its label7 bond with
root5. Replacing old root5 by d gives changing finite forest E6+A1.
Replacing old root4 by d gives D5+I2(7). Both are actual vertices by
elliptic projection. Together with the old vertex, their three common
rank28 edges close the polygon at old pair(4,5), so it is an actual
triangle which cannot be further cut by another chamber wall.

Its original weights are1,0,2. The rank-two weight2 corner is not a T
vertex and stays fixed. Thus the facewise refund removes the old E7
unit at(4,5). This conclusion applies to every essential quad under the
lemma's hypotheses, not just the previously constructed witness.

## 8. The universal T allowance and the remaining obstruction

If any designated triangle occurs, Section2 saves six units and gives
charge<=81. Otherwise every E7 factor's chosen face(3,5) is nontriangular.
For each of the three factors:

* a long face refunds its chosen unit;
* an essential quad forces a full refund at that factor's different
  angle(4,5), by Sections3--7;
* any other quad refunds at least1/2 at the chosen angle, by Section1.

The three counted angles are distinct: selected(3,5) and forced(4,5)
are different, and different factors have different root pairs. Thus
the total refund is at least3/2, proving the asserted upper171/2.

The precise unresolved case is now a facewise multiplicity problem.
A nonessential chosen quad whose other paid corners are unmodified
already gives a full refund. The remaining loss can arise only when its
positive corners are all eligible T angles; if there are exactly two,
the symmetric rule saves only1/2 at each. To reach84 at every T vertex,
one must rule out enough of these T-to-T quad configurations, direct
their surplus with a proved global incidence balance, or find additional
distinct refundable angles. An abstract at-most-one-active-factor
statement does not do this. Other extremal forest types must also be
treated before any unconditional compact bound improves.

The exact original maxima must not be obscured: besides the one modified
T type, there remain1447 typed forests of original allowance87,660 of
allowance86 and18 of allowance85. None is changed by this particular
refund selector. These are the2,125 other extremal types from the saved
savings frontier. Even an independently proved active density of1/3
among T vertices would control only the average of this T subclass:
`(1/3)*81+(2/3)*(171/2)=84`. No such density is proved, and it would
not by itself bound the remaining vertex types.

Inherited inputs: the original Vinberg weighting/face theorem and finite
classification; compact acute simplicity, actual-face projection and
affine exclusion; the completed designated-triangle contact theorem;
and the original triangle-refund note's six finite corner diagrams.
The new table, Schur reductions, outside-contact exclusion, essential-
quad classification and facewise refund composition are supplied above
as new paper arguments. They have no new executable certificate or
independent mathematical review yet. The unconditional compact ceiling
therefore remains29.
