# Complete extreme-E8 quadrilateral consumers, including D5 spectator contacts

18 October 2026. **Written proof and separately replayed exact coefficients.** At a
compact arithmetic Q(sqrt5) vertex D5+2E8, consider the two-face marked by
canonical E8 roots (0,6). If it is a quadrilateral and its other three
corners all have zero basic distance-at-most-two payment, its two new
walls are exactly F0 and F6. Both designated E8 triangle bundles are
therefore present. This treats every finite endpoint and every possible
D5 spectator contact, without an opposite-distance cutoff.

This is completeness for one particular E8 pair, not for all E8 long
pairs. It suffices for a complementary inactive-state charge saving
at this pair. The proof uses the earlier complete field and finite-angle
inputs. Independent review remains to be completed before this packet
is frozen.

## Endpoint domain and normalizations

The old E8 is the chain 0-1-2-3-4-5-6 with root 7 attached to 2.
Use old black roots A=0,B=6 and white D6={1,2,3,4,5,7}; the whole D5
and second E8 are additional white spectators. The latter E8 forces
both new sides to be orthogonal to it. Write the quadrilateral sides
in order A,B,c,d. If c,d have a nonzero mutual bond, their corner
already has a basic unit, so the zero-other-basic hypothesis forces
c dot d=0.

Deleting A leaves D7. A new c which is not orthogonal to that D7
can extend it only to D8 or E8: the spherical classification leaves
one label-3 attachment at a leaf. Attachment at B makes the marked
B,c pair adjacent, giving a basic unit. The zero-payment possibilities
are therefore S=attachment at white 1, or T=attachment at white 7,
each extending D7 to E8. Otherwise c is N, orthogonal to D7.

Deleting B leaves E7. Its only nonorthogonal spherical extension is
E8, obtained by attaching d at white 5 with label 3; call this L.
Otherwise d is N. An extension S,T,L is orthogonal to the whole
D5: joining D7 and D5 would produce a rank-thirteen spherical tree
with two branch vertices, while joining E7 and D5 would require a
spherical extension of E7 larger than E8. Both are impossible by
the complete finite classification.

Only an N root can contact D5. Its spherical D5 extension is D6
or E6, by a single label-3 contact at the long leaf or one of the
short leaves, respectively. After physical norm-two normalization
its D5 projection square is delta=0,1,5/4. If both N roots contact
D5 and c,d are orthogonal, their combined spherical extension
allows only one long-leaf and one short-leaf attachment. Two short
attachments give the affine three-arm tree (2,2,2); two at the same
leaf give two branch vertices. The permitted distinct attachments
give E7, with mutual D5 projection 1/2.

This proves that the complete endpoint domain is

  c in {N,S,T}, d in {N,L},

with the stated and initially unrestricted D5 contacts on each N.
There are no omitted high labels: all connected rank-eight
extensions of D7 or E7, and all rank-six extensions of D5, have
been covered by the spherical classification. A root cannot join
two of these old connected factors into an unlisted larger finite
component.

For S,T write alpha=A dot c; for L write beta=B dot d, with
norm-two field normalizations along the fixed label-three paths.
Their cyclic products along those paths give alpha,beta in O=Z[phi].
For an N root put its opposite product -sqrtX or -sqrtY. The
Cartan squares X,Y belong to O even when a norm-two normalization
is not field-valued. Independent old components need not have the
same field square class; no such assumption is used. If both N
roots contact D5, the simple cycle through D5, the two new roots
and old E8 gives t=sqrt(XY) in O. All physical opposite products
are nonpositive.

The inverse E8 entries needed below are

  J00=4, J01=7, J07=5, J11=14, J77=8,
  J55=6, J56=3, J66=2, J06=2.

Eliminating the white D6 gives old-black matrix [1/2,-1/2;-1/2,1].
The whole positive D5 is then projected off as well.

## Complete short unit intervals used in the branches

Every positive unit of O is phi^r, as proved in the reviewed Lambert
packet. If its conjugate is negative, r is odd. Both 2 and 3 are
inert, since t^2-t-1 has no root modulo either prime. Thus for
U,V in O with UV=p, physical U,V positive and conjugates negative,
one has U=p^a phi^r, a=0 or 1 and r odd. The conjugate intervals
below leave these exact possibilities:

|Product and conjugate interval for U|Possible U|
|---|---|
|UV=1, -1<=U'<=-1/2|phi|
|UV=1, -2<=U'<=-1|phi^(-1)|
|UV=2, -2<=U'<=-1|phi^(-1), 2phi|
|UV=3, -3<=U'<=-3/2|phi^(-1), 3phi|
|UV=3, -2<=U'<=-1|phi^(-1), 3phi|

Each row follows from p^a phi^(-r) lying in the indicated closed
interval. Its neighboring odd exponents differ by two, multiplying
the value by phi^2; the next possible exponent on either side
already lies outside the interval. This exhausts all integers r,
not a sampled range. Parity obstructions below take place in O/(2),
which is the field with four elements; in particular 1 and 3 are
not divisible by 2.

Conjugate positive semidefiniteness gives the following bounds.
For N at A with D5 projection delta, let a=2-delta; then
0<=X'<=a/4. For N at B let b=2-delta; then 0<=Y'<=b/2.
For S,T,L the conjugate height polynomials give, respectively,

  alpha' in [3/2,2], [1,3/2], and beta' in [1,2].

Indeed their heights are 12-14alpha+4alpha^2,
6-10alpha+4alpha^2, and 4-6beta+2beta^2. These bounds provide
the conjugate intervals in every unit reduction below.

## Neither new root is N

For S-L, the white residual splits into opposite symmetric and
antisymmetric pairs. Corank one gives (alpha-1)beta=1.
Put U=1-alpha,V=-beta. Then UV=1 and U' is in [-1,-1/2],
so U=phi. Thus alpha=1-phi is strictly between -1 and zero,
which is not a finite Coxeter product; it is not dotted either.
For T-L the equation is (1-2alpha)(-beta)=2. Here U=1-2alpha
has U' in [-2,-1], so U=phi^(-1) or 2phi. The first makes
alpha positive; the second makes alpha=1/2-phi nonintegral.
Thus both patterns are impossible.

## One new root is N

For N-L let a=2-delta on c. Eliminating c gives

  (2X-a)(-beta)=a.

For S-N and T-N let b=2-delta on d. Their equations are

  (1-alpha)(Y-b)=b/2,
  (1-2alpha)(Y-b)=b,

respectively. Opposite indefiniteness ensures every displayed
factor is positive; the factors involving alpha are at least one
by acuteness. The exact branch reductions are as follows.

|Pattern|delta|Reduction and obstruction|
|---|---:|---|
|N-L|0|(X-1)(-beta)=1 gives X=phi^2, beta=-phi^(-1), not a Coxeter product|
|N-L|1|(2X-1)(-beta)=1 gives X=phi^2/2, not integral|
|N-L|5/4|(8X-3)(-beta)=3 gives 8X-3=phi^(-1) or 3phi, both making X nonintegral|
|S-N|0|(1-alpha)(Y-2)=1 gives alpha=1-phi, not a Coxeter product|
|S-N|1|2(1-alpha)(Y-1)=1 contradicts integrality modulo 2|
|S-N|5/4|(1-alpha)(8Y-6)=3 contradicts divisibility by 2|
|T-N|0|(1-2alpha)(Y-2)=2 gives 1-2alpha=phi^(-1) or 2phi, hence alpha positive or nonintegral|
|T-N|1|(1-2alpha)(Y-1)=1 gives 1-2alpha=phi^(-1)<1, contradicting acuteness|
|T-N|5/4|(1-2alpha)(4Y-3)=3 gives 1-2alpha=phi^(-1) or 3phi, hence alpha positive or nonintegral|

For example X=(3+phi^(-1))/8 or (3+3phi)/8 in the third
row, whose coefficients in the integral basis 1,phi are not
integers. The other asserted nonintegral values have an odd
half-coefficient of phi. The short unit-interval table proves
all listed alternatives without any bound on opposite distances.
Every mixed endpoint branch is therefore excluded.

## N-N, including every D5 contact profile

Let delta_c,delta_d be the two D5 projection squares, put
a=2-delta_c,b=2-delta_d, and let C be their mutual D5 projection.
The heights are Hc=4X-a,Hd=2Y-b. The old-E8 mutual positive
projection is 2sqrt(XY), so c dot d=0 gives

  4XY-4bX-2aY+ab-4C sqrt(XY)-C^2=0.                 (6)

Opposite indefiniteness gives X>a/2 and Y>b. When at most one
root contacts D5 one has C=0, and (6) factors as

  (2X-a)(Y-b)=ab/2.

All one-contact cases exclude:

|delta_c,delta_d|Reduction and obstruction|
|---|---|
|1,0|(2X-1)(Y-2)=1 gives X=phi^2/2, nonintegral|
|5/4,0|(8X-3)(Y-2)=3 gives 8X-3=phi^(-1) or 3phi, both giving nonintegral X|
|0,1|2(X-1)(Y-1)=1, impossible modulo 2|
|0,5/4|(X-1)(8Y-6)=3, impossible modulo 2|

If both roots contact D5, the complete finite domain leaves only
long/short or short/long. In either case C=1/2 and
t=sqrt(XY) belongs to O by the closed Cartan cycle described above.
For delta_c=1,delta_d=5/4, twice (6) is

  8t^2-4t+1-6X-4Y=0,

impossible modulo 2. For delta_c=5/4,delta_d=1, it is

  8t^2-4t+1-8X-3Y=0.                                (7)

Here 0<=X'<=3/16 and 0<=Y'<=1/2. The complete finite-angle
Cartan-square list in k is 0,1,2,phi^2,3,(5+s)/2. The physical
bound X>3/8 excludes zero; every other finite candidate has
conjugate greater than 3/16. Thus the X pair is dotted. Both
pairs cannot be dotted or parallel: if X,Y>=4, the left side of
(6) is

  4XY-4X-3Y/2+1/2-2sqrt(XY)
  >=7XY/2-4X-3Y/2+1/2>=69/2>0.

So Y must be finite. The bounds Y>1,Y'<=1/2 leave only Y=phi^2.
Substituting X=t^2/Y into (6) yields a quadratic in t with
discriminant 20+8s, whose norm 80 is not a rational square.
It has no root in k, contradicting t in O. Thus every nonzero
D5 contact profile is excluded.

The remaining case delta_c=delta_d=0 has

  (X-1)(Y-2)=1.

Both factors are positive, their conjugates lie in [-1,-1/2]
and [-2,-1], and the unit lemma gives X-1=phi. Consequently
X=Y=phi^2. The new walls have the single old contacts -phi at
canonical roots 0 and 6 and are orthogonal to every other old
root. They are exactly F0 and F6. The earlier complete triangle
and product-face proof therefore forces both E8 triangle bundles.
This proves the stated complete consumer classification.

## Complementary charge rule and an actual linear-model successor

Let I_E indicate any designated E8 triangle at the D5+2E8 vertex.
In the inactive state I_E=0, remove the connected bonus at pair
(0,6) in each of its two E8 factors. Such a marked face cannot
be triangular, by inactivity. If it is a quad, the classification
above forces another basic unit at one of its other corners.
That unit is retained: the whole old D5 is white on this E8 face,
and cannot join the connected white D6 in a spherical corner: that
would contain two old fork vertices in a component of rank at least
twelve, whereas finite components of that rank are A, B or D.
Its component therefore has rank at most seven after adjoining the
two black roots, so it cannot become an E8. At a corner whose full
forest is D5+2E8 it must remain the entire D5 component, consisting
only of white nodes. Hence the marked black pair cannot be its
fork pair, and the active-state D5 subtraction cannot remove the
other basic unit. Long faces retain nonnegative weights. This proves
the deletion is safe simultaneously with the active-state D5
fork deletion from the earlier all-face corollary.

If b is the old bonus on this particular canonical E8 pair, the
two E8 components have the same whole-vertex context in the old
model, so the combined saving is

  1 if I_E=1, and 2b if I_E=0.

Therefore the pointwise uniform saving on this forest is at least
min(1,2b). No assumption b>0 is made. A concrete successor of the
old bonus model keeps all existing unpaid quad rows and replaces
the D5+2E8 forest budget by the two linear inequalities

  old_total-1<=z,    old_total-2b<=z.

Why may the old unpaid rows remain unchanged? They have no basic
unit at any corner. If one of their extreme E8 bonuses could be
deleted in the inactive state, the completeness theorem would
produce another basic unit, a contradiction. Basic-positive
quad rows affected by the D5 deletion have already been covered
by the eleven-pattern theorem and its
[strengthening without fixed half-units](fresh_compact21_d5_mixed_fork_quad_exclusion_20261018.md).
That strengthening is essential for the old golden-field baseline,
which does not include fixed D5 half-units. Every triangle route remains
unchanged and every long face remains nonnegative. Thus this is
a genuine face-dependent rule, not a guessed subtraction from
an old forest allowance.

This successor has not been optimized or checked against a saved
dual. It may still leave the best maximum above forty, and the
A1+D5+E7+E8 forest needs its own coupling analysis. A positive
or even large saving at this forest alone is not a dimension
theorem. The full written branch proof and exact coefficient
certificate must be independently reviewed before any such claim.

## Reproduction and resource record

The [small coefficient checker](../scripts/fresh_compact21_extreme_e8_quad_consumer_completeness_20261018.py)
and [certificate](fresh_compact21_extreme_e8_quad_consumer_completeness_20261018.json)
check the inverse entries, the stated finite endpoint cases, all
eleven spherical D5 two-root profiles, sparse determinant identities
and every short unit interval. The latter check includes exact
monotone tail inequalities excluding all omitted odd exponents.
The nonsquare field obstruction is verified exactly. Written
spherical and finite-angle completeness arguments are still needed;
a bounded numerical scan is not used in their place.

Generation passed in chunk `a1bdc1`, in 0.262 seconds, peak RSS
28,196,864 bytes. Separate-process replay passed in `b3fe23`, in
0.227 seconds, peak RSS 28,180,480 bytes. Both held the shared
nonblocking proof lock, used a 59-second alarm and checked RSS
below one GiB. Matrices have order at most eight. No optimizer,
distance queue, numerical linear algebra, temporary proof output,
bytecode cache or Lean process was used.

```sh
python3 -B scripts/fresh_compact21_extreme_e8_quad_consumer_completeness_20261018.py --verify
```
