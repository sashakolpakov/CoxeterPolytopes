# Cycle12 antipodal transfer with sealed finite coupons

**FINAL/FROZEN after generation, fresh saved replay and full source/proof review.** An
actual quadrilateral whose varying affine factor is a12-cycle, with old
marks at antipodal distance6, has a retained finite8/ideal5 payment or
receives one unit from these finite coupons:

* A marked whole A1 against the unique middle node of a whole A11:
  weight1.
* Either marked node of a whole A2 against the unique middle node of a
  whole A11: weight1/2.

Thus its six antipodal cycle pairs can be removed, retaining the60 pairs
at distance at most5. This is a local quadrilateral transfer. A finite
vertex budget, triangle-retention check and global dimension composition
are separate obligations; none is claimed here.

The [source](../scripts/fresh_cycle12_ideal5_transfer_20261017.py) and
[certificate](fresh_cycle12_ideal5_transfer_20261017.json) adapt the complete
[cycle13 proof](fresh_cycle13_ideal5_transfer_20261017.md), using the frozen
[energy-one classification](fresh_spectator_energy_one_boundary_20261017.md).
The literal checks use rational congruences and Fraction congruence
elimination; the formal determinant check uses SymPy. They corroborate
the written arbitrary-rank argument rather than enumerating spectators.

## Complete endpoint profiles

Use the actual factor-fiber face star W+a+b+c+d, with opposite pairs a,c
and b,d. At the old ideal corner its varying factor F is the12-cycle and
the remaining diagram Z is positive definite and orthogonal to F. The
white-plus-one-black edge diagrams are strictly positive definite. A
nonzero adjacent black contact or a shared spectator contact already
gives the retained short-path payment, so assume these are absent. Write
the nonnegative contact vectors to Z as x,y; their supports are disjoint.

Deleting a from F gives A11. The complete positive-semidefinite one-root
extension list is isolated, one3 or4 bond at either endpoint, one3 bond
at either second node, or3 bonds at both endpoints. This inherited
all-label classification still applies: no exceptional finite or affine
diagram contains A11. In each single-contact case the adjacent marked
distance is at most6. The opposite fixed simple leaf and its fixed
degree-two neighbor exclude every classical affine completion; spectators
can attach only at the new root and cannot alter those fixed degrees.
Thus the entire connected corner component is finite and pays finite8.
All other components are principal parts of Z and are positive definite.
The source verifies the six core positive-definite diagrams, distances
and fixed leaf/neighbor roles.

Only isolated and both-end3 profiles remain. Put alpha=1/6+1/6=1/3.
The normalized Schur determinant is separately concave in reciprocal
heights u=1/h,v=1/k. Spectators add

\[
\alpha\|x/h-y/k\|_{G_Z^{-1}}^2\ge0.
\]

The four ordered rectangle minima are−1/3 for two isolated profiles,
zero for either mixed order, and1/3 for two both-end3 profiles. All16
rectangle vertices are checked exactly. Two both-end3 profiles are
therefore impossible under the unpaid hypothesis.

For isolated/both-end3 the base determinant is

\[
1-\frac23u^2-\frac23v,\qquad0<u\le1,\quad0<v\le\frac12.
\]

It vanishes only at u=1,v=1/2. Corank forces equality in the nonnegative
spectator correction. Positive definiteness of Z gives x/h=y/k;
disjoint supports then give x=y=0. The adjacent corner containing the
isolated new root has that root as a whole A1 and the other old mark as
the middle of a whole A11. This corner is finite, so the unit coupon
pays. Reversing the ordered profiles is identical.

## The isolated/isolated boundary

The formally verified determinant identity is

\[
0=1-\frac{2}{3h^2}-\frac{2}{3k^2}
 +\frac13\|x/h-y/k\|_{G_Z^{-1}}^2.
\]

The positive native heights are1 or at least sqrt2. If both are at least
sqrt2, the base term is at least1/3. If one exceeds sqrt2, it is strictly
positive even when the other is1. If the heights are1 and sqrt2 in either
order, the base term is zero; disjoint-support equality again forces
x=y=0. Both adjacent corners have whole A1/A11 components, and one unit
coupon suffices.

The remaining case has h=k=1 and

\[
(x-y)^T G_Z^{-1}(x-y)=1.
\]

Both separate augmented diagrams Z+c,Z+d are strictly positive definite
principal edge diagrams. Consequently the complete energy-one theorem
applies, with no bound on the rank of Z.

If one support is empty, choose the adjacent corner containing that
spectator-free new root. It is a whole A1; the old cycle with the other
mark deleted is a whole A11, marked at its unique middle. The remaining
spectators form finite components. This corner is finite and its unit
coupon pays, regardless of which allowed one-sided energy-one profile
occurs on the excluded new root. In particular no bounded list of B/D
spectator ranks is needed for this argument.

If both supports are nonempty, the energy-one theorem gives unit contacts
to two distinct whole A1 components, to the two endpoints of a whole A3,
or to the two short leaves of a whole D_r with4<=r<=7. In the same-component
cases there is a c–spectator–spectator–spectator–d path of length4.
At the opposite actual corner its connected component is finite or affine;
all other components are positive-definite principal parts of the old
white diagram and Z. It therefore receives the retained finite8 or
ideal5 unit. The source illustrates A3→A5, D4→E6 and D5→affine E6, and
detects a negative direction for the proposed D6/D7 opposite corners.
The short-path argument already suffices without this illustrative split.

For two distinct whole A1 spectator components, each adjacent corner has
a whole A2 containing its new black root, against the whole A11 middle.
The two half-coupons sum to one. The energy classification seals the two
spectator components, and all other spectators are orthogonal to the core;
there is no unaccounted finite extension that changes these coupon types.

## Whole-component and verification scope

Whenever a coupon is used, spectators are right to all of F and the
selected new root has either no spectator contact or exactly its sealed
whole-A1 contact. Hence the displayed A11 and A1/A2 components are whole
components of the selected corner. Every other component is finite, so
the entire corner is finite. In A11's path order the marked node is index5
(zero-based), its unique middle. The source checks this role directly.

The checker constructs15 literal cases: two native-height equality cases,
two mixed-profile equalities, the two-A1 split, five one-sided examples,
and five same-component twin-leaf examples. Each white-plus-one-black
edge Gram is positive definite. Every admissible full local star has one
negative and one zero direction, its old corner is affine, and the saved
corner coupons or length-four path supply a unit. The two impossible
twin-leaf extensions are explicitly rejected by their corner inertia.

Generation and a fresh separate-process saved replay both terminated
exit0/PASS on the six single-profile checks,16 rectangle values, formal
determinant,15 literal cases and60/6 cycle-pair partition. Direct source
bindings, Python AST and local report links passed.

Root read the full source and proof, including the inherited value and
profile formulas: unique mixed equality, both-end positivity, the fixed
leaf/degree-two exclusion of affine completion with arbitrary spectators,
native-height reduction before energy one, disjoint-support equality,
untouched-side whole A1/A11 sealing, the split A2/A11 halves and the
same-component path-four alternative including nonsimple ideal corners.
The review also checked the rational congruences and unique-middle coupon
role, and found no dependence on an A17 endpoint payment. It did not
duplicate the15 literal inertia or16 rectangle calculations. No new
ambient bound, physical chamber, outside refereeing or proof-assistant
formalization is asserted.

```sh
python3 scripts/fresh_cycle12_ideal5_transfer_20261017.py
python3 scripts/fresh_cycle12_ideal5_transfer_20261017.py --verify
```
