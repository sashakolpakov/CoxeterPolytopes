# Radius-four transfers at classical affine tree corners

An actual quadrilateral with a classical affine B/C/D corner has the following
local replacement. Suppose the marked roots are at distance P in {5,6,7,8,9},
and both are farther than four steps from every feature of their varying
affine factor. A feature is a vertex of degree at least three or an endpoint
of a nonsimple bond. Then another corner pays by one of these rules:

* a finite marked path of length at most eight;
* an ideal marked path of length at most four;
* an ideal marked path of length five through nine with a marked root within
  four steps of a feature of that corner's varying affine factor.

The proof allows arbitrary supporting spectators and nonsimple cusps. It
uses no new finite coupon. This is a local tree theorem: the cycle case,
triangle payments, preferred small-factor obligations and a global dimension
comparison are separate. No new ambient bound follows from this package alone.

The [checker](../scripts/fresh_radius4_affine_tree_transfer_20261017.py) and
[certificate](fresh_radius4_affine_tree_transfer_20261017.json) record all
405 profile/distance cases, the finite E8 exception and the energy-one
path classification used below. Verification status appears at the end.

## Actual-star setup and the ninth profile

Use the actual-star and factor-fiber setup of the completed
[radius-five theorem](fresh_finite8_ideal5_local_replacement_20261017.md).
The chosen corner is F orthogonal to a positive-definite spectator diagram
Z, with the two old marks a,b in the same varying affine factor F. The two
opposite roots are c,d. The white diagram and each actual edge diagram
W+c and W+d are strictly positive definite. All adjacent black contacts
may be assumed zero, since otherwise another corner pays immediately.

The feature-distance condition puts a,b on the simple central spine.
Deleting either mark leaves an outer B component of rank at least six,
or an outer D component of rank at least seven. The other component has
rank at least eleven. Thus the complete eight one-root profiles in the
[weak-profile theorem](fresh_distance8_ideal_corner_replacement_20261017.md)
remain exhaustive except for an extension of an outer D7.

Here is the classification boundary explicitly. An outer B6 or larger
contains a doubled bond, so it cannot be a subdiagram of a simply-laced
exceptional diagram; finite F4 and affine F4 are too small. For an outer
D7, a finite exceptional extension by one new root can only be E8,
obtained by prolonging either short fork leaf with a simple bond. The
finite E8 arms have lengths 1,2,4, and deleting the outer vertex of the
length-two arm leaves the D7 arms 1,1,4. No affine exceptional diagram
of order eight has a D7 deletion: the affine E7 arms are 1,3,3. The outer
D8 possibility is the already treated affine E8 extension, contained
entirely in W+c and therefore forbidden by strict edge positivity.
Larger D ranks, the inner component, and an extension joining both old
components are too large for an exceptional addition. Noncrystallographic
finite diagrams have rank at most four, apart from the rank-two family.
The classical possibilities are the inherited eight profiles.

The extra profile has especially small cost. Number D7 by the chain
0--1--2--3--4 with short leaves 5 and 6 at 4. Its inverse column at 5 is

    (1/2, 1, 3/2, 2, 5/2, 7/4, 5/4).

The checker verifies the entire column by multiplication. Thus the
projection norm of the new simple contact is 7/4, and its pairing with
the long port 0 is 1/2. In the inherited Schur notation this gives

    l=0, c0=1/4, h=x+1/2 >= 1/2,
    g=2*l/h+c0/h^2 <= 1.

Here x is the nonnegative native opposite contact magnitude. The E8
extension is positive definite; its exact pivots and its three arm
lengths are recorded. This ninth profile cannot simply be discarded,
but its bound is better than the old maximum four.

## The complete nine-profile comparison

For each side the inherited profiles have scalar cost
g=2*l/h+c0/h^2. Their maxima, followed by the new profile, are:

|Profile|Maximum g|
|---|---:|
|isolated|2|
|outer endpoint3|1|
|outer endpoint4|0|
|outer second-node3|0|
|inner endpoint3|3|
|inner second-node3|4|
|both endpoints3|1|
|inner endpoint4|2|
|outer D7 short-leaf E8|1|

The actual normalized determinant, after eliminating Z, is

    1-(g_c+g_d)/P
      +(1/P)*(x/h-y/k)^T G_Z^(-1)*(x/h-y/k),

where x,y now denote the two spectator contact vectors. This is the
inherited formal identity, with the newly checked outer-only profile
having l=0. The last term is nonnegative. The actual star has corank one,
so the determinant must vanish. A common nonzero spectator coordinate
already supplies a paid path c--e--d of length two. Hence throughout the
remaining analysis the spectator supports are disjoint.

At P=5 through 8, an inner second-node3 profile gives an adjacent marked
path of length P-1. At P=5 this has ideal distance four. For the other
values the old attachment vertex retains three incident edges in the
actual new corner, so it remains a feature one step from the new mark.
If the corner is finite, the path has length at most seven and pays.
This accounts for 68 immediate cases. At P=9 even the two maximum costs
give a strictly positive gap 1/9.

After those cases are removed, 331 profile pairs have strictly positive
gaps. Five pairs have zero minimum, and one has negative minimum.
All of these are treated explicitly rather than inferred impossible
from a relaxed numerical bound.

## The five zero cases

At P=6 the only zero case is inner endpoint3 on both sides. Equality
forces both heights to be one and the spectator correction to vanish.
Disjoint supports then force x=y=0. The middle A5 with c,d at its two
ends is a whole finite A7 component, with marked distance six. This is
the sealed equality already proved in the radius-five theorem.

At P=5 the four zero cases are inner endpoint3 paired, in either order,
with isolated or inner endpoint4. Each scalar cost must attain its
maximum, and the spectator correction again vanishes. Positive
definiteness and disjoint supports force both spectator vectors to be
zero, even when the two heights differ. Take the adjacent corner on the
inner endpoint3 side. Its component is a finite B or D diagram obtained
by extending the inner long endpoint by a simple bond. The other outer
component and all spectators are orthogonal positive-definite blocks.
Thus the whole corner is finite, and its marked path has length five.

## The remaining energy-one case

Only P=5 with inner endpoint3 on both sides has negative relaxed minimum.
For this profile h=x+1, so the native alphabet gives h=1 or h>=2; the
same holds for k. If either height is at least two, that side's cost is
at most 5/4, and the normalized determinant is at least

    1-(5/4+3)/5 = 3/20 > 0.

Consequently h=k=1 before any unscaled contact classification is applied.
The zero determinant now requires precisely

    (x-y)^T G_Z^(-1)*(x-y)=1.

Both G_Z+x and G_Z+y are strictly positive definite as principal
subdiagrams of the actual edges. The complete
[energy-one classification](fresh_spectator_energy_one_boundary_20261017.md)
therefore applies. In the opposite corner, c,d are also joined by the
simple path through the old middle A4, of length five. These middle
vertices have no contact with Z. The outer B/D blocks are separate.

More generally, join c,d by a simple path of length L>=1, and give them
exactly the classified energy-one contacts with Z. Up to exchanging
the endpoints, the resulting marked component is as follows:

|Spectator contacts|Resulting component|
|---|---|
|one sqrt(2) contact on whole A1|finite B_(L+2)|
|two unit contacts on separate A1, on one side|finite D_(L+3)|
|one unit contact at the middle of whole A3|finite D_(L+4)|
|one unit contact at the long leaf of whole B_r|finite B_(r+L+1)|
|one unit contact at the long leaf of whole D_r|finite D_(r+L+1)|
|opposite unit contacts on separate A1 components|finite A_(L+3)|
|opposite unit contacts at whole A3 endpoints|affine cycle of order L+4|
|opposite unit contacts on D4 through D7 short leaves|impossible for a PSD corner|

These statements follow by literally concatenating the path with the
whole components from the classification, in arbitrary rank and length.
The last row contains the same simple cycle as the preceding row with
an attached nonorthogonal tail. The vector equal to one on the cycle
and zero elsewhere has zero quadratic form, but its Gram product is
minus one at the first tail vertex. A positive-semidefinite matrix
cannot have a zero quadratic vector outside its kernel, excluding it.
This also excludes any additional completion of that configuration.

At L=5 every finite outcome pays at marked distance five. The sole
affine outcome is the whole nine-cycle; its alternate path through
the A3 spectator has length four, so it pays by ideal distance four.
All uncontacted spectator components remain orthogonal and positive
definite. This proves the remaining case without assuming simple cusps
or bounding spectator rank.

## Verification scope

Primary generation passes all 405 tree cases, the full D7 inverse column
and E8 pivots, and 51 representative path concatenations and positive
semidefinite obstructions. The latter corroborate the explicit all-rank
diagram proof; a tested-rank cutoff is not its justification. Doubled
bonds use a rational congruence scaling of their terminal root, which
preserves positive definiteness and changes its diagonal from two to four.

Generation and a fresh separate-process saved replay both terminated
exit0/PASS. A separate complete source and mathematical review passed the
outer-rank classification boundary, ninth profile and cost, persistent
feature payments, zero-correction sealing, height-first energy-one case,
and all-rank finite/cycle/twin-tail concatenations. That review also read
the three inherited local proof reports. It did not repeat the 405-case
arithmetic or 51 representative diagram checks. Direct source bindings,
Python syntax and all local report links pass. Source, data and report
are final/frozen.

The inherited finite/affine classification, factor-fiber theorem,
formal determinant and spectator classification retain their existing
validation scope. This package does not re-enumerate those inputs and
does not establish a new global dimension theorem.

```
python3 scripts/fresh_radius4_affine_tree_transfer_20261017.py
python3 scripts/fresh_radius4_affine_tree_transfer_20261017.py --verify
```
