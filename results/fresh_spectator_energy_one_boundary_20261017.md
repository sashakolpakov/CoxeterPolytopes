# The strict spectator boundary at energy one

Let G be a positive-definite norm-two acute Coxeter Gram matrix, and let x,y
be nonnegative native contact vectors with disjoint supports. Both separate
augmented Grams G+x and G+y are required to be strictly positive definite.
The equality

    (x-y)^T G^(-1)(x-y)=1

has exactly the following possibilities, up to exchanging x and y and
diagram automorphism. Components not mentioned have no contacts.

* If only one support is nonempty: one sqrt2 contact on a whole A1; two
  unit contacts on two distinct whole A1 components; one unit contact at
  the middle of a whole A3; one unit contact at the remote leaf of a whole
  B_r, r>=2; or one unit contact at the long leaf of a whole D_r, r>=4.
  In B2 either endpoint is allowed, and in D4 any leaf is allowed.
* If both supports are nonempty: one unit contact from each side, on two
  distinct whole A1 components, on the two endpoints of a whole A3, or on
  the two short fork leaves of a whole D_r with4<=r<=7.

Together with the [strict spectrum below one](fresh_spectator_energy_below_one_20261017.md),
this classifies the complete interval0<=D<=1. The [source](../scripts/fresh_spectator_energy_one_boundary_20261017.py)
and [certificate](fresh_spectator_energy_one_boundary_20261017.json) check exact
identities. This local statement does not by itself prove a new corner
transfer or lower an ambient dimension bound.

## One nonempty support

As in the preceding proof, inverse entries are nonnegative and inverse
diagonals are at least1/2. Thus at most two contacts are possible. If there
are two, both are unit contacts, both inverse diagonals equal1/2 and their
inverse cross entry vanishes. An inverse diagonal equals1/2 only for an
isolated root: any incident bond gives a strictly positive Schur correction.
This gives the two whole A1 components.

For one contact of magnitude p, the bound p^2/2<=1 and the native alphabet
give p=1 or sqrt2. The latter forces a whole A1. For p=1 the contacted
inverse diagonal equals1. On A_r the formula i(r+1-i)/(r+1)=1 is equivalent
to (i-1)(r-i)=1, so only the middle of A3 occurs.

For a component that is not a simple path, take the principal B or D
obstruction from the preceding proof, beginning at the marked vertex.
The remote-leaf inverse diagonal is1, and its entire inverse column is
strictly positive. Adding any connected vertex therefore increases that
diagonal strictly: the Schur correction contains the nonzero product of
this positive column with a nonnegative nonzero contact vector. Increasing
any bond likewise increases it strictly, by the inverse derivative formula
and positivity of that column. Consequently the principal obstruction must
already be the whole component, with exactly its minimal bond labels.
This gives B_r or D_r at their stated leaves. If the marked vertex is itself
a branch, its principal A3 middle already has inverse diagonal1, and the
third neighbor makes it strictly greater, excluding this case.

All these one-sided cases have projection norm1<2, so their separate edge
extensions are strictly positive definite.

## Two nonempty supports

Suppose first that an x-support node i and a y-support node j are
nonadjacent. Gram-metric Cauchy gives

    1 >= (x_i+y_j)^2/4 >= 1.

Equality forces x_i=y_j=1 and

    G^(-1)(x-y)=(e_i-e_j)/2.

Hence every other coordinate of x-y is (G_ki-G_kj)/2. Since G is positive
definite, all off-diagonal entries are in(-2,0]. Their half-differences are
strictly between-1 and1. A nonzero native contact has magnitude at least1,
so every other coordinate vanishes. Disjoint supports now imply x=e_i,
y=e_j, and i,j have identical neighbor sets with identical bonds.

If there are no neighbors, they are two whole A1 components. Otherwise the
forest condition permits exactly one common neighbor k. The two equal bonds
must be simple: a three-node principal Gram with equal magnitude t has
determinant8-4t^2>0, hence t<sqrt2. Thus i,j are simple twin leaves.

The remaining connected component is A3 or D_r. Here is an elementary
all-rank argument. The common neighbor has at most one further neighbor,
since four simple leaves would contain an affine D4 subdiagram. Follow the
remaining simple path. A first doubled bond or second fork would contain an
affine B or D principal subdiagram, contradicting positive definiteness.
Explicitly, give the first two fork leaves coefficient1/2 and every stem
vertex coefficient1. At a terminal doubled bond give the terminal vertex
coefficient1/sqrt2; at a terminal fork give its two leaves coefficient1/2.
The resulting minimal diagram has positive null vector. Larger bond
magnitudes make its quadratic form negative. Hence the remaining path ends
at a leaf without a nonsimple bond or second branch.

For D_r with short leaves i,j, including D3=A3, the inverse column at i has
stem coordinates1/2,2/2,...,(r-2)/2, short-leaf coordinate r/4 at i and
(r-2)/4 at j. Multiplying by G verifies this in every rank. Thus the two
separate projection norms are r/4. Strict edge positivity requires r<8,
giving precisely A3 and D4,...,D7. Their inverse difference is
(e_i-e_j)/2, so the energy is1.

It remains to exclude the case in which every cross-support pair is
adjacent. Either support having two nodes would make a cycle in its
separate augmented Gram. Thus x=p*e_i,y=q*e_j for adjacent i,j. The
preceding below-one proof actually shows that any nonunit p or q gives
energy strictly greater than1. For p=q=1, a simple marked bond forces the
whole component to be A_r by the strict marked-edge obstruction, giving
energy r/(r+1)<1. A nonsimple marked bond forces a whole I2(4) or I2(5),
also with energy strictly below1. No adjacent case remains.

## Verification scope

The written identities and sealing arguments are in arbitrary rank. The
checker corroborates125 one-contact B/D identities through rank64,
62 twin-leaf identities through rank64, and the three isolated cases.
It verifies positive pivots, entire inverse contact columns, exact energy
and the strict rank7/rank8 boundary. The sqrt2/A1 case uses a rational
congruence scaling its isolated root by sqrt2.

Generation and a fresh separate-process saved replay both terminate
exit0/PASS. A separate full source/mathematical read passes the Cauchy
equality, native-gap argument, all-rank twin-leaf classification, strict
rank7/rank8 boundary, whole-component sealing and rational congruences.
That review does not duplicate the125/62 finite identity runs. Direct
bindings, syntax and local report links pass. Source, data and report are
final/frozen. No independent arithmetic reconstruction, outside refereeing
or proof-assistant formalization is claimed.

    python3 scripts/fresh_spectator_energy_one_boundary_20261017.py
    python3 scripts/fresh_spectator_energy_one_boundary_20261017.py --verify
