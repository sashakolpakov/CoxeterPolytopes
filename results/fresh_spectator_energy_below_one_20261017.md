# The complete strict spectator spectrum below one

The new all-rank assertion is as follows. Let G be a positive-definite
norm-two acute Coxeter Gram matrix, and let x,y be nonnegative native contact
vectors with disjoint supports. Require both separate augmented Grams
G+x and G+y to be strictly positive definite. For

    D=(x-y)^T G^(-1)(x-y)<1,

the complete set of possible values is

    0; r/(r+1) for integers r>=1; 2-sqrt2; 1-sqrt5/5.

At the rational nonzero value r/(r+1), the contacts lie on one whole simple
A_r component: either one unit contact at a leaf, or two unit contacts at
adjacent vertices. Both separate extensions must still be positive definite;
this excludes some interior adjacent positions in large A_r. Every other
component has zero contacts. The two irrational values have two unit contacts
at the endpoints of a whole I2(4) or I2(5), respectively. Conversely all the
listed energies occur: the one-contact A_r leaf realizes every rational value.

The [source](../scripts/fresh_spectator_energy_below_one_20261017.py) and
[certificate](fresh_spectator_energy_below_one_20261017.json) check the exact
identities below. This strengthens the [spectrum through three quarters](fresh_small_spectator_difference_gap_20261017.md).
It is a local lemma, without a new ambient dimension or construction claim.

## Two elementary obstructions in arbitrary rank

All positive-definite acute Coxeter diagrams are forests. A cycle's all-ones
vector would otherwise have nonpositive square. Their inverses are entrywise
nonnegative, and an inverse diagonal cannot decrease when a principal diagram
is enlarged or its bond magnitudes are increased. The principal comparison
is a positive Schur correction; bond monotonicity follows from the
nonnegative Neumann series, or differentiation of the inverse along the
positive-definite line segment.

Two explicit principal subtrees detect any branching or nonsimple bond.
On a simple path ending in a bond of magnitude sqrt2 (type B), the inverse
diagonal at its remote leaf is1. At the next vertex it is2 if that next
vertex precedes the doubled bond. On a simple path ending in a fork with
two short simple leaves (type D), the remote leaf and next vertex again
have inverse diagonals1 and2. This includes D4 with the fork center as the
second vertex.

Here is an all-rank calculation. For B, scale the terminal root by sqrt2.
The rational congruent Gram has last diagonal4 and last bond magnitude2.
Eliminating the last root makes the preceding diagonal1; eliminating
successive roots repeats the diagonal1. The first inverse column is
(1,...,1,1/2). The second is(1,2,...,2,1) when the path has at least3 nodes.
The first two diagonal entries are unaffected by the terminal congruence.
For D, eliminating the two short leaves changes the fork diagonal to1,
and the same recurrence applies. The first inverse column is
(1,...,1,1/2,1/2); the second is(1,2,...,2,1,1). These vector identities
and the positive elimination pivots prove the assertions for all ranks.

From any vertex of a tree that is not a simple A path, follow a path to
the first nonsimple bond or branching node. Keeping the terminal doubled
bond, or two outgoing fork leaves, gives one of these subtrees; larger
native bonds only increase inverse diagonals. Hence every inverse diagonal
of a component that is not a simple path is at least1.
If the starting vertex itself is a branch, retain two of its neighbors
instead: the resulting simple A3 has inverse diagonal1 at that middle
vertex, and larger bonds only increase it. Thus this zero-length path case
does not require an incoming edge for the fork construction.

There is a stronger statement for an adjacent simple marked edge i-j.
If its component has a nonsimple bond or branch anywhere, it lies in one
of the two sides of that edge. Retain a path from the opposite marked
endpoint to that feature, together with its doubled bond or fork leaves.
The opposite endpoint becomes the remote leaf and the other marked
endpoint becomes its next vertex. One marked inverse diagonal is therefore
at least2. This contradicts strict positive definiteness of that marked
root's separate unit-contact extension.

## One nonempty contact support

Nonnegative inverse entries and diagonal lower bound1/2 give
D>=sum_i x_i^2/2 when y=0. Under D<1, the support has only one coordinate,
with native magnitude1. A nonsimple magnitude is at least sqrt2 and would
already give D>=1.

The subtree obstruction forces the whole contacted component to be a
simple A_r path. Its inverse is, in one-based indices,

    (G_A_r^(-1))ij=min(i,j)*(r+1-max(i,j))/(r+1).

Every interior diagonal is at least1. Thus the contact is at a leaf and
D=r/(r+1). The same argument applies after interchanging x and y. If both
supports vanish then D=0; conversely positive definiteness and disjoint
supports force both to vanish when D=0.

## Two nonempty supports

For any i in the x-support and j in the y-support, Gram-metric Cauchy gives

    D >= (x_i+y_j)^2/(4-2G_ij).

Nonadjacent support nodes would give D>=1. Every opposite-support pair is
therefore adjacent. If one support had two vertices, adjoining its new
root would create a cycle through a vertex of the other support. Its
strictly positive-definite augmented diagram forbids this. Hence
x=p*e_i and y=q*e_j for adjacent i,j, with p,q native magnitudes.

Put t=-G_ij. If p>=sqrt2, the separate edge condition and the two-node
principal inverse give

    p^2*2/(4-t^2)<2,

so t^2<2 and the native alphabet forces t=1. Both p and q are then below
sqrt3<2. On the two-node principal Gram the difference energy is

    2*(p^2+q^2-p*q)/3.

For fixed p<2 this is increasing in q>=1. At q=1 it is increasing for
p>=sqrt2. Its minimum is therefore2-(2/3)sqrt2>1, contradicting D<1.
Interchanging the contacts gives the same contradiction when q>=sqrt2.
Thus p=q=1.

If t is nonsimple, the separate strict edge tests imply t<sqrt3.
Any further neighbor of either endpoint with magnitude u>=1 would give
inverse diagonal2/(4-t^2-u^2)>=2 in that three-node principal path,
again impossible. The whole component is consequently I2(m), with
sqrt2<=t<sqrt3 forcing m=4 or5. Its energy2/(2+t) is exactly one of
the two listed irrational values.

If t=1, the marked-edge subtree obstruction excludes every nonsimple
bond and branch in the whole component. It is A_r. For adjacent contacts
at k,k+1, its inverse difference vector has entries

    j/(r+1) for j<=k;
    -(r+1-j)/(r+1) for j>=k+1.

Every entry is nonzero and the difference energy is r/(r+1), independent
of k. The remaining strict extension conditions are the two explicit
diagonal inequalities(G^(-1))kk<2 and(G^(-1))k+1,k+1<2. This completes the
all-rank classification; no cutoff on spectator rank or native labels
is imposed.

## Consequences and verification

For a rational0<D<1, the forced rank is exactly D/(1-D). Thus4/7,
11/14,9/16 and8/17 are impossible;2/3,3/4 and4/5 force whole A2,A3,A4
contact components. The11/14 exclusion is a new input available for
shorter cycle transfers. Such transfers require separate native-height
and finite-coupon arguments and are not asserted here.

The checker verifies all A inverse columns through rank64, all2,016
adjacent-contact identities, and124 B/D obstruction templates with positive
elimination pivots and exact inverse columns. The formulas and subtree
arguments above prove arbitrary rank; these finite checks corroborate
them. The two irrational literal cases retain their frozen predecessor's
exact verification. Generation and a fresh separate-process saved replay
both terminated exit0/PASS. An independent full source/mathematical read
passes the arbitrary-rank obstruction, support reduction, strict nonunit
bound and exceptional rank-two cases. It identified the zero-length path
exposition case above, now made explicit. That review did not independently
rerun the finite identity checks. Source, data and report are final/frozen
after direct binding, syntax and local-link checks. All reviews are internal;
external refereeing and proof-assistant formalization are not claimed.
