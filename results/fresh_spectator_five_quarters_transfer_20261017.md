# A short-path alternative at spectator energy five quarters

Let G be a positive-definite norm-two acute Coxeter Gram matrix. Let x,y
be nonnegative native contact vectors with disjoint nonempty supports,
and require both separately augmented diagrams G+x and G+y to be strictly
positive definite. If

    (x-y)^T G^(-1)(x-y) = 5/4,

then either an x-contact and a y-contact have distance at most3 in one
connected spectator component, or the only contacted components are
whole A1 and whole A3, with one unit leaf contact in each and opposite
contact sides. The first alternative gives a path of length at most5
between the two new roots. No bound on spectator rank or native labels
is imposed.

The [source](../scripts/fresh_spectator_five_quarters_transfer_20261017.py)
and [certificate](fresh_spectator_five_quarters_transfer_20261017.json)
check the exact identities below. The finite Coxeter classification and
the [one-sided spectrum below1](fresh_spectator_energy_below_one_20261017.md)
are inherited inputs. This is a local spectator lemma, not a complete
face-transfer or dimension theorem.

## Contacts in the same component

Every finite Coxeter diagram is a forest. Thus a new root can contact at
most one node in each connected component: two contacts would make a
cycle in its separately positive-definite augmented diagram. Write two
opposite contacts in one component as p at i and q at j, with p,q>=1 and
i different from j. Put J=G^(-1), a=Jii, b=Jjj, c=Jij. Strict extension
positivity forces a<2 and b<2, even before accounting for other components.
The contribution to the difference energy is

    E(p,q)=a*p^2+b*q^2-2*c*p*q.

When a-c and b-c are nonnegative, expand around p=q=1. The linear terms
are nonnegative and the remaining quadratic form is positive semidefinite,
so E(p,q)>=a+b-2c for every p,q>=1. This does not assume coordinatewise
monotonicity at arbitrary p,q.

For A_r, put N=r+1 and use one-based indices. The all-rank inverse is

    Jij=min(i,j)*(N-max(i,j))/N.

For i<j and d=j-i, the row sums are a-c=i*d/N and
b-c=(N-j)*d/N, both positive. The unit-contact energy is d*(N-d)/N.
Since N>=d+2, for d>=4 we have

    d*(N-d)/N >= 2d/(d+2) >= 4/3.

The two differences are respectively d^2*(N-d-2)/(N*(d+2)) and
2*(d-4)/(3*(d+2)). These identities prove the assertion in arbitrary rank.
The certificate separately multiplies the proposed inverse by each Gram
and checks every far pair in ranks5 through32 as corroborating examples.

For B_r number the simple spine0,...,r-2 and the terminal doubled-bond
node r-1. Its native norm-two inverse is

    Jij=min(i+1,j+1)                      (i,j<r-1),
    Ji,r-1=(i+1)/sqrt2                    (i<r-1),
    Jr-1,r-1=r/2.

Thus a possible unit contact lies only at0, or at r-1 when r<=3. Two
distinct allowed positions never have distance at least4. The formula
follows by direct multiplication with the tridiagonal Gram, including
the two terminal equations, for every r>=2.

For D_r number the long spine0,...,r-3 and the two short leaves r-2,r-1.
Its inverse has Jij=min(i+1,j+1) on the spine, Ji,leaf=(i+1)/2,
diagonal r/4 at each short leaf and off-diagonal (r-2)/4 between the
short leaves. Direct multiplication at the spine and fork proves these
identities for every r>=4. The only positions of diagonal below2 are0
and, when r<=7, the two short leaves. Short leaves are distance2 apart.
A long/short pair of distance at least4 therefore has r=6 or7; its
row sums are1/2 and r/4-1/2, and its unit energy is r/4>=3/2.
The certificate verifies both inverse families through rank24, without
using that sample as a rank cutoff.

For E6,E7,E8 the entire finite case is checked exactly. E6 has two
possible contact nodes, at distance4, with a=b=4/3 and c=2/3; their
energy is at least4/3. E7 has only one possible contact node and E8 has
none. The complete Grams and inverses are saved. F4,H3,H4 and every
I2(m) have diagram diameter at most3 and need no energy estimate.

Every connected component contributes nonnegative energy. Consequently
one mixed component with its contact nodes at distance at least4 would
already contribute at least4/3>5/4. Every mixed component instead yields
the stated short-path alternative.

## No component has contacts on both sides

Suppose the short-path alternative does not occur. Every active component
has contacts from only one side. Its energy is at least1/2: a nonzero
native contact is at least1 and an inverse diagonal is at least1/2.
Since both supports are nonempty, there are at least two active components;
three would contribute at least3/2. There are exactly two, one per side,
and each has energy between1/2 and3/4.

The inherited one-sided spectrum below1 forces a unit leaf contact on
a whole A_r, with energy r/(r+1). Here r<=3. The only two such energies
summing to5/4 are1/2 and3/4, from A1 and A3. Every other spectator
component has zero contacts. This proves the second alternative.

## Verification status

Generation and a fresh separate-process saved replay both terminated
exit0/PASS on all exact identities and finite exceptional cases.
The arbitrary-rank argument is the written inverse calculation and
classification above; finite samples are supplementary checks. A separate
internal full source and written-proof review passed the support reduction,
positive-quadratic remainder argument, all-rank formulas, strict contact
position restrictions and final component-energy split. That review did
not duplicate the finite arithmetic tables. Source, certificate and report
are final/frozen. No external
refereeing, formalization or new finite-volume example is claimed.
