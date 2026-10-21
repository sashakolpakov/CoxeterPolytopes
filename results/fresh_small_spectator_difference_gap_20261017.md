# The discrete spectator-energy spectrum through three quarters

**FINAL/FROZEN: exact generation, fresh saved replay and a separate full
mathematical/source review PASS.**

Let G be a positive-definite norm-two acute Coxeter Gram matrix. Let x,y be
nonnegative vectors of native contact magnitudes with disjoint supports.
Assume both separate augmented Grams

    [[G,-x],[-x^T,2]],  [[G,-y],[-y^T,2]]

are **positive definite**. This is the strict actual-edge hypothesis, not
merely semidefiniteness at a corner. Put D=(x-y)^T G^(-1)(x-y). Then

    D<=3/4 implies
    D in {0,1/2,2/3,3/4,2-sqrt(2),1-sqrt(5)/5}.

In particular **D cannot equal4/7**. Equality D=3/4 has only two forms,
up to swapping x,y and diagram symmetries: one unit contact at an endpoint
of a whole simple A3 component, or one unit contact from each side at
adjacent vertices of a whole simple A3 component. Every other spectator
contact vanishes.

The [source](../scripts/fresh_small_spectator_difference_gap_20261017.py) and
[certificate](fresh_small_spectator_difference_gap_20261017.json) verify the
seven nonzero literal patterns, their separate strict edge extensions,
Schur identities and exact values. Completeness follows from the argument
below, which permits arbitrary spectator rank and native finite labels.
No ambient dimension or complete face-payment theorem is asserted.

## Preliminary positive-definite facts

Every nonzero native magnitude is at least1; if greater than1 it is at
least sqrt2. A positive-definite acute Coxeter graph is a forest: on any
cycle, the vector with entries1 on that cycle has square at most0, since
each cycle bond has magnitude at least1 and extra bonds only decrease
the square. This applies to G and to each of the two strict augmented
Grams, even if their graph is disconnected.

The inverse of G is entrywise nonnegative and has diagonal at least1/2.
One proof writes G=2I-A with A nonnegative and expands the convergent
Neumann series. Restricting to a principal diagram cannot increase the
corresponding inverse quadratic energy: the full inverse adds a
positive-semidefinite Schur correction. We also use the inherited
[rank bound](fresh_selective_distance9_dimension108_20261017.md): a connected
positive-definite norm-two Coxeter component of rank r has every inverse
diagonal at least r/(r+1). That argument covers all finite labels and
indirect tails, not a finite diagram sample.

## One nonempty support

If y=0, nonnegative inverse entries give D>=sum_i x_i^2/2. Thus D<=3/4
allows only one nonzero coordinate, and its native magnitude must be1.
The rank bound forces its whole component to have rank at most3.

Rank1 gives A1 and D=1/2. At rank2, the inverse diagonal is
2/(4-t^2), where t>=1 is its bond magnitude. Being at most3/4 forces
t^2<=4/3, hence t=1 and D=2/3. At rank3 the graph is a path. Its middle
inverse diagonal is at least1. At a leaf, the minimum is3/4, attained
only when both bonds are simple: increasing either bond strictly increases
that diagonal. This is the whole simple A3 endpoint case. Swapping the
two supports gives the same list. If both supports vanish, D=0.

## Two nonempty supports must be singletons

For any i in the x-support and j in the y-support, Cauchy in the Gram
metric gives

    D >= (x_i+y_j)^2/(4-2G_ij).

If i and j are not adjacent, this is at least1. Hence under D<=3/4,
every vertex of one support is adjacent to every vertex of the other.
If one support had two vertices i1,i2 and the other contained j, adjoining
the root with contacts to i1,i2 would create the cycle
i1-j-i2-new-i1 in its strict positive-definite edge diagram. This is
impossible. Thus x=p*e_i and y=q*e_j for adjacent i,j.

Put t=-G_ij. Positive definiteness of each separate extension gives
(G^(-1))_ii<2/p^2<=2 and (G^(-1))_jj<2/q^2<=2. The two-node principal
inverse is2/(4-t^2), so t<sqrt3. If either p or q exceeded1, the native
gap and the Cauchy bound would give

    D >= (1+sqrt2)^2/(4+2t) > 96/125 > 3/4.

The displayed rational lower bound uses sqrt2>7/5 and sqrt3<7/4.
Consequently p=q=1.

## A nonsimple spectator bond is isolated

Suppose t>=sqrt2. If either endpoint of i-j had another neighbor with
bond magnitude u>=1, its inverse diagonal in that three-node principal
path would be

    2/(4-t^2-u^2) >= 2/(3-t^2) >= 2,

contradicting the strict edge-extension inequality. The denominator is
positive because that principal path is positive definite. Therefore i-j
is a whole rank-two component. The native conditions sqrt2<=t<sqrt3
leave exactly labels4 and5. For the two unit contacts,

    D=2/(2+t),

giving2-sqrt2 and1-sqrt5/5 respectively. This excludes arbitrary high
labels by a proved inequality, not an imposed cutoff.

## A simple spectator bond has at most one additional root

Now t=1. If i-j is the whole component, it is A2 and D=2/3. Otherwise
choose a third node adjoining one of its endpoints, with bond magnitude
u>=1. It cannot adjoin both endpoints because G is a forest. On A2,
G_A2^(-1)(e_i-e_j)=(1/3,-1/3). Schur complementation therefore gives
the three-node principal energy

    2/3 + u^2/[6(3-u^2)],
    energy-3/4=(u^2-1)/[4(3-u^2)].

Positive definiteness ensures u^2<3. Hence the energy is at least3/4,
with equality only at u=1. The principal diagram is then simple A3,
and the two contacts occupy adjacent vertices. In one orientation its
inverse difference vector is

    G_A3^(-1)(e_0-e_1)=(1/4,-1/2,-1/4),

which has no zero coordinate. Any additional root attached to this A3
can meet only one of its nodes, since a second attachment would create
a cycle. Its product with the displayed inverse difference is nonzero.
The positive Schur correction is then strictly positive. Thus D=3/4
forces this A3 to be the whole connected component, including all tails.
Every other component has zero contact from both sides.

This proves the stated list and the sealed A3 equality classification.
The exact values2-sqrt2 and1-sqrt5/5 lie on opposite sides of4/7;
none of the other listed values equals4/7. The latter exclusion is useful
for the cycle14 resistance calculation, while the A3 equality description
is a possible input to a smaller affine-D transfer. Those face arguments
remain separate obligations.

## Verification scope

The literal matrices and separate edge extensions are checked exactly over
the relevant quadratic fields. The source also verifies the three-node
Schur identity, its positive derivative, the nonvanishing A3 inverse
difference, the rational nonunit-contact bound and the4/7 exclusion.
Generation and a separate fresh saved replay both terminated with exit0/PASS,
checking all seven literal patterns. A second agent read the complete proof
and source and found no gap in the one-support rank reduction, strict-edge
singleton argument, native-label exclusions or whole-A3 sealing. That review
does not claim a second implementation of the exact arithmetic. The root
also read the final proof and source. No external refereeing or formalization
is claimed.
