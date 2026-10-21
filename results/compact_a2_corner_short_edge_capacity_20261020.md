# Short directions at a whole A2^r finite corner

20 October 2026. Paper proof only: no execution, enumeration, optimizer,
build or formal verification. The argument applies to every compact
hyperbolic Coxeter polytope; arithmeticity is not used.

**Theorem.** Let p be an actual finite vertex whose whole diagram is
A2^r, and put n=2r. Write L_i for its n incident edge lengths. Then

    sum_i tanh(L_i/2) > (n-2)*173/500 + 11/25.             (1)

Consequently this vertex cannot satisfy the common-white drift inequality
A_f(p)<||rho_f(p)|| in dimension n>=14, for ANY single positive function
f on common-white diagrams. This excludes one complete vertex family
from the global drift witness, not all compact polytopes in those
dimensions. Mixed A1/A2 forests remain outside the theorem.

## 1. Coordinates and the threshold

Use p^2=-1 and outward unit old normals. Each A2 block has Gram matrix

    B=[[1,-1/2],[-1/2,1]],
    B^-1=[[4/3,2/3],[2/3,4/3]].

For an actual facet c not through p, set a=-p.c>0 and let beta_j=-c.n_j
be its nonnegative old contacts. If s_c is its spatial projection, then

    c=a p+s_c,       s_c^2=beta^T B_old^-1 beta=1+a^2.

Every nonzero contact is cos(pi/m), m>=3, or is greater than1. Contact1
is impossible in a compact chamber: the associated positive null sum
would be a feasible ideal point by all-wall acuteness.

Set

    sinh^2 L0=31/50,       tau0=tanh(L0/2)>173/500.        (2)

The last inequality follows by substituting 173/500 in
sinh^2 L=4tau^2/(1-tau^2)^2 and cross-multiplying positive integers.
If c closes an edge shorter than L0, eliminating its n-1 old white
normals gives

    delta=1-beta_W^T B_W^-1 beta_W>0,
    a^2=delta*sinh^2 L<31/50.                            (3)

Here delta<=1. For two DISTINCT walls c,d closing such short directions,

    c.d=beta^T B_old^-1 gamma-a b >-31/50.               (4)

In particular their pair is either orthogonal or has label3; label4,
larger finite labels and disjointness are excluded by (4). The more
restrictive orthogonality statements below use their actual heights.

We will use the elementary exact bounds, with q=cos(2pi/7),

    623/1000<q<5/8,       q>31/50,
    4sqrt2/9>31/50.                                    (5)

For the first bounds q lies in (1/2,1), satisfies
8q^3+4q^2-4q-1=0, and this polynomial is strictly increasing on
[1/2,1]. Its values at623/1000 and5/8 have signs negative and positive.
The other comparisons are immediate rational squarings.

## 2. Complete short-wall inventory

First a short closing wall cannot contact BOTH nodes of an A2 block.
If both contacts are3, that block contributes1 to the full spatial
square. With no other contact a=0; otherwise another block contributes
at least1/3. Releasing a node of the double-contact block then gives

    sinh^2 L >= (1/3)/(3/4-1/3)=4/5.

Releasing another block leaves the entire cost1 in the white support,
so cannot give a finite endpoint. If the two contacts are not both3,
their componentwise smallest possibility is (3,4). Even releasing the
larger contact gives

    a^2>=sqrt2/3,       delta<=3/4,
    sinh^2 L>=4sqrt2/9>31/50.

Larger contacts or further blocks only increase this ratio. Thus each
contacted A2 block supplies exactly one nonzero contact.

Write its full contribution as q_j=(4/3)beta_j^2>=1/3. A release in a
contacted block has delta=1-sum_(other blocks)q_j if it releases the
contacted node. Releasing its uncontacted mate leaves the additional
white cost beta_j^2, and therefore gives a strictly longer edge whenever
that endpoint is finite. A release in an untouched block is impossible,
since its white cost is already the full spatial square1+a^2>1.

There are at most three contacted blocks: four leave white cost>=1.
For three blocks, three3 contacts give a=0; the next native increment
already gives a^2>=1/3 and delta<=1/3, hence sinh^2 L>=1. Thus three
blocks also cannot occur here.

For two blocks let their contributions be q_i,q_j. Shortness of a
contacted-node release means

    q_i+(81/50)q_j<81/50,       q_i+q_j>1.              (6)

If both contacts are at least4, q_i,q_j>=2/3, and the corresponding
ratio is at least1. Otherwise one contact is3, with contribution1/3.
The other must exceed4: (3,3) and (3,4) have total contribution<=1.
It cannot be6, because its A2 block and c contain the affine3--6 path.
For label m>=7 the smaller possible edge has

    sinh^2 L = 2cos^2(pi/m)-1
             >=cos(2pi/7)>31/50.

The ONLY surviving two-block profile is therefore (3,5). Call it F.
Put phi=(1+sqrt5)/2. Its fixed squared height is

    h=a^2=(sqrt5-1)/6<5/24.

Only its5-contacted-node release is short: the other contacted-node
release has sinh^2 L=phi. The mate of its5 node has
delta=2/3-phi^2/4<1/12, giving sinh^2 L>12h>2, while the mate of
its3 node has delta=3/4-phi^2/3<0 and is not finite. Its short release
has sinh^2 L=(sqrt5-1)/4.

The remaining profile I has just ONE contact beta in ONE A2 block.
Its squared height and possible short directions are

    a^2=(4/3)beta^2-1,
    contacted node:       sinh^2 L=a^2,
    uncontacted mate:     sinh^2 L=a^2/(1-beta^2).       (7)

The mate is finite only for finite contact. Thus beta is a finite label
m>=7, or is dotted, and beta^2<243/200 by (3). Among finite labels only7
can give a second short direction: label8 already gives
sinh^2 L_mate=2sqrt2/3>31/50, and the ratio increases with beta^2.
The inequalities (5) show that the label7 mate is indeed below L0.
This completes the inventory, including zero contact at a released mate.

## 3. Pair compatibility and counting DIRECTIONS

Two F walls have height product h<1/2, so must be orthogonal. A common
block involving a5 contact contributes at least phi/6>h to their spatial
product, where phi=(1+sqrt5)/2. Two common blocks contribute at least1/3.
The only remaining shared-block possibility is two3 contacts at opposite
nodes, contributing1/6, which is not h. No common block contributes0,
also not h. Thus two distinct F walls cannot coexist.

An I wall and an F wall likewise must be orthogonal, since
a_I^2 h<(31/50)(5/24)<1/4. They must share the I block. Sharing F's
5-block, or its3 node, makes the spatial product larger than a_I*a_F.
Indeed the squared spatial product is then greater than1/6 or1/3,
respectively, while a_I^2 h<31/240<1/6.
The remaining possibility is the OPPOSITE node of F's3-block. Its
product is beta/3, so orthogonality requires

    h*((4/3)beta^2-1)=beta^2/9,
    beta^2=9h/(12h-1)>5/4>243/200.

Here h<5/24 and h>1/12 justify the displayed lower bound and positive
denominator. This contradicts (3). Hence I and F cannot coexist.

Two I walls in the SAME block cannot coexist. Even opposite nodes give
spatial product (2/3)beta*gamma, whereas

    a_I*a_J/(beta*gamma)
      =sqrt[(4/3-beta^-2)(4/3-gamma^-2)]
      <124/243<2/3.

For I walls in DIFFERENT blocks the spatial product is0, so (4) forces
c.d=-1/2, or a_I*a_J=1/2. Three such walls would form an affineA2
subdiagram, impossible in a compact chamber. Thus at most two I walls
occur. A finite I wall has a_I^2<1/3 and cannot coexist even with one
other I wall, because a_I*a_J<sqrt(31/150)<1/2. Therefore the two-wall
case consists solely of dotted I walls. Each closes only one finite
direction and has sinh^2 L>1/3.

There are consequently at most two short DIRECTIONS. Either they come
from one isolated7 wall, or from two dotted walls; the one-F and other
one-I cases supply at most one short direction. This argument does not
mistake a wall count for a direction count.

For the isolated7 wall write q=cos(2pi/7). Formula (7) gives

    sinh^2 L_contact=(2q-1)/3>41/500,
    sinh^2 L_mate=2(2q-1)/(3(1-q))>164/377.

Using tau=tanh(L/2), equivalently sinh^2 L=4tau^2/(1-tau^2)^2,
these imply tau_contact>7/50 and tau_mate>3/10. (The second rational
comparison is164*8281-3600*377=884>0.) Their sum is greater than11/25.
Each dotted short direction has tau>2-sqrt3>11/50, so the two-dotted
case has the same lower bound. Every one-wall/one-short-direction case
has tau>7/50. Since tau0>173/500>3/10, adding one further long direction
then also contributes more than11/25. With no short directions,
2tau0>11/25. Retaining all other n-2 long directions proves (1).

## 4. The common-white drift interface

Use the exact undirected identity in
`results/compact_finite_edge_white_weighted_drift_20261020.md`.
Every released white diagram here is A1+A2^(r-1), so ANY globally shared
function f takes one common value at these n releases. The inverse
diagonal is4/3, while each block contributes4 to the squared norm of the
sum of its two dual vectors. Therefore

    A_f=f*sqrt(4/3)*sum_i tanh(L_i/2),
    ||rho_f||=f*sqrt(2n).

For n14, the right side of (1), with tau0 replaced by its rational lower
bound, is574/125=4.592, whose square is greater than21=3n/2.
The ratio ((173n-126)/500)/sqrt(n) is increasing for n>=14. Hence
A_f>||rho_f|| for every possible whole A2^r vertex in dimension n>=14.

Only this local family is discharged. A dimension theorem still needs
all other finite diagrams, especially A1/A2 mixtures, and one choice of
shared white weights throughout the full compact polytope. No uniform
dimension ceiling is asserted.
