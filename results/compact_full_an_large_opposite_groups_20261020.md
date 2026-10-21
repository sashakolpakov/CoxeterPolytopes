# Large opposite-wall groups at a whole A_n corner

20 October 2026. Own paper proof; independent review pending.
No computation, enumeration, build, new dimension ceiling or global
chamber realization is asserted.

## 1. A necessary classification with no rank-by-rank census

Let P be a compact hyperbolic Coxeter polytope of dimension n>=14,
and let its WHOLE diagram at an actual vertex p be the simple chain
A_n, with outward unit normals n_1,...,n_n. For an actual outside
wall c put beta_j=-c.n_j>=0. Say that c closes i if it is the actual
opposite of the outgoing i-edge.

**Geometric classification.** If c closes at least three directions,
then, up to reversing the old chain, its COMPLETE old contact vector
has one of the following two forms:

    H: beta_1=cos(pi/5), all other contacts zero;
       c closes exactly {1,2,3,4}.

    B: beta_1=cos(pi/4), beta_3=1/2, all others zero;
       c closes exactly {1,2,3}.

This is a necessary actual-wall statement. It does not assert that
either local star extends to a compact polytope in arbitrary dimension.
Every other actual outside wall closes at most two directions.

The inputs from the standard finite Coxeter classification are only:
finite diagrams are forests; a connected finite diagram of rank>=9
is A, B or D; a finite connected diagram containing label5 is I2(5),
H3 or H4; and a finite connected diagram containing a label>=7 has
rank2. The affine exclusions used below have explicit positive null
combinations. These inputs concern ACTUAL finite endpoint Grams.

## 2. At least three old contacts cannot make a large group

Choose three supported old nodes j<k<l. The underlying graph on
c and the old path j,...,k is a cycle, as is the graph on c and
k,...,l. Deleting any one old node other than k retains at least one
cycle. Consequently its white-plus-c Gram cannot be positive definite.
Thus a wall with at least three contacts closes at most the middle
node k, and in particular never three directions.

For exactly two contacts at j<k, a direction outside [j,k] also retains
the cycle c,j,...,k, so cannot be closed by c. The contact indices j,k
can supply at most two directions. A group of size at least three
therefore requires some closed i with j<i<k.

## 3. Two contacts and an interior closed direction

For j<i<k, deleting old i separates the old chain into two nonempty
paths, and c joins them. The endpoint's finite diagram is CONNECTED
of rank n>=14, hence is A_n, B_n or D_n.

In A_n both new bonds are3. But c and the old path j,...,k then form
an all3 affine cycle in the full actual wall system, forbidden by
compactness. The same conclusion holds in D_n, whose bonds are all3.
One can also see its limited shapes directly: any fork is at an old
contact node; its branch through c has length at least two, so its two
other arms must both have length one. This does not avoid the all3
cycle.

In B_n there is exactly one4 bond, at an END of its path. Since the
old bonds are all3, this bond joins c to an old contact node. For it
to be an end bond, that entire old white path must be a singleton.
Up to reversal this forces i=2, j=1, beta_1=cos(pi/4). The other
contact has label3 and must be at an endpoint of the right white path:
either k=3 or k=n.

The k=n alternative contains a positive affine F4 path with consecutive
labels3,4,3,3. It is induced by old nodes2,1,c,n,n-1, respectively.
These are distinct and have no extra edges for n>=14. In unit-normal
coordinates the positive null vector is proportional to
(sqrt2,2sqrt2,3,2,1). Thus k=n is impossible in compact P.

The sole remaining two-contact pattern is B in Section1. Deleting1
gives the finite D_n diagram; deleting2 gives B_n; deleting3 gives
B3+A_(n-3). Deleting any other direction retains its nonfinite
four-cycle. Hence an ACTUAL wall of this pattern closes exactly
{1,2,3}. Actualness follows either from acute positive-definite
projection or the exact reciprocal-shadow criterion. In particular
the finite B_n receiver at direction2 must not be omitted.

## 4. One old contact

Suppose c contacts only old node k, with native value beta.

For beta=1/2, the full star is finite A_(n+1) or D_(n+1) if k is at
distance at most one from an end. Such a star cannot contain the old
rank-n positive definite vertex together with an actual outside wall:
the full normals span Lorentz space, so their Gram is indefinite.
At the remaining positions the full star contains one of the all3
affine trees with arm lengths (1,2,5) or (1,3,3). Their positive null
vectors have central coefficients6 and4, with arm coefficients

    (3); (4,2); (5,4,3,2,1),
    (2); (3,2,1); (3,2,1),

respectively. Thus an isolated3 contact is impossible here.

For beta=cos(pi/4), an end contact gives finite B_(n+1); an internal
contact gives the affine three-leaf star with squared bond cosines
1/2,1/4,1/4. Both are excluded. For beta=cos(pi/6), any old neighbor
gives the positive affine3--6 path, also excluded. Value beta=1 is
forbidden by the positive null sum c+n_k.

For beta=cos(pi/5), retaining k in a finite replacement requires the
connected component containing c to be I2(5), H3 or H4. Thus old k
must be an end of its retained old path and that path has at most
three nodes. If k is an old chain end, precisely the first four
directions are closed: this is H. If k is internal, only direction k
and possibly one neighboring direction can be closed. The two old
tails cannot both be short when n>=14.

For a finite native label>=7 the same argument uses I2(m) only, so
there can be at most two closed directions. For a dotted contact
beta>1, any replacement retaining old k has a non-positive-definite
two-wall submatrix. Only direction k can be closed. These cases
complete the one-contact inventory and prove Section1.

## 5. Arithmetic removes the three-direction pattern from n14 onward

Assume now the usual arithmetic admissibility property: the ground
field is totally real and every nonphysical embedding gives positive
definite wall Gram matrices after positive diagonal normalization.
This holds for the compact arithmetic Coxeter polytopes under study.
No finite list of ground fields is used.

For B, its four-cycle product is a nonzero rational multiple of sqrt2,
so sqrt2 lies in the cyclic ground field. An embedding sending sqrt2
to -sqrt2 extends to a nonphysical embedding of that field. Apply it
to the full old A_n plus c star. The old Gram is unchanged. Put N=n+1.
The Schur projection cost of c at that embedding is

    q'=[5N-11-2sqrt2*(N-3)]/(2N).                    (1)

Indeed C_11=2(N-1)/N, C_33=6(N-3)/N, C_13=2(N-3)/N,
and the conjugated contacts are -sqrt2/2 and1/2. Positive definiteness
requires q'<1. But q' increases with N, and at N=15,

    q'=32/15-4sqrt2/5>1,

because17>12sqrt2, equivalently289>288. Thus B is impossible for
arithmetic n>=14. This exact threshold matters: at N=14 the cost is
below one, so it cannot be discarded by finite classification alone.

**Arithmetic consequence.** At a whole A_n vertex in arithmetic
dimension n>=14, every opposite closes at most two directions, except
possibly ONE isolated5 wall at one old chain end, closing four.

There cannot be two such end walls even geometrically. Their common
height square is cos^2(pi/5)*2(N-1)/N-1, so their unit product is

    1-phi^2/2+phi^2/N,    phi=(1+sqrt5)/2.

For N>=15 it lies strictly between -1/2 and0, forbidden for distinct
actual outward Coxeter normals. Hence at most one four-direction
group is possible.

This additional native same-wall information excludes the two
four-end-group interpretation of the formal A14 relaxation in
`compact_an_pair_shadow_energy_relaxation_barrier_20261020.md`.
It does NOT yet prove that every remaining exact reciprocal profile
has farthest-vertex energy at least one. That simultaneous energy and
interior-halfspace obligation remains open; no new dimension ceiling
or universal facet-descent claim follows from this classification.
