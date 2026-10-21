# Exact adaptive corridor at finite distance eight and feature depth six

The retained-path calculation has99 candidate contact vectors,13 accepted
profiles and576 adaptive Schur branches. Every branch has negative Schur
diagonals and a determinant polynomial with nonnegative exact Q(sqrt2)
coefficients and positive constant. The minimum constant is3/44. There is
no exceptional or equality branch.

This is a diagnostic for a proposed stronger local rule. It does **not**
by itself prove that the complete payment rule remains valid or improve
the completed unrestricted bound58. Short-seven replacements, the changed
cycle boundary cases and simultaneous composition are separate obligations.
The [producer](../scripts/fresh_depth6_finite8_adaptive_pilot_20261018.py)
and [certificate](fresh_depth6_finite8_adaptive_pilot_20261018.json) retain
that limited scope explicitly.

## Complete retained domain

Use an actual quadrilateral star W+a+b+c+d, with opposite pairs a,c and
b,d and roots of square two. Each W plus one black is positive definite;
adjacent corners are positive semidefinite; opposite pairs with W are
indefinite. Suppose the selected whole finite A/B/D component has marked
distance eight and both marks have distance at least six from every leaf,
branch and endpoint of a nonsimple bond. It supplies a simple A21 path
numbered0 through20, with a=6,b=14. Additional supporting whites remain
arbitrary. Features at the retained endpoints are allowed.

Assume no adjacent black pair has a white-interior path of length at
most seven. Then c.b=d.a=c.d=0. The retained whites split as A6+A7+A6.
Their inverse entries are nonnegative with minimum diagonal6/7. Three
native contacts give projection at least18/7>2. A label at least five
has squared contact magnitude greater than5/2 and gives projection
greater than15/7>2. Thus c has at most two retained white contacts, each
one or sqrt2. Its permitted positions are0,1,2,3,4,5,7. All others give
a short adjacent path to b.

Strict positivity of retained W+c and semidefiniteness of the adjacent
corner (A21 minus a)+c leave these13 profiles from all99 candidates:

    empty; {7}; {7:sqrt2}; {5}; {5,7}; {5:sqrt2};
    {4}; {3}; {2}; {1}; {0}; {0,7}; {0:sqrt2}.

Unwritten contacts are one. Reflect i to20-i for d. The two opposite
magnitudes t=-a.c and u=-b.d are independently zero or at least one.
They are not bounded above.

## Forced additional supporting whites

Let q be the projection of c onto A21, including its a-contact. If q<=2,
then t=0 because J21[6,6]=105/22>2. Only empty,{1},{0},{0:sqrt2} can
have q<=2. None has a retained contact within four steps of a.

Choose a minimal indefinite diagram K in W+a+c. The inherited all-label
classification bounds its order by ten; both marked deletions are
positive definite. The component through a in K minus c is a subtree of
the selected finite component. Every branch or nonsimple-bond feature
is at distance at least six from a. Its distance from a to c is at least
six: a closer contact would be inside the retained path, while leaving
that path already requires seven white steps. Diagrams of order at most
five cannot occur.

The producer screens all768 ordered distinct finite-deletion pairs of
the inherited26-diagram higher-order catalogue, allowing disconnected
K minus c. Only diagram25 with (a,c)=(8,1) or(8,2) survives. These are
the two E10 positions from the
[distance-seven theorem](fresh_depth6_finite7_adaptive_transfer_20261018.md).
The forced c-contact is five steps from a. It cannot point toward b,
where it gives a c-b path of length four. It is therefore retained
position1, so the complete low profile must be {1}.

Writing c=21 and the additional white as z=23, the maps from E10 nodes
0 through9 are

    extension: 1,21,0,23,2,3,4,5,6,7;
    separate:  1,0,21,23,2,3,4,5,6,7.

In the first case z extends A21 simply at0. The selected finite tree
excludes any second retained contact. In the second z joins c simply
and is orthogonal to retained positions0 through7. Every remaining
interior position is within five steps of a or b, so another white
contact there would create a forbidden branch feature. A possible
contact at20 is ruled out by strict W+c positivity: the left A6 costs
10/7, and the far A6 with z costs at least7/8, totaling129/56>2.
For a native far-end magnitude p>=1 the latter term is
1/(2-p^2*6/7); a nonpositive denominator already violates supporting-white
positivity. Thus this z is orthogonal to the whole A21.

In both cases d.z=0, since otherwise there is an adjacent c-d path of
length at most four. If reflected witnesses are both needed, they are
distinct and orthogonal. A common white or mutual contact gives a path
of length at most five whenever a separate white is involved. Two
extensions instead lie in the selected finite tree and cannot form
a cycle. These restrictions justify the entire adaptive positive block;
unspecified supporting contacts are not assigned zero by convention.

## Unbounded exact Schur certificate

Enlarge the positive block by the forced simple endpoint extensions and
separate A1 whites. A separate white contributes1/2 to its own black
projection and zero to the cross projection. Endpoint extensions use
the A22 or A23 inverse. All enlarged projections q,r exceed two.

For every pair of retained profiles, first discard any already supplying
a c-d path of length at most seven. Split t into0 or1+T and u into0 or1+U,
with T,U>=0. Low-projection branches other than the forced singleton
are excluded by the catalogue argument; each singleton tests both actual
white alternatives. All576 remaining branches have

    q-2>0, r-2>0, (q-2)(r-2)-s^2>0.

Every coefficient is signed exactly in Q(sqrt2), with no numerical
positivity test. The minimum polynomial constant is3/44. Eliminating
the positive block consequently gives a negative-definite two-by-two
Schur block, impossible in an actual Lorentz Gram. The central singleton
pair {7},{13} has determinant12/11 and needs no separate treatment.
Under the stated actual-star and classification inputs, the no-short-seven
branch is therefore excluded. Turning that shorter path into a retained
payment is deliberately outside this packet.

## Counts and verification

The candidate depth-six distance-eight deletion counts for A_n,B_n,D_n
are max(0,n-20),max(0,n-21),max(0,n-22). The current depth-eight counts
are max(0,n-24),max(0,n-25),max(0,n-26). The eligible marks form a simple
spine interval, so these formulas count its possible starting positions.
Direct graph enumeration checks236 diagrams through rank80. These are
hypothetical payment savings until the composition obligation is proved.

Generation8074/chunkd22e0f and separate fresh replay6472/chunka3690c both
terminated exit0/PASS. Each reconstructs the complete retained domain,
768 catalogue pair checks, the two exact E10 maps,576 Schur polynomial
branches and236 graph counts. They do not regenerate the historical
catalogue or formally verify the supporting-white reasoning. The minimal
core classification, actual-star facts and finite/affine classification
remain explicit inherited inputs. No complete local payment theorem,
global dimension improvement, finite-volume example, external refereeing
or proof-assistant verification is claimed here.

    python3 scripts/fresh_depth6_finite8_adaptive_pilot_20261018.py --verify
