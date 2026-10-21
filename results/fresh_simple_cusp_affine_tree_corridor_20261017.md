# Long central corridors at affine tree corners

There is an all-rank paid branch for quadrilaterals with an affine
B/C/D corner. It is narrower than the full tree case. Suppose its two
marked roots a,b lie on the central simple-bond path, at distance P>=4,
and the two end diagrams beyond a,b are finite B or D diagrams of ranks
r,s>=9, attached at their long endpoints. Then the opposite corner has
adjacent marked roots and pays the quadrilateral. Here n=r+s+P+1>=23.
Marked roots in the bounded decorated-end regions remain open. This
lemma gives **no new ambient dimension bound**.

The [checker](../scripts/fresh_simple_cusp_affine_tree_corridor_20261017.py)
and [data](fresh_simple_cusp_affine_tree_corridor_20261017.json) bind the
[affine-A proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) and
its [independent review](review_simple_cusp_affine_a_quadrilateral_20261017.md).
The ambient assumptions are again finite volume and simplicity of all
vertices. Thus the [global exclusion](fresh_simple_cusp_affine_exclusion_20261017.md)
forbids every affine subdiagram of order below n.

## The complete four-profile list

Deleting a from the old affine corner gives two finite components, each
of type B or D and rank>=9. The same holds after deleting b. A new root
c replacing a can therefore have only these contacts with the remaining
old diagram:

1. no contacts;
2. one bond3 at the exposed long endpoint of the outer component;
3. one bond3 at the exposed long endpoint of the inner component;
4. bonds3 at both of those long endpoints.

To see completeness, every old component remains whole in any connected
new component. A connected finite extension of one B_r or D_r, r>=9,
can only be B_(r+1) or D_(r+1), extending its long endpoint by a simple
bond: exceptional types are too small, and the retained double bond or
fork rules out the other classical finite types. A connected affine
extension of only one old component would leave the other component
disconnected. This is neither a spherical corner nor the irreducible
affine diagram of a simple ideal corner, so it is forbidden.

If c touches both old components, the new connected diagram has all n
nodes. Two double-bond ends, a double-bond end and a fork, or two forks
can occur only in the corresponding affine C, B, or D diagram. In each
case the joining root has simple bonds to exactly the two long endpoints.
An extra attachment would produce a cycle or additional branch absent
from those diagrams. This proves the four-profile list using the full
finite/affine classification, in every rank. It does not extrapolate a
finite enumeration. In particular there are no free unbounded-label
contacts to old white roots in this list.

For c the outer component lies beyond a, and the inner component contains
b. The nearest contact of c on that inner component is the first white
node after a. Since P>=4, c.b=0. Analogously d.a=0 for the root d opposite
b. If c.d<0 the opposite marked pair is adjacent and already paid; assume
c.d=0 for a contradiction.

## Explicit low-contact exclusions

Normalize all roots to norm2. Put x=-a.c. For the isolated profile, x=0
would make c orthogonal to the entire simple cusp, hence null. The value
x=sqrt(2) creates the four-node affine B3 star centered at a. The value
x=1 makes c a simple leaf on the central corridor.

Every such simple leaf creates a proper affine E8 with arms(1,2,5).
For a leaf at a, take the length5 arm towards the long outer block and
the length2 arm towards b. For a leaf at the exposed outer port, take
five steps into that block and two through a. For a leaf at the first
inner white node, take two steps towards a and the outer port, and five
towards b and, if necessary, into the other long block. Ranks r,s>=9
and P>=4 provide these simple-bond arms; the long port of D9 still has
six simple steps before its fork, and B9 has seven before its double
bond. The original diagram is a tree, so no extra edge closes these
arms. The new leaf has no other old contacts in each case.

Thus isolated x=1 is excluded. A single-bond profile with x=0 is the
same forbidden leaf construction; with x=1 it creates an affine A2
triangle. For the both-ends profile, x=0 creates an affine A3 four-cycle,
and x=1 creates an affine A2 triangle. All these affine diagrams have
order<=9<n. The literal positive kernels are checked in the inherited
affine-A package.

The normalized cusp kernel has coefficient1 at a,b and their neighboring
long ports. Consequently the height h=-c.nu equals x, x+1, or x+2 for
the isolated, single, or both profiles. The same complete Coxeter
alphabet argument as before gives conservative reciprocal bounds

    u=1/h <= 1/sqrt(2), sqrt(2)-1, or 1-sqrt(2)/2,

respectively. All heights are positive. The identical bounds hold for
v=1/k of d.

## Schur reduction and uniform gap

The inverse Gram at the long port of either finite B_r or D_r equals1.
This holds in every rank: its inverse column is1 along the main path,
with the last B-node coefficient1/sqrt(2), or the two D-fork leaves
each1/2. Direct multiplication gives the unit vector at the port.
The checker confirms these literal columns in184 examples of ranks9..100;
the displayed multiplication identity proves every rank.

Eliminating W leaves the residual a,b block alpha[[1,-1],[-1,1]], now
with alpha=1/P. The intervening A_(P-1) path supplies the only cross
term between c and d. Introduce moment l and constant c0 for c:

| Profile | l | c0 | upper bound on 2*l*u+c0*u^2 |
|---|---:|---:|---|
| isolated |0|2|1|
| outer single3 |0|1|3-2sqrt(2)|
| inner single3 |1|1|1|
| both3 |1|0|2-sqrt(2)|

Then B=l/P and C=c0+l^2/P in the notation of the affine-A Schur
identity. For d use m,d0, giving A=m/P, D=d0+m^2/P and E=l*m/P.
Substituting into that formally checked determinant identity cancels the
quadratic cross term and gives

    detS/(h^2*k^2)
      =1-[2*(l*u+m*v)+c0*u^2+d0*v^2]/P.

Each single-root expression is increasing on its permitted reciprocal
height interval and is at most1 by the exact table. Therefore the right
side is at least 1-2/P>=1/2. The actual quadrilateral star has corank1,
so detS=0, a contradiction. Hence c.d<0, as claimed.

## Status

The four exact profile costs and184 inverse-column checks are generated
and freshly replayed by the checker. Syntax, all ten direct source hashes
and all five local links pass. A separate internal mathematical/source
read confirms the full B/D extension list, the three leaf-E8 embeddings
including P=4, the all-rank port columns and the scalar Schur cancellation.
That read does not claim a duplicate184-case arithmetic replay.
The completeness arguments above
are classification and all-rank identities, not finite cutoff claims.
Only this long central-corridor branch is closed. The full affine-tree
quadrilateral obligation, including all decorated-end positions, remains
open. The mixed rule and earlier dimension bounds are unchanged.
