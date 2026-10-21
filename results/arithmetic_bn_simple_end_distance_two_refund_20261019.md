# A full-B19/B20 simple-end distance-two refund in the eight high-label fields

19 October 2026. **New written consumer theorem, awaiting independent
review.** No numerical job, proof script, optimizer, Lean process or build
was run. Existing pinned evidence and weight sources are unchanged.

This note proves a different proposed deletion from the adjacent-end
unit obstructed by the Q(sqrt(3)) Lambert model. It deletes the mark at
the FIRST and THIRD roots at B_n's simple end. All simultaneous deletions
of this same role on a face are included. The theorem gives one saving
at a full B20 vertex, reducing its cost from 37 to 36; it does not reach
the allowance 35 and does not lower any dimension ceiling.

## 1. Field scope, charge and statement

Let P be a compact arithmetic simplest-type Coxeter polytope of dimension
n=19 or 20 whose ground field is one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Use the charge already proved in
`results/fresh_arithmetic_high_label_20261016.md`:

* one unit for a marked pair at graph distance at most two;
* an additional unit for an entire marked I2(m) component with m>=6.

At a vertex whose ENTIRE diagram is B_n, start from the simple end
and call its first three roots a,b,c. Thus a-b-c-t1-...-t_(n-3)
is the old path, every displayed early bond has label 3, and the only
4-bond is at the remote end. Delete the one baseline unit at (a,c).
Perform this deletion at EVERY vertex of this whole-B_n role.

**Consumer theorem.** The resulting nonnegative charge still pays at
least two on every triangle, at least one on every quadrilateral,
and the required nonpositive amount 5-k on every k-gon with k>=5.

The theorem is restricted to whole B_n vertices. It does not extend
the deletion to B_(n-1)+A1 or to arbitrary outside factors.

Inputs are compactness, acute-wall/finite-face correspondence, finite
Coxeter classification, arithmetic cyclic products and nonphysical
positive definiteness, and the existing eight-field triangle/quadrilateral
charge theorem. Each specified field has degree at most five. The new
argument does not require a new enumeration or extend an unproved
finite-rank computation to arbitrary rank.

## 2. Triangles are protected before analyzing the new quads

In these eight fields the inherited triangle proof finds a distinguished
black bond of label at least six. At its corner it is an entire I2
component, and the baseline plus high-label unit pay two there. That
corner cannot have an entire large connected B_n diagram. This payment
is therefore unchanged by EVERY deletion considered here, including
simultaneous deletions at several other triangle corners.

All weights remain nonnegative, since the deleted unit is an existing
baseline one. Hence polygons with at least five sides need no further
argument. Only a quadrilateral newly reduced to zero could defeat
the proposed rule.

## 3. Complete finite-corner reduction of a newly unpaid quad

Fix such a quad at its old selected corner (a,c), and let E,D be the
other side roots in cyclic order a,c,E,D. Its fixed white support is

    W={b}+{t1,...,t_(n-3)} = A1+B_(n-3).

All roots below are unit roots; all contacts not explicitly permitted
must still be determined from finite classification or the zero-charge
condition. The opposite products are

    (a,E)=-x,       (c,D)=-y,       x,y>=0.

Since any remaining positive unit would pay this quad, each of its other
three corners has modified charge zero. In particular each adjacent
black product is zero: a nonright adjacency would retain a baseline
unit, whereas a deleted full-B_n mark is an ORTHOGONAL distance-two pair.
Thus (c,E)=(E,D)=(D,a)=0.

### The E possibilities

The actual finite corner W+c+E contains the entire B_(n-1) path
b-c-t1-...-t_(n-3). A finite one-root continuation at these ranks is
either isolated, or a simple extension at its simple endpoint b.
Indeed its connected component, already containing a 4-bond and having
rank at least 18, can only be a B path. The terminal 4-bond must remain
terminal; a branch, a second 4-bond, or attachment at the short endpoint
cannot be a finite B path. Exceptional finite diagrams have insufficient
rank.

Thus precisely two E types can have modified charge zero:

    I: E is perpendicular to W+c;
    L: E.b=-1/2, and E is perpendicular to the other W roots and c.

For I the marked corner (c,E) had original charge zero. For L it is
the same full-B_n deleted role, with path E-b-c-t1-... . This includes
the possible adjacent simultaneous deletion rather than assuming that
the other corner was originally unpaid.

### The D possibilities

The corner W+a+D initially contains A2(a,b)+B_(n-3). Zero modified
charge implies a.D=0. If b.D is nonzero, the original distance-two
unit survives unless this corner itself is the selected full B_n.
In that selected case D must join b to t1, by two simple bonds; a
finite B path leaves no other placement or label.

If b.D=0, D has no contact with A2(a,b). Its finite continuation of
the remote B_(n-3) is isolated or a simple extension at t1. These are
the only extensions for the same large-B-path reason used above.

Consequently precisely three D types remain:

    I: D is perpendicular to W+a;
    T: D.t1=-1/2, and D is perpendicular to the other W roots and a;
    J: D.b=D.t1=-1/2, and D is perpendicular to the other W roots and a.

Types I,T give original zero at (a,D). Type J gives a simultaneous
deleted full-B_n mark, with path a-b-D-t1-... . No unspecified contact
is being set to zero merely for convenience.

Finally E.D=0 is forced by the zero-charge opposite corner. In the
combination L,J, W+E+D is again the whole B_n path E-b-D-t1-...,
and its marked (E,D) unit is also deleted. In the other five
combinations that marked baseline was already zero. Thus these SIX
combinations include every newly unpaid quad and all simultaneous
deletion patterns. There is no unexamined three-selected-corner branch:
the L,J combination automatically has all four corners selected.

## 4. The exact projected Gram data

Eliminate W, whose B_(n-3) simple-end inverse entry is 2. The old
projected norms and cross product are

    a_bar^2=3/4,       c_bar^2=1/4,       a_bar.c_bar=-1/4.

The new projected data are as follows. An entry in the cross columns
is a projected product, not a claim about an ambient dihedral angle.

| Type | Norm | Product with a_bar | Product with c_bar |
|---|---:|---:|---:|
| E=I | 1 | -x | 0 |
| E=L | 3/4 | -x-1/4 | -1/4 |
| D=I | 1 | 0 | -y |
| D=T | 1/2 | 0 | -y-1/2 |
| D=J | 1/4 | -1/4 | -y-3/4 |

The projected E,D product is zero except for L,J, when it is -1/4.
These terms follow directly from the shared white roots b,t1; in
particular the products created by eliminating b must not be dropped.

Opposite sides do not meet in this actual compact two-face. A positive
definite opposite extension would produce an actual intersection by
acute projection, and a positive semidefinite singular one would be
affine in the compact chamber. Therefore its two-root projected block
is indefinite. At the a,E pair this gives

    E=I: x^2>3/4;         E=L: x>1/2.                 (1)

The case x=1 is impossible for actual distinct walls of a compact
acute chamber: their positive null sum would give a feasible ideal
point. The same ambient exclusion applies to y=1 when needed below.

The old n roots together with E already have Lorentz rank n+1, as
also follows from the projection costs in Section 5. Hence the full
n+2 local root system has rank exactly n+1. Its four-by-four white
Schur matrix has corank one. Put

    z=4x^2,       w=4y^2,       u=2x,       v=2y.

Its determinant gives these exact identities:

| E,D types | Corank identity |
|---|---|
| I,I | (z-3)(w-1)=1 |
| I,T | (z-3)[2(v+1)^2-1]=1 |
| I,J | (v+1)(z-3)=1 |
| L,I | (u-1)(w-1)=1 |
| L,T | (u-1)[2(v+1)^2-1]=1 |
| L,J | (x-1/2)(y+1/2)=1/4 |

For example I,I follows by eliminating E,D and setting
(3/4-x^2)(1/4-y^2)-1/16=0. For I,J use the symmetric and antisymmetric
combinations of c,D; for L,I and L,T use those of a,E. In L,J both
pairs split: the symmetric two-by-two block has diagonal 1/2-x,
-1/2-y and cross product -1/2, giving the last identity. These
calculations use only the displayed projected matrix.

## 5. The finite-angle alternatives for E are very short

This section is needed only when x is a finite ambient dihedral
cosine. Dotted products x>1 will be excluded directly from the
corank identities. Products equal to 1 are forbidden by compactness.

### I: only labels 12 and 18 survive

An E of type I contacts only the old endpoint a. The old whole-B_n
projection cost is q=2x^2. At every nonphysical embedding it satisfies
q'<1, and the squared contact remains positive. The negative physical
determinant ratio of the independent old-B_n-plus-E block also shows

    K=Q(x^2),       0<(x^2)'<1/2.                    (2)

For if K properly extended Q(x^2), a nonidentity embedding fixing x^2
would preserve that negative physical determinant, contradicting
nonphysical positive definiteness.

If x=cos(pi/m), (1) gives m>=7. Degree at most five and
[Q(cos^2(pi/m)):Q]=phi(m)/2 give exactly

    m=7,8,9,10,11,12,14,15,16,18,20,22,24,30.

Labels 7,14 generate K7; 8 generates Q(sqrt(2)); 10 generates
Q(sqrt(5)). None is one of the eight allowed ground fields. For
each remaining label except 12,18 the table exhibits a nonphysical
conjugate cos(2j*pi/m)>0; the corresponding x^2 conjugate then
exceeds 1/2, contradicting (2).

| m | j |
|---:|---:|
|9|2|
|11|2|
|15|2|
|16|3|
|20|3|
|22|3|
|24|5|
|30|7|

Every displayed j is coprime to m, distinct from the physical sign
class, and 0<2j*pi/m<pi/2. This proves, without a numerical test,

    E=I, x finite: m=12 or18.                         (3)

### L: only labels 6 and 9 survive

Here E contacts a by -x and b by -1/2. The relevant old-B_n inverse
entries at a,b are 2,2,4, so its projection cost is

    q=2x^2+2x+1.

The triangle a-b-E has cyclic product -2x. Thus x belongs to K,
2x is integral, and the same independent-block argument gives

    K=Q(x),       -1<x'<0 at every nonphysical embedding.  (4)

For finite x=cos(pi/m), (1) gives m>=4. Degree at most five gives
exactly m=4,5,6,7,8,9,10,11,12,15. Labels 4,5 generate the excluded
fields Q(sqrt(2)), Q(sqrt(5)). Each of 7,8,10,11 has the positive
nonphysical conjugate cos(3pi/m), while 12 has cos(5pi/12)>0 and
15 has cos(7pi/15)>0. These violate (4). Hence

    E=L, x finite: m=6 or9.                          (5)

Only the stated eight-field scope is used in (3),(5). These conclusions
must not be transported to Q(sqrt(2)), K7 or the golden fields.

## 6. Exclusion of all six quad types

Ambient actual products are either zero, finite cosines at least 1/2,
or at least 1. Thus if y is nonzero then y>=1/2. A finite value of
y strictly between cos(pi/4) and cos(pi/5) is impossible, as is a
dotted value in that interval.

### I,I and L,I

Both identities have the form (h-1)(w-1)=1, where h=z-2 for I,I
and h=u for L,I. Positivity in (1) gives w>1, hence y>1/2 and w>=2.
Consequently h<=2. If x>=1 this forces x=1 or is already contradictory,
and x=1 is forbidden. So x is finite.

In I,I the possibilities (3) give

    z-3=sqrt(3)-1, or 2cos(pi/9)-1.

In L,I the possibilities (5) give exactly the same two possible
values of u-1. Since sqrt(3)<2cos(pi/9)<2, in both cases

    2<w<=1+1/(sqrt(3)-1)=(3+sqrt(3))/2
       <(3+sqrt(5))/2=4cos^2(pi/5).

Thus cos^2(pi/4)<y^2<cos^2(pi/5), an impossible ambient contact.
This excludes both cases, including any dotted possibility for x.

### I,T and I,J

If y=0, either identity gives z=4 and x=1, forbidden. If y>0,
then v>=1. The I,T identity gives z<=3+1/7; I,J gives z<=3+1/2.
Both rule out x>=1. But the finite alternatives (3) satisfy

    z>=2+sqrt(3)>7/2,

so they are also impossible.

### L,T

If y=0 the identity gives u=2 and x=1, forbidden. Otherwise v>=1
and u<=1+1/7, so x<=4/7. This contradicts x>1/2 and ambient
cosine discreteness: the first permitted finite value above 1/2 is
cos(pi/4)>4/7. Dotted x is larger still. Thus this case is excluded
even without (5).

### L,J, including four simultaneous deletions

If y=0 the last identity gives x=1, forbidden. Otherwise y>=1/2
and 1/2<x<=3/4. The only ambient cosine in that interval is
x=cos(pi/4)=sqrt(2)/2; a dotted product is impossible. Substitution
gives y=sqrt(2)/2 as well. Equation (4) then forces K=Q(sqrt(2)),
outside the eight-field scope. This excludes the final case.

Every newly zero quad therefore leads to a contradiction. With
Section 2 this proves the consumer theorem of Section 1, conditional
on independent review of the written argument.

## 7. Composition with the D_n whisker-pair deletion

The separate notes
`results/arithmetic_dn_whisker_pair_triangle_retention_20261019.md` and
`results/arithmetic_dn_whisker_pair_quad_retention_20261019.md` address
deletion of the short-whisker-pair unit at whole D19/D20 vertices.
Their independent written review is
`results/review_dn_whisker_pair_retention_finite57_20261019.md`.

The two deletions cannot interfere on a quadrilateral. Any face
containing the new B_n deleted role has fixed white support containing
the remote 4-bond. Every corner of that face contains that same bond.
It cannot therefore have a whole D_n corner, whose bonds are all simple.
This remains true regardless of the number of other selected B_n
corners on the quad, which were included above.

In the eight fields, all triangles remain paid by their unchanged
entire high-label I2 corner after both kinds of deletion. Larger
polygons have nonnegative charges. Consequently the two consumer
theorems compose to a valid partial charge improvement within this
explicit field scope.

At dimension20 the resulting full connected costs are

    A20: 37;       B20: 36;       D20: 37.

The D20 entry starts from 38 and uses the separate whisker theorem.
No deletion is being asserted at B19+A1, D19+A1, or any other forest
with outside factors. Thus the combined rule does not supply a
vertex allowance 35, a dimension-20 exclusion, or a new universal
arithmetic dimension ceiling. At least the remaining A20 and D20
costs, the second B20 saving, and the large multifactor cases still
need complete consumer arguments.

The earlier adjacent-end B_n Lambert survivor is preserved in
`results/arithmetic_bn_simple_end_quad_survivor_20261019.md`. The new
distance-two theorem does not delete its sole paid adjacent mark,
so there is no contradiction between the two results.
