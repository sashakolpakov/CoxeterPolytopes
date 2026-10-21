# A spectator-aware half-unit end refund at whole B19+A1 vertices

19 October2026. **New written consumer theorem candidate, awaiting
independent review.** Source and paper work only. No scientific job,
enumeration, optimization, Lean process or build was run. All previous
source-bound A/B/D reports are unchanged.

## 1. Statement and precise charge scope

Let P be a compact arithmetic simplest-type hyperbolic Coxeter
polytope of dimension20, with ground field one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Start with the eight-field distance-at-most-two charge and its extra
unit on an entire I2(m>=6) component. At a vertex whose ENTIRE
diagram is B19+A1, select the first/third pair at the UNIQUE SIMPLE
end of the B19 component. Subtract1/2 at this angle, changing its
baseline1 to1/2. Do this at every such whole-forest vertex.

**Consumer claim.** This half-unit deletion preserves the triangle
and quadrilateral payments. It is also compatible with the earlier
whole-A20 and whole-D20 deletions, whole-B20 simple-end distance-two
deletion, and balanced simple-end quarter rule, because their changed
angles have different complete white graphs. No extension to other
spectator forests is asserted.

The existing balanced quarter rule adds1/4 at a DIFFERENT marked
pair in a whole B19+A1 forest: the isolated root and the simple
endpoint of B19. Thus the whole-forest cost would decrease from
141/4 to139/4, rather than merely to35. This does not resolve B20
or the remaining D/E two-component families and yields no dimension
ceiling19 by itself.

The geometric inputs are compact Coxeter simplicity, finite diagram
classification, acute positive-definite projection, absence of actual
positive-null principal diagrams, arithmetic cyclic integrality and
strict conjugate positivity. No enumerated continuation catalogue or
new field list is assumed. The exact eight-field charge and its
triangle payment are inherited from
`fresh_arithmetic_high_label_20261016.md`.

## 2. The white support and why only one selected corner matters

At a selected vertex write the simple end of B19 as

    a--b--c--S-- ... --(the far 4-bond end),

all displayed bonds3. Let h be the A1 spectator. The white graph of
the selected first/third pair (a,c) is

    W={b,h}+B16.

Here b,h are isolated white roots, and S is the unique simple
endpoint of the remote B16 component. This is the COMPLETE white
support, of rank18. The unit-normal B16 inverse has endpoint entry
(G_B16^-1)_(SS)=2. It follows, by eliminating W, that

    a_bar^2=3/4,  c_bar^2=1/4,  a_bar.c_bar=-1/4.

A quadrilateral with at least two selected B19+A1 angles retains
at least1 from those angles alone. Consequently a newly unpaid quad
would have exactly one selected source, and its other three corner
payments would be original zeros. All selected source angles have
weight1/2 after deletion; every unselected original positive angle
in these faces still has integral weight at least1.

For completeness, no earlier modified charge invalidates that last
assertion. Their white graphs are

* A20 end and D20 fork: A1+A17;
* D20 short-leaf pair: A18;
* B20 simple-end distance-two: A1+B17;
* B20 simple-end adjacent and its positive quarter recipients: B18.

None is B16+2A1. Thus every other positive charge here is unchanged,
and this proposed spectator rule is not being silently combined by
separate consumer arguments on a common white graph.

Use cyclic sides a,c,E,D, so the opposite pairs are (a,E),(c,D).
Assume the other three original weights are zero. In particular

    c.E=a.D=E.D=0.                                  (1)

Every other contact below is derived from this zero-weight condition
and complete finite continuation, rather than being set to zero by
convenience.

## 3. Complete finite continuation with the spectator retained

At the corner (c,E), a contact of E to b would make marked distance2
through b. Hence E.b=0. More generally W+c has type B18+A1: the
simple end of B18 is b. A finite one-root continuation touching that
large component can only extend its simple endpoint by a simple
bond. An interior attachment makes a branch, an additional 4-bond
makes a nonfinite large path, and exceptional finite diagrams are
too small. The only permitted contact would therefore be to b,
which was just excluded. Consequently E has NO remote B16 contact.

The complete possibilities for E against W are thus

    E isolated;   or   E.h=-gamma, gamma=cos(pi/ell), ell>=3.

The spectator dihedral label is initially unrestricted. Put
t=1-gamma^2, taking gamma=0 for the isolated case.

At (a,D), the zero mark similarly gives D.b=0. Finite classification
of W+a+D then gives exactly these possibilities for D:

    H: only a possible finite contact D.h=-delta (including delta=0);
    T: a simple contact D.S=-1/2, with delta=0 or1/2.

Indeed a remote contact can only be a simple endpoint extension of
B16. If it is combined with h, the new large finite path permits
only a further simple bond to h. Two endpoint 4-bonds, an interior
branch, or any other remote contact are nonfinite. D never contacts
b or a. These statements include arbitrary contacts to the old
spectator, rather than assuming h remains a spectator afterwards.

If gamma is nonzero, E.h and D.h cannot both be nonzero: that would
give marked distance2 at (E,D). Hence

    gamma!=0  implies  delta=0.                     (2)

Write a.E=-x, c.D=-y. After eliminating W, E and D are orthogonal:
their only possible shared white contact was h, forbidden by (2).
Set

    H: u=1-delta^2,       Y=y;
    T: u=1/2-delta^2,     Y=y+1/2.

The four-black Schur Gram is therefore

    [3/4   -1/4    -x     0 ]
    [-1/4   1/4     0    -Y ]
    [-x      0      t     0 ]
    [ 0     -Y      0     u ].

Every t,u is strictly positive because its old white extension is
finite. Opposite projected sides of an actual compact acute polygon
have indefinite two-wall Gram: a positive-definite pair gives an
actual intersection by acute projection, and a positive-null pair
gives a feasible ideal point. Thus

    x^2>3t/4,       Y^2>u/4.

The n+2 actual normals have corank at least1 in Lorentz dimension
n+1. Their Schur determinant gives the exact relation

    P*(Y^2/u-1/4)=1/16,
    P=x^2/t-3/4>0.                                  (3)

No native opposite product has been assumed dotted.

## 4. Independent old-vertex cost and its field consequence

The old whole B19+A1 vertex is a positive20-dimensional span, and
the new wall E does not contain that simple vertex. Its old-span
projection cost is

    q_E=2x^2+gamma^2>1.                              (4)

For example the coefficient2 is the a,a entry of the old inverse;
it also follows by inverting the displayed old a,c Schur block,
whose inverse is [2,2;2,6]. At every nonphysical embedding the
independent old-plus-E block is positive definite, giving q_E'<1.
The quantities x^2,gamma^2 lie in K by the length-two cyclic
products, and their determinant-ratio combination (4) is valid
under embeddings with the usual positive root scaling.

In particular

    K=Q(x^2,gamma^2).                                (5)

If this were a proper subfield of K, a nonidentity embedding fixing
both squares would preserve q_E>1, contradicting conjugate positivity.
We do not assume the individual unit-normal square roots lie in K.

## 5. E contacts the spectator: a finite four-label reduction

Suppose gamma!=0. By (2), D.h=0. In case H, u=1 and projected
opposition gives y>1/2; the next native value is cos(pi/4), so
y^2-1/4>=1/4. In case T, u=1/2 and y>=0 gives
(y+1/2)^2/u-1/4>=1/4. Thus (3) forces

    P<=1/4,       x^2<=t=1-gamma^2<=3/4.             (6)

The nonzero native x is consequently cos(pi/m), m=3,4,5 or6.
Combine (4), (6), and gamma's finite native values:

* m=3: gamma^2>1/2 and gamma^2<=3/4, so ell=5 or6.
  Equation (5) gives respectively Q(sqrt(5)) or Q.
* m=4: gamma^2<=1/2, so ell=3 or4. Both squares are rational,
  and (5) gives Q.
* m=5: gamma^2<1/2, so ell=3 only. Equation (5) gives Q(sqrt(5)).
* m=6: gamma^2<=1/4, so ell=3 only. Equation (5) gives Q.

None is one of the eight fields. This excludes EVERY spectator
label, not merely a bounded list chosen before the geometry.

## 6. E is spectator-orthogonal and x is finite

Now gamma=0, so t=1 and (4) gives K=Q(x^2). Equation (3) requires
x^2>3/4. At each nonphysical embedding,

    0<(x^2)'<1/2.

The same elementary native-label step as in the independently
reviewed whole-B simple-end distance-two source applies; it uses
only these squared-contact and field-generation statements.
Here is a full finite reduction to make the inheritance explicit.
For x=cos(pi/m), [Q(x^2):Q]=phi(m)/2<=5. The labels with m>6 are

    7,8,9,10,11,12,14,15,16,18,20,22,24,30.

Labels7,14 generate K7, absent from the eight; 8 generates Q(sqrt(2));
10 generates Q(sqrt(5)). The following remaining labels have a
nonphysical conjugate cos^2(a*pi/m)>1/2, since a/m<1/4:

    m: 9,11,15,16,20,22,24,30;
    a: 2, 2, 2, 3, 3, 3, 5, 7.

Each a is coprime to m and gives a nonidentity cosine-square
embedding. Thus only m=12 or18 remains, with actual fields
Q(sqrt(3)) or K9. Put

    A=1/4+1/(16P),       P=x^2-3/4.

Since cos(pi/18)>cos(pi/12), the larger A occurs at m=12 and is
(3+sqrt(3))/8<3/5. Also P<1/4, so

    1/2<A<3/5,       Y^2=u*A.                        (7)

For H with delta=0, (7) puts y^2 strictly between1/2 and3/5,
inside the native gap from cos^2(pi/4)=1/2 to
cos^2(pi/5)>3/5. For H with delta!=0, u<=3/4, so y^2<9/20<1/2.
Projected opposition excludes y=0. Hence y=1/2. Equation (7) then
forces u in (5/12,1/2), or delta^2 in (1/2,7/12), another native
gap. Thus neither H possibility occurs.

For T with delta=0, u=1/2 and (7) give
(y+1/2)^2 in (1/4,3/10). Therefore 0<y<1/2, impossible for a
native contact. For T with delta=1/2, u=1/4 and (7) give
(y+1/2)^2<3/20<1/4, impossible for y>=0. The finite-x branch is
fully excluded with the spectator contacts retained.

## 7. E is spectator-orthogonal and x is dotted

Native x=1 gives a positive-null pair of actual walls and is
forbidden by compactness. Suppose x>1. Then P=x^2-3/4>1/4.
Case T has u<=1/2 and Y=y+1/2>=1/2, so its second factor in
(3) is at least1/4. Their product is strictly greater than1/16,
a contradiction.

For H, equation (3) gives y^2/u<1/2, so native y can only be1/2.
The option delta=0 would give second factor0 and is impossible.
For nonzero delta, (3) gives

    P=u/(4(1-u)).

Since P>1/4, u>1/2. Native delta therefore has delta^2<1/2 and
must be1/2. Thus u=3/4, P=3/4 and x^2=3/2. But (4) is then
the rational physical value q_E=3>1. Any nonidentity ground-field
embedding would preserve that negative determinant ratio, forcing
K=Q, outside the eight. No dotted x survives either.

## 8. Payment conclusion and exact remaining frontier

The single-source quad with all three other original weights zero
has been excluded in every case. A quad with two or more half-weight
sources retains at least1. A quad with an unchanged positive angle
retains at least1 there. These exhaust the consumers. Triangles keep
their untouched entire high-label I2 weight2, and larger polygons
retain nonnegative charge.

This is therefore a complete written HALF-unit refund candidate
at the one selected simple-end distance-two role in whole B19+A1,
subject to independent review. It includes arbitrary permitted new
spectator contacts and does not infer spectator safety from the
whole-B19 proof alone.

The cost at this whole forest becomes 35+1/4-1/2=139/4.
The whole B20 cost remains141/4 under the reviewed simple-end rules.
Other two-component D/E/high-I2 forests remain to be lowered. No
global allowance35, compact arithmetic bound19, new global example
or formal verification is claimed here.
