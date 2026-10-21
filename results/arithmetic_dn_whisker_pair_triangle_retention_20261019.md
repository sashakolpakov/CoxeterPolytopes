# Full D19/D20 whisker-pair triangles retain both other baseline units

19 October 2026. **New written geometric lemma, awaiting independent
review.** No computation, enumeration or build was run. The source-bound
opposed-fork model notes are unchanged. This note addresses triangles
only; it proves no universal refund or dimension bound.

## 1. Statement and weight convention

Let an actual finite vertex of a compact hyperbolic Coxeter polytope
have whole diagram D_n, with n=19 or 20. Let a be its fork, b,e
its short whiskers, and a-l1-...-lL its long arm, L=n-3. All old
edges have label 3. Work with square-one outward roots.

Let B be the ordinary distance-two baseline on finite marked diagrams:

    B(C;i,j)=1 if i,j lie in the same component at distance at most 2,
              0 otherwise,

with any nonnegative existing extras left in place. There are no H4
components in the finite corners considered below.

**Triangle-retention lemma.** If an actual triangle occurs at the old
marked pair (b,e), then both its other corners have baseline 1.
Neither of those two angles is any of the three distinguished
full-D_n angles (fork with either whisker, or the two whiskers).

Consequently deleting the whisker-pair unit at every whole D_n vertex
would leave the two required triangle units intact, even if both
fork/whisker units were also deleted at whole D_n vertices. This is
only a triangle-consumer conclusion. Quads are a separate obligation.
No arithmetic input is used.

## 2. A zero other corner would force a very short continuation list

Let c be the actual third wall of the triangle. Its white support is

    W={a,l1,...,lL},     W=A_(n-2),

whose two ends are a,lL. The other actual finite corner W+b+c
contains the long path

    b--a--l1-- ... --lL

of type A_(n-1) before c is added.

Suppose its marked pair (b,c) has baseline zero. A nonzero product
b.c would make the marked roots adjacent, so b.c=0. By finite
classification, a one-root extension of this path of rank at least 18
is isolated, an endpoint extension of label 3 or 4, or a simple
next-to-end extension producing a D diagram. A second contact
creates a cycle and cannot be finite. Exceptional finite diagrams
are too small.

The near endpoint b cannot be contacted because b.c=0. The near
next-to-end node a cannot be contacted because the marked path
b-a-c would then have distance 2 and baseline 1. Thus exactly four
types remain for c against W+b:

1. c is perpendicular to W+b;
2. c has a simple contact with lL and no other contact;
3. c has a 4-bond contact with lL and no other contact;
4. c has a simple contact with l_(L-1) and no other contact.

The product c.e is still unrestricted at this stage. It occurs in
the separate actual finite corner W+e+c and must be checked there.

## 3. Every continuation in that list is impossible

In the isolated type, c.e is zero or a finite dihedral cosine. The
corner W+e is the path e-a-l1-...-lL of rank n-1. If c.e is nonzero,
finite classification of its endpoint extension restricts its label
to 3 or 4. If c.e=0, c is perpendicular to the whole positive old
D_n space; its perpendicular is a timelike line and cannot contain
the unit spacelike normal c.

If c.e has label 3, the actual roots a,b,e,c,l1,...,l5 form an
induced affine E8 diagram, with arms of lengths 1,2,5 at a.
If c.e has label 4, the roots {a,b,l1,e,c} form the affine B4
diagram with its terminal 4-bond at e-c. Both are impossible in a
compact chamber.

For each of the three far types, a nonzero c.e would create a cycle
in the actual finite diagram W+e+c. Hence c.e=0.

* A simple far-end contact makes all old D_n roots plus c a positive
  definite D_(n+1). A positive (n+1)-space cannot fit in Lorentz
  signature (n,1).
* A far-end 4-bond makes an induced affine B diagram on all old
  roots plus c. Its positive null coefficients are 2 at a and all
  old long-arm nodes, 1 at b,e, and sqrt(2) at c.
* A simple far next-to-end contact gives two forks a,l_(L-1),
  joined by the old simple spine, with leaves b,e at the first and
  lL,c at the second. Its positive null coefficients are 2 along
  the whole spine and 1 at all four leaves: an affine D diagram.

For the two small affine witnesses above, the square-one coefficients
are: affine E8 center a coefficient 6, arm b coefficient 3, arm e,c
coefficients 4,2, and arm l1,...,l5 coefficients 5,4,3,2,1;
affine B4 coefficients a=e=2, b=l1=1, c=sqrt(2).
Each kernel equation is read directly from the displayed bonds.

A positive-null combination of actual acute outward normals gives a
nonzero future null vector feasible against every wall. Nonzeroness
follows from its negative pairing with an interior timelike point;
feasibility follows from acuteness. Compactness forbids the resulting
ideal ray. Thus all listed affine cases are impossible.

This contradicts the assumed baseline zero at (b,c). Exchanging b,e
gives the same contradiction for (e,c). Both other triangle corners
therefore have baseline 1.

## 4. Neither remaining unit can be one of the proposed D_n deletions

At the corner W+b+c, the old root b already has exactly one white
neighbor, namely a. Adding the marked root c can raise its degree
by at most one. It cannot become a fork of degree three. If c were
the fork and b its short whisker, their bond would be nonzero, but
b would then have two neighbors a,c and would not be a whisker.
Thus (b,c) cannot be a full-D_n fork/whisker pair.

Suppose instead that W+b+c were whole D_n with b,c its two short
whiskers. Since b has its fixed neighbor a, the fork must be a.
Consequently c has only its simple contact with a in this corner.
The actual third corner W+e+c forces c.e=0, because a,e,c would
otherwise create a cycle in its finite diagram. Also b.c=0 by the
two-whisker diagram. Then {a,b,e,c,l1} is the induced affine D4
star with center a and four simple leaves. Its coefficients 2 at a
and 1 at the leaves contradict compactness.

Therefore (b,c) cannot be a full-D_n whisker pair either. The same
argument applies to (e,c). This proves the retention statement even
under simultaneous deletion of all three distinguished units, solely
as far as these triangle consumers are concerned.

## 5. The separate quad problem remains

For a quadrilateral at (b,e), an unpaid adjacent corner similarly
forces its new wall to be isolated, or to have a far-end 3/4 contact,
or a simple far next-to-end contact. The finite opposite corner then
couples the two new profiles. If both are isolated, compactness and
the quad determinant force one ambient opposite product finite and
the other dotted; the previously reviewed double-whisker arithmetic
lemma excludes that special branch. This observation does not cover
every joint continuation or every simultaneous new zero corner.

The safe conclusion of this note is restricted to the fully proved
triangle-retention lemma. A whisker-pair unit refund still needs
all quadrilaterals and all interactions with the other selected
angles. The compact arithmetic ceiling remains 20.
