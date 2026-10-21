# A golden opposed-fork quadrilateral survives the D19/D20 local tests

19 October 2026. **New paper-checked local model and capacity argument,
awaiting independent review.** No program, enumeration, optimizer,
proof process or build was run. Existing reports and manuscripts are
unchanged. This note supplies neither a compact polytope in dimension
19 or 20 nor a new dimension bound.

## 1. What this model changes

A new whole D_n corner can put its fork at either end of the retained
long chain. Restricting it to the old near end loses a genuine finite
continuation. Keeping the far end gives the following explicit local
model at n=19 and 20:

* four adjacent corner diagrams are finite;
* the white perpendicular space contains a genuine compact hyperbolic
  quadrilateral, with those four actual sides and finite vertices;
* the complete root Gram matrix has Lorentz signature and integral
  entries over Q(sqrt(5)), and its nonidentity conjugate is positive
  semidefinite of the required rank;
* no proper principal Gram submatrix is singular;
* the distance-at-most-two unit baseline pays the quadrilateral only
  at its two opposite whole-D_n fork/whisker corners.

Thus these local conditions do not justify deleting both fork/whisker
units at every whole D_n vertex. The model is not a globally realized
compact n-polytope: closing its other faces remains unproved. A
universal refund could still use additional geometric exclusion,
an existing bonus, or an explicit compensating charge. It cannot
discard this model merely as a nonfinite diagram or an arithmetic
conjugate failure.

There is also a useful actual-wall capacity statement: the two fork/
whisker angles of a given whole D_n vertex cannot both carry this
specific golden profile. This is proved in Section 7.

## 2. The complete root matrix

Let n be 19 or 20, put N=n-2 and L=n-3=N-1. The old D_n roots are

    a,b,e,l1,...,lL,

with fork a, whiskers b,e and long arm a-l1-...-lL. All old bonds
are simple. Adjoin roots c,d. Set

    x=(1+sqrt(5))/4=cos(pi/5),
    beta=(3+sqrt(5))/4>1.

Use square-one roots. The additional nonzero off-diagonal products are

    c.e=c.lL=c.d=-1/2,
    a.c=-x,
    b.d=-beta.

Every other new product is zero. In particular

    a.d=b.c=0.

The far-end contact c.lL is essential. The old fork a and the new
fork c share e as a short whisker, but their long chains run in
opposite directions. The opposite side pairs of the marked two-face
are (a,c) and (b,d), with cyclic order a,b,c,d.

The two golden parameters satisfy

    x*(beta-1)=1/4.                              (1)

Let W={e,l1,...,lL}; its Gram matrix is A1+A_L and is positive
definite. Direct inversion gives

    s=a_bar^2=c_bar^2=(N+2)/(4*N),
    a_W*G_W^-1*c_W^T=s.

The bar denotes orthogonal projection off W. Both b,d are unchanged
by this projection. The projected Gram matrix in order a,b,c,d is

    [  s     -1/2   -(x+s)    0   ]
    [ -1/2    1       0     -beta ]
    [-(x+s)   0       s     -1/2  ]
    [  0    -beta   -1/2     1   ].              (2)

## 3. Physical signature and all four finite corners

Taking the sums and differences of a,c and of b,d splits (2) into
two symmetric two-by-two blocks:

    P_+ = [ -x     -1/2  ],
          [-1/2   1-beta ],

    P_- = [ x+2*s  -1/2  ].
          [ -1/2   1+beta ]

By (1), P_+ has determinant zero and one strictly negative eigenvalue.
The block P_- is positive definite: its first diagonal is positive
and its determinant is

    (11+5*sqrt(5))/8 + (7+sqrt(5))/(4*N)>0.

Hence the projected matrix has two positive eigenvalues, one negative
eigenvalue and one zero. Restoring W gives the full n+2-root matrix
signature (n positive, 1 negative, 1 zero). These roots therefore
exist in an actual Lorentz space of dimension n+1.

The four adjacent n-root corner diagrams are

|Corner roots|Whole finite diagram|
|---|---|
|W+a+b|D_n, with fork a and whiskers b,e|
|W+b+c|A_(n-1)+A1, with b isolated|
|W+c+d|D_n, with fork c and whiskers d,e|
|W+d+a|A_(n-1)+A1, with d isolated|

In the second row the long path is e-c-lL-...-l1. In the fourth it
is e-a-l1-...-lL. Thus every adjacent corner really is finite.

## 4. The four projected walls bound a compact quadrilateral

The normal relation in the three-dimensional white perpendicular space is

    a_bar+c_bar=2*x*(b+d).

The vector b+d is timelike, since its square is 2*(1-beta)<0.
Choose its future direction. Its pairings with b,d are 1-beta<0,
and its pairings with a_bar,c_bar are -1/2. After normalizing it we
obtain a future timelike point strictly inside all four halfspaces.

Each adjacent pair in (2) has positive-definite Gram matrix: for
the nonorthogonal pairs the determinant is s-1/4=1/(2*N)>0, and
the other two pairs are orthogonal. Projection of the interior point
onto the perpendicular of an adjacent pair remains timelike. The
two-root inverse is entrywise nonnegative, so the correction is a
nonnegative combination of the two outward normals. Acuteness keeps
the other two wall inequalities strict. Thus all four adjacent finite
vertices are actual vertices of the four-halfspace intersection.

One can also see compactness directly in the Klein model centered at
that interior point. The spatial parts of a_bar,c_bar are opposite,
as are those of b,d. They span the two-dimensional spatial plane.
The four inequalities therefore bound a Euclidean parallelogram.
Its four vertices are the finite vertices just constructed, lying
strictly inside the ball. Their convex hull is compact inside the
ball. This is a genuine compact hyperbolic quadrilateral.

Its angles are theta,pi/2,theta,pi/2, where

    cos(theta)=1/(2*sqrt(s))=sqrt(N/(N+2)).

This realization concerns the white perpendicular two-face and its
finite corner root data. It does not supply all other walls of a
compact n-dimensional Coxeter polytope.

## 5. Arithmetic conjugate check and absence of affine principal blocks

Write r=sqrt(5). Under the nonidentity embedding put

    x'=(1-r)/4<0,     beta'=(3-r)/4 in (0,1).

The conjugate P_+ is positive semidefinite of rank one: both its
diagonals are positive and (1) remains an identity. For the conjugate
P_-, its first diagonal is

    x'+2*s=beta'+1/N>0.

Its determinant is

    [11*N+14-(5*N+2)*r]/(8*N)>0.

At N=17 the numerator is 201-87*r, positive because
201^2-5*87^2=2556>0. At N=18 it is 212-92*r, positive because
212^2-5*92^2=2624>0. Both rational terms being compared are positive.
Thus the full nonphysical Gram matrix is positive semidefinite with
n+1 positive eigenvalues and a one-dimensional kernel.

For norm-two roots the diagonal entries are 2 and all old simple
products are -1. The two exceptional entries are

    -2*x=-(1+r)/2,
    -2*beta=-(3+r)/2.

Both belong to the ring of integers of Q(r). Thus every Gram entry,
and consequently every cyclic product, is integral. The entry and
cyclic fields are Q(r). This checks the local integrality and
conjugate-form conditions without treating arbitrary unit-normal
entries as ground-field elements.

The full kernel relation, before white projection, is

    a+c+e+l1+...+lL-2*x*(b+d)=0.                 (3)

Every coefficient is nonzero under both real embeddings. Because
the conjugate Gram matrix is positive semidefinite with this
one-dimensional full-support kernel, every proper principal
conjugate submatrix is positive definite. Indeed a kernel vector of
a singular principal submatrix, extended by zeros, would be a kernel
vector of the full positive semidefinite matrix, contradicting the
full support in (3). Every proper principal determinant is therefore
nonzero, physically as well as at the conjugate. In particular there
is no affine principal subdiagram to remove the model.

Choosing any independent n+1 roots gives a quadratic space over
Q(r) of Lorentz signature at the physical embedding and positive
definite at the other embedding. This space is anisotropic over
Q(r): a nonzero rational-field isotropic vector would remain a
nonzero isotropic vector at the positive-definite embedding. These
facts describe admissible arithmetic ambient data. They do not
prove that this finite root set extends to a finite-volume chamber.

If such a profile occurs in an actual arithmetic polytope, its actual
ground field must be exactly Q(r). The old D_n plus c block is
physically independent, with determinant ratio 1-q_c<0 and

    q_c=1-s+2*N*(x+s)^2 in Q(r).

Any nonidentity embedding of a larger totally real field fixing Q(r)
would preserve that negative ratio, contrary to arithmetic conjugate
positive definiteness. This explains why unrelated field-specific
exceptional profiles cannot automatically be combined with this one.

## 6. What a facewise refund would have to pay

For the baseline used in the compact arithmetic manuscript,

    B(C;i,j)=1[distance_C(i,j)<=2]+the stated H4 extras,

the four corner values in this model are exactly (1,0,1,0). The
two units are at the opposite full-D_n fork/whisker angles. The
other corners mark different whole factors A1 and A_(n-1), so they
have zero baseline. No H4 component is present. Deleting both units
would leave zero instead of the required quadrilateral payment 1.

This is a baseline statement. The complete existing rule may have
additional bonuses at the two cross-component angles; their actual
coefficient/context must be checked before alleging a failure of
that complete rule.

For a useful quantitative formulation, suppose a diagram-invariant
proposal refunds delta from every full-D_n fork/whisker unit and
adds a cross weight t at the marked second-from-end node of
A_(n-1)+A1. The two zero corners in the model belong to that same
rooted type, so its quad payment would be

    2-2*delta+2*t.

Any facewise proof retaining this local model therefore needs

    delta<=1/2+t.                                (4)

This is a necessary charge row for the local-proof inventory, or for
an actual polytope containing this face. It is not a statement that
every abstract local model occurs in a global compact polytope.

At an entire A_(n-1)+A1 vertex the relevant path orbit has two nodes,
so this compensation adds total cost 2*t. Its bare baseline is
2*n-5, equal to 35 at n=20. Thus t>=1/2, needed by the naive
delta=1 proposal, would raise that vertex cost to at least 36 unless
another charge is refunded there. Moving payment to this cross orbit
alone does not establish the desired uniform target 35.

## 7. At most one golden exceptional face at a given D_n vertex

This section is a genuine actual-wall exclusion. It assumes a
compact polytope containing the specified wall profiles; it does not
assume their global existence.

At the old whole D_n vertex, the golden face at (a,b) requires a new
wall c_b with old contacts

    a.c_b=-x,     e.c_b=lL.c_b=-1/2,     b.c_b=0,

and all other old contacts zero. The analogous golden face at (a,e)
requires a wall c_e with b,e interchanged. These are distinct walls.

Let G be the unit old D_n Gram matrix, let v_b,v_e be their old
contact rows, and put q=v_b*G^-1*v_b^T. Swapping the old whiskers
is a diagram automorphism, so v_e*G^-1*v_e^T=q. Their difference is

    v_b-v_e=(e_b-e_e)/2.

The old inverse entries are

    (G^-1)_bb=(G^-1)_ee=n/2,
    (G^-1)_be=(n-2)/2.

Therefore

    (v_b-v_e)*G^-1*(v_b-v_e)^T=1/2,
    v_b*G^-1*v_e^T=q-1/4.                       (5)

Both walls lie strictly outside the old finite vertex. In the
one-dimensional timelike perpendicular of the old positive space,
their components thus have the same time orientation and square
-(q-1). Their actual scalar product is forced to be

    c_b.c_e=(q-1/4)-(q-1)=3/4>0,

contrary to acuteness of distinct outward walls. The two profiles
cannot coexist at one vertex.

More generally, if two proposed actual walls have old contact rows
related by the whisker swap, agree on all other old roots, and their
b,e contact difference is h, the same argument gives product
1-h^2. Any h with |h|<1 excludes coexistence. This is a local
capacity principle for exact wall profiles, independent of arithmetic.

The capacity result does not yet provide a complete activation rule.
Other possible faces at the two angles, and every triangle consumer,
must still be classified. It does show that simply summing the two
golden exceptional costs at a single vertex would overcount.

## 8. Status and next precise obligation

The far-end rooted continuation must remain in the full-D_n
classification. The golden model is a verified local survivor of
finite-corner, Lorentz, compact-two-face, arithmetic-integrality,
conjugate-positivity and affine-principal tests. It blocks an argument
that deletes both opposite fork/whisker payments using only those
tests. It does not disprove a global theorem about actual compact
polytopes, because global completion is not supplied.

A next viable refund argument must either exclude this profile by
additional global geometry, use an already present cross payment,
or balance an activation rule with the exact capacity in Section 7.
Any such rule must then cover all other selected triangle and
quadrilateral consumers and all whole-forest vertex costs. No
dimension improvement or new executable/formal verification is
claimed here.
