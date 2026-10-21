# An integral Q(sqrt(3)) quadrilateral with four full B19/B20 corners

19 October 2026. **New exact local countermodel, awaiting independent
review.** No scientific job, proof script, optimizer, Lean process or build
was run. Existing evidence is unchanged. The model is not a compact
finite-volume hyperbolic n-polytope construction.

The earlier simple-end adjacent B_n countermodel had two corners of
type B_(n-2)+2A1. It suggested paying for an adjacent-end deletion by
adding half a unit at their marked isolated pair. The present model
shows that this proposed recipient does not cover all consumers:
there can instead be FOUR full B_n corners and no isolated-pair corner.

## 1. An explicit reflector Gram

Let n=19 or20, put r=sqrt(3), and use K=Q(r). The construction works
at every n>=4. Let W=(R1,...,R_(n-2)) be a raw integral B_(n-2) path:
all squared norms are 2 except R_(n-2)^2=1, and every consecutive
product is -1. All its other products are zero.

Add four roots A,B,C,D of squared norm 2. Their cyclic products are

    (A,B)=(B,C)=(C,D)=(D,A)=-1.

Their opposite products are

    (A,C)=-(1+r),          (B,D)=-(2+2r).

The white contacts are

    (B,R1)=(D,R1)=-1,

and every other black-white contact is zero. All entries are in
O_K=Z[r], every root norm is totally positive, and every Cartan ratio
2(G_i,G_j)/G_i^2 is integral. For black-black pairs it is just the
displayed inner product; on the old far B bond it is -1 or -2.
Thus the reflections preserve the O_K root module in the nondegenerate
quotient described below.

The complete old corner W+A+B is B_n, with A its simple endpoint
and B its neighbor. C,D add no finite bonds other than labels 3 or 4.
After unit normalization the opposite products are -x,-y, where

    x=(1+r)/2>1,             y=1+r>1.                (1)

In particular neither opposite product is an unnoticed finite angle.

## 2. Rank, inertia, radical and conjugates

The inverse entry at R1 of the raw B_(n-2) Gram is 1. Eliminate W.
The residual norms are A^2=C^2=2 and B^2=D^2=1. The projected
B,D product is -(3+2r), while all four cyclic products remain -1.
In unit-normal notation this is equivalently

    A_bar^2=C_bar^2=1,       B_bar^2=D_bar^2=1/2,
    A_bar.C_bar=-x,          B_bar.D_bar=-y-1/2,
    each cyclic cross product =-1/2.

Split into the symmetric and antisymmetric combinations of A,C and
B,D. The two raw antisymmetric eigenvalues are

    2+2x,       2+2y.

They are positive physically and after conjugation. The raw symmetric
two-by-two matrix is

    [ 2(1-x)     -2  ]
    [   -2       -2y ].                              (2)

The identity y(x-1)=1 makes its determinant zero. Its trace is
negative at r>0 and its diagonal becomes positive after r -> -r.
Consequently the full physical Gram has inertia

    (n positive, 1 negative, 1 zero),

and its conjugate is positive semidefinite of rank n+1.

An exact radical vector has coefficients

    A=C=2,
    B=D=1-r,
    each white root R_i=2(1-r).                       (3)

Multiplication by the Gram verifies every row. This vector has full
support, and all its conjugate coefficients are positive. Since the
conjugate Gram is positive semidefinite with one-dimensional radical,
every proper principal conjugate Gram is positive definite. Hence
no physical proper principal determinant vanishes either: its
entries lie in K, and conjugation would preserve a zero determinant.

The nondegenerate K-quotient by the radical is Lorentz at the physical
embedding and positive definite at the other embedding. The integral
Cartan ratios prove that the reflections preserve its root module.
The cyclic products generate exactly K; for example the triangle
A,B,C has a cyclic product proportional to 1+r. Thus this is an
explicit admissible integral local reflector model over one of the
eight high-label fields, rather than only a real Gram matrix.

No assertion about finite index, finite covolume or completion of this
reflection subgroup is made.

## 3. All four corners are actual finite B_n corners

Every adjacent black pair together with W is a positive definite
full B_n path:

    A-B-R1-...,
    C-B-R1-...,
    C-D-R1-...,
    A-D-R1-... .

The far white 4-bond is unchanged. Thus all four marked cyclic pairs
are the SIMPLE-END ADJACENT role of a whole B_n vertex. Each receives
exactly one unit from the distance-two baseline, and none is an
entire high-label I2 component. The four ordinary payments are

    (1,1,1,1).                                       (4)

For completeness the compact quad is actually realized in an acute
halfspace system. The physical Gram is connected with nonpositive
off-diagonal entries and has a negative eigenvalue. A strictly positive
eigenvector for its smallest eigenvalue supplies a timelike root
combination whose products with every normal are negative. Hence
the n+2 halfspaces have a strict hyperbolic interior point.

Each of the four positive-definite n-root corners gives an actual
finite vertex by the all-wall acute projection lemma, with every
other inequality strict. Every root is therefore an actual facet.
The white W intersection is an actual two-face. Its four consecutive
vertices are the above finite corners, while the opposite pairs are
ultraparallel by (1). The four remaining halfplanes consequently
bound the actual compact quadrilateral in that white Lorentz plane.

At each corner the projected raw norms are 2 and 1 with product -1,
so every intrinsic angle is pi/4. It is a compact four-acute-angle
quadrilateral. Strictly speaking it is not a Lambert quadrilateral,
which would have three right angles; the filename retains the
research lane's name, while this distinction is explicit here.

The surrounding n-dimensional chamber's volume and global combinatorics
are unproved. In particular this is not a new finite-volume example
in dimensions 19 or20.

## 4. Consequences for proposed balanced refunds

Deleting the simple-end adjacent baseline unit at every full B_n
vertex would delete all four payments in (4). Adding a half-unit
only at marked B_(n-2)+2A1 isolated-pair corners cannot repair this
face: it has no such corner.

The separately proved simple-end DISTANCE-TWO deletion does not
affect (4). The quad's fixed white graph is connected B_(n-2);
the white graph for that distance-two mark is A1+B_(n-3). The two
roles cannot occur on the same face. Nor can a full D_n whisker
deletion occur on this face, since every corner contains the white
4-bond and D_n is simply laced.

Any universal rule retaining the same adjacent weight at these four
whole-B_n marked roles therefore needs weight at least 1/4 there,
unless it adds another payment on these same four actual corners.
This is a constraint on a local charge argument, not a claim that
all global compact chambers must realize the pattern.

For example, weight 1/4 on the adjacent mark, together with a
3/8 bonus at the isolated-pair role of B_(n-2)+2A1, would pay this
four-B_n face and the earlier one-B_n Lambert face. Their B20
vertex cost after the separate distance-two deletion would be
37-1-3/4=35.25, and the recipient forest cost would be 33.375.
This calculation does NOT validate that proposed rule: quads with
two alternating white-attached black roots and only some adjacent
simple bonds must still be classified. In particular the two
displayed countermodels are not an exhaustive consumer list.

No new weight source, complete balanced-refund theorem, allowance35,
or lower dimension ceiling is asserted by this report.
