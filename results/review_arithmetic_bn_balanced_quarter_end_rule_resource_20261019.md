# Independent full review of the balanced B_n quarter-unit rule

19 October 2026. **Written PASS; no flaw found in the specified whole-
forest and eight-field theorem.** Reviewed source:
`results/arithmetic_bn_balanced_quarter_end_rule_20261019.md`, SHA256
`514e41e5d5f17c61350c9031af47c72c4671b115228d5e98e496d009bcad9b6d`.
I read its complete source, the full-B_n distance-two refund and its
baseline charge input, and the supporting four-corner countermodel
`results/arithmetic_bn_four_corner_lambert_survivor_20261019.md`, SHA256
`1e6154b98ef0dd6292c3c119a74b27472199f20fea9de480b1c4507769972c3b`.
All sources were preserved unchanged. No numerical job, finite census,
optimizer, proof script, Lean process or build ran in this review.

The new quarter rule has complete face-consumer coverage within its
stated field scope. Combined full B20 and B19+A1 costs are141/4,
not35. Neither a new global allowance nor a dimension ceiling follows.

## 1. Fixed white support seals every affected corner

Only a full-B_n simple-end adjacent unit decreases. For any quad
carrying that role, the actual n-2 white normals form the connected
B_(n-2) path, including its remote4-bond. At every corner these same
white normals plus the two marked black normals are the entire rank-n
finite diagram: compact Coxeter vertices are simple. There are no
unmentioned spectator factors in this inference.

Each black root must individually extend that fixed W to a finite
diagram. The connected component containing this rank-at-least17
B path cannot be exceptional, simply laced, or dihedral. It is
another B path. Its terminal4-bond cannot acquire a neighbor at
the far end, another double bond, a branch or a cycle. Thus its
only possible contact is a simple bond at the simple endpoint of W,
or no contact. These are exactly N and I.

This is a complete finite-extension argument, not an assumption
that remote white contacts vanish. The finite B/C Coxeter paths
are the same diagram type under unit normalization.

## 2. Every finite corner on this quad pays one quarter or more

For I,N, a nonright black bond makes the entire finite corner the
B_n path. It must be simple, with the isolated-from-W root as the
new simple endpoint. The marked pair is exactly the source role,
whose new weight is1/4. If the black bond is right, the entire
forest is B_(n-1)+A1, and the marked connected root is precisely
its simple endpoint. That is exactly the cross recipient of1/4.

For I,I, a right black pair gives the entire B_(n-2)+2A1 recipient,
paying1/4. A nonright pair is an entire separate I2(m), retaining
its original marked unit, plus the existing high-label unit when
applicable. It never qualifies as a reduced whole-B_n role.

For N,N, two roots attach at the same W endpoint. If they are
orthogonal this gives a branch in the large component carrying the
old4-bond; if not, it gives a cycle. Neither is a finite large B
path, and no exceptional finite component has the required rank.
Thus this case is impossible.

The list covers every corner independently of its original payment,
of which other corners are selected, and of the opposite products.
Summing four lower bounds1/4 proves the quad payment1. In particular
the alternating I,N patterns, which defeated an isolated-pair-only
repair, are covered by the cross recipient.

## 3. Other faces and simultaneous compositions

Triangles still receive two at the inherited entire high-label I2
marked corner. A whole large B_n source role cannot be that corner,
and the two new recipients only increase weights. A quad with no
new reduced role has no decrease under this three-part rule; a
quad with one is covered in Section2. All weights remain nonnegative,
so polygons with at least five sides are also covered.

The separate B_n first/third deletion has white graph
A1+B_(n-3), whereas the present adjacent reduction has connected
white graph B_(n-2). A two-face has one fixed set of white facet
normals and one fixed induced Gram graph, so the two source roles
cannot occur at different corners of the same quad. The positive
coupons can only help the separate distance-two proof.

Any quad affected by either B source contains a white4-bond at
every corner. No such corner can be a whole simply-laced D_n,
so the separate D_n whisker deletion cannot remove another payment
on these quads. Triangles after all modifications retain their
untouched high-label I2 payment. This composition inherits the
separately reviewed D_n theorem, rather than extending its scope.

## 4. Exact per-vertex costs, including recipients

The source role and each recipient role are unique unordered marked
pairs at the specified whole forest. There is one B simple endpoint,
one isolated root in B_(n-1)+A1, or one unordered isolated pair in
B_(n-2)+2A1. No multiplicity or spectator cost is omitted.

At rank20 the correct costs are

    B20:       37-1-3/4 =141/4;
    B19+A1:    35+1/4   =141/4;
    B18+2A1:   33+1/4   =133/4.

The full A20 and D20 costs remain37 under the listed rules. The
remaining quarter at BOTH B20 and B19+A1 is therefore a real
obstruction to allowance35, not something absorbed by strictness
of the global average. No outside-factor extension is inferred.

## 5. Independent check of the supporting four-corner model

The cited raw Q(sqrt(3)) Gram has all four black norms2, cyclic
products-1, opposite products-(1+r),-(2+2r), and B,D white-end
contacts-1. Eliminating raw B_(n-2) subtracts1 from each of B,D's
norms and1 from their mutual product. The residual norms are2,1,2,1.
The two antisymmetric eigenvalues are2+2x,2+2y, both positive at
both embeddings. The symmetric block is exactly

    [2(1-x)  -2 ]
    [  -2   -2y],  x=(1+r)/2, y=1+r.

The identity y(x-1)=1 gives determinant0. Its physical nonzero
eigenvalue is negative; its conjugate block is positive semidefinite
of rank1. With the white block this gives the stated physical
inertia and positive conjugate quotient.

I checked the full radical with A=C=2, B=D=1-r and every white
coefficient2(1-r). For example the B row reduces to
-4-(2+2r)(1-r)=0, and the first white row also vanishes. Its
conjugate is strictly positive and full support. The integral
Cartan ratios and a triangle cyclic product generating K are correct.

Each adjacent black pair with W is a finite full B_n path whose
marked pair is its simple-end adjacent role. The connected acute
Lorentz Gram supplies a strict timelike center by its negative
Perron eigenvector; all-wall projection makes the four corners
actual finite vertices and W an actual compact quad face. Each
intrinsic angle is pi/4. It is correctly distinguished from a
three-right-angle Lambert quadrilateral despite the historical filename.

This validates the local uniform1/4 lower constraint when no other
payments are added on those same four corners. It does not show
that a globally compact arithmetic chamber must contain this pattern,
nor that its surrounding halfspace system has finite volume.

## Review limits

The balanced rule is a complete written partial consumer theorem
using the inherited eight-field charge, finite classification,
compact acute-face correspondence and separately reviewed deletions.
It has not been implemented or executed as a new global weight
certificate here, nor checked in Lean. Its source makes no improper
dimension-bound or finite-volume promotion.
