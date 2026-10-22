# Glued golden E7/A2 host: a complete ambient D6 edge and its first mixed wall

22 October2026. **New complete SOURCE-ONLY CANDIDATE, pending full
different-author mathematical review and root acceptance.** This proves
one actual ambient edge at the written inputs below. It does not prove
finite chamber volume, reflectivity, nonreflectivity, a chamber period,
or a compact H9 example. No mathematical worker, import, enumeration,
symbolic process, test, Lean, TeX, compiler, plan or Git operation ran.
No frozen source, receipt, paper, rule, guide or bound was changed.

I authored the first two construction sources and this continuation.
This is new mathematics, not an independent review of those inputs.
The two earlier sources are now root accepted. Their compact triangle
controls an E7 two-face; it does not settle any ambient edge obtained
by releasing an E7 wall. The original module is kept throughout.

## 1. Exact inputs, original module and actual starting germ

The following direct bindings were checked by byte reads only.

| Source | SHA256 | Exact interface used |
|---|---|---|
| results/compact_h9_glued_e7_a2_ramified_five_module_and_first_edge_candidate_20261022.md |33281068c43fda60f3356dc2cfcc226fcc72aae03d818ce49598e14382deb17c|Admissible ORIGINAL module, primitive norm pool, full initial germ and actual first edge to g|
| results/compact_h9_glued_e7_a2_second_edge_and_full_triangle_face_candidate_20261022.md |a7e4bceb77715cbeabd8858fcd8146c2786736bae37c9cbc6505f10221db7bc3|The same actual chamber and its three known ternary facet walls; no ambient-volume conclusion|
| results/review_glued_e7_a2_original_module_first_two_edges_and_triangle_root_acceptance_20261022.md |27eaff5edc42bae509302052667ddc842abdfd52b062596a87f8792a649b8968|Acceptance of precisely those original-module and actual-germ interfaces after the complete outside reviews|

Standard finite Coxeter classification, simplicity of finite vertices
of a reflection chamber, and local finiteness of the admissible integral
mirror arrangement remain explicit structural inputs. No named arithmetic
field list, earlier551 isometry, or finite-prefix closure is imported.

Write phi^2=phi+1, K=Q(phi), O=Z[phi], and b=2+phi. On square2 roots
f1,...,f7 the E7 edges, with product-1, are

    12,23,34,45,56,37.

The ORIGINAL module is the orthogonal sum of this E7(O) and the ternary
module on e1,e2,g with Gram

    [[2,-1,1],[-1,2,0],[1,0,-2phi]].                 (1)

Put n=(phi-1)g+2e1+e2. The three vectors g,e2,n are orthogonal, with

    g^2=-2phi, e2^2=2, n^2=2b,
    g.e1=1, e2.e1=-1, n.e1=b.                     (2)

The accepted actual chamber C has at g the FULL outward germ

    -f1,...,-f7, -e2, n.                          (3)

The old wall-e1 is also an ACTUAL facet of this same C. Its product
with g is-1. The complete necessary normalized primitive norm pool is

    r^2=2d, d in{1,2,b,2b}.                      (4)

This is the previously proved exponent bound d|2(1+3phi) and golden
unit normalization for this exact module. The original module is even
and O is a PID. For a primitive norm2d reflector, its original covector
divided by d must be integral. All four classes are retained below.

## 2. Release f1: an original integral finite upper

The highest E7 root in this labeling is

    theta=2f1+3f2+4f3+3f4+2f5+f6+2f7,
    theta^2=2, theta.fi=delta(i,1).                (5)

These identities follow directly from the displayed E7 edges. Delete
f1; the retained f2,...,f7 form D6, with fork leaves f2,f7, center f3,
and long arm f4,f5,f6. Its perpendicular in E7 is K*theta.
Release-f1 at g while retaining this D6 and-e2,n. The outgoing line is

    Q(tau)=g+tau*theta, 0<tau<sqrt(phi).           (6)

Its released product is-tau; its square is-2phi+2tau^2. Every retained
wall has product0 and the other known actual wall-e1 has product-1.

Consider the original vector

    R=phi*theta+g, R^2=2.                         (7)

It is integral and primitive, since its original g coefficient is1.
A norm2 integral vector in an integral module has an integral reflection:
the reflection subtracts (v.R)R. Its products with the known actual walls
are

    R.(-fi)=-phi*delta(i,1), R.(-e1)=-1,
    R.(-e2)=R.n=0.                               (8)

Thus it is acute to every known facet, with native5 against-f1 and
native3 against-e1. It cuts(6) at the finite CLOSED parameter tau=1:

    Q(tau).R=2phi*(tau-1),
    q=g+theta, q^2=2(1-phi)<0, q.g=-2phi<0.       (9)

At this stage R is an integral mirror supplying a genuine finite upper;
it is not yet assumed to be an actual facet. The full initial germ
supplies an initial actual edge interval. Local finiteness on the compact
segment up to q implies that an actual first wall occurs at0<tau<=1.
Its endpoint is finite, so adding that simple wall to the retained
eight walls gives a positive definite finite simple diagram.

## 3. Every first-simple wall contacting the retained D6

Let r be any such actual first simple wall, in normalization(4).
A finite irreducible diagram containing D6 and one further node is D7
or E7. All its bonds are native3; a branched diagram of this rank cannot
be B/C, F4 or H3/H4. It therefore requires d to be a K-square. Among
the four representatives in(4), only d1 is a square:2 is inert and b
has odd valuation at the ramified prime5. A new node can contact only
one D6 vertex, since two contacts create a circuit. The complete
possibilities are

| Helper in D6 | New connected type | D6 projection square |
|---|---|---|
|f2 or f7, the fork leaves|E7|3/2|
|f6, the long leaf|D7|1|

An internal helper either creates degree4 or a second branch point,
and is not a finite diagram. A simultaneous contact with the retained
old-e2 can extend the long-leaf D7 to D8, with one further native3
bond. For a fork helper it instead gives arms(1,3,3), the affine E7
diagram, so is impossible. There is no E8 option: the D6 has arms
(1,1,3), and the two added nodes, one of which is an old A1, cannot
produce E8's arms(1,2,4). Contact with both old A1 nodes creates another
branch point. A contact with n for d1 would require native10, the only
golden field-compatible finite angle with its norm2b; no finite diagram
containing D6 allows this bond. This exhausts all attachments of either
retained A1, as well as the pure D7/E7 cases.

Let lambda_j be the D6 fundamental weight at the indicated helper.
The E7 projection of r is a*theta+lambda_j. Write the remaining
ternary part x*g+y*e2+z*n. Since r is not incident at g, g.r<0 and
x>0. A collision with(6) gives a>0 and tau=phi*x/a<=1, hence
a>=phi*x. Its norm2 now implies

    2phi*x^2=2a^2+lambda_j^2+2y^2+2b*z^2-2,
    a^2<=phi^2*(1-lambda_j^2/2).                 (10)

For the fork helpers this gives a<=phi/2<1. But the old ACTUAL wall
-f1 must also be acute to r. The D6 fundamental weights obey

    f1=theta/2-lambda_2,
    lambda_2^2=3/2, lambda_7.lambda_2=1.

Thus r.f1=a-3/2 for helperf2 and a-1 for helperf7, both negative.
This contradicts r.(-f1)<=0. These fork cases are excluded using the
original released wall, not a conjectural germ symmetry.

For the long helper an even stronger original-coset obstruction applies.
Directly in the retained simple-root basis,

    lambda_6=f6+f5+f4+f3+(f2+f7)/2.              (11)

The f6 coordinate of a*theta+lambda_6 is a+1, forcing a in O if the
whole r is original integral. Its f7 coordinate would then be2a+1/2,
which is not in O. Ternary coefficients cannot clear an E7 denominator:
the ORIGINAL module is a direct sum at E7. This excludes both the pure
D7 and mixed D8 cases. Consequently EVERY possible first wall is
orthogonal to the retained D6. No first-wall D6 contact is omitted.

## 4. Complete twelve profiles orthogonal to D6

Its E7 part is now a*theta with a in O, since theta's f6 coefficient
is1. Use the orthogonal ternary coordinates(2) to write

    r=a*theta+x*g+(k/2)e2-(j/2)n,
    k=r.e2>=0, j=-r.n/b>=0.                     (12)

The j integrality follows because the ACTUAL norm2b reflection n is
integral and primitive: r.n/b is in O. Again x>0,a>0 and the collision
is tau=phi*x/a. The released-wall signs are retained: r.f1=a>=0 and
r.e1=x-k/2-j*b/2>=0. The second inequality can only shrink the domain
below, so it need not be used to discard any additional candidate.

The complete golden finite cyclic alphabet is

    4cos^2(pi/m) in{0,1,2,phi^2,3,b},
    m=2,3,4,5,6,10.                             (13)

Here the field-degree/native-angle argument is the same elementary one
used in the accepted first source; native8/12 lie in different quadratic
fields. The two finite pair invariants in(12) are k^2/d and b*j^2/d.
Exact square membership in K gives the following COMPLETE twelve
profiles before original membership:

| d | (k,j) |
|---|---|
|1|(0,0),(1,0),(phi,0),(0,1)|
|2|(0,0),(2,0)|
|b|(0,0),(b,0),(0,1),(0,phi)|
|2b|(0,0),(0,2)|

For example the d1,j1 and db,kb cases are native10, not orthogonal.
The finite pair alphabets individually also allow simultaneous nonzero
contacts when d1 or db. Each would give a three-node chain containing
native10 and at least native3. Its Gram determinant is nonpositive:
the sum of the two cosines squared is at least(b+1)/4>1. Thus these
simultaneous rows are correctly excluded. Parallel/dotted pairs cannot
occur in a positive definite incident Gram. These observations prove
completeness without assuming r is orthogonal to either old A1.

Expanding(12) on the ORIGINAL ternary basis gives

    r_B=(x-j*(phi-1)/2)g-j*e1+(k-j)e2/2.        (14)

The condition (k-j)/2 in O eliminates all odd/phi rows in the table.
Exactly six necessary profiles survive original membership:

| d | (k,j) | D=d-(k^2+b*j^2)/4 |
|---|---|---|
|1|(0,0)|1|
|2|(0,0)|2|
|b|(0,0)|b|
|2b|(0,0)|2b|
|2|(2,0)|1|
|2b|(0,2)|b|

In all six x is in O: for j2 its original g coefficient is x-(phi-1).
Their norm equation is

    a^2-phi*x^2=D.                              (15)

This is WHOLE original-module membership, not merely separate integral
coordinates in the nonintegral orthogonal basis g,e2,n.

## 5. Complete seven-value coefficient strip and finite hand table

The closed finite upper tau<=1 gives a>=phi*x. Thus(15), with x>0,
and its conjugate imply

    D<a^2<=phi^2*D, a'^2<=D'.                   (16)

Across all six rows D<=2b and D'<=2b'. Hence

    0<a<9/2, |a'|<7/4.                         (17)

Indeed phi^2*2b=6+8phi<19<81/4, using phi<13/8;
2b'=6-2phi<3<49/16, using phi>3/2. Write a=A+Bphi with
integers A,B. Since a-a'=B*sqrt5, (17) gives -1<B<3, using
25/4<3sqrt5 and7/4<sqrt5. Thus B is0,1 or2. The exact remaining
integer possibilities are

    B0: A1;
    B1: A=-1,0,1,2;
    B2: A0,1.

For B2, A2 already violates a<9/2. The ENTIRE seven-value list is

    a in{phi-1,1,phi,phi^2,2phi,b,phi^3}.        (18)

No unbounded Pell family or unexecuted coefficient census is hidden.
Substitution in the exact inequalities(16), including the conjugate
inequality, gives this complete hand table:

| D | Surviving a | Required x^2=(a^2-D)/phi |
|---|---|---|
|1|phi|1|
|2|phi|2-phi=(phi-1)^2|
|b|phi^2|2|
|2b|2phi|2|
|2b|b|b|
|2b|phi^3|5+phi|

The a1,D1 zero value x0 is expressly excluded because g.r<0 requires
x>0; it is not a new opposite. The other discarded a values violate
(16). The last four rows have no K-square x^2:2 is not a square in
Q(sqrt5), and N(b)=5 and N(5+phi)=29 are not rational squares.
These norm tests cover both possible signs of a proposed square root's
rational norm. None of the surviving positive-time rows has a zero root.

The remaining D1/D2 rows require a=phi. For D1 the orthogonal norm2
profile is precisely R, with positive x1. But the mixed d2,k2 profile
would instead be phi*theta+g+e2, of norm4. For D2 the independent d2
profile would be phi*theta+(phi-1)g, also of norm4. Both are integral
vectors, but neither is an integral reflector: the f1 covector entry
divided by d is a/d=phi/2, which is not in O. Primitivity is not silently
assumed to fail; reflection integrality itself fails. Thus all six
original-membership profiles have been resolved, and R is the UNIQUE
possible first actual simple wall before or at tau1.

## 6. Actual firstness, full new germ and precise remaining problem

The finite upper and local finiteness guaranteed a first simple wall.
The complete finite-white and original-module analysis proves that this
wall must be R at tau1. In particular R IS an ACTUAL facet of C.
Feasibility of the entire closed segment against all unknown original
mirrors follows by firstness, not merely from the known-wall products(8).

At the future finite point q=g+theta the outward incident roots are

    -f2,...,-f7, -e2, n, R.

Their positive definite Gram is

    D6 perpendicular[2]perpendicular[2b]perpendicular[2],
    determinant=32b>0.                          (19)

These are nine independent actual facets, so finite-vertex simplicity
gives the FULL actual germ D6+3A1. The released-f1 has product-1 at q;
old-e1 also has product-1. No hidden incident facet or E8 recognition is
asserted. The new mixed wall R has an E7 component, so the earlier
triangle calculation never implied all ambient walls were E7-orthogonal.

The achieved ambient walk is therefore an actual finite edge

    E7+[2]+[2b] at g -> D6+[2]+[2b]+[2] at q.

This advances the full original chamber outside its known E7 triangle.
It supplies no integral map between equivalent full outward germs, no
period, no upper bound on the number of other ambient facets and no
finite-volume conclusion. Further actual edges from q remain open.
No global bound, example, commensurability-wide result or execution
credit is claimed. Full mathematical outside review and root acceptance
are still required for this NEW edge theorem.
