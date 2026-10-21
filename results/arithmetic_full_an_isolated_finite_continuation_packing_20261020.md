# Isolated finite continuations of a whole A_n corner are three-spaced

20 October 2026. Paper proof only; no execution, census, optimizer, build
or formal verification. This is one uniform restriction for a complete
geometric class of continuations, not a bound on every incident edge.

Let P be a compact arithmetic hyperbolic Coxeter n-polytope, n>=13, with
an actual whole A_n finite vertex q. Normalize its incident simple roots
alpha_1,...,alpha_n to norm two, so their Gram matrix is the Cartan matrix
C of A_n. Write N=n+1 and

    d_i=(C^-1)_ii=i(N-i)/N.

This normalization is available over the arithmetic ground field: along
each3-bond, field-valued root rescalings make the two incident roots have
the same norm, and the path has no cycle obstruction. If their common
norm is b, replace the entire quadratic form by the totally positive
scalar2/b. No square root of b is required in the ground field.

An opposite endpoint wall c_i of the outgoing i-th edge is called
**isolated finite** here when its only nonzero product with these n old
roots is a finite Coxeter contact with alpha_i. Its index is called bulk
when d_i>2. A shared wall may also close another outgoing direction;
only its isolated CONTACT INDEX is counted in the theorem below, not
every direction that it closes.

**Theorem.** All bulk isolated finite continuations present at this one
vertex have the same label, in {5,8,12}. Their indices are at least three
apart. In particular, if M is the size of the FULL contiguous interval
defined by d_i>2, their number is at
most ceil(M/3). For arbitrary nonnegative index weights, their total weight
is at most the maximum weight of a subset with pairwise index separation
at least three. No integrality of equally normalized raw root pairings
is assumed.

The precise admissibility windows are

    label5:  d_i<2phi^2,
    label8:  d_i<2+sqrt(2),
    label12: d_i<4+2sqrt(3),       phi=(1+sqrt(5))/2.       (1)

## 1. Correct field normalization and the finite-label list

Put x=4cos^2(pi/m) for an isolated contact of label m. The root on that
wall can be normalized OVER THE GROUND FIELD to

    c_i^2=2x,       <c_i,alpha_i>=-x,
    <c_i,alpha_j>=0 for j!=i.

Indeed, if an initial field-valued root has norm b and product g with
alpha_i, then x=2g^2/b. Multiplying it by -x/g gives exactly this
normalization. Norm two would instead require a square root outside the
ground field at labels8 and12; it is not used for the field argument.

The full-rank Lorentzian star has Gram matrix over k=Q(x). Its Schur
complement after C is x(2-x*d_i)<0. Its cyclic field is exactly k, since
the only nontrivial cycle generator is x. The actual arithmetic ground
field cannot properly extend k: a nonidentity real embedding fixing k
would preserve this full Lorentzian Gram matrix, contradicting the
required positive definiteness at every nonphysical embedding. Thus the
ground field of P is k. This uses the ordinary arithmetic conjugate-form
criterion, not a claim that arbitrary Coxeter groups are arithmetic.

For every nonphysical conjugate x' of x, positive definiteness gives

    x'(2-x'*d_i)>0,       hence 0<x'<2/d_i<1.             (2)

Here is an elementary complete cyclotomic restriction. If m is odd and
m>=7, the coprime exponent a=2 lies in (1,m/3). If m is even, 3 does not
divide m, and m>=10, use a=3. If m=6k>=18, use a=2k-1 unless k=2 mod3;
in that case use a=2k-3. In both cases 1<a<m/3 and gcd(a,m)=1: the
possible common factors divide3 or9 respectively, and the chosen
congruence excludes them. Such an exponent produces a nonphysical
conjugate x_a=2+2cos(2pi*a/m)>1, contradicting (2).

The only surviving labels are therefore 3,4,5,6,8,12 (label2 has x=0
and cannot be an opposite endpoint wall). The conjugates of x at the
three irrational survivors are

    5:  x=phi^2,       x'=phi^-2,
    8:  x=2+sqrt(2),   x'=2-sqrt(2),
    12: x=2+sqrt(3),   x'=2-sqrt(3).

Equation (2) gives (1). These are three structural alternatives, not a
restriction of the theorem to a preselected list of arithmetic fields.

Labels3,4,6 are excluded by compactness. Since d_i>2, both old arms at
alpha_i have at least two nodes. For label4 the center alpha_i with
its two old neighbors and the new4-neighbor is the affine star with
positive null coefficients (2;1,1,sqrt(2)). For label6, the new root,
alpha_i and either old neighbor form the affine6--3 path, with positive
null coefficients (sqrt(3),2,1). For label3, if both old arms have at
least three nodes, the new leaf and three nodes of each arm form the
affine simply-laced star with arm lengths (1,3,3), whose positive null
coefficients are central4, short arm2, long arms (3,2,1). Otherwise one
old arm has length two and the other has at least five; use arm lengths
(1,2,5), with central6 and arm coefficients (3),(4,2),(5,4,3,2,1).

These affine null coefficients use the usual equal-norm real Gram
convention (each root has square two), obtained by positive real
rescaling solely for this geometric affine-subdiagram check.

For clarity, any such affine subdiagram is impossible in a compact
acute Coxeter chamber: its positive null combination of outward normals
is a nonzero future null vector feasible against every other wall by
acuteness. It would be an ideal point of P. Thus only5,8,12 remain.
Their fields are distinct. The full-rank argument above also proves
that different surviving labels cannot coexist at this same corner.

## 2. The common timelike residual fixes every pair product

Fix the common surviving x. Let omega_i be the C-dual spatial weights
at q, and write c_i=H_i*q-x*omega_i. The positive physical residual is

    H_i=sqrt(x*(x*d_i-2)).

The sign is positive because q lies strictly inside every opposite wall.
Hence for i<j, with C_ij^*=i(N-j)/N, the forced twice-unit-normal product is

    R_ij=<c_i,c_j>/x
         =x*C_ij^*-sqrt((x*d_i-2)*(x*d_j-2)).            (3)

Distinct isolated supports give distinct actual walls, so R_ij<=0 by
all-wall acuteness. An identical wall cannot be counted twice. Since
both roots have norm2x over the ground field, R_ij belongs to that field.

If the pair meets, R_ij=-2cos(pi/l) for an integer l>=2. The real
cyclotomic degree of 2cos(pi/l), together with the three quadratic
fields, leaves the following finite possibilities:

    Q(sqrt(5)): {0,-1,-phi},
    Q(sqrt(2)): {0,-1,-sqrt(2)},
    Q(sqrt(3)): {0,-1,-sqrt(3)}.                         (4)

If the pair is disjoint, compactness gives R_ij<-2. No arithmetic
integrality assumption is needed for either assertion.

## 3. Adjacent indices and distance-two indices are impossible

Arithmetic-geometric mean in (3) gives

    R_ij>=2-x*((d_i+d_j)/2-C_ij^*).                     (5)

For j=i+1 the parenthesis is (N-1)/(2N). Since x<4, the right side
is positive. This contradicts acuteness.

For j=i+2, put C_*=i(N-i-2)/N and delta=2-4/N. The parenthesis in
(5) is (N-2)/N, so R_ij>2-x. For label5 this leaves only R_ij=0
among (4); for labels8,12 it leaves only0 or-1. Disjointness is also
excluded, since 2-x>-2.

If a remaining value is R_ij=r, squaring (3), necessarily, gives

    2C_*[x^2-(2-r)x]-2x*delta+4-r^2=0.                 (6)

Both C_* and delta are rational. Reduce x by its quadratic equation:

* Label5 has x^2=3x-1. For r=0, the constant and x coefficients of
  (6) force C_*=delta=2, impossible since delta<2.
* Label8 has x^2=4x-2. For r=0 they force C_*=1, delta=2,
  again impossible. For r=-1 they force C_*=delta=3/4,
  giving N=16/5, not an integer.
* Label12 has x^2=4x-1. For r=0 they force C_*=2, delta=4,
  impossible. For r=-1 they force C_*=delta=3/2 and N=8,
  excluded already by n>=13. Even without that rank hypothesis,
  C_*=i(6-i)/8<=9/8 when N=8, so it cannot equal3/2.

This excludes distance two and proves the asserted three-spacing.

## 4. What this does and does not count

The indices with d_i>2 form a contiguous interval of the A_n path; a
three-spaced subset of its M consecutive positions has size at most
ceil(M/3). This M is not the number of positions left after additional
label-window restrictions.
For weighted use, the exact combinatorial envelope is the maximum of
sum_i w_i over such subsets. It can equivalently be defined by
F_i=max(F_(i-1),w_i+F_(i-3)) at an available index, and F_i=F_(i-1) at an
unavailable index, with the initial values zero. The indices retain their
original path positions; unavailable positions are not compressed. This
is a definition, not an executed optimizer.

This theorem does not cover one-sided contacts to white components,
two-component bridges, exceptional small-component contacts, isolated
dotted continuations, or indices d_i<=2. A single wall may close two old
edges and must be counted as such in a full incident-edge argument,
even if its continuation at one of those indices is isolated. The known
[6,3,3] short-end family is another shared wall example, with two nonzero
old contacts. No
uniform dimension exclusion follows from this isolated subtheorem alone.

The geometric interfaces used here are standard all-wall acuteness,
the arithmetic conjugate-positivity/cyclic-field criterion and finite
Coxeter angle labels. The An inverse formula, cyclotomic exclusion,
affine null combinations, common-residual pair formula and spacing
equations are derived explicitly above. Existing direct-descent and
trace-lift barriers remain untouched:
`results/first_type_projection_and_trace_lift_obstructions_20261020.md`.
