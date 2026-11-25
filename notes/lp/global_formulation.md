# A global necessary-data formulation

25 November 2026. This replaces the discarded diagram-by-diagram search.
The established ceilings remain57/17/15. The construction below is a
research formulation, conditional on the stated geometric and literature
inputs, and has not received a separate author's source-bound review.
No numerical result is asserted by this specification.

## Target and inputs

For a proposed N, the geometric hypothesis is that a finite-volume
hyperbolic Coxeter polytope exists in some integer dimension n>=N.
Write t=1/n. Odd and even dimensions have different finite allowances;
each parity is represented by one problem on0<t<=1/N_parity. Artificial
mixtures of dimensions and closure points at t=0 are permitted by the
outer relaxation. Every actual polytope still maps into it.

The exceptional indicator e=[n=58] satisfies e(58t-1)=0. The
dimension-gap rows admit every integer dimension in the stated parity.
In an even tail starting below 58, the additional gap uses
(1-e)(1-56t)(1-60t)>=0; the simpler (1-e)(1-60t)>=0 is used only
when the even endpoint is already at least 58.

The completed exclusion for n>=58 remains available. An optional anchor
restricts the new problem to N<=n<=57, or1/57<=t<=1/N. This excludes no
hypothetical polytope left by the established theorem. It is one uniform
band formulation, not a collection of selected diagram cases. At N=57
the odd band has t=1/57 and the complementary dimensions are already
excluded. The new constraints themselves are valid throughout the tail.

The inputs are:

* The complete spherical and affine Coxeter classifications and vertex
  figures, including products of affine simplex factors at nonsimple
  cusps. Their total affine rank is n-1; the whole cusp Gram matrix is
  not assumed to have corank one.
* The intrinsic26 corner functions and the incidence-excess theorem of
  the [LP manuscript](../../arxiv/linear_programming_dimension_bounds/main.tex).
  The table source is pinned in `data.py`. Its old displayed certificate
  coefficients are not inputs to the new search.
* Vinberg's positive-face and parabolic-vertex facts, and the common
  Lorentz geometry of the supporting walls of an actual face.
* The uniform Prokhorov and face/link rows recorded in
  [global_sources.md](global_sources.md). Their exact hypotheses are
  retained; a projected three-face is not assumed to be Coxeter.

All26 intrinsic functions remain the same functions on both sides of an
incidence marginal. The extra n=58 cusp correction is retained as a
separate feature. Two uniform indicators are also used: ideal corner,
and a finite corner whose marks belong to distinct spherical components.
The latter is precisely a right angle in the intrinsic two-face.
Their aggregate formulas are stated below. Coefficients of all these
features are free dual variables, found together with the multipliers
for the new constraints.

One further common observable is the square of the intrinsic planar
angle cosine. It couples the complete vertex Gram data directly to the
shared face Gram variables, and has its own free dual potential.

## Complete vertex populations with finitely many coordinates

For finite vertices, use fixed connected-type bins A_r for1<=r<=20,
B_r for2<=r<=20, D_r for4<=r<=21, D28, and all exceptional spherical
types. The three remaining families are A_r,B_r for r>=21 and D_r for
r>=22. Assign D28 to its fixed bin. Allowing it additionally in the
relaxed D tail only enlarges the domain. For ideal vertices, use all
connected affine types of rank below25, and the four classical affine
tails of ranks at least25. These are complete classification schemas,
not a list of chosen forests. There is no maximum component rank.

For a fixed type T let u_T=m_T/n. For each growing family retain

    U=m/n,  R1=sum r/n,  R2=sum r^2/n^2,  R3=sum r^3/n^3.

All are bounded and nonnegative. The same state contains every fixed
counter, growing-family moment and contextual guard. Whole-vertex rank
is exactly1 for a finite vertex and1-t for an ideal vertex. If a tail
starts at a, valid conditions include

    R1>=aU,  U>=tR1,
    a*t*R1<=R2<=R1,  a*t*R2<=R3<=R2,
    U*R2>=t*R1^2,
    [[R1,R2],[R2,R3]]>=0,
    [[U,tR1],[tR1,tR2]]>=0.

The matrix statements follow from Cauchy or from a component measure
on the relative ranks r/n. Fixed-bin upper bounds are u_T<=1/rank(T).
Integer populations give u_T(u_T-t)>=0, and the same factorial-moment
condition holds for every tail count U. These are uniform population
constraints. Higher integer factorial moments can strengthen later levels.

A guard for an integer count m>=k is represented by a Boolean g with

    g*(u-k*t)>=0,  (1-g)*((k-1)*t-u)>=0.

Apply this also to the stated sums of counters. Products and Boolean
unions give precisely the old contextual guards. No vertex budget is
computed from an independent maximum for each guard or each feature.

There is also a uniform rank-support inequality for each old contextual
feature. If its nonzero support needs at least R_j total rank, then

    (1-R_j*t)*Phi_j>=0.

When n<R_j the feature is zero; otherwise both factors are nonnegative.
The required ranks add the ranks of the marked components and their
disjoint witnesses, with the specified multiplicities. This gives
R_j=58 for j=16,17,18,24, so these coordinates vanish at dimension 57.
The code applies this construction to every contextual coordinate.

Set sigma_j=sum r_i^j/n^j, including fixed and growing components.
The rank identities and nonnegative elementary symmetric sums give
further polynomial inequalities. For example, with total relative rank
q, sigma2<=q^2 and q^3-3q*sigma2+2sigma3>=0.

## Uniform intrinsic sums

Let S_T be the sum of the manuscript's ordinary pair values in one
connected finite component. Its literal large-rank definitions give

    S_A(r)=5r+23 (r>=21),
    S_B(r)=5r+26 (r>=21),
    S_D(r)=5r+30 (r>=22).

Here is a direct all-rank count for these ranges. Distances1 through5
contribute5r-15 in a path, and5r-14 in D. At each distance6,7,8 the
eligible boundary positions number12,13,14 in A,B,D respectively.
At distance9 the total is2 in each family; the two short D leaves
each contribute one-half. The boundary regions do not overlap in the
stated ranges. All other entries vanish. This proves the formulas
without extrapolating from a finite table of ranks.

For an affine factor of order m=r+1 in the rank>=25 ranges,

    S_aff,A=5m, S_aff,B=5m+46,
    S_aff,C=5m+41, S_aff,D=5m+51.

The cycle has five paid pairs per node. The three tree families have
short-distance totals5m-14,5m-15,5m-13, and respectively15,14,16
eligible pairs at each distance6 through9. Their two bounded end
regions are disjoint for m>=26. These counts follow directly from the
literal feature-node definition, including the nodes on label-four
edges. Lower ranks and exceptional components use the fixed intrinsic
tables, not these formulas.

For an affine tail also retain Vminus=sum(1/(r+1))/n. It satisfies

    (a+1)*Vminus<=U, (1+t)*Vminus>=tU,
    [[R1+U,U],[U,Vminus]]>=0.

The last matrix is Cauchy for the sequences sqrt(r+1) and its reciprocal.
Thus both count and reciprocal contributions remain coupled to rank.
The ordinary affine sums per order are5,5+46/(r+1),5+41/(r+1),
5+51/(r+1), so their complete aggregate is polynomial in these fields.

The other old finite aggregate features are the manuscript's exact
counter polynomials. Same-component context features are m_T*g; cross
features are their orbit size times m_T*m_U*g. The repeated A2 feature
is2m_A2(m_A2-1)*g, not a square without the subtraction. All fixed
types used in these polynomials are present in the complete binning.

## Vertex columns and normalizations

Write E_v for the number of corners and beta_v for its allowance.
For a finite vertex use the column scale n^4; for an ideal vertex with
factor ranks r_i and V=product(r_i+1), use V*n^2. The scale multiplies
the entire occupation column, including its allowance and every
incidence statistic. This is a positive change of measure, not a change
of corner weights.

The finite normalized columns include

    E_v/n^4 = t^2(1-t)/2,
    Phi1/n^4 = t^3 * (sum S_T*m_T/n),
    beta_v/n^4 = t^2(1-6t+5*[n odd]*t^2)/8,
    right-corners/n^4 = t^2(1-sigma2)/2,
    edge-endpoints/n^4 = t^3.

The right-corner identity counts pairs in different spherical components.
Let Y be the sum of independent uniforms on {-r_i,-r_i+2,...,r_i} at an
ideal vertex, and b=E|Y|/n. Then

    E_v/(V*n^2)=t(1-t)/2,
    Phi2/(V*n^2)=t * (sum S_aff,T/(r_T+1))/n,
    beta_v/(V*n^2)=t*(3(1-t)-5b)/4,
    edge-endpoints/(V*n^2)=t^2.

All other old intrinsic coordinates vanish at an ideal corner. The
exceptional cusp-correction column is its corner-count column times a
fraction in[0,1], supported at n=58. The ideal-corner column equals the
whole ideal corner count; its right-corner column is zero.

The additional angle-square feature is described below. Its ideal column
also equals the whole ideal corner count, since an ideal planar angle
has cosine one.

Use3b^2<=sigma2+2t(1-t), from the exact variance of Y and Cauchy, and
b<=1/2. For the latter, convolution preserves symmetric unimodality;
a symmetric unimodal distribution on its lattice interval is a mixture
of uniforms on central intervals. Each has no greater absolute first
moment than the uniform on the full interval. Dividing by n gives the
stated bound. These are valid envelopes; the full product distribution
is not claimed to be characterized by these moments.

The same mixture argument improves the variance bound to

    4b^2<=sigma2+2t(1-t)+[n even]*t^2.

For a centered uniform on an even lattice, 4(E|Y|)^2<=3 Var(Y).
For a centered uniform on an odd lattice, equality holds after adding
one to the right side. A symmetric unimodal law is a mixture of these
central uniforms on its own lattice. Jensen's inequality for the
absolute first moments gives the displayed bound, since their means
are zero and variance averages under this mixture. The lattice is odd
exactly when n is even. This is a uniform condition on all cusp products.

The code also uses positive constant column scaling and the coordinate
N_parity*t to improve numerical conditioning. Every row uses the same
scaling. Its exact factors are recorded in model metadata.

## One shared Gram schema for all small faces

An actual two-face has a positive definite supporting block B of rank
n-2. Project all its side normals away from that same white block.
Normalize each projected side normal to square one. If H is their Gram
matrix, H has diagonal1, nonpositive off-diagonal entries and lies in
a Lorentz space of dimension three. No projected angle is assumed to
be a Coxeter angle.

For s=3 or4 set adjacent entries H_i,i+1=-c_i, with0<=c_i<=1.
The two opposite entries of a quadrilateral are additional nonnegative
variables with a minus sign. Choose a unit future point p in the
relative interior of the actual face, and put h_i=-<p,u_i>>0. Then

    M=H+h*h^T is PSD and has rank at most2.

Indeed, M is the Gram matrix of the projections of the side normals to
the Euclidean plane p-perp. This uses the same p and the same walls at
every corner. Strict positivity is encoded polynomially by h_i*z_i=1,
h_i,z_i>=0. For PSD M, rank at most2 is equivalently the vanishing of
all principal order-three minors. Also impose det H=0 for a
quadrilateral and the necessary nonpositive order-three principal
minors of H. Rank equations are polynomial equalities; the program
does not mistake them for convex PSD constraints.

An ideal corner has c_i=1. A finite corner has c_i<1; an auxiliary
nonnegative inverse encodes that strict gap. A finite right corner has
c_i=0. Hyperbolic area gives at most one right angle on a triangle and
at most three on a quadrilateral. These necessary conditions hold for
every actual face, including ideal vertices. The common interior point
prevents admitting an opposite pair of equal, oppositely oriented walls
as a genuine polygon merely because some principal Gram tests pass.

For a quadrilateral also impose -H_02>=1 and -H_13>=1. To prove
this, suppose an opposite pair has a positive definite two-wall Gram
block B. Project an interior point orthogonally away from that block.
The correction has nonnegative coefficients B^-1*b, where b_i>0.
All other wall entries are nonpositive, so the projected point still
strictly satisfies their inequalities. Its square is the original
negative square minus b^T*B^-1*b, hence it lies in hyperbolic space.
It belongs to both opposite sides inside the polygon, contradicting
their incidence. Thus their two-wall Gram cannot be positive definite.
The weak lower bound one is sufficient here; no upper bound is imposed.

The rank test adds information to separate-corner and signature tests.
For example, the formal matrix with diagonal 1, cyclic adjacent entries
-1/2 and both opposite entries -3/2 has eigenvalues
-3/2,1/2,5/2,5/2. Its adjacent pairs are positive definite and it has
one negative direction, but its determinant is -75/16. Its rank four
excludes it from the shared two-face Lorentz space. This is a matrix
separation witness, not a claimed realizable star.

## Binding intrinsic features to the shared geometry

For a positive spherical Coxeter tree, let marked nodes have graph
distance d. Its projected two-wall angle satisfies

    cos(alpha)>=1/(d+1).

To prove this, eliminate branches off the unique marked path. The
remaining positive tridiagonal matrix has diagonal at most1 and
negative edge magnitudes at least1/2. The endpoint inverse correlation
is the product of these edge magnitudes divided by the square root of
the two endpoint-deleted determinants. Increasing a diagonal increases
the determinant; decreasing a negative edge magnitude increases it,
since the inverse of a positive Stieltjes matrix is nonnegative. The
interpolation remains a positive Stieltjes matrix. Each deleted
determinant is consequently at most(d+1)/2^d, its value on the simple
A_d path, and the numerator is at least2^-d. The correlation, which
equals the projected angle cosine, is at least1/(d+1).

For disconnected marked components that cosine is zero. Conversely,
the inverse of an irreducible positive Stieltjes block is strictly
positive, so a finite right corner has disconnected marks. In
particular every connected finite corner has c_i>=1/n. Every positive
ordinary finite table entry has d<=16: ranks up to17 have that
diameter bound, while the literal higher-rank branch has support only
through distance9. Thus its supporting intrinsic atom has c_i>=1/17.
The contextual same-component atoms have the precise distances listed
in `data.py`, giving stronger bounds1/(d+1).

The26 atoms of the manuscript's corner nonnegativity cone have their
usual literal feature vectors. Every actual corner vector is a convex
combination of them. This containment proof uses the feature definitions
alone; it does not require every atom to be realizable at the dimension
being searched. Same-component corners have at most one contextual
indicator; cross-component contexts are disjoint except that23 implies22;
their ordinary value vanishes on the contextual support. These cases,
and the ideal segment, give the complete convex containment.

Disintegrate each corner's mass into those atom weights. The four
geometric groups are neutral, finite connected, ideal, and finite
cross-component. The atom submeasures keep their separate t,e moments,
while each group shares the face's actual c_i and ideal/right flags.
Every group/atom and group/face pushforward is matched through degree2.
The group's first angle moment also bounds the sum of its atom masses
times their respective inverse-Gram lower bounds. This preserves all26
feature channels without introducing a catalogue of marked stars.

## A common inverse-Gram angle observable

At a finite vertex with positive Gram G, put
c_ab=(G^-1)_ab/sqrt((G^-1)_aa*(G^-1)_bb). This is the cosine of the
planar angle after projecting away every other incident wall. It is
zero between different components. Let C_T=sum_(a<b)c_ab^2 in one
connected component. With R=r+1 and H_r^(2)=sum_(k=1)^r 1/k^2,

    C_A(r)=R^2*H_r^(2)-(3R^2-R-2)/2,
    C_B(r)=r(r-1)/4,
    C_D(r)=r(r-1)/4-1/2-2/r+4/r^2.

These formulas cover the whole families. In A, numbering by 1,...,r,
c_ij^2=i(R-j)/(j(R-i)) for i<j. Sum first over i to obtain

    sum_j (R/j-1)*(R*(H_(R-1)-H_(R-j))-j+1).

Use sum_(j=1)^(R-1) H_(R-j)/j=H_R^2-H_R^(2), together with
sum_(k=1)^(R-1)H_k=R*H_(R-1)-(R-1). Substitution gives C_A above.
The first identity follows by induction: both sides increase by
2H_(R-1)/R when R increases by one.
For B, the root basis e_i-e_(i+1), e_r has inverse Gram entries
min(i,j), so c_ij^2=i/j and the sum is r(r-1)/4.
For D, chain pairs have c_ij^2=i/j, each chain-leaf pair has i/r,
and the leaf pair has (r-2)^2/r^2. Summing gives C_D.

All fixed crystallographic components use exact rational inverse-Gram
calculations or these formulas. H3,H4 use the already proved path
lower bound and the pair-count upper bound. The symbolic I2(m), m>=5,
uses 5/8<=cos^2(pi/m)<=1, since cos^2(pi/5)=(3+sqrt(5))/8>5/8.
These are outer bounds for the complete noncrystallographic families,
rather than sampled labels.

For the growing A family, use rational a_-=1644933/10^6 and
a_+=1644935/10^6 with

    a_-<=sum_(k>=1)1/k^2<=a_+.

The exact finite sum through 1000 and the integral tail bounds
1/1001<=sum_(k>1000)1/k^2<=1/1000 prove these inequalities;
the checker verifies both rational comparisons. At every r>=1,

    (a_--3/2)R^2-R/2-1/r <= C_A(r)
        <= (a_+-3/2)R^2-R/2+1.

In D, its correction to r(r-1)/4 lies between -1/2-2/a and -1/2
when r>=a. The normalized whole-vertex sum is therefore enclosed by
polynomials in the existing U,R1,R2,t fields. For example,
sum m(r+1)^2/n^2=R2+2tR1+tU. These identities give a finite uniform
tail representation, without imposing a maximum component rank.

Let z=sum_T m_T*C_T/n^2. Keep its proved lower and upper polynomial
bounds and 0<=z<=(sigma2-t)/2. The finite normalized feature column
is t^2*z; the ideal column is its corner-count column. Match their
global sum to sum_i c_i^2 on every small face and to c^2 on residual
corners. The same observable is used on both sides, so its coefficient
is a new unrestricted dual weight. The matrix formulas and tail
identities above provide the actual-polytope lift for this feature.

## Retaining the uniform triangle theorem

The triangle proof in the retained
[finite-volume article](../../manuscripts/finite_volume57_consolidated.tex),
Appendix `app:routing`, selects two distinct pure ordinary-unit corners.
Its inputs include the minimal-indefinite-core classification and the
retained 1025 marked-core test. No new execution of that test is claimed.
The proof treats its four exceptional marked triples and the three
reduced small affine pairs explicitly.

Its extension to this dimension range uses the content of that proof.
A minimal indefinite core contains all three side roots and has order
at most ten. The selected paths and their persistent features survive
arbitrary additional supporting roots. The exceptional-case kernel and
inverse-Gram arguments likewise allow arbitrary spectator rank. The
small affine exclusions use only two positive deleted-factor blocks
and the actual cusp vector. None of these steps uses n=58.
The selected finite paths have length at most seven; all contextual
same-component features have longer marked paths, and cross-component
contexts cannot occur on those connected paths. Each selected corner
therefore has exactly atom E1 or E2 in the present coordinates. Hence

    sum_(three corners)(mass(atom1)+mass(atom2))>=2*mass(mu3),

also multiplied by t and e, is a uniform necessary row. The source
classification remains an explicit input, and this extension has not
received a separate author's review. Its multiplier is free alongside
the feature potentials.

Preserving all 26 functions does not automatically preserve every old
face theorem. Some quadrilateral transfers use equalities forced by the
dimension-58 vertex LP. Those implications are not imported at dimension
57 without their premises. This initial model retains the corner cone
and the uniform triangle cover, but does not yet encode all the old
quadrilateral transfers. Common supporting-diagram marginals are still
needed beyond the projected-angle matching already present.

## Global incidence and the positive excess

Let mu3,mu4 be the face occupation measures, muR the corners on faces
with at least five sides, mu5 their face-count measure, and muD their
sum(sides-5) measure. Let nuF,nuI be the scaled whole-vertex measures.
For every feature j and potential multiplier q in{1,t,[n=58]}, impose

    integral q*Phi_j dnuF + integral q*Phi_j dnuI
      = integral q*sum_corners phi_j dmu3
        + integral q*sum_corners phi_j dmu4
        + integral q*phi_j dmuR.

Here Phi_j denotes the complete column with its stated scaling, not
an independently chosen upper envelope for that coordinate. All old
features, the cusp correction and both extra indicators use these
same measures. Match also muR=5mu5+muD and

    f2=mu3+mu4+mu5,
    2f1=integral edge_endpoints dnuF+dnuI.

Prokhorov's P,Q,A,S ports are integrals over nuF,nuI themselves. With
sigma_j as above their finite P,A columns are

    P=t(1-t)(1-2t)/6,
    A=t(1-3sigma2+2sigma3)/6,

and Q=0. The ideal columns are

    P=(sigma2-t(1-t))/6,
    Q=((1-t)^2-sigma2)/8, A=0.

The distance-charge column S is(nW-J)/2 at a finite vertex and
(V/2)*sum_T H_T/(r_T+1) at an ideal vertex, where W is Prokhorov's
fixed distance charge, J its(d+1)-weighted version and
H_T=sum_deleted-node W(T minus node). This is computed from the same
complete vertex population. In the growing affine tails the deletion
average has the form

    H_T/(r+1)=(29r+b_T+c_T/(r+1))/3,
    (b_A,c_A)=(-176,0), (b_B,c_B)=(-349,1687),
    (b_C,c_C)=(-352,1696), (b_D,c_D)=(-346,1678).

For the tree families H_T=mW(T)-J(T). Their distance histograms through15
are those of a path, with one extra distance-two pair in B and two in D.
For an affine A cycle every deletion is A_r. This proves the displayed
tail formulas; all low-rank and exceptional deletion averages are fixed
input calculations on their connected graphs.

The exact edge/face/link marginal rows and Prokhorov capacities are in
`incidence.py`; none of P,Q,A,S is exported as an unbound free port.
Their coefficients are polynomials in n. For a row of degree d multiply
by t^d before taking moments. This is positive for every actual dimension
and gives polynomial coefficients on the entire parity tail. The large
four-face bin has no maximum facet count.

Finally impose the normalization

    2*mass(mu3)+mass(mu4)-mass(muD)
      -integral beta_F dnuF-integral beta_I dnuI=1.

Its unnormalized left side is D5(P)-sum_v beta_v, strictly positive by
the retained incidence-excess theorem.

## Actual-polytope lift

Given any hypothetical P in the specified tail, set t=1/n and the
parity indicator to its actual value. If the retained anchor is used,
P lies in the stated band. Divide all actual incidence counts by
L=D5(P)-sum_v beta_v>0. At each vertex record its whole connected
population, all its contextual guards and moments, and, at a cusp,
its actual b and exceptional-corner fraction. Use mass n^4/L for a
finite vertex and V*n^2/L for an ideal one, followed by the recorded
positive constant numerical column rescaling.

At every small face record its projected common Gram matrix and an
actual interior point. At every corner split the literal feature
vector into its nonnegative atom coordinates, retaining its geometric
group. Use actual residual corners and the actual counts of all the
remaining face flags. The rank, inverse-Gram, moment and incidence
arguments above show that these data satisfy every pointwise domain
condition and every shared marginal. The normalization is the definition
of L. Thus P supplies nonnegative occupation measures satisfying the
entire necessary-data problem. Truncating to the prescribed finite
moments preserves all equations and PSD/localizing conditions used by
the compiler. This proves the required necessary direction for every
polytope, without a selected list of face or vertex types.

## Finite LP/SDP and its dual

For each occupation measure keep the explicitly selected moments and
moment/localizing matrices. Polynomial identities become linear moment
equations. Polynomial nonnegativity and matrix positivity give scalar
and PSD rows. All variables are nonnegative state coordinates, so the
compiler also uses valid monomial positivity and bound-product rows.
Only variables with proved upper bounds are boxed; opposite Gram
entries, interior distances and inverse strictness witnesses remain
unbounded. Sparse moments need not have representing measures. That
enlarges the feasible set and is compatible with the lift above.

For each non-Boolean nonnegative coordinate x also impose the support
localizer [[L(x),L(x^2)],[L(x^2),L(x^3)]] PSD. If x<=u, impose
[[L(u-x),L((u-x)x)],[L((u-x)x),L((u-x)x^2)]] PSD as well. These follow
by integrating x vv^T and (u-x)vv^T for v=(1,x). Unbounded positive
geometry coordinates therefore have more than nonnegative individual
moments; no artificial upper bound is introduced.

The resulting finite model is

    F_b(x)=F_b0+sum_i x_i F_bi >=0,   E(x)=E0+sum_i x_i Ei=0.

Its explicit certificate dual seeks Z_b PSD and free y with

    sum_b <F_bi,Z_b>+y^T Ei=0 for every i,
    sum_b <F_b0,Z_b>+y^T E0=-1.

A primal point would then imply-1=sum_b <F_b(x),Z_b>>=0. The free
multipliers of the feature marginals are the unknown corner weights;
the other multipliers supply the controlled global compensation missing
from independent vertex budgets. Positive row scalings are recorded,
so the displayed weight coefficient is its multiplier divided by that
row's recorded scale. Off-diagonal trace pairings have factor two.
Solving this dual seeks feasibility. Infeasibility of the dual search
does not establish geometric nonexistence. Weak SDP infeasibility can
also prevent obtaining a separator of this normalized form.

The exploratory LP uses a declared outer approximation of PSD blocks
by e_i and e_i+-e_j test directions. The complete SDP remains separately
exported. A numerical result is only a candidate. Exact verification
checks every original equality, the complete adjoint and cone membership.
A new accepted global bound additionally requires the stated source
review and root acceptance.

Rational PSD separation can strengthen that LP without changing the
global schemas. A violated matrix supplies a numerical eigenvector,
which is rounded to a rational v; the exact row v^T M v>=0 is necessary.
An LP certificate with nonnegative multipliers lambda lifts by
Z_b=sum lambda*v*v^T in each original block. The checker verifies the
complete original adjoint and PSD matrices. Failure to find a negative
eigenvalue at a floating tolerance is not an exact feasibility proof.

The initial hierarchy does not enforce every integer population relation,
the full ideal-link product law, or all overlaps between different face
stars along edges. Those are possible uniform refinements of this fixed
global model. They are omitted constraints, not omitted classes of
polytopes; the actual-polytope lift still applies.
