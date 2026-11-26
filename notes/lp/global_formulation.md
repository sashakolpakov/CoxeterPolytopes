# A global necessary-data formulation

25 November 2026. This replaces the discarded diagram-by-diagram search.
The established ceilings remain57/17/15. The construction below is a
research formulation, conditional on the stated geometric and literature
inputs, and has not received a separate author's source-bound review.
The separate computational sections identify the frozen inputs to which
each exact feasibility result applies. Newly added levels require their
own feasibility decisions.

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
constraints.

The second level applies the following additional families to **every**
fixed component counter and every growing-family count. Write u=m*t,
where m is a nonnegative integer, let r be the component's rank or the
tail's minimum rank, and let q=1 at finite vertices and q=1-t at ideal
vertices. For d=2,3,4 put

    f_d(u,t)=product_(j=0)^(d-1)(u-j*t).

Then f_d>=0. If m<d this product is zero; if m>=d the total rank gives
q>=d*r*t. Therefore also

    (q-d*r*t)*f_d(u,t)>=0.

These conditions are valid throughout the dimension tail. They impose
no fixed upper bound on component multiplicity. In a fixed dimension,
when d*r exceeds the available rank, the two inequalities force the
corresponding falling factorial to vanish.

For each k=0,1,2,3 the adjacent-integer gap

    g_k=(u-k*t)*(u-(k+1)*t)

is nonnegative: (m-k)(m-k-1)>=0 for every integer m. The second level
adds these scalar inequalities and the matrix localizers

    [[L(g_k),L(u*g_k)],[L(u*g_k),L(u^2*g_k)]] >= 0.

Their validity follows by integrating the pointwise PSD matrix
g_k*(1,u)^T*(1,u). No independence between counters is assumed. Existing
fixed-counter moment blocks already contain their fourth moments;
only the seven growing-family count fourth moments are new coordinates
in the dimension-57 export.

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

### Second-level shared Gram constraints

The second level keeps every first-level constraint and strengthens
the shared Gram schema uniformly for both face sizes. Write C_3(A) for
the third compound matrix: its rows and columns are indexed by triples,
and its (I,J) entry is det A[I,J]. The pointwise rank bound for M gives

    C_3(M)=0.

The first level already imposed its diagonal entries. The second level
also imposes the mixed minors; these give six additional equations on
quadrilaterals. Although the full rank condition follows from principal
minors for an actual PSD matrix, the additional identities need not
follow at a truncated moment level. Each identity is integrated before
taking moments; no rank bound on the averaged matrix is asserted.

There is also a stronger signature localizer. Write H=U^T*J*U, with
U having three rows and J=diag(1,1,-1). Cauchy--Binet gives

    C_3(H)=-v*v^T,   v_I=det U[:,I].

Consequently -C_3(H) is PSD. Its diagonal entries are the old principal
three-minor tests; keeping the full matrix couples those tests. This is
valid for singular configurations as well.

Finally the pointwise positive semidefiniteness of M and 0<=c_i<=1 give

    M tensor M >= 0,   c_i*M >= 0,   (1-c_i)*M >= 0.

The tensor matrix has entries M_ik*M_jl, indexed by pairs (i,j),(k,l).
It equals (W tensor W)^T*(W tensor W) when M=W^T*W. Integrating these
polynomial matrices produces necessary PSD moment blocks. In particular
they involve products of entries of the same actual face, rather than
products of independently averaged Gram entries. All entries use
existing Gram moments of degree at most four. Thus these geometric
families add no moment coordinates in the current export.

All these deductions hold in the three-dimensional residual Lorentz
space of a two-face, independently of the ambient dimension. The
integer-population families above are also uniform in dimension. Their
current anchored export tests dimension 57, with the established n>=58
exclusion covering the complementary tail.

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

## Supporting forests and attachment responses: levels 3--5

These levels retain every constraint of level 2. Level 3 adds a common
supporting-forest response model and deletion-count marginals. Level 4
adds joint closed-face conditioning and a shared edge statistic. Level 5
adds sparse polynomial-ideal closure. They remain necessary-data moment
relaxations; their feasibility is not a geometric realization claim.

### A component representation of the same white block

Use the original unit facet normals before projection. For an s-sided
two-face, write its star Gram matrix as

    G=[[B,C],[C^T,D]],  B=direct_sum B_alpha,
    T=C^T*B^-1*C=sum_alpha Q_alpha,  S=D-T.

Here B is the positive supporting Gram block of rank n-2. Adding one
sidewall to B gives the positive diagram of an edge. Its connected
components are spherical Coxeter diagrams, which are trees of maximum
degree three. Consequently each sidewall meets at most three white
nodes, and at most one in any one connected white component: two
connections to the same tree would make a cycle. At most 3s white
components have any attachment to these sidewalls. All remaining white
components are retained as spectators, with their rank and count.

For each white component, rescale its simple roots to the rational
Cartan normalization Bbar=L*B_alpha*L, with L diagonal and L_pp^2
equal to 2, except at short B/F roots, where it equals 1. If side i
attaches at node p_i, put b_i=-C_(p_i,i)*L_(p_i,p_i); otherwise b_i=0.
Then, simultaneously for all side pairs,

    (Q_alpha)_ij=b_i*b_j*(Bbar^-1)_(p_i,p_j).

These are the same marked nodes and amplitudes in every entry of the
response matrix. Every Q_alpha is PSD and entrywise nonnegative. Each
entry is at most one: its diagonal is bounded by the corresponding
diagonal of T, and positivity of the edge block gives T_ii<1.

All connected spherical families occur in the component domain:
A_r, B_r, D_r, E6/E7/E8, F4, H3/H4 and I2(m). There is no rank cutoff.
The final family marginals use the canonical convention I2(m), m>=5;
labels three and four are represented by A2 and B2. Ranks and marked positions are continuous outer
variables in this finite level, with the stated family bounds. No
exhaustiveness claim is made for a finite sample of ranks or labels.

For the growing families the Cartan inverse entries are exact:

    A_r: K_ij=min(i,j)*(r+1-max(i,j))/(r+1),  1<=i,j<=r;
    B_r: K_ij=min(i,j), with the short terminal root indexed by r;
    D_r: K_ij=min(i,j) on the chain 1,...,r-2,
         K_(i,leaf)=i/2,
         K_(leaf,leaf)=r/4,
         K_(different leaves)=(r-2)/4.

For A, the second difference 2K_ij-K_(i-1,j)-K_(i+1,j) is delta_ij,
with zero boundary values at 0 and r+1. For B, the same interior
identity and K_rj-K_(r-1,j)=delta_rj verify the inverse. For D, the
interior second differences, the fork row and the two leaf rows verify
the displayed matrix. This proves the formulas at every rank. In the
program use R=r/n, x_i=p_i/n, t=1/n and explicit ordered minima;
multiplication by the positive denominators gives polynomial identities.

The E/F inverses are exact rational matrices. For H3/H4, invert the
fixed Cartan matrix over Q(phi), phi^2=phi+1, then enclose each inverse
entry by rational bounds. The bracket

    1618033988749894/10^15 < phi < 1618033988749895/10^15

follows by evaluating x^2-x-1 at both endpoints, where it is strictly
increasing. Each entry a+b*phi consequently has checked rational lower
and upper bounds. Both bounds, together with PSD, are used in the model.
This is an enclosure of the fixed H matrices, not an exact rational
replacement of their entries.

For I2, with c=cos(pi/m), K=[[1,c],[c,1]]/[2(1-c^2)]. Its polynomial
response equation is used with 1/2<=c<=1 and c^2>=5/8 in that family.
The latter is a necessary envelope for m>=5. The endpoint is an outer
closure point. For each nonzero attachment, its original squared
cosine u satisfies 1/4<=u<=1 and

    (u-1/4)*(u-1/2)>=0,
    (u-1/2)*(u-5/8)>=0.

The second gap uses cos^2(pi/5)>5/8. Inactive attachments use the
homogenized versions with their Boolean activity flag. These are
necessary gaps, not an exact discrete-label encoding.

Each face retains T, positive projection lengths a_i=sqrt(1-T_ii),
and the original sidewall Gram entries D_ij=-d_ij. Thus

    T_ii+a_i^2=1,
    d_ij+T_ij=c_corner*a_i*a_j       (adjacent sides),
    d_ij+T_ij=opposite*a_i*a_j       (opposite sides).

The a_i have positive inverse witnesses. Adjacent d_ij is zero or in
[1/2,1]; nonadjacent original entries have no artificial upper bound.
The existing projected Gram matrix uses these same a_i and T entries.

### Component and vertex occupation marginals

For each triangle/quadrilateral, let k be its white component count,
kap=k/n, d_i the number of attachments of side i, and w_i the number
of white components meeting both sides at corner i. The side degrees
belong to {0,1,2,3}. If e_i records an original edge between those two
sidewall nodes, then

    e_i+w_i<=1+eta_i,
    (1-eta_i)*(e_i+w_i-1+right_i)=0.

At a finite corner the enlarged diagram is a forest. At an ideal
corner the selected n-node block has one affine component and otherwise
spherical components: it is obtained by taking an edge in one simplex
factor of the cusp link and vertices in the others. It has at most one
cycle, giving the displayed bound. The full cusp diagram can still have
arbitrarily many affine factors; no corank-one assertion is made for it.

Adding the two corner nodes to the white forest gives, at a finite
corner, its full normalized component count

    kap_v=kap+t*(2-d_i-d_(i+1)-e_i).

The shared edge along side i has normalized component count

    kap_E=kap+t*(1-d_i).

A component occupation measure records each connected white component
once per actual face. Its count, rank, degree and response pushforwards
are matched to the common face measure. For every retained common-state
test polynomial p,

    integral t*p dxi             = integral p*kap dmu,
    integral R*p dxi             = integral p*(1-2t) dmu,
    integral attach_i*p dxi      = integral p*d_i dmu,
    integral (Q_alpha)_ij*p dxi  = integral p*T_ij dmu,
    integral attach_i*attach_j*p dxi = integral p*w_corner dmu.

Spectators have every attachment flag and response entry zero, but
their positive rank and component count remain in the first two rows.
There is also a two-port residual-corner measure for larger faces. It
retains the same necessary corner/supporting data, without claiming to
encode an entire face with five or more sides.

The new observables are connected to the original whole-vertex
populations by exact deletion identities. For a spherical forest on n
nodes with k components, summing the component count after deleting
two nodes gives

    Psi=choose(n,2)*(k-2)+(2n-3)*(n-k).

Indeed the count for deleted nodes u,v is k-2+deg(u)+deg(v)-[uv is
an edge]. Sum this identity over all pairs, using sum deg=2(n-k).
The sum after deleting one node is (n-2)k+n. Thus, with kap_v=k/n,
the finite normalized columns used by the program are

    white: t^2*(1-2t)*((1-3t)*kap_v/2+t),
    finite component count: t^2*(1-t)*kap_v/2,
    edge component count: t^3*((1-2t)*kap_v+t).

For a cusp, let its affine factors have orders m_i=r_i+1 and unweighted
edge counts E_i, and put V=product m_i. A two-face through it selects
two deleted nodes in one factor and one in every other factor. Counting
the surviving edges gives

    Psi=(V/2)*((n-1)*(n-2)-(n-3)*sum_i E_i*(m_i-2)/m_i).

A link vertex deletes one node in each factor. Its mean spherical
component count is the sum of 1 for a cycle factor and 2r_i/(r_i+1)
for a tree factor. This includes the two-node affine A1 tree. These
expressions use exactly the existing ideal component counts, ranks and
reciprocal tail moments. They do not replace the product cusp by an
independent upper estimate for each statistic.

The final profile formulation also matches the count and rank sum of
**each** connected white family. This is needed because a total count
alone could match different supporting families to the vertex data.
Let A_(S,T,j) and B_(S,T,j) sum the T-family statistic after deleting
one or two nodes of a connected component S; j=0 counts components,
and j=1 sums their ranks. A finite component of rank r contributes

    choose(n-r,2)*[family(S)=T]*r^j
       +(n-r)*A_(S,T,j)+B_(S,T,j)

to the whole-vertex two-deletion sum. The three terms correspond to
zero, one and two deleted nodes in this component. They count the
resulting white components, so deleting nodes in different components
introduces no spurious factor two.

For a product cusp the corresponding sum divided by V is

    sum_S ((n-1-r_S)*A_(S,T,j)/2+B_(S,T,j))/(r_S+1).

The normalization is the same as above: divide the finite sums by
n^5, and the cusp sums divided by V by n^3. Their global totals match
the same component profile measures: t*type_T for j=0 and
R*type_T for j=1, with s copies per s-sided face. Thus the family
and rank data used in inverse-Gram responses are tied to the existing
whole-vertex spherical and affine populations.

For a path A_r the one-deletion count/rank sums are
(2(r-1),r(r-1)); the two-deletion sums are
(3(r-1)(r-2)/2,r(r-1)(r-2)/2), all of family A. For B_r, a B
component survives precisely when the distinguished label-four edge
survives. Its count sums are r-2 and choose(r-2,2); its rank sums
are r(r-1)/2-1 and (r-2)(r-3)(r+2)/6. For D_r, retain the fork
and its three neighbors. The D count sums are r-4 and choose(r-4,2),
and its rank sums are r(r-1)/2-6 and
(r-4)(r-5)(r+6)/6. The other pieces are A components; subtract
these contributions from the total forest count/rank identities.
These arguments prove all-rank formulas, including the zero boundary
cases, without interpolation.

The affine A family is a cycle, so every proper deletion piece is A.
Affine C is a path with label-four edges at both ends: retain either
end edge to count B pieces. Affine B has one such edge and one fork;
retain these respectively to count B and D pieces. Affine D has two
forks; after a proper deletion, surviving forks belong to distinct
spherical components. A surviving fork contributes a D piece exactly
when its four-node core survives. Counts are binomial counts of
deleted nodes outside these cores. Rank sums additionally sum the
number of choices avoiding the path from the retained edge/core to
each surviving node. These are arithmetic/binomial sums along the
classical chains. Their resulting polynomials are implemented in
`supporting_counts.py`. The affine tails begin at rank 25, so their
two fork cores/end regions are distinct. All smaller and exceptional
connected inputs instead use their complete exact labelled deletion
histograms, including I2/G2/H2 identifications.

Finite tail expressions are polynomials of degree at most three in
component rank and use the existing U,R1,R2,R3 fields. For affine
tails, divide the polynomial numerator by r+1; its quotient uses the
rank moments and its remainder uses the existing reciprocal moment.
No extra rank cutoff or independently adjustable family budget is
introduced by these twenty global family-statistic equations.

At levels 4--5 these marginals also hold separately at ideal corners.
For each family and j=0,1, the ideal-vertex deletion sum alone equals
the component-profile moment of

    (sum_(corners i) eta_i)*type_T*(t if j=0 else R).

Subtracting this equation from the combined marginal gives the finite
corner equation. The ideal indicators belong to the same common face
state as the attachments; the joint moment matrix includes the family
indicators together with the face coordinates. Thus a family distribution at ideal corners
cannot be funded by unrelated finite vertices. Its source is deletion
of complete classified affine cusp factors. For example, their deletion
histograms contain no H3/H4 component, so those components cannot have
positive ideal-corner mass. This is a uniform classified-cusp marginal,
not a selected prohibition of individual vertex forests. Full
conditional cusp completions require further moments beyond these
count and rank marginals.

### Closed-face conditioning and the edge marginal

Level 3 retains p=1. Level 4 also retains every linear common-state
coordinate and every full corner-kind pattern. A corner has three
categories eta_i, right_i and 1-eta_i-right_i. Products of one category
per corner give 27 triangle patterns and 81 quadrilateral patterns.
All component pushforwards above are matched with these same products.
Their nonnegative masses and pattern-weighted component, total-white
and spatial Gram PSD matrices are retained. The same response matrix
and attachment marks therefore participate around the whole cycle.
Finite matching of real moments still does not assert equality of the
entire conditional distributions.

Level 4 also retains the strict incidence of opposite sides. Their
projected Gram entry is strictly less than -1. Indeed, equality would
make the sum of their two unit normals a nonzero null vector v.
An interior point pairs negatively with v, fixing its future direction.
Every other side normal pairs nonpositively with v, because all
off-diagonal Gram entries are nonpositive. Thus its ideal point is in
the polygon's closure and on both supposed opposite sides. They would
be incident at an ideal vertex. This contradicts their opposite
positions in a genuine quadrilateral. Encode the strict gap with
(opposite-1)*z=1, z>=0; no uniform positive epsilon is assumed.
The gap also has the PSD localizer gap*(1,z)^T*(1,z). Together with
the inverse identity this prevents a zero expected gap from being
hidden by unconstrained inverse moments.

An additional edge measure carries kap_E. Its count matches f1, and
twice its first component moment equals the whole-vertex sum of edge
endpoint component moments. An edge belongs to n-1 two-faces, so

    t*sum_(two-face sides) kap_E = (1-t)*sum_edges kap_E.

The residual-corner measure selects its outgoing edge once per corner.
This is a common edge statistic linking faces and both endpoints. It
does not yet match every statistic of the full edge diagram, and no
complete three-face supporting-diagram measure is claimed.

### Sparse pointwise equation closure

Level 5 distinguishes pointwise domain equations q=0 from aggregate
incidence equations between measures. Only pointwise equations can be
multiplied by an arbitrary state polynomial. It adds L(q*x)=0 for
every coordinate x whose product is already in the retained registry.

For reduced q of degree at most four, at most eight variables and at
most 24 terms, it also adds L(q^2)=0 and a PSD moment clique containing
1, those variables and the monomials of q. An invertible constant
basis change replaces the leading monomial by q itself. Its zero
diagonal exposes the corresponding kernel explicitly. All entries use
the same measure's monomial registry, so cliques agree on overlaps.
The kernel equations L(q*b)=0 for the clique basis polynomials b are
also written explicitly. They already follow from the PSD zero row,
and give the LP outer approximation the same known linear identities.
The resulting moments have degree at most eight. The finite sparse
budgets and omitted multipliers are recorded for every measure; this
is not a claim to have closed the entire polynomial ideal.

Every addition in these three levels is valid in all dimensions of
the declared domain. Their current exports are anchored at dimension
57, with the previous n>=58 exclusion covering the complementary tail.

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

At levels 3--5 additionally record the actual white components and
their attachment nodes in the component measures, using their common
face state at every occurrence. Record the actual original-to-projected
scales and sidewall entries. The inverse identities, positive edge
blocks and forest/link counts above establish the new pointwise rows
and additive marginals. Use the actual corner-kind patterns and edge
supporting counts for level 4. Pointwise identities remain zero after
multiplication and squaring, proving every level-5 closure row and
moment clique. Integer ranks, discrete labels and exact H matrices
satisfy the declared outer envelopes. Thus no actual polytope is
excluded by replacing these exact data by the stated relaxations.

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

A bounded objective supplies another solver formulation after exact
equality elimination. Choose a free coordinate x_g with an existing pure
nonnegative scalar row. Write h=1/(1+x_g), y_i=h*x_i for i!=g, and
y_g=1-h. The equivalent search for a positive-scale point is

    maximize h,
    F_bg + sum_(i!=g) y_i F_bi + h(F_b0-F_bg) >=0 for every b,
    h>=0.

The old nonnegative x_g row becomes a positive multiple of 1-h>=0,
so 0<=h<=1. Each displayed matrix is exactly h*F_b(x) when h>0,
with x_i=y_i/h and x_g=(1-h)/h. Every original feasible point maps to
a feasible point with h>0, and every such positive-scale point maps back.
This changes neither the global constraint families nor their feature
coefficients. The implemented gauge is the nonnegative shared-edge
measure mass in the prepared coordinates. Its finite positive changes
of units are included in the recorded substitution.

An h=0 point is recession data and cannot certify original feasibility.
A numerical zero optimum cannot certify original infeasibility either;
it still needs a checked bound or another exact infeasibility argument.
Positive-scale numerical answers are rounded in the projective
coordinates, divided by their exact positive h, and tested in every
original equality and cone. The mapping identities are also checked
coefficient by coefficient. This bounded optimization is an alternative
to the zero-objective and normalized Farkas searches, not an extra
geometric assumption.

For numerical face discovery one can instead solve the auxiliary problem

    maximize a,   a<=1,   F_b(x)-a I>=0 for every b.

It has a known strictly feasible starting point: take x=0 and a=-rho,
where rho>=1 and every constant diagonal plus rho exceeds its row's
sum of absolute off-diagonal entries by at least one. Each shifted
matrix is then strictly diagonally dominant with positive diagonal.
The upper-bound cone is strictly positive as well.
The implementation verifies this rational inequality for every row.

The auxiliary dual has Z_b>=0, lambda>=0 and

    sum_b <F_bi,Z_b>=0,
    sum_b trace(Z_b)+lambda=1,

and gives the upper bound a<=sum_b <F_b0,Z_b>+lambda. A positive primal
margin gives original feasibility, subject to checking the actual original
matrices. Exact cone weights with zero constant pairing and lambda=0
expose a face: every positive summand in sum_b <F_b(x),Z_b>=0 must vanish.
For a positive rational rank-one summand v v^T this forces F_b(x)v=0.
These forced linear equations can be eliminated before another solve.
Exact negative constant pairing gives the usual normalized Farkas
certificate. Numerical zero margins or approximate kernel directions
alone establish none of these conclusions. This formulation differs from
minimizing a nonnegative common violation, whose zero objective can have
only the trivial cone multipliers and therefore need not expose a face.

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

The current dimension-57 level-5 SDP is exactly feasible. The
[rational witness](../../results/global_shared_cusp_scalar_circuits01/original_feasible_point.exact.json)
passes all 9,388 original equalities, 3,643 PSD blocks and 134,833 scalar
inequalities. A separate [integer Bareiss verification](../../results/global_shared_cusp_independent/result.json)
also checks all 193,440 adjoint equations in a certificate ruling out
the normalized Farkas dual. The original input is unchanged. This
computational conclusion concerns the stated outer relaxation and
does not lower the accepted dimension bound.

## Complete rooted edge hosts and shared edge populations: levels 6--7

These levels modify the compatibility retained by the global schemas.
They preserve the frozen level-5 constraints and all free feature
potentials. No new list of chosen forests or face diagrams is used.
The geometric premises are Vinberg's face theorem and cusp theorem,
Theorems 3.1--3.2 of [Hyperbolic reflection groups (1985)](https://www.mathnet.ru/eng/rm2140),
and the complete spherical/affine classifications already in use.
The deductions and their normalization are given here; separate author's
review is still pending.

### The same root controls deletion and projection

For a two-face F let B be its supporting spherical forest, and let i
be a sidewall. The edge E_i has spherical diagram B+i, of rank n-1.
Let H_i be its connected component containing i. Deleting i from H_i
gives exactly the white components attached to i. Every other white
component is an unchanged spectator. This is an exact graph operation,
simultaneously at all sides of this SAME face.

Let a_i be the existing projected length. A block inverse of the unit
Gram matrix G(H_i) gives

    a_i^2 = 1/(G(H_i)^(-1))_(i,i).

Indeed, the inverse entry is the reciprocal of the Schur complement of
the white principal block; the complement is 1-T_ii=a_i^2. Spectator
blocks are orthogonal and do not change it. Introduce v_i=a_i^2, but
continue to use the original a_i in every face Gram equation.

The complete rooted classification is finite with symbolic ranks:

* For A_r with root p, deletion gives A_(p-1) and A_(r-p), omitting
  zero ranks; v_i=(r+1)/(2p(r+1-p)).
* For B_r with the short node numbered r, a long root p deletes to
  A_(p-1)+B_(r-p). B_1 means A_1, and B_0 is empty. At the short
  root the deletion is A_(r-1). The inverse equations are 2p*v_i=1
  at a long root and r*v_i=1 at the short root.
* In D_r number the chain 1,...,r-2 and then the two leaves. A chain
  root p<=r-4 leaves A_(p-1)+D_(r-p). At p=r-3 the right component
  is A_3; at the fork p=r-2 it is 2A_1. A leaf deletion leaves A_(r-1).
  The inverse equations are 2p*v_i=1 on the chain and r*v_i=2 at a leaf.
* E6/E7/E8, F4 and H3/H4 use every root of their complete connected
  diagrams. Nodes are grouped only when both their deletion statistics
  and inverse-diagonal data coincide. Rational inverse entries are exact;
  H entries use the previously justified golden-field brackets.
* I2(m), m>=5, deletes to A1, and v_i=sin^2(pi/m). The retained
  envelope 0<=v_i<=3/8 covers every such label. A2 and B2 belong
  to their canonical classical families.

The classical inverse equations follow from the inverse matrices proved
above. The deletion descriptions follow by cutting a path or fork at
the root, with the low-rank boundaries explicitly assigned as stated.
With R=r/n and x=p/n, for example the A equation becomes
2*v_i*x*(R+t-x)=t*(R+t). Boolean root-case selectors and nonempty-arm
selectors give polynomial descriptions covering all ranks. Their
continuous rank/position relaxation is intentional.

Write W_(T,0)=t*number_of_white_T_components and
W_(T,1)=t*sum_of_their_ranks, in the common face state. All sides use
these same variables. Let delta_(i,T,j) be the rooted deletion statistic,
and let epsilon_(i,T) indicate the host family. Then

    W_(T,j)-delta_(i,T,j)>=0,
    E_(i,T,j)=W_(T,j)-delta_(i,T,j)
                 +epsilon_(i,T)*(t if j=0 else R_i).

The first quantity counts or sums the ranks of spectators of family T.
The second is the statistic of the WHOLE edge diagram B+i. It retains
the distinction between touched white components and spectators.
The white component occupation measure is matched to W, and its
attach_i-weighted family statistics are matched to delta_i, with the
declared finite common-state tests. Thus the root's classification,
projection length and actual white attachment families are coupled.

Two useful consequences use only old face coordinates:

    1-a_i^2-degree_i/4>=0,
    (degree_i-1)(degree_i-2)(degree_i-3)(1-a_i^2)=0.

For the inequality, each nonzero original Coxeter attachment has squared
cosine at least 1/4. The inverse of a positive unit-diagonal block has
diagonal at least one, by its scalar Schur complement. Distinct
attachments meet distinct white components, so their diagonal costs
add and T_ii>=degree_i/4. For the equality, degree_i is one of 0,1,2,3;
when it is zero there is no response and a_i=1. A 2-by-2 PSD localizer
of the inequality is retained as well. These are pointwise conditions,
not products of unrelated averaged constraints.

### One spherical population for every edge

The existing shared-edge measure now contains a whole spherical forest
population of rank n-1, with all fixed connected inputs and the same
symbolic A/B/D tails. Its counters have the integer factorial, rank and
PSD conditions used for vertex populations. Let C_(T,j)(E) denote its
normalized family count or rank sum, both normalized by t as for W.
Both endpoints and every incident face must produce this SAME measure.

For a connected spherical S of rank r write
A_(S,T,j)=sum_(p in S) statistic_T,j(S-p), with unnormalized statistics.
Its contribution to the endpoint sum at a finite vertex is

    (n-r)*[family(S)=T]*r^j + A_(S,T,j).

At a cusp with affine factors S of orders r+1 and V=product(r+1),
the contribution is V*A_(S,T,j)/(r+1). An edge at that cusp selects
one deleted node in each factor, so this is an exact product-link sum.
After normalizing the statistics by t and using the existing vertex
measure scales, the finite column is t^5 times its unnormalized sum;
the ideal column is t^3*sum A/(r+1). Hence, for every family and j=0,1,

    2*integral C_(T,j) dmu_E
       = integral endpoint_F_(T,j) dnu_F
         + integral endpoint_I_(T,j) dnu_I.

Each edge lies in exactly n-1 two-faces: deleting one of the n-1
supporting nodes gives those faces, by Vinberg's positive-face theorem.
Consequently another exact marginal is

    sum_small_faces integral t*sum_sides E_(i,T,j) dmu_F
      + integral t*E_(outgoing,T,j) dmu_R
        = integral (1-t)*C_(T,j) dmu_E.

The residual corner measure selects one outgoing side per corner; this
visits every side of every larger face once. It does not presume a
three-dimensional face to be Coxeter.

Deleting a node of the whole edge gives the same white supporting
forest as on the corresponding face flag. In a connected edge component
S of rank r its deletion sum is

    (n-1-r)*[family(S)=T]*r^j + A_(S,T,j).

Multiplying by t^2 gives the common-edge column matched to t times
the sum of the incident face white statistics. This provides a second
route through the same edge population. The all-rank classical A
polynomials and the fixed exceptional deletion tables are shared with
the existing two-deletion implementation, rather than fitted anew.

### Stronger pointwise equation closure

Level 7 raises the variance-clique identity-degree budget from four to
six and its variable budget from eight to twelve, retaining the term
budget. In particular this admits mixed A/B inverse-response identities
that previously received only sparse first-moment equations. Eligible
identities also give L(q^2)=0 and PSD clique kernel equations.
It also imposes L(q*m)=0 for every retained multiplier monomial m of
degree two for which all product moments are present. All products use
one measure's existing registry; aggregate equations relating different
occupation measures are never multiplied by local state variables.
These operations preserve every actual-polytope Dirac/counting lift.

The finite common-state marginals still do not assert a discrete global
realization. Additional edge state moments and genuinely joint three-face
states remain possible extensions. Feasibility of levels 6 or 7 must be
decided on their own frozen inputs; the level-5 feasible point proves
nothing about their feasibility.

## Dual certificate after exact equality reduction

This is an algebraic change of representation of the frozen research model.
It adds no geometric condition. All matrices below are restricted to the
previously proved face of the certificate cone; extending a point of the
corresponding compressed primal to the full original primal is a separate task.

Write the prepared source blocks as
`F(x)=F_0+sum_j F_j*x_j`. Pack symmetric matrices by their actual upper-triangle
entries. Let `S` multiply off-diagonal entries by `sqrt(2)`, so `W=S^T*S` has
diagonal weights one and off-diagonal weights two. In these coordinates the
source certificate equalities are `E_j q=0`, where `E_j=vec(F_j)^T W`, and
`1+E_0 q=0`.

The verified rational elimination gives a homogeneous map `q=Tz` with
`E_j T=0` for all ordinary coordinates. Its free coordinates occur as identity
rows of `T`. The positive common unit is `tau=62,500,000`. The balanced
certificate feasibility problem is

    tau*E_0*T*z = -1,       S*T*z in K.

Its canonical matrix is the vertical concatenation of `tau*E_0*T` and `-S*T`.
For each ordinary source coordinate `j`, the exact elimination witness selects
a scalar diagonal `p_j` such that `F_0[p_j]=0`, `F_j[p_j]=c_j != 0`, and
`F_i[p_j]=0` for `i != j`. The checked source scan verifies these statements.
Define

    G(u) = tau*F_0 + sum_j (F_j/c_j)*u_j.

Every `u` for which all blocks `G(u)` are PSD yields a normalized dual
infeasibility certificate for the balanced problem:

    y_equal = 1,             y_cone = S*vec(G(u)).

Indeed its right-side pairing is `-1`, and its adjoint is identically zero:

    tau*(E_0*T)^T - T^T*W*vec(G(u))
      = tau*T^T*E_0^T - T^T*(tau*E_0^T + sum_j E_j^T*u_j/c_j)
      = 0.

Thus all dual linear identities hold over any real field before choosing `u`.
Only the PSD conditions remain. The chosen diagonal of `G(u)` is exactly `u_j`,
which provides a numerical coordinate inverse. This parametrization uses 61,513
coordinates, compared with the balanced certificate problem's 300,139 free
variables. Preparation and a separate replay check every source coefficient,
cone partition, pivot and numerical export. The field-certificate checker also
recomputes all 300,139 adjoint identities directly from `T`, with exact trace
weights; it does not accept floating residuals as identities.

Kernel relations inferred from a numerical candidate restrict reconstruction
only. After exact affine projection over `Q` or `Q(sqrt(d))`, every scalar and
PSD block must be checked at the specified real embedding. For quadratic
fields this checker selects `sqrt(d)>0` and decides signs by rational norm
comparisons. An accepted candidate still requires a separate full replay.
Such a certificate proves the balanced certificate search infeasible. It does
not by itself prove full original-primal feasibility or geometric realizability.
