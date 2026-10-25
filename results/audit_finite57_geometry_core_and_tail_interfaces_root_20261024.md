# Finite57: fresh review of the common geometry and tail interfaces

24 October 2026. This is a manual, source-bound audit by the coordinating
author. It covers the common actual-face geometry, the counting comparison,
the elementary high-core reduction, and the geometric reduction used by the
unbounded tail. It does **not** certify all distributed local receipts,
every dimension-specific marked table, or the complete scientific codebase.
Those remaining obligations are listed below. No numerical bound changes.

The reviewed source is `manuscripts/finite_volume57.tex`, SHA256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
I reconstructed the arguments from the current TeX before comparing the
older review verdicts. I read Sections 1--7, the core/extension and routing
appendix, the common bounded-tail reduction, and the complete unbounded-tail
argument. This is not a claim that every numerical table or certificate in
the other appendices was independently regenerated. The original source,
PDF and historical reviews are unchanged.

## 1. Actual faces, rather than formal principal diagrams

Let the outward roots have square two. A positive-definite principal Gram
matrix B with nonpositive off-diagonal entries has nonnegative inverse.
For a strict interior timelike point p, put b_i=-(p,e_i)>0 and
q=p+sum_i(B^{-1}b)_i e_i. Then

    (q,e_i)=0,   q^2=p^2-b^T B^{-1}b<0,   (p,q)<0.

The last inequality fixes the future sheet. Every unselected original wall
r satisfies (q,r)=(p,r)+sum_i(B^{-1}b)_i(e_i,r)<0. Thus the construction
produces a finite point of precisely the selected intersection, with all
other wall inequalities checked. It does not assume an arbitrary face is
Coxeter or that an abstract Gram matrix has a chamber realization.

For a connected affine principal set with positive kernel coefficients
u_i, the represented vector nu=sum_i u_i e_i is null. Its pairing with p
is strictly negative, so it is nonzero and future. Every other wall has
nonpositive product with nu. The ideal ray is feasible. Finite volume and
the finite-sided hypothesis identify it with a cusp of the original
chamber. At that cusp, a further incident root has zero pairing with nu;
positivity of the u_i and acuteness force each individual pairing with the
affine set to vanish. This proves the whole-factor conclusion directly,
including at a nonsimple ideal vertex.

An indefinite represented core has a different kind of kernel. Its span
contains a timelike vector. A represented Gram-kernel combination is null
and orthogonal to that timelike vector, whose perpendicular space is
positive definite. It therefore represents the zero vector. The resulting
linear relation genuinely constrains every additional original wall. It
must not be confused with the nonzero affine null vector above.

These arguments use the all-wall acuteness theorem and the standard
finite-volume cusp structure as external geometric foundations. No field,
integrality or nonphysical embedding condition enters.

## 2. Minimal cores and the finite contact alphabets

For opposite sides of a quadrilateral, a positive-definite completion by
the supporting whites would give an actual finite intersection by the
projection argument. A singular semidefinite completion would give an
actual ideal intersection by an affine component. Both contradict the
disjointness of the opposite closed sides. A minimal indefinite subset
must contain both marked sides, because deleting either leaves a subset
of an actual edge diagram.

For minimal indefinite order at least six, the compact Lanner case is
excluded by its order bound. A proper affine connected block is available.
Adjoining one adjacent root already makes it indefinite: its nonzero
contact with the positive affine kernel gives a negative direction.
Minimality consequently makes that block plus one root the entire core.
Its affine rank is between four and eight. The four classical families in
those ranks and the four exceptional affine sources give the stated 24
source types before finite attachment enumeration.

A finite label at least five would occur in a connected proper five-root
subdiagram of such a core. The finite/affine classification forbids that
connected five-root diagram. A dotted pair is already indefinite, and a
parallel pair with a contacting third root is already indefinite. Thus
the high-core attachment alphabet is a consequence of classification;
it is not a numerical truncation of arbitrary angles.

I checked the relevant original statement in [Chein's article, printed
page 4](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf): its principal
deletion criterion and order ceiling have the scope used here. Dotted
order-two and continuous order-three cases still require their separate
analytic treatment; this citation does not discretize them.

One useful explicit check of the union alphabet is worth retaining in
the next editorial pass. A new white--white bond lying in a connected
spherical component of rank at least three has finite label at most five.
The only spherical components admitting unbounded finite labels are the
rank-two dihedral components. If a new bond joined two formerly isolated
white roots from the two opposite cores, each is adjacent to a black mark
in its own connected core. Their new bond would give an adjacent-black
path of length three, already handled by the earlier short-path rule.
Thus that dihedral exception cannot silently introduce an unbounded
alphabet in the residual high/high branch. Cross bonds involving a black
mark similarly give the earlier short path. This observation supports
the finite-alphabet claim; it is not a substitute for checking the
generator's complete matching and cross-bond domains.

## 3. Whole-component extension estimates

For a connected positive-definite norm-two Coxeter Gram of rank k, delete
the selected node and let the remaining component ranks be k_j. Every
component has a native contact of magnitude at least one. Inverse
positivity and induction give Schur cost at least sum_j k_j/(k_j+1),
which is at least (k-1)/k. The selected inverse diagonal is therefore at
least k/(k+1). The rank-one base is 1/2.

For an attached forest of total rank M, contacted in each component, this
implies u^T B^{-1}u>=M/(M+1). A nonsingular indefinite core with inverse
diagonal r>1 gives the residual condition B-r uu^T positive semidefinite,
hence u^T B^{-1}u<=1/r. Combining the two inequalities yields
M<=1/(r-1). The forest includes all indirect tails; bounding just the
first neighboring roots would not prove this assertion.

The manuscript's inverse and residual-vector arguments use the positive
orthogonal complement of a nondegenerate Lorentzian core. Its radical
arguments use the actual zero relation proved above. The text correctly
separates nonnegative-vector monotonicity from arbitrary changes to the
negative index. I found no general sign or indirect-extension error in
these common lemmas. The applicability of each particular seal, including
the claimed permitted contact positions, remains a local certificate
obligation.

## 4. Corner fibers and the strict counting correction

At a cusp with affine ranks r_i, the horospherical link is a product of
simplices and sum r_i=n-1. An edge of the link varies in one simplex
factor. A marked pair in that factor has one free vertex choice in every
other factor, giving exactly product_{j!=i}(r_j+1) actual corner fibers.
This proves the multiplicity V/(r_i+1), rather than an unweighted sum of
factor costs. A proper path inside a larger factor cannot acquire the
whole-factor weight by this argument.

The shape comparison can also be checked directly. Choose a generic
linear function on the closed Klein polytope. Each two-face has exactly
one minimum and one maximum. Its corner at a vertex corresponds to an
edge of that vertex figure; that edge is uncut precisely at those two
extrema. Therefore I=sum_v E_v and 2f_2=sum_v E_v^0, giving
D5=sum_v(5E_v^0/2-E_v).

Orient the simple vertex figure by the induced generic function. Assign
an uncut edge below the section to its upper endpoint, and one above to
its lower endpoint. A vertex with j descending edges receives j or d-j
assignments, at least min(j,d-j). The index distribution is its h-vector.
The original [Khovanskii paper, Theorem 9 and its proof, printed page
49](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf) supplies
the incidence reduction in exactly this vertex-figure setting. The index
interpretation is the separate standard simple-polytope input cited in
the manuscript; it is not a script output.

In ambient dimension 58, the figure dimension is 57. For a simplex,
V=58 and G=sum_{j=0}^{57}min(j,57-j)=812, so its curvature is
(5G-57V)/2=377. At either global extremum the whole figure is uncut.
Symmetry of its h-vector gives

    E-G=sum_{j=29}^{57}h_j(2j-57)>=1+3+...+57=841.

Each relevant h-polynomial is a simplex polynomial or a product of them,
so every coefficient in its support is at least one. The global minimum
and maximum are distinct, even when ideal. Their two extra contributions
are (5/2)*2*841=4205, counted once. This proves the claimed strict
correction without assuming that either extremum is finite. The payment
comparison consequently handles the all-finite case as well.

## 5. Actual framing walls in the unbounded tail

I directly checked [Nikulin's corrected paper, Section 4, printed page
10](https://arxiv.org/pdf/0708.3991v2). The narrow-face sentence explicitly
extends its framing bound to general hyperbolic closed convex polyhedra.
Its arithmetic surrounding discussion does not impose arithmeticity on
that sentence. Reversing the form's sign convention yields -(x,y)<14
for outward roots of square two.

At a cusp in the selected facet, deleting one node in each affine factor
gives an actual outgoing edge. The chosen omissions avoid the facet
normal. Repeating an omission is allowed in the affine rank-one factor,
whose modified target is zero. At a finite other endpoint a terminating
wall makes a positive-definite pair with the facet normal. At an ideal
other endpoint, choose the missing wall in the facet's factor when its
rank is at least two; in the rank-one case choose a missing wall in
another factor. The dimension range guarantees another factor exists.
The projection argument proves that such a positive-definite pair really
defines an actual ridge, so the terminating wall is in the framing set
to which narrowness applies.

The terminating wall cannot be an old omitted wall, which remains strict
away from the starting cusp along the edge. Its old contact row is
nonnegative and every nonzero native magnitude is at least one. Kernel
normalization makes the affine heights positive. A zero height in one
factor would force every contact to vanish, leaving a square-two vector
orthogonal to the full cusp span, an impossibility.

The two endpoint Schur budgets are separately at most two. They are not
replaced by residual budgets of 3/2 before selecting a failing factor.
This order is essential to the later packing argument.

## 6. Scale cancellation and the two height ranges

Let the positive kernel vectors of the factors be c_a=k_a c. In a null
coordinate decomposition, put b=-(x,c), b'=-(y,c). The squared Euclidean
quotient difference gives

    -(x,y)=1/2 sum_a sigma_a tau_a E_a-q-q^{-1},
    sigma_a=k_a b, tau_a=k_a b', q=b/b'.

The k_a factors cancel in each Euclidean summand. This works for a
product cusp without rationality or commensurability. With
M=max(q,q^{-1}) and W_a=min(sigma_a,tau_a)^2 E_a, it gives
-(x,y)>=M(sum W_a/2-2).

The ordinary and modified targets sum to (2N-2)/9, N=n-1. The supplied
near-range loss bound four yields M(N-37)/9. The far-range loss bound
16/3 yields M(N-43)/9. At N=163+t, the first is at least 14+t/9.
For M>=21/20, the second is at least 14+7t/60. Both contradict the
strict framing bound, including the boundary n=164.

I independently checked the classical B/C/D contact table used at the
ORIGINAL budgets. Its point-potential polarization gives W>=E and
unit inverse-diagonal costs bound the weighted linear losses. Every
nontrivial contact has loss at most its unit cost; the marked trivial
ones already cost more than the target. The omission distances r/4,
(r-1)/4, r/2 and (r-1)/2 exceed the respective targets. This coefficient-
one estimate is needed before the exceptional first-failure dispatch.
The later coefficient-4/3 estimate alone would not establish it.

The small-A and exceptional support/KKT propositions and the ordered
packing maxima are still explicit computed inputs in this audit. I
checked how those inputs imply the stated two loss bounds, including
the two original budgets and the single modified-factor allowance. I
did not regenerate every support, KKT witness or packing state.

## 7. The compact-facet alternative

If the chosen facet has no ideal vertex, its closed Klein face is
compact. It need not be a Coxeter chamber. All its two-faces and their
corner roots are nevertheless actual faces and roots of the ambient
chamber. Minimal cores give two distinct paths of length at most eight
for triangles. For quadrilaterals, two opposite minimal cores have at
most eight white roots each. They overlap or have a joining bond,
because two orthogonal indefinite spans would give two negative
directions. A shortest cross path has length at most 17 and joins
adjacent black sides.

Unit ambient payment at spherical distance at most 17 therefore pays
each required face. Its finite-vertex cost is at most 17n: small
components use their complete-pair bound, and large paths or D forks
have at most their rank many pairs at any fixed distance. Simplicity
comes from the ambient finite vertices, not from an assumed Coxeter
structure on the facet. The symmetric h-vector comparison gives the
stated strict face-incidence lower bound. The arithmetic differences
at facet dimensions 143+2t and 144+2t are respectively
t^2/2+36t+3/2 and t^2/2+73t/2+19. They are positive for t>=0.
This branch excludes n>=144, as needed in the unbounded-tail split.

## 8. Remaining audit obligations and manuscript consequence

No defect was found in the common geometric arguments audited above.
That verdict does not certify the complete theorem independently of its
other inputs. In particular, the following still need their separate
source-bound audits:

* every matching, cross-bond and whole-component extension domain in
  the distributed local routing atlas, with current chosen outputs;
* the complete typed context split, finite domination and all residual
  certificate application premises;
* the identification of each printed finite/affine table with one
  simultaneous rule, retaining every earlier consumer and debt;
* all dimension-specific rules and terminal domains for 59--163;
* the exceptional/small-A KKT and packing propositions in the all-rank
  tail, including the precise checker implementations.

The separate fresh algorithm audit discharges a coherent low-eight and
abstract Bellman subdomain, not that entire list. It also records real
reproduction defects that must be fixed before the successor entry is
advertised as verified. No scientific, Lean or TeX process ran for this
report, and no old successful module was rerun.

For consolidation, retain the theorem 57 and the finite/ideal split.
Later dimension-57 financing research toward 56 is not a missing part of
this proof and should remain outside its TeX. The present source needs
an eventual guarded-command/reproducibility delta and any clarifications
arising from the remaining audits; it does not need an appendix of the
unfinished newer case search.
