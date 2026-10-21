# Exact strict feasibility over \(\mathbf Q(\phi)\)

## 1. Result

Let \(K=\mathbf Q(\phi)\) with its distinguished real ordering, let
\(q(x)=x^TQx\) have inertia \((1,n-1)\), and let \(p\in K^n\) satisfy
\(q(p)>0\).  Suppose finitely many oriented covectors \(r_i\in(K^n)^*\)
satisfy \(r_i(p)>0\).  For any selected covector \(r_0\), the repository now
decides exactly whether

\[
 r_0(x)=0,\qquad r_i(x)>0\ (i\ne0),\qquad q(x)>0,
 \qquad B(p,x)>0.                                      \tag{1}
\]

If (1) is feasible, it returns a point in \(K^n\) and verifies every equality,
sign, and quadratic value exactly.  If it is infeasible, it returns either an
exact phase-I linear infeasibility certificate or an exact KKT certificate
that the global Lorentzian maximum on the normalized face is nonpositive.
Thus linearly redundant walls and linearly exposed but causally inactive
walls are not conflated.

The implementation is in
`src/coxeter_polytopes/geometric/golden_halfspaces.py`:

- `find_nonnegative_span_phase_one` decides exact conic membership;
- `decide_lorentzian_strict_face` decides (1);
- `audit_golden_halfspace_certificate` independently checks facet witnesses,
  positive dependences, and causal-inactivity decisions.

The Shimada reducer in
`scripts/explore_shimada_disc11_global_eigen_nef.py` uses these routines and no
floating-point linear program.

## 2. Exact phase I

For a finite homogeneous system

\[
r_0(x)=0,\qquad r_i(x)>0,\qquad B(p,x)>0,               \tag{2}
\]

strict feasibility is equivalent to feasibility after replacing every strict
inequality by a lower bound \(\ge1\).  Indeed, a solution of (2) can be
multiplied by the reciprocal of the smallest positive pairing.  Write each
free coordinate as \(x_j=x_j^+-x_j^-\), with
\(x_j^\pm\ge0\), and introduce slack, surplus, and artificial variables.

The resulting phase-I problem is solved over \(K\), using exact sign
comparison and Bland's entering and leaving rules.  Bland's rule prevents
cycling, so this is a finite decision algorithm rather than a denominator
search.  The terminal primal and dual vectors are checked for feasibility and
equal objective.  Objective zero is equivalent to feasibility of (2); a
positive optimum is a Farkas certificate of infeasibility.

The same kernel decides

\[
       v=\sum_j\lambda_jv_j,\qquad \lambda_j\ge0.       \tag{3}
\]

Zero coefficients are removed from a feasible basic solution, leaving a
strictly positive exact dependence that the half-space verifier checks
coordinate by coordinate.  This replaces both the former numerical support
proposal and the exponential support fallback in the Shimada traversal.

## 3. Why the timelike condition is also decidable

A basic solution of (2) need not be timelike.  This is a genuine issue, not a
numerical nuisance.  Put

\[
t(x)=B(p,x),\qquad c=q(p)>0,
\]

and normalize the strict linear witness to \(t(x)=c\).  Every point on this
affine hyperplane has a unique expression \(x=p+z\), with
\(z\in p^\perp\), and

\[
q(x)=c+q(z).
\]

Because \(Q\) has inertia \((1,n-1)\), its restriction to \(p^\perp\) is
negative definite.  Therefore \(q\) is strictly concave and coercive on
\(t=c\).  It has a unique global maximum on the closed polyhedron

\[
r_0(x)=0,\qquad r_i(x)\ge0,\qquad t(x)=c.               \tag{4}
\]

The strict phase-I point gives Slater's condition.  Hence the KKT conditions
are necessary and sufficient.  A minimal normal-cone expression has an
independent support of at most \(n-2\) active inequalities, modulo the two
equalities in (4).  The code enumerates those finitely many supports, solves
the bordered KKT systems over \(K\), and checks

- primal feasibility;
- nonnegative active multipliers;
- both equality constraints;
- exact stationarity; and
- the exact value of \(q\).

If the maximum is nonpositive, (1) is impossible.  If it is positive, mixing
the maximizer with the strict phase-I point makes every inequality strict
while preserving positive quadratic value.  The implementation uses the
unbounded exact sequence \(2^{-k}\); positivity at \(k=\infty\) proves that a
finite term succeeds.

The tests include both indispensable edge cases for
\(Q=\operatorname{diag}(1,-1,-1)\): a strict future linear face whose exact
maximum is \(-27/4\), and a spacelike phase-I basic point whose exact maximum
is \(27/16\) and whose mixed strict witness has square \(63/64\).

## 4. Complete half-space reduction

The reducer first partitions primitive wall restrictions into exact positive
proportionality classes.  It then uses the following certificate-preserving
fixed point:

1. try the cached exact Lorentzian projection of the known interior point;
2. certify unresolved redundancies by (3) against the facets already found;
3. apply the complete phase-I/KKT decision to one remaining class;
4. after finding a new facet, retry (3), since the enlarged facet set can make
   several other classes immediately redundant; and
5. record a genuine causal-inactivity certificate separately when a
   linearly exposed face has nonpositive Lorentzian maximum.

The fast paths affect runtime only.  Every unresolved class reaches the finite
decision in step 3, so the reduction is complete.  On the physical Shimada
\(D_{68}\) state, exact projection finds all eleven facets and phase I proves
all eight remaining proportionality classes redundant; there are no causal
inactive classes.  The reduction itself takes under one second on the current
development machine.

## 5. No crossing or gallery cutoff

The earlier implementation tried only the denominators
\(2,4,\ldots,1024\) when crossing an internal facet and stopped a chamber
gallery after 100 steps.  Both bounds have been removed.

For a strict relative-interior facet point, local finiteness supplies an open
transverse interval that crosses the selected wall and no unrelated wall.
The unbounded dyadic sequence therefore reaches a valid opposite point after
finitely many trials.  To locate that point, the atlas walk crosses any wall
separating the current chamber from it.  A compact segment meets only finitely
many walls of the locally finite induced arrangement, so this gallery is
finite.  A repeated physical state is treated as an exact inconsistency, not
as an undecided timeout.

Consequently a *specified* Shimada face and a *specified* transverse crossing
no longer have a numerical-support, denominator, or gallery-length caveat.
This does not by itself make the decorated breadth-first traversal complete.
A restricted facet can be subdivided by the auxiliary Shimada arrangement;
different relative-interior points of the same restricted facet can therefore
enter different opposite induced chambers.  Completeness of the global graph
requires an exact enumeration of every cell in that finite facet subdivision,
not one witness per restricted facet.  The current traversal is being repaired
at precisely this separate star-enumeration gate, so no finite-radius output is
promoted to a complete centralizer ball on the strength of this feasibility
theorem alone.

## 6. Verification

The focused tests are in `tests/test_golden_halfspaces.py`.  They cover:

- feasible and infeasible phase-I systems with primal-dual equality;
- exact golden conic dependences;
- linear infeasibility of a future face;
- positive and nonpositive Lorentzian KKT maxima;
- rejection of a timelike witness on the opposite Lorentz sheet;
- separation of causal inactivity from conic redundancy; and
- the pre-existing independent witness/dependence certificate audit.

The current focused result is `11 passed`.  The external 45 MB Shimada atlas
is not vendored; the centralizer production test remains conditional on its
SHA-pinned supplement.
