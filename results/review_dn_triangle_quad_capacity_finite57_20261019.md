# Independent triangle, quadrilateral and exact-source review

19 October 2026. **Written mathematical and source review: PASS within the
scope below.** No checker, matrix program, enumeration, optimizer, Lean job
or build was executed in this review. Saved literal orbit records were
streamed and reformatted for reading. No global compact completion, complete
face classification, valid full fork refund or dimension improvement follows.

## Reviewed bytes

* `results/arithmetic_dn_golden_triangle_quad_capacity_20261019.md`, SHA256
  `74f07fee8ca7240f53e9ed8f8ffe33a4539c60c2de6bd2f98bb2f49c7b8dcf46`.
* `results/arithmetic_dn_golden_opposed_fork_quad_survivor_20261019.md`, SHA256
  `06900f90735d59dfd1b2a85e17d964438749124f5d2e7a5c710b2b2eda4f06a8`.
* `results/arithmetic_dn_silver_opposed_fork_quad_survivor_20261019.md`, SHA256
  `5d82def07c8543f45a00fffd17301f57c58f36065a5e3c26920d5ad6d2cf6f42`.
* `scripts/check_arithmetic_dn_opposed_fork_models_20261019.py`, SHA256
  `185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9`.

The first report is the main independently reviewed capacity theorem. The
checker is approved here as source. After this source read, the coordinator
reported separate successful generation and replay, retained as
`results/arithmetic_dn_opposed_fork_models_20261019.json`,
`results/arithmetic_dn_opposed_fork_models_generation_20261019.json` and
`results/arithmetic_dn_opposed_fork_models_replay_20261019.json` (about 0.91s
and 55.8MB per bounded run). Those executions are separate evidence, not jobs
performed by this reviewer. The golden and silver reports supply the local model
definitions; their global nonexistence/existence is not asserted by this
review.

## 1. Actual-wall capacity is correctly normalized

For unit old whole-D_n roots, put N=n-2 and write the fork as a, its
whiskers as b,e, and its long arm as l1,...,lL, where L=n-3. The inverse
entries in Section 2 of the capacity note are correct. In particular the
two short dual weights have squares n/2 and mutual product (n-2)/2.

All actual new outward walls have strictly negative product with the old
finite vertex. Compact Coxeter chambers are simple, so a distinct wall
cannot also pass through that vertex. In the one-dimensional timelike
perpendicular of the old positive span, the new walls therefore have the
same orientation. Their product is the cross projection MINUS the product
of their positive timelike magnitudes. One cannot independently reverse
that sign for a second wall.

This proves the note's formula (1). The equal-profile difference calculation
then gives 1-h^2: the two singleton triangle walls have product
1-x^2=(5-sqrt(5))/8>0; the two opposed-fork quad walls have product 3/4>0.
The walls are distinct because their old contact rows differ. Both are
actual acuteness exclusions, conditional on the walls belonging to one
polytope.

I separately derived the mixed case before reading its written expansion.
With x=(1+sqrt(5))/4 and 4x^2=2x+1, the three quantities are

    q_triangle-1 = (N+2)x^2/2-1,
    q_quad-1     = 2Nx^2+(N+2)x+(N+2)/8,
    cross        = Nx^2+(N+4)x/4.

Their exact difference of squares is

    cross^2-(q_triangle-1)(q_quad-1)
       = [22(N+1)x+7N+3]/16 > 0.

The cross term and both residual energies are positive at N=17,18. Thus
the forced actual wall product is positive. This independently confirms
the union capacity: at most one old fork/whisker angle can support any of
the two specified golden exceptional profiles. Field separation from the
silver model is correctly based on the independent old-D_n-plus-new-wall
determinant, not just containment of a quadratic subfield.

## 2. The singleton golden triangle is a valid local obstruction

This is useful additional detail for the triangle audit. Let x=cos(pi/5)
and adjoin a unit wall t with b.t=-x and all other old contacts zero.
Let W={e,l1,...,lL}. Its three finite corners are

    W+a+b: D_n;
    W+a+t: A_(n-1)+A1;
    W+b+t: A1+A_(n-3)+I2(5).

The middle corner marks the isolated t and the SECOND node a of the path
e-a-l1-...-lL. In the white perpendicular Lorentz space the projected
three-wall Gram, in order a,b,t, is

    [ s,  -1/2,  0 ],       s=n/[4(n-2)],
    [-1/2,  1,  -x ],
    [ 0,   -x,   1 ].

Every two-root corner block is positive definite. The vector
a_bar+(5/4)b+t has strictly negative product with all three walls:
the products are s-5/8, 3/4-x and 1-5x/4. Its square is consequently
negative. Projection of that interior timelike point perpendicular to any
adjacent pair remains timelike and feasible, since the pair inverse is
entrywise nonnegative. Thus all three finite vertices are actual vertices
of a compact triangle in this white perpendicular space.

For norm-two roots the complete Gram has diagonal 2, old simple entries
-1, and sole new entry b.t=-(1+sqrt(5))/2. It is integral over the ring of
integers of Q(sqrt(5)). The old-D_n projection cost is

    q=n(3+sqrt(5))/16 > 1,
    q'=n(3-sqrt(5))/16 < 1       (n=19,20).

At n=20 the latter strict inequality is equivalent to 121<125; n=19 is
smaller. Hence the complete physical Gram has signature (n,1), and its
nonidentity conjugate is positive definite. The two-cycle generates
Q(sqrt(5)); an embedding of a larger totally real ground field fixing this
subfield would preserve the negative physical determinant ratio, so an
actual ambient arithmetic chamber would have exactly that field. These
are admissible local arithmetic data, not a completed finite-volume or
compact n-dimensional chamber.

The original triangle charge is exactly (1,0,1). Its I2(5) pair receives
one intrinsic unit. There is no H4 component or relevant same-hub pair.
The middle second-node A1 cross bonus is absent from the saved golden
registry. Directly read literal bonds in BOTH the A18 and A19 records give
the near-end path

    marked 1 -- 5 -- 6 -- 4(endpoint).

The other neighbor of 1 is 11. Thus their sole saved large-path rooted
orbit is three edges from the endpoint, whereas the required orbit is one
edge from it. The actual router keys the complete rooted component and
returns zero on an absent key. The unchanged saved bonuses therefore do
not pay the missing unit if this particular D_n triangle unit is deleted.
This is an obstruction to that proposed local charge proof; it is not a
counterexample to the existing compact20 theorem, which retains the unit.

## 3. The proposed sqrt(3) triangle uses a different current rule

For comparison, if c is perpendicular to W and a.c=-1/2, b.c=-beta, its
old-D_n projection cost is

    q=[n beta^2+2(n-2)beta+(n-2)]/2,
    q-1=(beta+1)(n beta+n-4)/2.

At beta=sqrt(3)/2 this becomes

    q =(7n-8)/8+(n-2)sqrt(3)/2,
    q'=(7n-8)/8-(n-2)sqrt(3)/2.

For n=20, q'=33/2-9sqrt(3) is strictly between zero and one; for n=19 it
is 125/8-17sqrt(3)/2 and has the same property. Norm-two entries -1 and
-sqrt(3) are integral. The two full-D_n corners and the whole-I2(6)
corner are locally arithmetically possible.

Nevertheless `results/fresh_arithmetic_high_label_20261016.md`, Section 4,
gives a further unit to a whole I2(m) component with m>=6 in the actual
Q(sqrt(3)) branch. That corner alone pays two. Substituting the golden
distance-two baseline would manufacture a spurious obstruction to the
complete field-specific proof.

There is also a short arithmetic classification for this exact two-D_n
triangle profile. For beta=cos(pi/m), m>=3, its black three-cycle makes
A=2(beta+1) an algebraic integer. At any nonphysical conjugate, finite
corner positivity and q'<1 force 0<A'<8/n. The physical A lies in [3,4).
Degree at least three gives a nonzero positive integer norm less than
4(8/19)^2<1. In degree two its norm is 1 and its trace is 4, giving
A=2+sqrt(3), hence m=6. Degree one gives m=3, whose three black walls
form an affine triangle and violate compactness. The right m=2 case has
an induced affine D4 on a,b,e,c,l1. Thus this exact profile leaves only
the already paid m=6 alternative. This argument does not classify every
other contact pattern of a third triangle wall.

## 4. Exact checker source and remaining interface

The reviewed pure-Fraction quadratic field operations, sign test and
integer-ring test are correct for discriminants 2 and 5. Symmetric Schur
elimination computes inertia and rejects an unsupported zero pivot with
nonzero residual row; its fixed old-positive-first ordering is appropriate
for these explicit models. Full-support kernel checks at both embeddings,
the white projection, complete finite corner types, timelike center and
the twin/mixed projection identities agree with the reports.

The checker records the corner baseline tuple as a literal. The tuple is
justified here by the checked corner diagrams and marked ordering; it is
not independently recounted by a separate weight-router routine. The
checker also does not prove universal wall feasibility, the field list,
all proper-submatrix statements for arbitrary models, complete consumer
coverage, or global completion. In the four fixed models its conjugate
positive semidefinite corank-one matrix and full-support kernel imply
positive definiteness of every proper principal submatrix by elementary
linear algebra. This written implication should not be confused with an
exhaustive model census.

The capacity theorem gives a possible conditional activation interface:
retain a fork unit on any angle supporting one of the specified exceptional
profiles. Their union uses at most one of the two angles at each whole
D_n vertex. A full valid refund still requires proving that EVERY affected
triangle and quadrilateral is paid, including all other actual contact
patterns, arbitrary spectator roles, and simultaneous changes at other
corners. No such completeness is claimed here.

Even a complete argument saving one unit per whole D20 vertex changes its
bare cost only from 38 to 37. Excluding dimension20 by this charge method
needs a uniform vertex bound at most35 and the other rank-20 forests too.
The new capacity does not close that remaining gap.
