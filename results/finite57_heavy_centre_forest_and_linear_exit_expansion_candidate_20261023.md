# Finite57: heavy centre edges form a forest with linear exit expansion

23 October 2026. **New complete written candidate by root, pending full
different-author review.** No mathematical worker, enumeration, symbolic
process, compiler, Lean or TeX ran. The existing rule, allocations, papers
and dimension ceilings remain unchanged. This strengthens the accepted
heavy-component escape theorem; it does not supply the missing receiver
capacity or prove dimension56.

## 1. Exact inputs and statement

Work in the ORIGINAL finite-volume hyperbolic Coxeter57-polytope P.
Actual outward wall normals are independently scaled to square2. A
heavy vertex has ENTIRE finite forest 2D17+D23. There are three distinguished
actual edges, obtained by releasing one fork centre in each component.
Let G have all heavy vertices as vertices and the actual centre edges
whose two endpoints are heavy as edges. It is a finite simple graph.

The following literal inputs are used at their stated interfaces:

| Input | SHA256 | Restricted use |
|---|---|---|
| manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 | All-wall acuteness, finite simplicity, actual vertex/edge interpretation and norm-two finite diagrams |
| results/finite57_complete_finite_centre_recipients_and_heavy_component_escape_candidate_20261023.md | 994d3b9ed36a32e8698610140ee656b40a88507c7f0a7b6b1978faef21b00186 | Complete finite-centre classification: row9 is the ONLY heavy completion; exact W support; complete finite reverse-port bounds |
| results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md | 1c2e887918d4277efffd9c568aed5d9e203a61dce2ca7c98eca96665831792f9 | Accepted written scope of the preceding classification and allocation |

Root derived the metric forest proof below after reading the entire
finite-centre source. The pinned acceptance records Replay as author of
that inherited classification and Foundations as its outside reviewer.
Composition checked the proposed long-chain mechanism in discussion.
This is not yet an outside review of the complete new theorem.

**Theorem.** G is a forest. If S is any nonempty set of heavy vertices
and c(S) is the number of components of the induced forest G[S], then
the number b(S) of ACTUAL centre-edge incidences leaving S is exactly

    b(S)=|S|+2c(S).                                      (1)

Leaving S includes an endpoint at another heavy vertex outside S or a
nonheavy finite or ideal vertex. In particular a whole heavy component
of H vertices has exactly H+2 edges to nonheavy vertices. If all those
endpoints are finite, they include at least ceil((H+2)/3) distinct vertices.
The final assertion uses the complete finite-centre union's reverse-port
bound, including collisions between the two rank cases.

## 2. Actual heavy-to-heavy edge lengths

Take an actual edge p--q of G. Delete its old fork centre a. The FULL
56-root support W is

    A_s+2A1+D_j+D_k,
    (s,j,k)=(14,17,23) or (20,17,17).

At p, a contacts one endpoint of A_s and both isolated A1 roots by
simple3, and is orthogonal to the inert factors. Since q is heavy,
row9 of the complete finite-centre table forces its new wall a' to
have EXACTLY the same type of contacts. It may select either endpoint
of the retained A_s. The product a.a' remains unrestricted apart from
actual all-wall acuteness a.a'<=0. In particular it is not set to zero
and no arithmetic restriction on it is used.

Let B be the positive definite norm-two Gram of ALL W, and let h,h'
be the positive contact columns of a,a' to W. Orthogonally project the
actual normals a,a' to the Lorentz two-plane W-perpendicular. Their
projected vectors u,u' have

    u^2=(u')^2=2-h^T B^-1 h=1/(s+1).                    (2)

Indeed the A_s endpoint inverse diagonal is s/(s+1), and each isolated
A1 contributes 1/2. The two contact columns agree on the isolated A1s.
Their cross contribution on A_s is s/(s+1) if the endpoints agree or
1/(s+1) if they differ. Consequently

    h^T B^-1 h' >= 1+1/(s+1),
    u.u'=a.a'-h^T B^-1 h' <= -(s+2)/(s+1).              (3)

The rays p and q lie on the ACTUAL hyperbolic line in W-perpendicular;
u and u' are its two outward endpoint normals. For endpoints at distance
L, their normalized product is -cosh L. This follows directly by taking
p=(1,0), q=(cosh L,sinh L) in the unit hyperboloid: the two outward unit
normals can be (0,-1) and (sinh L,cosh L). Equations(2)--(3) therefore give

    cosh L=-(s+1)*(u.u') >= s+2 >= 16.                 (4)

This is an estimate for the length of an ACTUAL edge of P. It does not
infer a new endpoint from an abstract completion or replace the endpoint
by a vertex assumed farthest from an interior point.

At a heavy vertex, any two of the three centre-edge directions are
orthogonal. They lie in different orthogonal D components of its full
spherical diagram; the corresponding dual directions are orthogonal
as well. Thus a nonbacktracking path in G turns through angle pi/2 at
every interior vertex, in the tangent space of the ORIGINAL ambient
hyperbolic space. No lower-dimensional face is assumed Coxeter.

## 3. Long orthogonal broken geodesics cannot close

Here is the geometric argument in any real hyperbolic dimension. Suppose
successive geodesic segments have cosh length at least16, and consecutive
segments meet at angle pi/2. Such a nonbacktracking chain cannot close.

Fix its first point z0. Inductively, at a later point z_k the direction
back to z0 makes an angle less than pi/4 with the direction back to
z_(k-1). This is zero at k=1. Orthogonality of the next segment and
the triangle inequality for tangent-sphere angles imply that the angle
B between z_k--z0 and z_k--z_(k+1) lies strictly between pi/4 and
3pi/4 (or equals pi/2 at the first step). In particular the next point
cannot be z0, since that would make B zero.

Let R=d(z0,z_k), L=d(z_k,z_(k+1)), and D=d(z0,z_(k+1)). Put
C_L=cosh L, S_L=sinh L. The hyperbolic cosine rule gives

    cosh D=cosh R*C_L-sinh R*S_L*cos B
           > cosh R*(C_L-S_L/sqrt2) > 4*cosh R.          (5)

For the last step, S_L<C_L, 1/sqrt2<3/4 and C_L>=16.
Thus the distance from z0 strictly increases.

If C is the angle at z_(k+1) in this triangle, the cosine and sine
rules give the exact identity

    cot C=(coth R*S_L-C_L*cos B)/sin B.                 (6)

Since R>0, cos B<=1/sqrt2<3/4 and 0<sin B<=1, its numerator is
greater than S_L-3C_L/4. Moreover C_L>=16 implies
S_L/C_L>=sqrt255/16>15/16. Therefore

    cot C > 3C_L/16 >= 3 > 1.                           (7)

It follows that C<pi/4, establishing the induction. Each triangle
lies in its own hyperbolic two-plane; no common plane for the chain
is needed. In particular the proof applies to skew paths in H57.

After m edges, the quantitative bound is

    cosh d(z0,z_m) >= 16*4^(m-1),

with strict inequality for m>=2. A first repeated point would give
a closed subchain with the same right-angle and length hypotheses;
applying the argument from that point excludes the repetition.

## 4. Forest and exit count

A simple cycle in G would be a closed chain of the kind just excluded,
by(4) and the orthogonality of its distinct centre directions. Hence G
has no cycle and is a forest.

Every heavy vertex has precisely three actual centre edges in P. Finite
volume and the already classified actual-edge setting supply their
other endpoints, finite or ideal. For any finite nonempty S, each edge
internal to G[S] contributes two incidences and each leaving centre
edge contributes one. An induced forest with c(S) components has
|S|-c(S) edges. Thus

    3|S|=2*(|S|-c(S))+b(S),

which proves(1). In a full connected heavy component the remaining
endpoints are nonheavy, giving H+2 actual exits. If these are all
finite, each has at most three reverse ports in the COMPLETE finite
table, so the number of distinct receivers is at least ceil((H+2)/3).
No bound on ideal reverse multiplicity is needed or claimed.

## 5. What this changes, and what remains unproved

The accepted escape argument only forced at least one nonheavy exit
per heavy component. Equation(1) forces linearly many actual exits,
even though all inside moves may restore heavy diagrams. This is a
new geometric restriction on the residual graph and its subsets.

It does NOT give eight units of available capacity per exit. Some
finite exits have held charge above364, others have little surplus,
several sources may share a receiver, and some receivers already carry
older debt. If S omits previously discharged heavy vertices, its
boundary count includes edges to those vertices; these are not new
unencumbered receiver capacity. No source is discharged by this note.

The original C0 rule and every accepted earlier allocation are retained.
This theorem is a possible input to a future quantitative capacity or
actual-angle compatibility argument. It is not a Hall completion,
a universal finite cap364, dimension56, a new example, a Lean theorem,
or evidence of a mathematical execution. Full outside/root acceptance
is required before promoting this new candidate.
