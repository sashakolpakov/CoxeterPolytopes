# Independent adversarial review: three maximum-height descents

20 October 2026. Full mathematical/source review, including independent
algebra and primary-source reading. No scientific job, checker, optimizer,
Lean or build was run in this lane.

Reviewed in full:
`results/right_angled_dimension8_descending_height_proof_20261020.md`,
SHA256 `9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe`.

**Verdict: PASS.** The proposed proof excludes every finite-volume
right-angled nine-polytope and hence gives the universal dimension ceiling
eight, subject to its explicit standard geometric facts and published
six-dimensional facet minimum. I found no omitted vertex-type case,
projection-sign error, intrinsic-dimension mismatch or spherical equality
gap. Existence of dimension-eight examples is a separate external input
for sharpness, not part of this nonexistence argument.

## 1. Actual projection and the role of the maximum

The source uses a strictly interior finite unit future point q. Thus
every outward unit normal has positive finite height `t_m=-<q,m>`.
There are finitely many actual facets, so the largest height t is attained
at an actual facet n. The displayed projection has square minus one,
is orthogonal to n, and is future: its numerator has product
`-1-t^2<0` with q.

All-wall acuteness applies to every different actual wall m, not merely
to ridge neighbors. It makes

    <q+t*n,m>=-t_m+t*<n,m><0.

Consequently the projected point lies in the relative interior of the
actual facet. This is stronger than a formal projection to its supporting
hyperplane and is enough to repeat the argument without degeneracy.

Each intrinsic facet of the selected facet comes from a genuine
orthogonal ridge neighbor m. Its ambient normal already lies in n-perp
and remains unit, so its intrinsic height is exactly
`t_m/sqrt(1+t^2)`, with no extra normalization or angle factor. The maximum
condition bounds it by `t/sqrt(1+t^2)`.

The argument repeats using the new intrinsic *maximum*, not the old
ambient maximum or a maximum only over visible walls. If the first height
is A0, induction gives the finite upper bound

    A_j<=A0/sqrt(1+j*A0^2)<1/sqrt(j).

Thus three actual descents from dimension nine to six give strict
`0<t_i<1/sqrt3`. Ties in any maximum do not remove the strict inequality:
it comes from the finite initial A0, not uniqueness of the maximizing
wall. All faces used have dimension at least six, so there is no mistaken
inheritance of right angles at ideal polygon vertices.

## 2. Actual intrinsic walls give the correct spherical code

The resulting six-face is finite-volume and right-angled by the standard
Coxeter-face inheritance fact. In its Lorentz seven-space, its actual
intrinsic wall normals decompose as

    m_i=(t_i,sqrt(1+t_i^2)*u_i),   u_i in S^5.

All-wall acuteness and the strict height bound give pairwise products
strictly below one quarter. The spatial directions are distinct: equal
directions would give a strictly positive Lorentz product between their
unit normals. Thus the number of code points is exactly the number of
actual facets of this six-polytope, with no multiplicities or candidate
walls counted as actual ones. The Euclidean code dimension is six,
not five or seven.

## 3. Independently expanded tensor certificate

I rederived the two tensor products in the source. For unit u,v in R6
and x=u.v, the traceless quadratic tensor has product `x^2-1/6`.
For `S(u)=sym(u tensor I)`, the two contractions against a cubic tensor
are both `3x`, while `S(u).S(v)=24x`. Consequently

    (u^3-S(u)/8).(v^3-S(v)/8)=x^3-3x/8.

Combining these products with the positive coefficients in the source
gives exactly

    (x+1/2)^2*(x-1/4)
       =1/16+(3/8)x+(3/4)(x^2-1/6)+(x^3-3x/8).

The ordered-pair sum is therefore N^2/16 plus three nonnegative squared
norms, exactly as written. This is a finite algebraic identity valid for
arbitrary real vectors. It does not assume the addition theorem, a
numerical spherical-code table, classification of an equality code, or
positivity merely observed at sample vectors.

The polynomial is nonpositive for every distinct-pair product below one
quarter and has value 27/16 on the diagonal. Thus N<=27. If N=27, the
two inequalities force every off-diagonal term individually to vanish.
Strictness excludes the root one quarter, leaving only minus one half.
Four such unit vectors would have squared sum `4-6=-2`, impossible.
Hence the exact integer upper bound is 26. The endpoint argument is
essential and is explicitly present; simply reporting the weak bound
27 would not suffice for the theorem.

## 4. Published lower bound checked and the final scope

I independently checked the actual text of Alexandrov's
[*On ideal vertices of right-angled hyperbolic polyhedra*](https://arxiv.org/html/2303.09533v1#S4).
Section 4 explicitly states the 27-facet lower bound for six-dimensional
finite-volume right-angled polyhedra, citing Dufour's Lemma 1 and
Proposition 4. Section 6 records the underlying two statements as
Propositions 6.4 and 6.5 and displays their dimension-six consequence.
Thus the needed 27 is a facet count in dimension six; it is not a cusp
count or a bound restricted to entirely ideal polytopes. Like the author,
I did not retrieve the original Dufour full text in this pass. This is an
explicit bibliographic limitation, not a replacement of the input by an
unverified numerical database.

The contradiction is now between at most 26 actual facets and at least
27 actual facets of the same finite-volume right-angled six-face. Higher
dimensions have actual nine-faces, to which the same argument applies.
Neither the new ideal-four theorem, the previous mixed-ten theorem,
arithmetic restrictions nor any opposition-graph estimate is needed.

The scope is therefore a full written geometric exclusion, relative to
the stated classical acuteness/face-inheritance inputs and published
facet minimum. This review makes no claim that Lorentz geometry or the
published theorem has been formalized in Lean, and no execution is being
used to infer universal geometric validity.
