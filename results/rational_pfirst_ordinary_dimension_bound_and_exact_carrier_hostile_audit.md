# Hostile audit: ordinary P-first dimension bound and exact carrier

**Audit status (2026-08-11): PASS after repairing the positive-sheet order
in the timelike branch.**  The result is scoped to the standard marked
canonical vector on \(I_{1,n}\), \(n\geq3\), and facets of physical square
\(-1\) or \(-2\).  No claim is made for a twice-primitive \((-4)\) wall.

This audit checks
[the combined theorem](rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md)
and its two pinned positive-branch dependencies.

## 1. The spacelike escape is genuinely homogeneous

For the finite set of oriented walls \(r_j\), put

\[
                 \mathcal W=\{x:x.r_j\geq0\}.
\]

Let \(A\) be an interior positive vector and \(w\in\mathcal W\) have
\(w^2<0\).  The exact two facts used are

\[
 (A+tw).r_j=A.r_j+t(w.r_j)>0
\]

for every \(t\geq0\), and

\[
 (A+tw)^2=A^2+2tA.w+t^2w^2\longrightarrow-\infty.
\]

Thus the first positive light-cone intersection is strict for every wall.
Finiteness of the wall list promotes strictness to an open ideal-sphere
patch.  A radial collar over such a patch has infinite hyperbolic volume
because the ball-model density is proportional to
\((1-|x|^2)^{-n}\).

No assertion that \(w\) is nef, future-directed, or itself a projective
point of \(P\) is used.  In particular the argument remains valid when
\(w=F=-K\) is spacelike and therefore lies on neither positive sheet.

## 2. The canonical sign gives the sharp cutoff

For every ordinary physical facet,

\[
 F.C=-K.C\in\{1,0\}.
\]

Therefore \(F\in\mathcal W\) solely from the raw oriented colors, while

\[
                         F^2=9-n.
\]

At \(n\geq10\), the wall cone contains the spacelike vector \(F\), so finite
volume is impossible.  Reflection integrality, maximality, Riemann--Roch,
effectivity, and a carrier are not being smuggled into this step.

The cutoff is sharp at the level of this argument:

- \(n\leq8\): \(F\) is timelike;
- \(n=9\): \(F\) is null;
- \(n\geq10\): \(F\) is spacelike.

## 3. Null orientation and cusp extremality

The rank-ten composition needs \(F\) to be a cusp on the selected future
sheet, not merely a null vector satisfying the homogeneous inequalities.
This is proved rather than assumed.

If a nonzero causal \(w\in\mathcal W\) were past-directed, then the
Lorentzian plane spanned by \(A\) and \(w\) would contain a spacelike
positive combination, unless \(w\) were proportional to \(-A\).  The
proportional case contradicts \(A.r_j>0\) and \(w.r_j\geq0\).  Hence every
causal vector of \(\mathcal W\), and in particular the null \(F\), is
future-directed.

If a future-null \(F\) were not extreme, write \(F=u+v\) with
nonproportional \(u,v\in\mathcal W\).  The spacelike exclusion makes
\(u,v\) causal, and the preceding orientation argument makes them future.
Reverse Lorentz Cauchy--Schwarz forces two future causal summands of a null
vector to be proportional null vectors.  This contradiction makes \(F\)
an extreme future-null ray, hence an ideal vertex.  The rank-ten cusp
theorem is therefore invoked with all its premises proved.

## 4. The timelike orientation was repaired

The earlier timelike report originally read the signs \(F.C\geq0\) as
\(F\in P\) before checking the chosen positive component.  That implication
was too quick.  The repaired argument runs in the correct order:

1. relative maximality and the all-reflections descent first give an
   oriented reflection isometry \(w:P\to P_n\);
2. the transformed characteristic vector \(G=wF\) satisfies every
   homogeneous standard inequality and \(G^2=9-n\);
3. those inequalities themselves give
   \(a\geq m_1+m_2+m_3\) and
   \(m_1\geq\cdots\geq m_n\geq0\);
4. the exact reduced-characteristic gap
   \[
   (9-n)k(k+1)+4k(A+B)+2AB+2A+2B
   \]
   forces \(G=F\);
5. only then is the standard finite vertex transported back to the chosen
   sheet.

Thus the conjugating reflection fixes \(K=-F\), and the carrier marking is
canonical.  No separate positive-sheet hypothesis is required in the
timelike range.

## 5. Positive carrier composition

For \(3\leq n\leq8\), the repaired timelike theorem constructs the first
\(n\) blowups of the flex-cubic base cluster.  The residual base ideal is
\((x^{9-n},w)\), so the anticanonical pencil has no fixed curve and
independently certifies \(-K\) nef.  Adjunction plus reflection saturation
then returns the complete nef cone.

For \(n=9\), the now-proved future characteristic cusp enters the rank-ten
cusp-rigidity theorem.  The ninefold flex-pencil resolution supplies the
unique section, the affine nine-component fibre, and exact full-nef return.

These are explicit constructions, not an appeal to weak-del-Pezzo or
rational-elliptic-surface classification.

## 6. Scope boundary

The complete proved implication is

\[
\begin{gathered}
L=I_{1,n},\ K=(-3,1,\ldots,1),\ n\geq3,\\
\text{finite-volume integral Coxeter chamber with ordinary colors},\\
\text{relative reflection maximality}
\end{gathered}
\Longrightarrow
\begin{gathered}
n\leq9,\\
\text{automatic canonical escape }K\notin P,\\
\text{explicit marked rational carrier},\\
P=\mathbf P(\operatorname {Nef}S).
\end{gathered}
\]

For a twice-primitive physical boundary \(B\),

\[
                         F.B=-2.
\]

Consequently \(F\notin\mathcal W\), so neither the spacelike cutoff nor the
anticanonical carrier composition applies.  The exponent-two adapted-host
\((-4)\) lane remains open at the independent weak-log
nef/semiampleness/effectivity gate.
