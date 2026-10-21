# Hostile audit: weak-log reflection maximality and exact nef return

**Audit status (2026-08-10): PASS, conditionally on the five geometric
premises stated in the theorem.  Universal carrier existence remains OPEN.**

This audit checks only
results/reflection_maximal_weak_log_exact_nef_return.md.  It does not turn a
raw Coxeter packet into a surface, and it does not infer nefness from the
proposed chamber.

## 1. Adversarial checklist

### 1.1 Arbitrary irreducible negative curves

Let \(C\) be an irreducible negative curve which is not a component of the
disjoint boundary \(B\).  Adjunction and
\(D=-(K_S+B/2)\) give

\[
 C^2=2p_a(C)-2+D.C+\frac12B.C.                         \tag{A.1}
\]

Both correction terms are nonnegative: the first by actual nefness of
\(D\), the second because distinct irreducible effective curves have
nonnegative intersection.  If \(p_a(C)\ge1\), (A.1) contradicts
\(C^2<0\).  Hence \(p_a(C)=0\).  For an integral projective curve,

\[
 p_a(C)=g(\widetilde C)+\sum_p\delta_p,
\]

so \(C\) is smooth rational, not merely geometrically rational.  Since its
square is an integer, (A.1) leaves only:

- \(C^2=-2\), with \(D.C=B.C=0\); or
- \(C^2=-1\), with \(D.C+B.C/2=1\).

Thus Lemma 2.1 covers singular curves, positive-genus curves, and curves
meeting several boundary components.  No Mori-cone generation premise is
hidden in this step.

### 1.2 Integrality on the adapted host

Put

\[
 M=\langle L,B_1/2,\ldots,B_b/2\rangle,\qquad
 (x,y)_M=2(x,y)_L.
\]

The form is integral because
\((x,B_i/2)_M=x.B_i\), each half-boundary has square \(-2\), and distinct
half-boundaries are orthogonal.  The three reflection formulas are:

\[
\begin{array}{ll}
C^2=-2:&s_C(x)=x+(x.C)C,\quad s_C(B_i/2)=B_i/2;\\[2mm]
C^2=-1:&s_C(x)=x+2(x.C)C,\quad
        s_C(B_i/2)=B_i/2+(B_i.C)C;\\[2mm]
C=B_i:&s_{B_i/2}(x)=x+(B_i.x)B_i/2.
\end{array}                                             \tag{A.2}
\]

Every coefficient in (A.2) is integral in the displayed generators of
\(M\).  The first line uses the already proved equality \(B_i.C=0\).
Consequently every possible negative curve supplies an integral
hyperbolic reflection of \(M\).

### 1.3 The exact maximality mode

The proof is valid under any one of:

1. no proper discrete reflection overgroup of \(W(P)\) inside \(O^+(M)\);
2. absolute maximality among discrete hyperbolic reflection groups; or
3. maximality among arithmetic reflection overgroups in the ambient real
   hyperbolic group, allowing commensurable integral hosts.

Indeed, \(O^+(M)\) is discrete.  Adjoining any integral reflection remains
a discrete reflection group.  Since \(W(P)\) has finite covolume, every
such discrete overgroup has finite index.

The following weaker statements do not suffice:

- maximality only in \(O^+(L)\) when \(M\ne L\);
- saturation of the root span or of a displayed diagram;
- absence of one more facet root in a bounded search; or
- maximality among chambers in a fixed enumeration.

These distinctions are now explicit in the source theorem.

### 1.4 Positive-cone orientation and ideal closure

The positive component is the one containing an ample class \(A\).
Because all displayed facet curves are effective, \(A.C_j>0\), so \(A\)
lies strictly inside the chosen chamber.  A reflection mirror of the
Coxeter group cannot cross the interior of its fundamental chamber; this
fixes the nonnegative side for every hidden negative curve.

For an irreducible curve with \(C^2\ge0\), the inequalities
\(A.C>0\) and \(C^2\ge0\) put \(C\) in the closed future Lorentz cone.
The closed cone over \(P\) is in the same component, and reverse
Cauchy--Schwarz gives \(x.C\ge0\).  This includes ideal rays by closure.
Conversely, a nonzero nef class lies in the same closed component and is
nonnegative on every displayed effective curve.  Hence both containments
hold for the complete ideal closures.

### 1.5 The one genuinely independent geometric gate

The inequalities \(D.C_j\ge0\) on the proposed facets prove only that
\(D\) lies in the proposed dual cone.  They cannot prove that \(D\) is nef,
because exact return of that cone is the conclusion.  A valid finite
certificate must instead provide at least one of:

1. an integer \(m>0\) such that \(mD\) is Cartier and \(|mD|\) is
   base-point-free;
2. a morphism \(g:S\to T\) and a nef divisor \(H\) with \(mD=g^*H\); or
3. an independently complete curve-cone calculation.

Without this gate, (A.1) has no sign control and the reflection argument
cannot exclude hidden curves.

## 2. Exact non-Coble calibration

Let \(S\) be the blow-up of a torus fixed point on the positive section of
\(\mathbb F_4\), away from the negative section \(B\).  In the basis
\((B,F,E)\),

\[
 \begin{pmatrix}-4&1&0\\1&0&0\\0&0&-1\end{pmatrix},
 \qquad K=(-2,-6,1),\qquad G=F-E.
\]

The complete toric Mori cone is

\[
 \overline{NE}(S)=\operatorname{Cone}(B,E,G).
\]

For

\[
 D=-\left(K+\frac B2\right)=\frac32B+6F-E
\]

one obtains

\[
 D^2=8,\qquad (D.B,D.E,D.G)=(0,1,1/2).
\]

Thus nefness is certified independently of the proposed reflection
chamber.  In the adapted basis \((B/2,F,E)\), the doubled form and facet
root Gram are

\[
 G_M=\begin{pmatrix}-2&1&0\\1&0&0\\0&0&-2\end{pmatrix},
\qquad
 G_{\rm root}=
\begin{pmatrix}-2&0&1\\0&-2&2\\1&2&-2\end{pmatrix}.
\]

The Coxeter triangle has labels \((2,3,\infty)\) and area \(\pi/6\).
Every noncompact hyperbolic Coxeter polygon has area at least \(\pi/6\);
any proper discrete reflection overgroup would have strictly smaller area.
It is therefore reflection-maximal.  This gives an exact-return example
with \(D\ne0\) and \(B\ne-2K\).

The companion exact certifier is
scripts/certify_reflection_maximal_weak_log_exact_nef_return.py.

## 3. Verdict

The implication

\[
\text{constructed weak-log surface + independently nef \(D\)
      + adapted-host reflection maximality}
\Longrightarrow
P=\mathbf P(\operatorname{Nef}S)
\]

is sound.  It is a closure lemma, not an existence theorem.  In particular,
the following remain OPEN:

- a universal raw-polytope-to-surface construction;
- a universal finite semiampleness certificate for \(D\);
- exhaustiveness of any compiler or atlas; and
- arithmetic dimension bounds without the separate boundary-contraction
  and orbifold hypotheses.
