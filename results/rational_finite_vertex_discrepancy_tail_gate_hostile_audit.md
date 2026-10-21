# Hostile audit: finite-vertex discrepancy barrier and exceptional tails

**Audit date:** 2026-08-10  
**Verdict:** **PASS after six scope/rigor repairs.** The M-matrix argument, the
three displayed discrepancy solutions, and the arbitrary-tail formulas are
correct. The result is a local, carrier-relative obstruction theorem. It is
not an absolute Picard-rank theorem and does not say that the complete nef
cone in the tail examples is Coxeter, finite-volume, or reflection-maximal.

The audited owner artifacts are:

| artifact | terminal SHA-256 |
|---|---|
| `results/rational_finite_vertex_discrepancy_tail_gate.md` | `817d30726dad66f7512c9eaa115823f34c2f38d7e86224876b4d50613b9a4548` |
| `scripts/certify_rational_finite_vertex_discrepancy_tail_gate.py` | `e3660d360ea0a98f706a7e543a2326e2851a67998e80baaebdd94af0712b1480` |
| `data/enumerations/rational_finite_vertex_discrepancy_tail_gate.json` | `2567341c8e717d93fb969aa8c6c47db615a60dd87cebb69333bf828f7ad37c81` |
| `tests/test_rational_finite_vertex_discrepancy_tail_gate.py` | `c4c879736fca9c1845b2a00d4ef42f3a9ae645fb91e091d362e24213eebad7d6` |
| `tests/test_rational_finite_vertex_discrepancy_tail_gate_hostile.py` | `9c9cd87d5ca4c5b5be6c109224797d5e7b942b011a039a6c179a84f9fab5398c` |

No Nikulin--Vinberg, ABSW, Prasad, finite-diagram classification, or
arithmetic-volume input enters the proof or this audit.

## 1. The M-matrix sign is correct

Let $Q=(C_i.C_j)$. Negative definiteness of $Q$ makes $A=-Q$
positive definite, while nonnegative intersection of distinct effective
curves makes every off-diagonal entry of $A$ nonpositive. Thus $A$ is a
symmetric nonsingular M-matrix.

The inverse sign can be proved without importing an M-matrix theorem. If
$Ax=y\geq0$, write $x=w-z$ with $w,z\geq0$, disjoint supports, and
$z$ the negative part. Were $z\ne0$, then

\[
 0\leq z^Ty=z^TAw-z^TAz<0.
\]

Indeed, $z^TAw\leq0$ by the off-diagonal sign and disjoint supports, while
$z^TAz>0$ by positive definiteness. Hence $x\geq0$ whenever $y\geq0$,
which is exactly $A^{-1}\geq0$ entrywise.

The discrepancy convention is also consistent:

\[
 Qa=k \quad\Longleftrightarrow\quad Aa=-k.
\]

For $u_B=1/2$ and $u_R=u_E=0$, the weak-incidence rows are

\[
 (Au-k)_B=4/2-2=0,
 \qquad
 (Au-k)_R=0,
 \qquad
 (Au-k)_E=1-\frac12\sum_iB_i.E\geq0.
\]

Consequently

\[
 a+u=A^{-1}(Au-k)\geq0,
\]

so $a_B\geq-1/2$ and $a_R,a_E\geq0$. The sign would have failed had the
proof used $A(a-u)$, $Au+k$, or a nonpositive inverse; none of those
errors occurs in the repaired artifact.

## 2. Exact contraction and klt scope

The linear system alone defines numerical coefficients. Calling them
discrepancies and concluding klt requires more. The final theorem now assumes
a projective birational morphism

\[
 f:S\longrightarrow Y
\]

to a normal complex projective $\mathbb Q$-Gorenstein surface, with
exceptional locus exactly the displayed SNC packet and with $f$ an
isomorphism away from it. Under those hypotheses $S$ is a log resolution
of the contracted points and

\[
 K_S=f^*K_Y+\sum_i a_iC_i
\]

has its usual discrepancy meaning. Without the $\mathbb Q$-Gorenstein and
full-exceptional-locus hypotheses, the original global klt wording was not
justified.

The displayed curves have log discrepancies $A_i^{\log}=1+a_i\geq1/2$.
This really controls all higher divisorial valuations, but SNC is what makes
the assertion true so directly. Blow up a point on $s=0,1,2$ exceptional
components. The new divisor has

\[
 a_F=1+\sum_{p\in C_i}a_i,
 \qquad
 A_F^{\log}=2+\sum_{p\in C_i}a_i.
\]

The three possibilities for $A_F^{\log}$ are $2$,
$A_i^{\log}+1$, and $A_i^{\log}+A_j^{\log}$. Starting from log
discrepancies at least $1/2$, none lowers the minimum. Point blowups retain
SNC, and every divisorial valuation over a smooth surface appears after a
finite sequence of point blowups. Therefore the all-valuations claim is
valid.

With the standard strict convention, “$\epsilon$-klt” means every log
discrepancy is strictly greater than $\epsilon$. The class is uniformly
$\epsilon$-klt for each fixed $0<\epsilon<1/2$. It is not uniformly
$1/2$-klt: contraction of the negative section of $\mathbb F_4$ has
discrepancy $-1/2$, hence log discrepancy exactly $1/2$.

## 3. Hostile packets and the SNC poison test

Exact row reduction gives:

\[
\begin{array}{c|c}
\text{packet}&a\\ \hline
B_1.B_2=2&(-1,-1)\\
B.R=2&(-1,-1)\\
E\text{ joined to three disjoint }B_i&(-2,-1,-1,-1).
\end{array}
\]

All three displayed matrices are negative definite. In the stated SNC and
$\mathbb Q$-Gorenstein contraction scope, H1 and H2 are lc but not klt,
while H3 is not lc.

There is an important geometric poison test. Intersection number two for H1
or H2 can also arise from an ordinary tangency, but then the displayed
surface is not a log resolution. Starting with the numerical coefficients
$(-1,-1)$, the first blowup creates a divisor of discrepancy
$1-1-1=-1$ and leaves a triple point. Blowing up the triple point creates
a divisor of discrepancy

\[
 1-1-1-1=-2.
\]

That realization is not lc. Thus “H1/H2 are lc” is correct only with the
SNC qualification now printed in the theorem; the intersection matrix by
itself is insufficient.

## 4. Arbitrary tails are geometric finite nef vertices

Let $X_m\to\mathbb P^2$ be the sequence of blowups in the report. In the
ordered basis of strict exceptional components, the exceptional matrix is
the chain with diagonal

\[
 -2,\ldots,-2,-1
\]

and adjacent entries one. For every proper leading block,

\[
 D_j=(-1)^j(j+1),
\]

while the full determinant is $D_m=(-1)^m$. Sylvester's criterion gives
negative definiteness, and the full determinant gives unimodularity. The
recurrence $Qa=k$, with $k=(0,\ldots,0,-1)$, has the symbolic solution

\[
 a_j=j\qquad(1\leq j\leq m).
\]

The canonical correction is not merely a declared integer:

\[
 -a^TQa=-a^Tk=m.
\]

The certifier now computes the left side from the exact matrix and exact
discrepancy vector before comparing it with $m$.

The geometric vertex claim also checks. For
$L=\pi^*\mathcal O_{\mathbb P^2}(1)$, one has $L^2=1$, and the only
irreducible curves orthogonal to $L$ are the exceptional components. They
are independent, $\rho(X_m)=m+1$, and their hyperplanes therefore meet in
the single positive ray $\mathbb R_{>0}L$. Each negative component gives a
genuine nef wall: for an ample $H$ and such a component $C$,

\[
 H+\frac{H.C}{-C^2}C
\]

is nef, vanishes on $C$, and is positive on every other irreducible curve.
It is also big, because its square is
$H^2+(H.C)^2/(-C^2)>0$. After denominators are cleared, its restriction to
$C\simeq\mathbb P^1$ is the trivial degree-zero line bundle. The standard
nef-big surface contraction criterion therefore makes it semiample and
contracts exactly $C$. Thus every component supplies a genuine facet, and
$[L]$ is an actual finite vertex of the projectivized nef cone in arbitrarily
large Picard rank.

In the doubled lattice, every reflection is integral. For a square-minus-two
component the reflection coefficient is $-x.R$; for the square-minus-one
end it is $-2x.E$. The Cartan products are one on an $R$--$R$ edge and
two on the terminal $R$--$E$ edge, giving orders three and four. This is
a finite crystallographic local packet. It does **not** prove that all other
nef walls form a Coxeter system or that the resulting global chamber is
finite-volume or reflection-maximal.

## 5. Repairs made during the audit

Six defects were real and were repaired:

1. The contraction was not explicitly normal, $\mathbb Q$-Gorenstein,
   projective, and exhaustive of the exceptional locus.
2. The assertion about all further blowups lacked the SNC blowup formula and
   the strict-$\epsilon$ qualification.
3. The lc labels for H1/H2 did not warn that a tangency with the same total
   intersection number is non-SNC and can be non-lc.
4. The arbitrary chain was presented principally as a matrix type; it is now
   realized as an actual finite nef vertex on iterated blowups of
   $\mathbb P^2$, with local reflection integrality proved.
5. The executable set `discrepancy_correction = length` instead of computing
   $-a^TQa$.
6. The final section claimed an H1/H2 log-MMP core bound for which this
   artifact supplied no proof and no supporting result. That claim was
   removed. The hostile lanes remain open here.

No change was needed to the central M-matrix inequality, the three exact
hostile vectors, the determinant recurrence, or the discrepancy formula
$(1,2,\ldots,m)$.

## 6. Exact residual boundary

What is proved is:

- every specified weak-incidence SNC packet with a normal
  $\mathbb Q$-Gorenstein contraction has log discrepancies at least $1/2$;
- local finite-nef-vertex singularity and reflection data alone cannot bound
  rank, because the $R\!\cdots R-E$ tail occurs in every length.

What is not proved is:

- that an arbitrary negative-definite packet admits the required algebraic
  $\mathbb Q$-Gorenstein contraction;
- a bound on the exceptional tail in a complete literal full-nef carrier;
- finite volume, reflection maximality, or compiler surjectivity for the
  blowup examples;
- any global reduction or rank bound in H1, H2, or H3.

Thus the exact surviving H-free gate is global: finite volume and reflection
maximality must either bound the tail or force a finite numerical orbit of a
geometric contraction. The local discrepancy theorem cannot supply that
bridge.

## 7. Reproducibility

The independent hostile suite checks the M-matrix sign convention, exact
barrier identity, negative definiteness and discrepancies of H1--H3, the SNC
higher-blowup formula, the non-SNC tangency poison, arbitrary-tail determinant
and quadratic-correction formulas, doubled-lattice reflection integrality,
and the final scope language.

Terminal replay at the pinned bytes:

```text
$ python3 scripts/certify_rational_finite_vertex_discrepancy_tail_gate.py --check-archive
PASS: weak finite-vertex discrepancy barrier and arbitrary tail

$ python3 -m pytest -q tests/test_rational_finite_vertex_discrepancy_tail_gate.py tests/test_rational_finite_vertex_discrepancy_tail_gate_hostile.py
..........                                                               [100%]
10 passed
```
