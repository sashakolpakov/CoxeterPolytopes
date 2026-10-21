# A finite arithmetic factory for golden K3 eigen-nef chambers

**Status.**  The theorem below is an exact existence and realization theorem.
It upgrades the Bugaenko calculation from one lattice to every primitively
K3-embeddable even golden lattice with the same trace normalization.  It does
not assert that the resulting chamber is finite-sided, finite-volume, or
compact.  Those three properties remain separate reflection-arrangement
calculations.

## 1. The finite input

Put

\[
 F=\mathbf Q(\sqrt5),\qquad
 \mathcal O=\mathbf Z[\phi],\qquad
 \phi^2=\phi+1,
 \]

and let

\[
 \alpha=-\frac{\phi}{\sqrt5},\qquad u=\phi-1.
\]

Both real embeddings of \(\alpha\) are negative.  A **golden K3 carrier
datum of dimension \(d\)** is the following finite arithmetic information.

1. A free \(\mathcal O\)-module \(M\) of rank \(d+1\), with an even symmetric
   \(\mathcal O\)-valued form \(h\).  Even means
   \(h(x,x)\in2\mathcal O\) for every \(x\in M\).
2. At the identity embedding, \(h\) has signature \((d,1)\); at the conjugate
   embedding it is positive definite.
3. The trace lattice
   \[
   N=(M,q),\qquad q(x,y)=\operatorname{Tr}_{F/\mathbf Q}
       \bigl(\alpha h(x,y)\bigr),
   \]
   together with a primitive embedding \(j:N\hookrightarrow\Lambda_{K3}\).
4. A vector \(p\in M\) which is time-like at the identity place and regular:
   \[
   h^{(0)}(p,p)<0,qquad h(p,r)\ne0
   \quad\text{for every }r\in M\text{ with }h(r,r)=2.
   \]

The codifferent condition on \(\alpha\) and evenness of \(h\) make \(N\) an
even integral lattice.  Its signature is \((1,2d+1)\).  Multiplication by
\(\phi\) gives a self-adjoint integral endomorphism

\[
 \tau^2=\tau+1,qquad q(\tau x,y)=q(x,\tau y).
\]

The regularity test is itself finite.  The orthogonal complement of a
time-like \(p\) is positive definite at the identity place and is positive
definite at the conjugate place.  Hence it contains only finitely many
vectors of field norm two, which can be enumerated exactly.

Let \(E\subset N_{\mathbf R}\) be the \(\phi\)-eigenspace belonging to the
identity embedding.  Let \(\mathcal H_2\) be the mirrors in \(E\) of all
vectors \(r\in M\) with \(h(r,r)=2\), and let \(C_p\) be the component
containing \([p]\) of

\[
 \mathbb H(E)\setminus\bigcup_{H\in\mathcal H_2}H.
\]

The vector \(p\) is finite data selecting a chamber; an infinite list of
signs is not part of the input.

## 2. Universal shell theorem

### Theorem 2.1 — shell collapse

For every even golden K3 carrier datum, the restriction to \(E\) of the
complete \((-2)\)-root arrangement of \(N\) is exactly \(\mathcal H_2\).
Moreover every mirror in \(\mathcal H_2\) has exactly four signed
trace-root representatives

\[
 \pm r,\qquad \pm ur,
\]

where \(h(r,r)=2\).  The two unsigned roots \(r,ur\) are orthogonal for
the trace form.

#### Proof

Write \(t=h(\delta,\delta)=a+b\phi\).  If \(q(\delta,\delta)=-2\) and the
hyperplane of \(\delta\) meets the positive cone of \(E\), then \(t\) is
positive at both real embeddings.  The trace equation is

\[
 \operatorname{Tr}(\alpha t)=-(a+b)=-2.
\]

Thus \(a=2-b\).  Total positivity gives \(-3\le b\le1\), so the complete
finite shell list is

\[
 5-3\phi,\quad4-2\phi,\quad3-\phi,\quad2,\quad1+\phi.
\]

Evenness requires \(t\in2\mathcal O\).  Only

\[
 4-2\phi\quad\text{and}\quad2
\]

survive.  Since \(u^2=2-\phi\),

\[
 2u^2=4-2\phi,
\]

and multiplication by the unit \(u\), or by \(\phi=u^{-1}\), identifies
the two shells without changing the restricted hyperplane.  Conversely, a
norm-two vector has trace square \(-2\), and its field reflection is integral:

\[
 s_r(x)=x-h(x,r)r.
\]

This proves equality of the arrangements.

If two integral vectors define the same \(F\)-rational hyperplane in \(E\),
their normals are proportional by an element of \(F\).  Comparing their
norms shows that the proportionality factor is one of
\(\pm1,\pm u\).  Finally

\[
 q(r,ur)=\operatorname{Tr}(2\alpha u)=0.
\]

This proves the last two assertions. \(\square\)

The important point is that this proof uses only five scalar shells.  It is
not a Vinberg enumeration and it is independent of the particular golden
lattice.  The Bugaenko lattice is one calibration of this general fact.

## 3. K3 realization theorem

### Theorem 3.1 — Weyl-framed eigen-nef factory

Let \((M,h,j,p)\) be a golden K3 carrier datum.  Then there is a projective
marked K3 surface \(X\), and a self-adjoint integral endomorphism

\[
 \tau_C\in\operatorname{End}_{\mathbf Z}\operatorname{NS}(X),
 \qquad \tau_C^2=\tau_C+1,
\]

such that

\[
 \boxed{
 \mathbb P\bigl(\operatorname{Nef}(X)\cap E_C^+\bigr)
 \;=\; \overline{C_p},
 }
\]

up to the marking isometry.  Here \(E_C\) is the distinguished eigenspace of
\(\tau_C\).  Consequently the output is automatically a chamber of a
discrete reflection arrangement, hence a Coxeter chamber, possibly with
infinitely many facets.

#### Proof and the division of labour

The primitive embedding already included in the carrier datum and the K3
period theorem produce a projective K3 surface with
\(\operatorname{NS}(X)\simeq N\); there is no additional period-existence
hypothesis.  Choose \(p\) in the interior
of \(C_p\).  Theorem 2.1 says that \(p\) misses every K3 root hyperplane.
It therefore lies in a unique chamber \(D\) of the full \((-2)\)-arrangement
in the positive cone of \(N\), and

\[
 D\cap E=C_p.
\]

The K3 Weyl group is transitive on the root chambers.  If \(A_X\) is the
ample chamber, choose \(w\in W(N)\) with \(wD=A_X\), and put

\[
 \tau_C=w\tau w^{-1},\qquad E_C=wE.
\]

Then \(A_X\cap E_C=wC_p\), and the same equality holds for the closed nef
cones.  Every K3 root reflection is trivial on the discriminant group of
\(N\), so \(w\) extends by the identity on \(N^\perp\) to the K3 lattice and
does not alter the period.  This is a change of Weyl frame, not an assertion
that \(w\) is induced by an automorphism of \(X\). \(\square\)

Rationally, \(\tau_C\) is induced on \(\operatorname{NS}(X)\) by an algebraic
correspondence in \(CH^2(X\times X)_{\mathbf Q}\), supported on products of
divisors.  Thus the output can be recorded as an enhanced algebraic object
\((X,Z_{\tau_C})\), rather than as an abstract real slice.  No natural real
multiplication on the transcendental Hodge structure is claimed.

### What is genuinely automatic

The finite carrier datum, period theorem, Riemann--Roch, and the K3 Weyl
theorem supply all of the following without designing an equation for the
curves:

- a projective surface with the prescribed full Neron--Severi lattice;
- the self-adjoint golden correspondence on its divisor lattice;
- the complete restricted divisor-root arrangement;
- a nonempty eigen-nef chamber selected by \(p\); and
- irreducible rational curves supporting every actual chamber facet.

They do **not** supply:

- termination of the norm-two Vinberg algorithm;
- a finite list of simple walls;
- finite covolume;
- compactness; or
- an explicit projective equation or explicit Weyl word.

Those are not period-theoretic questions.  They are exact arithmetic
properties of \((M,h,C_p)\).

## 4. The complete two-curve lift

Every restricted facet has two ambient K3 walls.  If \(r_i\) is its
norm-two representative, the two trace roots are

\[
 r_i,\qquad ur_i=(\tau-1)r_i.
\]

At a generic point of that facet no other restricted mirror passes through
the point.  The local ambient root system is therefore \(A_1\oplus A_1\), so
after the Weyl framing these two roots are classes of two disjoint smooth
rational curves.

Put

\[
 Q_{ij}=q(r_i,r_j),\qquad B_{ij}=q(\tau r_i,r_j).
\]

The full intersection matrix on
\((r_1,\ldots,r_n,ur_1,\ldots,ur_n)\) is

\[
 \boxed{
 \widetilde Q=
 \begin{pmatrix}
 Q&B-Q\\
 B-Q&2Q-B
 \end{pmatrix}. }
\]

This is the curve matrix, whereas
\(\left(\begin{smallmatrix}Q&B\\B&Q+B\end{smallmatrix}\right)\) is the
root-presentation matrix on \((r_i,\tau r_i)\); the vectors \(\tau r_i\) are
generally not roots.

Distinct irreducible curves have nonnegative intersection.  Hence every
completed eigen-nef facet system necessarily satisfies

\[
 Q_{ij}\ge0,\qquad B_{ij}-Q_{ij}\ge0,
 \qquad2Q_{ij}-B_{ij}\ge0\quad(i\ne j),
\]

or equivalently

\[
 \boxed{Q_{ij}\le B_{ij}\le2Q_{ij}.}
\]

If \(h(r_i,r_j)=a+b\phi\), this is simply \(a\le0\) and \(b\le0\).
It also recovers the lost distinction between labels 3 and 5:

| restricted edge | \((Q_{ij},B_{ij})\) | complete curve lift |
|---|---:|---|
| \(m=2\) | \((0,0)\) | \(A_1^4\) |
| \(m=3\) | \((1,1)\) | \(A_2\oplus A_2\) |
| \(m=5\) | \((1,2)\) | \(A_4\) |

The earlier statement that a 3-edge and a 5-edge both look like two curves
meeting once is correct only after choosing one root above each wall.  The
complete paired lift distinguishes them algebro-geometrically.

The exact construction is implemented as
`paired_trace_curve_gram` in
`src/coxeter_polytopes/geometric/k3_design.py` and is calibrated on all eleven
Bugaenko walls.

## 5. A finite arithmetic gluing object in dimension nine

For \(d=9\), the trace lattice has rank 20.  Fix \((M,h,p)\), hence \((N,\tau)\).
Define \(\mathfrak K(N,\tau,p)\) to be the set of triples

\[
 (T,\gamma,o),
\]

where

- \(T\) is an even positive-definite rank-two lattice;
- \(o\) is an orientation of \(T\); and
- \(\gamma:(A_N,q_N)\xrightarrow{\sim}(A_T,-q_T)\) is an anti-isometry such
  that gluing \(N\oplus T\) along the graph of \(\gamma\) is unimodular.

Hold the enhanced datum $(N,\tau,p)$ pointwise fixed and impose

\[
 (T,\gamma,o)\sim(T',\gamma',o')
\]

when there is an orientation-preserving isometry $u:T\to T'$ satisfying
$\bar u\circ\gamma=\gamma'$.  One may subsequently quotient by automorphisms
of $N$ preserving $\tau$ and the selected chamber if an unmarked object is
wanted.  With the pointwise convention above,
$\mathfrak K(N,\tau,p)$ is a **finite, zero-dimensional arithmetic gluing
parameter set**.  Finiteness follows because $|\det T|=|\det N|$, there are
only finitely many reduced positive binary even forms of a fixed determinant,
and the discriminant groups are finite.

Every point of \(\mathfrak K(N,\tau,p)\) gives an oriented primitive embedding

\[
 N\hookrightarrow \Lambda_{K3},\qquad N^\perp\simeq T.
\]

The Shioda--Inose theorem then gives the corresponding singular K3 surface,
and Theorem 3.1 equips it with a Weyl-framed correspondence whose closed
eigen-nef section is exactly $\overline{C_p}$.  The anti-isometry $\gamma$
retains embedding/marking information, so distinct parameters may give the
same unmarked underlying surface.  Thus this is not silently called a coarse
surface moduli space.  Every parameter nevertheless produces the stated
Coxeter chamber.

### Exact nonemptiness criterion

The object \(\mathfrak K(N,\tau,p)\) is nonempty if and only if the following
finite computation succeeds.

1. Check that \(N\) is even of signature \((1,19)\).
2. Enumerate reduced matrices
   \[
   T=\begin{pmatrix}2a&b\\b&2c\end{pmatrix}>0,
   \qquad \det T=|\det N|.
   \]
3. Compute their finite discriminant quadratic forms.
4. Find an anti-isometry \(q_N\simeq-q_T\).
5. Verify that the graph gluing is even and unimodular.

The resulting lattice has signature \((3,19)\), so uniqueness of the even
unimodular lattice of that signature identifies it with \(\Lambda_{K3}\).
Conversely every primitive rank-20 K3 embedding has exactly such a rank-two
complement and discriminant anti-isometry.  This proves necessity and
sufficiency, not merely a stable-range filter.

Useful immediate rejection tests are

\[
 \ell(A_N)\le2
 \quad\text{and}\quad
 |\det N|=\det T
\]

for some positive even binary form \(T\).

Because \(\rho=20\) is maximal, there is no generic-period caveat in this
dimension: the constructed singular K3 has
\(\operatorname{NS}(X)=N\) exactly.  If desired, Shioda--Inose turns the
oriented binary form \(T\) into an explicit Kummer/Inose model using CM
elliptic curves.

## 6. The finite-wall and compactness gates

For a proposed finite set \(R=\{r_1,\ldots,r_n\}\) of norm-two roots, the
surface-existence problem and the polytope-completeness problem should be
kept separate.

### Carrier gate

The two-Gram data \((Q,B)\), their saturated root presentation, and the
rank-two discriminant gluing decide whether a singular K3 carrier exists.
Once this gate passes, no further period or curve-existence condition is
needed.

### Chamber gate

One must still prove that \(R\) is the complete simple-root set of \(C_p\).
This can be done by an exact Vinberg enumeration or another exhaustive
root-closure certificate.  The shell theorem guarantees that this one
norm-two calculation simultaneously accounts for every divisor wall; there
is no second hidden algebro-geometric wall search.

### Volume and compactness gates

After completeness:

- finite volume is the elliptic/parabolic face criterion, or equivalently
  finite covolume of the norm-two reflection group;
- compactness is the absence of parabolic and hyperideal extreme rays; and
- on the surface it says that the eigen-nef section has no visible isotropic
  ray.  It does not say that the full K3 nef cone has no isotropic class.

This division is the useful primal--dual conclusion.  Algebraic geometry
constructs and rigidifies the carrier and proves that the arithmetic
arrangement is the complete divisor arrangement.  Reflection arithmetic
alone decides whether the selected chamber terminates and is compact.

## 7. Search consequence

The dimension-nine example track can now search a finite arithmetic gluing
parameter set, not a catalogue of surfaces:

\[
 (Q,B)
 \longrightarrow (N,\tau)
 \longrightarrow (T,\gamma,o)
 \longrightarrow X_T
 \longrightarrow C_p.
\]

In the reverse, design-first direction:

1. prescribe a golden Coxeter Gram matrix or its two integral shadows
   \((Q,B)\);
2. impose the full paired-curve matrix \(\widetilde Q\);
3. take the saturated rank-20 root-presentation quotient and enumerate its
   finite \(\tau\)-stable even overlattices;
4. apply the exact binary-complement gluing criterion;
5. only for surviving singular K3 carriers, run norm-two Vinberg completion
   and compactness.

The surface is therefore constructed *from* the combinatorial constraints.
The difficult unresolved boundary is no longer whether the divisor classes
can live on a K3 surface; it is whether their norm-two arithmetic chamber is
finite and compact.

## References

- I. I. Piatetski-Shapiro and I. R. Shafarevich, *A Torelli theorem for
  algebraic surfaces of type K3* (1971).
- I. V. Dolgachev, *Mirror symmetry for lattice polarized K3 surfaces*,
  [arXiv:alg-geom/9502005](https://arxiv.org/abs/alg-geom/9502005).
- T. Shioda and H. Inose, *On singular K3 surfaces* (1977).
- E. B. Vinberg, *Discrete linear groups generated by reflections* (1971).
