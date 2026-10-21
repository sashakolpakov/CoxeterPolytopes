# The rank-two golden OG10 route collapses, and the first reflective escape box is empty

## Result

Let

\[
 \mathcal O=\mathbf Z[\phi],\qquad \phi^2=\phi+1,
 \qquad
 F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle .
\]

There are two exact conclusions.

1. A genuine extension of \(F\) by two \(\mathcal O\)-lines cannot be an
   OG10 Picard lattice.  It has field rank \(12\), trace rank \(24\), and,
   with the required two real signatures, trace signature \((1,23)\).  The
   OG10 Beauville lattice has signature \((3,21)\), and a projective OG10
   manifold has Picard rank at most \(22\).
2. Even abstractly, a rank-two complement does not enlarge the projected
   norm-two wall alphabet under the natural fixed-slice involution.  If the
   involution fixes the saturated copy of \(F\) pointwise and has no further
   fixed direction, every projected norm-two root already spans a one-root
   carrier of the previously searched form
   
   \[
     \begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix}.
   \]

Thus the proposed rank-two norm-two construction is structurally redundant,
not merely unsuccessful in a coefficient search.

The first genuinely different *reflective* wall type available on OG10 is a
primitive Beauville class of square \(-6\) and ambient divisibility \(3\).
Its wall-first field compiler is

\[
 H_6(k)=
 \begin{pmatrix}F&3k\\3k^t&6\end{pmatrix}.
\]

It has an integral norm-six field reflection and can supply an OG10
stably-prime-exceptional wall if a compatible primitive embedding sends a
selected trace class to the unique ambient order-three discriminant class.
This embedding has **not** been constructed here.

An independent exact replay of the first three normalized covector boxes
finds no \(H_6(k)\) wall cutting the certified \(v_{61}\) escape:

| normalized family | signature-admissible wall tests | cuts of \(v_{61}\) |
|---|---:|---:|
| zero \(E_8\), coefficient bound \(30\) | 18 | 0 |
| scalar-unit \(E_8\), \(n\leq8, v^2\leq8\) | 696,000 | 0 |
| first mixed \(E_8\) shell | 5,760 | 0 |

The rank-two obstruction is global under its stated hypotheses.  The
norm-six conclusion is a finite-box no-go only.

## 1. The exact rank-two compiler

Write a two-line field extension as

\[
 H=
 \begin{pmatrix}
 F&L\\L^t&R
 \end{pmatrix},qquad
 L=(\ell_1,\ell_2),qquad
 P=F^{-1}L=(p_1,p_2).
\]

The Schur matrix is

\[
 S=R-L^tF^{-1}L
   =R-\bigl(h(p_i,p_j)\bigr)_{i,j=1}^2.
\tag{1}
\]

Since \(F\) has signatures \((9,1)\) and \((10,0)\), the desired signatures
\((11,1)\) and \((12,0)\) for \(H\) are equivalent to

\[
 S^{(0)}>0,qquad S^{(1)}>0.
\tag{2}
\]

For a \(2\times2\) matrix these are the exact principal-minor tests

\[
 S_{11}^{(j)}>0,qquad \det S^{(j)}>0,qquad j=0,1.
\tag{3}
\]

If both adjoined basis vectors have norm \(2\) and mutual product \(c\), then

\[
 S=
 \begin{pmatrix}
 2-p_1^2&c-(p_1,p_2)\\
 c-(p_1,p_2)&2-p_2^2
 \end{pmatrix}.
\tag{4}
\]

The natural central-inversion involution is

\[
 J|_F=1,qquad J(e_i)=2p_i-e_i.
\tag{5}
\]

It is integral precisely when every \(2p_i\) belongs to the fixed lattice
(for a basis extension), and then (1) immediately shows that it is an
isometric involution.  These are the requested exact two-place,
integrality, and involution gates.  The problem is that such a full
rank-two field carrier cannot pass the OG10 rank/signature gate.

## 2. The OG10 trace-rank obstruction

Put

\[
 \alpha=-\frac{\phi}{\sqrt5},\qquad
 B(x,y)=\operatorname{Tr}_{\mathbf Q(\sqrt5)/\mathbf Q}
          \bigl(\alpha H(x,y)\bigr).
\]

Both real values of \(\alpha\) are negative.  At the identity place the
signature \((11,1)\) therefore contributes \((1,11)\) to \(B\), while the
positive-definite conjugate place contributes \((0,12)\).  Hence

\[
 \operatorname{sig} B=(1,23),\qquad \operatorname{rk}B=24.
\tag{6}
\]

On the other hand

\[
 \Lambda_{\mathrm{OG10}}
 =3U\oplus2E_8(-1)\oplus A_2(-1)
\]

has signature \((3,21)\).  A real subspace of it cannot have negative index
\(23\).  Equivalently, the period of a projective irreducible symplectic
manifold consumes a positive real two-plane, so its Picard rank is at most
\(24-2=22\).  Thus no nondegenerate field-rank-12 trace lattice of the
required Hodge signatures embeds in \(\Lambda_{\mathrm{OG10}}\), primitively
or otherwise.

The certificate checks the exact model

\[
 p_1=\tfrac12u,qquad p_2=\tfrac\phi2u,qquad R=2I_2
\]

(\(u\) is the last basis vector of \(F\)).  Its two field signatures are
\((11,1)\), \((12,0)\), its even trace lattice has signature \((1,23)\), and
the integral involution (5) has fixed dimension \(10\) and anti-fixed
dimension \(2\).  This model is only an exact identity check; equation (6)
proves the general obstruction.

## 3. Central-inversion reduction to one root

The second obstruction is independent of the rank count.

Let \(M\) be an integral \(\mathcal O\)-lattice containing \(F\) as a
saturated sublattice.  Assume that at the identity embedding \(M\) is
Lorentzian with its unique negative direction already contained in \(F\),
and that at the conjugate embedding \(M\) is positive definite.  Suppose
\(J\in O(M)\) is an integral involution,

\[
 J|_F=1,qquad \operatorname{Fix}(J)_{K}=F_K,
 \qquad K=\mathbf Q(\sqrt5).
\tag{7}
\]

The \(+1\) and \(-1\) eigenspaces of an isometric involution are orthogonal.
Condition (7) therefore forces

\[
 J=-1\quad\text{on }F_K^\perp.
\tag{8}
\]

Let \(r\in M\) have \(r^2=2\), and put

\[
 p=\frac{r+Jr}{2},qquad u=\frac{r-Jr}{2}.
\tag{9}
\]

Then \(p\in F_K\), \(u\in F_K^\perp\), and \(p\) is exactly the orthogonal
projection of \(r\) to \(F_K\).  Integrality of \(J\) gives

\[
 2p=r+Jr\in M\cap F_K=F,
\tag{10}
\]

where the last equality is saturation.  If
\(\ell=(h(r,f_i))_i\), then \(\ell\in\mathcal O^{10}\) and
\(p=F^{-1}\ell\).  Moreover

\[
 \operatorname{Gram}(F+\mathcal Or)
 =\begin{pmatrix}F&\ell\\\ell^t&2\end{pmatrix},qquad
 r^2=p^2+u^2.
\tag{11}

The two ambient signature hypotheses make \(F_K^\perp\) positive definite
at both real places.  A mirror meets the old hyperbolic space exactly when
\(p^2>0\), so (11) gives

\[
 0<p^2<2,qquad (p')^2<2.
\tag{12}

These are exactly the old one-root Schur inequalities, and (10) is exactly
the old involution-integrality condition.  Thus every projected root of the
rank-two carrier belongs to the one-root feasibility universe.

There is no cooperative two-wall escape from this conclusion.  Orient every
new half-space to contain the archived interior point.  If their intersection
removes the \(v_{61}\) ray, at least one half-space removes it, hence at least
one covector satisfies

\[
 \ell(x_{\mathrm{int}})\ell(v_{61})<0.
\tag{13}

Cross-coupling in (4) cannot make two individually noncutting walls jointly
cut the ray.

If \(F\) is not saturated, (10) lands only in the saturation of \(F\).  That
changes the fixed divisor lattice and requires a new fixed-arrangement
calculation; it is an explicit escape from the theorem's hypotheses, not a
rank-two realization of the same carrier.

## 4. The first reflective escape: norm six and divisibility three

OG10 has a stably-prime-exceptional wall type

\[
 (D^2,\operatorname{div}_\Lambda D)=(-6,3).
\]

For a field vector \(r\) of norm \(6\), the field reflection is

\[
 s_r(x)=x-\frac{h(x,r)}3r.
\tag{14}

It is integral if

\[
 h(r,M)\subset3\mathcal O.
\tag{15}

This supplies a different reflective root length, rather than trying to get
new geometry from a second norm-two line.

The wall-first rank-11 compiler is

\[
 H_6(k)=
 \begin{pmatrix}F&3k\\3k^t&6\end{pmatrix},qquad
 p_0=F^{-1}k,qquad p=3p_0.
\tag{16}

The adjoined vector \(e\) is primitive of norm \(6\), condition (15) is
automatic, and its fixed-slice hyperplane is

\[
 p^\perp=p_0^\perp.
\tag{17}

Its Schur norm is

\[
 6-p^2=6-9p_0^2.
\]

Writing

\[
 k=(A,E_8w,C),qquad
 q=2p_0^2=-\frac{A^2}{\phi}+C^2+2w^2,
\tag{18}

the exact two-place and cutting gates become

\[
 0<q<\frac43,qquad q'<\frac43,qquad
 k(x_{\mathrm{int}})k(v_{61})<0.
\tag{19}

The factor \(3\) disappears from the hyperplane and from the sign test, but
it tightens the Schur budget from \(4\) to \(4/3\).

### The trace pair and its intrinsic reflection

Let \(u=\phi-1\).  For

\[
 D=r,qquad E=ur,
\]

the trace pairing gives

\[
 B(D,D)=B(E,E)=-6,qquad B(D,E)=0.
\tag{20}

If \(h(x,r)=a+b\phi\), then

\[
 B(x,D)=-(a+b),\qquad B(x,E)=-b.
\]

The two integral trace-lattice reflections therefore commute and their
product is

\[
 R_D R_E(x)
 =x-\frac{a+b\phi}{3}r
 =s_r(x).
\tag{21}

Thus (14) is an intrinsic integral Coxeter reflection of the golden lattice.
There is, however, an important ambient OG10 qualification.

## 5. The mod-\(p\) folding bound

The following elementary lemma is reusable.

> **Mod-\(p\) folding lemma.**  Let \(N\hookrightarrow\Lambda\) be a primitive
> embedding of integral lattices.  If \(D_1,\ldots,D_d\in N\) are independent
> in \(N/pN\) and each has ambient divisibility divisible by \(p\), then
> \[
> d\leq\dim_{\mathbf F_p}\operatorname{Rad}(\Lambda/p\Lambda).
> \]

Indeed, primitivity makes \(\Lambda/N\) torsion-free.  Tensoring

\[
 0\longrightarrow N\longrightarrow\Lambda
   \longrightarrow\Lambda/N\longrightarrow0
\]

with \(\mathbf F_p\) therefore remains left-exact: the Tor term vanishes and
\(N/pN\to\Lambda/p\Lambda\) is injective.  Ambient divisibility by \(p\)
puts every image \(D_i\bmod p\) in the radical of the ambient form.  Their
independence is preserved, proving the inequality.

For OG10, exact row reduction gives

\[
 \dim_{\mathbf F_3}\operatorname{Rad}
 (\Lambda_{\mathrm{OG10}}/3\Lambda_{\mathrm{OG10}})=1.
\tag{22}

The classes \(D=e\) and \(E=(\phi-1)e\) in (20) are independent modulo \(3\)
and both have divisibility \(3\) *inside the trace lattice*.  Hence they
cannot both have ambient OG10 divisibility \(3\) under a primitive embedding.
The intersection of their two-plane with the ambient radical may be zero or
may be any line; it need not be the line of \(D\) or of \(E\).

Consequently (21) is **not** automatically a product of two OG10 monodromy
wall reflections.  A viable embedding must separately arrange that one
selected class, say \(D/3\), maps to the nonzero class of
\(A_{\Lambda_{\mathrm{OG10}}}\cong\mathbf Z/3\).  Its discriminant value is

\[
 (D/3)^2=-\frac23\equiv\frac43\pmod{2\mathbf Z},
\tag{23}

as required.  Only after an explicit gluing proves this does \(D\) become an
actual OG10 \((-6,3)\) wall.  The other trace lift cannot simultaneously be
such a wall.  Existence of a compatible primitive embedding is an open gate
in this note.

This distinction is useful: the field reflection remains integral, and a
single selected trace lift can in principle provide the algebro-geometric
wall, but the especially clean two-monodromy-reflection lift available for
norm-two roots is impossible here.

## 6. Exact finite-box screen

The certificate independently replays the normalized covectors \(k\), not
merely the old norm-two counts, using the sharper exact test \(0<3q<4\) at
both embeddings.

The declared families are:

1. \(w=0\), with all four coefficients of \(A,C\) in \([-30,30]\);
2. \(w=\phi^n v\), \(1\leq n\leq8\), with
   \(v\in E_8\), \(v^2\in\{2,4,6,8\}\), and the same coefficient box;
3. \(w=a+\phi b\), \(b^2=2\), \(0<(w')^2<2\), again with coefficient
   bound \(30\).

The exact census is:

| family | signature parameter records | expanded wall tests | cutters |
|---|---:|---:|---:|
| zero \(E_8\) | 18 | 18 | 0 |
| scalar-unit \(E_8\) | 136 | 696,000 | 0 |
| first mixed shell (13,920 eligible \(w\)'s) | 24 | 5,760 | 0 |

All comparisons use exact arithmetic in \(\mathbf Z[\phi]\).  The \(E_8\)
shell counts are pinned by the theta series, and the mixed shell uses the
exact rational Schur recursion of the earlier box theorem.  Canonical hashes
of every accepted wall record are frozen in the JSON certificate.

This result could also be inferred as a subset of the old \(q<4\) no-cut
theorem, but the independent replay supplies the exact norm-six counts and
digests.

## 7. Remaining carrier and arrangement gates

A norm-six covector outside the excluded boxes is only the beginning.  A
genuine example must pass all of the following.

1. **Novelty and truncation:** \(p_0^\perp\) is not proportional to an old
   \(F\) norm-two mirror and satisfies (19).
2. **Exact Coxeter angles:** after orienting normals toward one chamber,
   every retained pair of root norms \(n_i,n_j\in\{2,6\}\) satisfies
   
   \[
   h(r_i,r_j)=-\sqrt{n_i n_j}\cos(\pi/m_{ij})
   \]
   
   for a declared integral label, or the exact parallel/ultraparallel
   alternative.  Condition (15) is part of this arithmetic gate; for
   example a norm-six/norm-two product \(-3\) gives label \(6\).
3. **Trace lattice:** the trace form is even of signature \((1,21)\), and
   the proposed root remains primitive.
4. **Primitive OG10 gluing:** find a positive rank-two complement and an
   explicit primitive overlattice isometric to
   \(3U\oplus2E_8(-1)\oplus A_2(-1)\), with the selected class in (23)
   surviving as the unique ambient order-three radical line.
5. **Wall completeness:** enumerate every \((-2,1)\) and \((-6,3)\)
   stably-prime-exceptional wall meeting the golden eigenslice.  If a Kähler
   rather than movable cone is claimed, the \((-4,1)\) and \((-24,3)\)
   small-contraction types must also be screened.
6. **Geometric chamber selection:** period realization, projectivity, and a
   marked birational model must select the intended chamber.
7. **Volume:** only a completed arrangement can be certified finite-volume
   or compact.

The simple exact test model \(k=(0^9,1)\) demonstrates that the arithmetic
ansatz is nonempty: its field determinant is \(-6\phi\), its trace lattice is
even of signature \((1,21)\) and determinant \(-36\), and \(e\) has trace
square \(-6\) and trace-lattice divisibility \(3\).  It does not cut
\(v_{61}\), and no primitive OG10 embedding is claimed for it.

The external geometric input identifying the reflective escape is
Mongardi--Onorati, [Proposition 3.1](https://arxiv.org/abs/2010.12511): the
stably prime-exceptional OG10 types are `(-2,1)` and `(-6,3)`.
Their Proposition 5.4 supplies the additional small-contraction Kähler-wall
types `(-4,1)` and `(-24,3)`.  The exact lattice calculations here do not
replace the separate MBM/effectivity and chamber theorems.

## Reproduction

The frozen artifacts are:

- `scripts/certify_og10_golden_rank_two_reduction.py`;
- `data/candidates/og10_golden_rank_two_reduction_no_go.json`;
- `tests/test_og10_golden_rank_two_reduction.py`.

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_og10_golden_rank_two_reduction.py \
  --verify-certificate
PYTHONPATH=src:. pytest -q tests/test_og10_golden_rank_two_reduction.py
```

The verifier recomputes the field and trace signatures, the two involution
identities, the norm-six reflection, both mod-\(3\) radicals, and every wall
test in the three declared boxes before comparing the frozen payload and
source hashes.
