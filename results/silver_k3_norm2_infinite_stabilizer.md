# The complete silver K3 norm-two chamber has infinite volume

## Result

The intrinsic norm-two reflection chamber of the rank-ten silver K3 carrier
does **not** have finite hyperbolic volume.  This is a theorem about the
complete chamber, not a diagnostic on the archived 1,291-wall prefix.

More precisely, let

\[
F=\mathbf Q(\sqrt2),\qquad O=\mathbf Z[\sqrt2],\qquad L=O^{10},
\]

and let \(h\) be the Alexeev-aligned silver form used by the exact Vinberg
run.  At the identity real place it has signature \((9,1)\); at the other
real place it is positive definite.  Let

\[
R_2=\{r\in L:h(r,r)=2\},\qquad
W=\langle s_r:r\in R_2\rangle,
\]

where \(s_r(z)=z-h(z,r)r\).  If \(C\subset\mathbb H^9\) is the complete
arrangement chamber selected by the Vinberg control point and its
height-zero Weyl sector, then

\[
\boxed{\operatorname{Stab}_{O^+(L,h)}(C)\text{ is infinite},\qquad
       [O^+(L,h):W]=\infty,\qquad
       \operatorname{vol}(C)=\infty.}
\]

Thus this carrier is now closed negatively for intrinsic norm-two
reflectivity.  It supplies no new finite-volume hyperbolic Coxeter polytope.
It does not rule out different reflective norms, a different algebraic
carrier, or a general finite-volume polytope in \(\mathbb H^9\).

The production certificate is
`data/candidates/silver_k3_norm2_infinite_stabilizer.json`; its SHA-256 field
is

```text
9b6fb3a41affed6e6b43b967c89d124daeab8ba5f4d0fe05a0c7ff6ccf2b9b12
```

The standard-library-only hostile replay is
`scripts/independent_silver_k3_norm2_infinite_stabilizer.py`, with frozen
summary SHA-256

```text
edbd1429942e9baa6d0f7e76b345cc8ee086bb7c6baaa3b210c0c31995d6b497
```

## 1. The carrier and anisotropy are recomputed directly

The certificate computes an exact \(LDL^T\) decomposition of the actual
rank-ten field Gram matrix.  Its pivots are

\[
2,\ {3\over2},\ {2\over3},\ {1\over2},\ -1-\sqrt2,
2,\ {3\over2},\ {2\over3},\ {1\over2},\ 2.
\]

Consequently

\[
\det h=-2-2\sqrt2,\qquad
\operatorname{sig}_{\mathrm{id}}h=(9,1),\qquad
\operatorname{sig}_{\sigma}h=(10,0).
\]

This does not rely on the inertia of the 1,291-by-1,291 wall Gram matrix.
It also proves anisotropy over \(F\): if a nonzero \(z\in F^{10}\) satisfied
\(h(z,z)=0\), applying the nonidentity embedding \(\sigma\) would give a
nonzero isotropic vector for the positive-definite form \(h^\sigma\), which
is impossible.

The control vector is

\[
p=(\sqrt2-1,\sqrt2-1,2-\sqrt2,1,1,0,0,0,0,0),
\]

with

\[
p^2=1-\sqrt2,\qquad Gp=(-1,0,\ldots,0).
\]

For an oriented positive-height root \(r\), therefore,

\[
(p,r)=-x,
\]

where \(x>0\) is the first coordinate of \(r\).

## 2. The equality-inclusive separator lemma

The global step is a short Lorentzian estimate.  Put

\[
P=p^2<0,\qquad V=v^2<0,\qquad C=(p,v)<0.
\]

Let \(r^2=2\), orient \(r\) by \((p,r)=-x<0\), and suppose

\[
(v,r)\geq 0.
\]

The weak inequality is important: it includes an unseen mirror containing
\(v\), not only one lying strictly between \(p\) and \(v\).  Write
\(A=-P>0\), \(D=-V>0\), and decompose orthogonally to \(p\).  Then

\[
r_\perp^2=2+{x^2\over A},\qquad
v_\perp^2={C^2\over A}-D,
\]

and \((v,r)\geq0\) implies
\((v_\perp,r_\perp)\geq-Cx/A\).  Cauchy--Schwarz in the
positive-definite space \(p^\perp\) gives

\[
\boxed{
x^2\leq {2\bigl(C^2+(-P)V\bigr)\over -V}.
}
\tag{1}
\]

This is a future-proof bound: it rules out every later shell at once once
the right-hand side lies below an already exhausted Vinberg height.

## 3. Why the control point's height-zero mirrors cause no gap

The point \(p\) is not regular: it lies on nine norm-two mirrors.  The proof
therefore explicitly redoes the height-zero calculation.

The positive trace blocks have exact short-vector counts

\[
(240,240,4),
\]

and their field norm-two counts are

\[
(32,48,2).
\]

Thus all 82 height-zero roots are exhausted.  Reflection closure of the nine
archived simple roots has exactly the same 82 elements, of finite type

\[
B_4+F_4+A_1.
\]

Both fixed witnesses below are strictly negative on those nine simple
roots.  Hence the segment leaving \(p\) toward either witness immediately
enters the chosen finite Weyl sector and crosses no height-zero mirror.

Now suppose a witness were not regular in the complete chamber.  Local
finiteness of the norm-two arrangement supplies a first mirror met by the
segment (allowing the witness itself as the first meeting point).  After the
height-zero exhaustion, an oriented root of such a mirror has \(x>0\), and
the mirror is a wall of the current arrangement chamber: equivalently it is
a Vinberg simple wall.  It also satisfies \((v,r)\geq0\), so (1) applies.

Every admissible positive \(x<4\) is

\[
1,\ \sqrt2,\ 2,\ 1+\sqrt2,\ 2+\sqrt2,\ 1+2\sqrt2.
\]

All these levels were exhaustively processed in the frozen Vinberg chain.
Therefore, if (1) gives \(x<4\) (or \(x<3\)), every possible first wall is
already among the archived accepted walls.  Strict negativity on those
walls then gives the contradiction.  This is the precise reason a bounded
1,291-wall computation can certify membership in the *complete* chamber.

For completeness, the arrangement is locally finite.  If a root mirror
meets a fixed compact subset of \(\mathbb H^9\), the identity-place vector
\(r\) is uniformly bounded because it has norm two in the positive-definite
orthogonal complement of a point on that compact set.  At the conjugate
place, \((r^\sigma)^2=2\) in a fixed positive-definite form, so \(r^\sigma\)
is also bounded.  The Minkowski embedding of \(O^{10}\) is discrete in
\(\mathbb R^{20}\), leaving only finitely many such roots.

## 4. Two exact complete-chamber involutions

### The seed reversal \(g\)

The first ten archived walls are a basis.  Sending them by

```text
(0 1)(2 9), with 3,4,5,6,7,8 fixed
```

produces an \(O\)-integral isometric involution \(g\).  Its matrix digest is

```text
bf8b56c5d8faeec170f2bb863a3526a1d525353a9c278fd3ee661ec7f5282a37
```

It fixes the exact vector

\[
\begin{aligned}
v_g={}&(1,
{147-30\sqrt2\over71},
{5+112\sqrt2\over71},
1+2\sqrt2,
{141+77\sqrt2\over71},\\
&{54+31\sqrt2\over71},
{96+63\sqrt2\over71},
{96+63\sqrt2\over71},
{54+31\sqrt2\over71},
{12-\sqrt2\over142}).
\end{aligned}
\]

The exact checks give

\[
v_g^2={3614-3011\sqrt2\over5041}<0,qquad
(p,v_g)=-1.
\]

The largest pairing with any of the 1,291 archived walls is

\[
{-12+\sqrt2\over71}<0.
\]

Formula (1) gives

\[
x^2\leq {4620+2005\sqrt2\over503}<16,
\]

so any first wall would have \(x<4\), an exhausted range.  Thus \(v_g\) is
a regular point of the complete chamber.  Since \(g\) fixes this future
timelike point, it preserves the chosen component of the negative cone; as
an arrangement isometry fixing \(v_g\), it stabilizes \(C\).

### The wall swap \(H\)

Archived walls 4 and 18 satisfy

\[
(r_4,r_{18})=-\sqrt2.
\]

Put \(d=r_4-r_{18}\).  Then

\[
d^2=4+2\sqrt2,qquad
H=I-{2\over d^2}d(d^TG)
  =I-\left(1-{\sqrt2\over2}\right)d(d^TG).
\]

The resulting matrix is \(O\)-integral, isometric, involutive, and swaps
\(r_4,r_{18}\).  Its digest is

```text
c754c65f8ea3f64d76163bdb9d077c292292329ecdf05b0e56f59990c9e419cb
```

The advertised boundary minimizer is independently checked to be fixed by
\(H\), to have norm \(-\sqrt2/3\), and to vanish precisely on archived walls

```text
0, 1, 2, 3, 5, 6, 7.
```

The strict midpoint witness \(v_H\) recorded literally in the certificate
satisfies

\[
v_H^2={29766\over508369}
      -{297893\over1016738}\sqrt2<0,qquad
(p,v_H)=-1.
\]

Its largest archived-wall pairing is

\[
-{16\over713}-{35\over1426}\sqrt2<0,
\]

and (1) gives

\[
x^2\leq {461210+253640\sqrt2\over171073}<9.
\]

Thus any first wall has \(x<3\), again an exhausted range.  The same
regular-point and future-cone argument proves \(H\in\operatorname{Stab}(C)\).

## 5. The stabilizer is infinite

For the product \(gH\), exact Faddeev--LeVerrier arithmetic gives

\[
\chi_{gH}(t)
=(t-1)^7(t+1)\bigl(t^2-(4+4\sqrt2)t+1\bigr)
\]

and

\[
\operatorname{tr}(gH)=10+4\sqrt2>10.
\]

Its matrix digest is

```text
f8672080597bbb10daf78bc893ffe4e17a763cb0e04df18d114ddce22c22db43
```

A finite-order real \(10\times10\) matrix has roots of unity as eigenvalues,
so the absolute value of its trace is at most 10.  Hence \(gH\) has infinite
order, and \(\operatorname{Stab}_{O^+(L,h)}(C)\) is infinite.

## 6. Infinite stabilizer implies infinite chamber volume

Set \(\Gamma=O^+(L,h)\).  The following chain separates direct deductions
from standard theorem inputs.

1. The positive-definite conjugate form proves that \(h\) is anisotropic
   over \(F\), as shown in Section 1.
2. By Borel--Harish-Chandra and the anisotropy/Godement compactness
   criterion, \(\Gamma\) is a uniform arithmetic lattice in
   \(\operatorname{Isom}(\mathbb H^9)\).  Therefore
   \(0<\operatorname{vol}(\Gamma\backslash\mathbb H^9)<\infty\).
3. Every \(\gamma\in\Gamma\) permutes \(R_2\), and
   \(\gamma s_r\gamma^{-1}=s_{\gamma r}\).  Hence \(W\triangleleft\Gamma\).
4. For the locally finite complete mirror arrangement, the standard
   reflection-group chamber theorem says that \(W\) acts simply transitively
   on its chambers.  Consequently

   \[
   \Gamma=W\rtimes\operatorname{Stab}_\Gamma(C),\qquad
   \Gamma/W\cong\operatorname{Stab}_\Gamma(C).
   \]

5. Since the stabilizer contains the infinite-order element \(gH\),
   \([\Gamma:W]=\infty\).  The orbifold covering
   \(W\backslash\mathbb H^9\to\Gamma\backslash\mathbb H^9\) therefore has
   infinite degree.  Volume multiplication gives

   \[
   \operatorname{vol}(C)
   =\operatorname{vol}(W\backslash\mathbb H^9)
   =[\Gamma:W]\operatorname{vol}(\Gamma\backslash\mathbb H^9)
   =\infty.
   \]

The cited inputs in steps 2 and 4 are not re-proved by the computational
certificate.  The arithmetic-lattice and compactness input is
[Borel--Harish-Chandra, *Arithmetic Subgroups of Algebraic Groups*, Annals of
Mathematics 75 (1962), 485--535](https://doi.org/10.2307/1970210).  The
reflection-chamber input is the standard framework of
[Vinberg, *Discrete Linear Groups Generated by Reflections*, Math. USSR-Izv.
5 (1971), 1083--1119](https://www.mathnet.ru/eng/im2148).

## 7. Replay commands

```bash
PYTHONPATH=src python3 scripts/certify_silver_k3_norm2_infinite_stabilizer.py --check
python3 scripts/independent_silver_k3_norm2_infinite_stabilizer.py
PYTHONPATH=src python3 -m pytest -q tests/test_silver_k3_norm2_infinite_stabilizer.py
```

The hostile replay imports no project search or exact-arithmetic module.  It
hard-codes the field Gram matrix, implements its own
\(\mathbf Q(\sqrt2)\) arithmetic, independently recomputes the full-form
signatures, exhausts the height-zero roots through an exact trace-lattice
short-vector search, reconstructs both involutions, checks integrality and
isometry, checks both fixed witnesses and separator bounds, and recomputes
the characteristic polynomial.
