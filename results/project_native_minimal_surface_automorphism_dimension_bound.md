# Minimal-surface canonical-isometry bound for intrinsic Coxeter trace modules

**Status (2026-08-10): project-native theorem, nonempty calibration, and
universality obstruction.** No Nikulin--Vinberg, ABSW, Prasad, or Mostow
input is used.

If the full restriction-of-scalars Coxeter representation is intrinsic in
the numerical Neron--Severi lattice of one minimal smooth projective complex
surface and every simple reflection extends to a pairing isometry fixing the
canonical class, then

\[
                    \boxed{[k:\mathbf Q](n+1)\leq20.}             \tag{0.1}
\]

Thus \(n\leq\lfloor20/[k:\mathbf Q]\rfloor-1\). Actual automorphisms are an
important automatic subcase, but are not needed for this rank theorem. In
positive characteristic the unconditional numerical ceiling is \(22\), not
\(20\). This is a genuine AG dimension--field inequality on a nonempty lane.
It is not yet a theorem about every arithmetic Coxeter polytope: there is no
universal polytope-to-minimal canonical carrier compiler, and actual
automorphism reflections cannot be the usual effective nef-wall
reflections.

## 1. Exact hypothesis and canonical orthogonality

Let \(k\) be totally real of degree \(d\). Let \((V,q)\) have dimension
\(n+1\) over \(k\), signature \((1,n)\) at one real place, and be negative
definite at the other real places. Put

\[
 T=\operatorname {Res}_{k/\mathbf Q}V,\qquad
 Q_T=\operatorname {Tr}_{k/\mathbf Q}q.
\]

Then

\[
 \operatorname {rank}T=d(n+1),\qquad
 \operatorname {sig}(T)=(1,d(n+1)-1).                            \tag{1.1}
\]

Let the facet roots \(r_i\) span \(V\) over \(k\). An **intrinsic minimal
canonical carrier** is a minimal smooth projective surface \(X\), an
isometric injection

\[
          \iota:T_{\mathbf Q}\hookrightarrow\operatorname {Num}(X)_{\mathbf Q},
                                                                    \tag{1.2}
\]

and pairing isometries

\[
 u_i\in O(\operatorname {Num}(X)_{\mathbf Q}),\qquad
 u_i(K_X)=K_X,\qquad
 u_i\iota=\iota\operatorname {Res}_{k/\mathbf Q}(s_{r_i}).       \tag{1.2a}
\]

An **intrinsic minimal automorphism carrier** is the subcase
\(u_i=g_i^*\) for actual automorphisms \(g_i\in\operatorname {Aut}(X)\).
An integral primitive injection is the arithmetic version; the rank argument
only needs (1.2)--(1.2a). Using only the distinguished real summand would
lose the factor \(d\).

### Lemma 1.1

\[
                              K_X\perp\iota(T).                    \tag{1.3}
\]

The \((-1)\)-eigenspace of
\(\operatorname {Res}(s_{r_i})\) is the rational space underlying \(kr_i\).
For \(a\in k\),

\[
 (K_X,\iota(ar_i))=(u_iK_X,u_i\iota(ar_i))
                  =(K_X,-\iota(ar_i))=0.
\]

The root lines span \(T\), proving (1.3). Hodge index gives
\(\operatorname {sig}\operatorname {Num}(X)=(1,\rho-1)\). Since \(T\)
already consumes the unique positive direction,

\[
                             \iota(T)^\perp\text{ is negative definite}.
                                                                    \tag{1.4}
\]

## 2. Dimension theorem

### Theorem 2.1 (complex minimal canonical carriers)

For \(n\ge2\), an intrinsic minimal canonical carrier over \(\mathbf C\)
has Kodaira dimension zero and

\[
                            d(n+1)\le\rho(X)\le20.                 \tag{2.1}
\]

| minimal surface | Picard ceiling | trace bound |
|---|---:|---:|
| K3 | \(20\) | \(d(n+1)\le20\) |
| Enriques | \(10\) | \(d(n+1)\le10\) |
| abelian | \(4\) | \(d(n+1)\le4\) |
| bielliptic | \(2\) | no \(n\ge2\) carrier |

**Proof.** If \(\kappa(X)=-\infty\), minimal classification makes \(X\)
rational or ruled with \(\rho\le2\), too small for \(d(n+1)\ge3\).
If \(\kappa(X)\ge0\), minimality makes \(K_X\) nef, hence \(K_X^2\ge0\).
Equations (1.3)--(1.4) force \(K_X\equiv0\). Minimal surface classification
then gives \(\kappa=0\). A complex K3 has \(h^{1,1}=20\), an Enriques
surface has numerical rank \(10\), an abelian surface has \(h^{1,1}=4\),
and a bielliptic surface has numerical rank \(2\). Injection (1.2) finishes
the proof. \(\square\)

This independently eliminates properly elliptic and general-type surfaces.
For \(\kappa=1\), the fixed canonical/Iitaka-fibre direction is nonzero nef
isotropic. For general type, the fixed canonical class has positive square.
Neither can lie in the negative-definite complement (1.4).

### Corollary 2.2 (cohomological and correspondence lane)

The same bound holds when the \(u_i\) arise from monodromy operators,
algebraic correspondences, or other cohomological operators, provided that
each operator actually preserves \(\operatorname {Num}(X)_{\mathbf Q}\), its
intersection pairing, and \(K_X\), and its restriction to \(\iota(T)\) is
the indicated restriction-of-scalars reflection. These conditions are
hypotheses, not automatic properties of an arbitrary correspondence. No
effective-cone conclusion follows in this broader lane.

### Positive-characteristic correction

Over an algebraically closed field of characteristic \(p>0\), with the same
trace-module hypothesis, the corresponding conclusion is

\[
                              d(n+1)\le22.                         \tag{2.2}
\]

K3 surfaces have \(\rho\le b_2=22\), and supersingular K3 surfaces can have
\(\rho=22\); using \(20\) would be false. Enriques surfaces retain numerical
rank \(10\), abelian surfaces have \(\rho\le6\), and the other
Kodaira-zero types are smaller. Recovering \(20\) by specialization requires
the marked module and the canonical-preserving operators to lift.

## 3. Pseudo-automorphisms add nothing to the automorphism lane

Every pseudo-automorphism of a smooth projective surface is an automorphism.
Indeed factorization of birational maps of smooth surfaces into point
blowups and blowdowns shows that any non-biregular map or its inverse has an
exceptional curve. Equivalently, resolve it by a minimal graph

\[
                         X\xleftarrow{p}Z\xrightarrow{q}X.
\]

If a curve \(E\) is \(p\)-exceptional and \(q(E)\) is a curve, that curve is
contracted by the inverse map. If the whole \(p\)-exceptional locus maps to
points under \(q\), the common exceptional part is removed in the minimal
factorization. Thus a minimal factorization has no blowups or blowdowns, and
the map is biregular. A Cremona transformation, Picard--Manin Weyl action,
or isometry between differently marked blowups does not qualify.

## 4. Effective-wall obstruction

If an automorphism sends a nonzero negative rational divisor class \(r\) to
\(-r\), neither sign can be rationally effective. Automorphisms preserve the
effective cone, while an ample class has positive pairing with every nonzero
effective class and therefore prevents both signs from being effective. In
particular, for an integral K3 \((-2)\)-root one sign is effective by
Riemann--Roch, so its Weyl reflection is never induced by a surface
automorphism. It changes nef chambers. This obstruction does not apply to
the canonical-isometry lane of Corollary 2.2.

There is a global form for the actual-automorphism lane, but it must be
stated at the distinguished real place. Decompose
\(T_{\mathbf R}=V_{\sigma_0}\oplus\bigoplus_{\sigma\ne\sigma_0}V_\sigma\).
Let \(\pi_0\) be orthogonal projection to \(V_{\sigma_0}\). Since all other
summands and \(T^\perp\) are negative definite, \(h^2>0\) implies
\((\pi_0h)^2>0\). Hence \(\pi_0(\operatorname {Amp}(X))\) is a nonempty
invariant open convex cone in one component of the positive cone. When the
reflections generate the finite-covolume Coxeter group of \(P\), its limit
set is the whole ideal boundary. In the Klein model, a proper invariant
convex projection would admit a separating hyperplane, whereas the convex
hull of any orbit has closure containing the convex hull of that full limit
set, namely all of hyperbolic space. Therefore

\[
 \pi_0(\operatorname {Amp}(X))=
 \operatorname {Pos}(V_{\sigma_0}).                              \tag{4.1}
\]

For \(d>1\), the analogous claim for the whole rational trace hyperbolic space
would be false: the restriction-of-scalars group lies in a proper real
subgroup and has infinite covolume there.

    actual automorphisms: P is a fundamental chamber in projected Amp(X)
    effective divisors:  P is a nef chamber; reflections are Weyl isometries

This is why strengthening every existing K3 nef carrier to an automorphism
carrier would be wrong.

## 5. Exact nonempty calibration: Wehler's ideal triangle

For a very general \((2,2,2)\) hypersurface \(X\subset(\mathbf P^1)^3\),
\(\operatorname {NS}(X)\) has basis \(H_1,H_2,H_3\) and Gram matrix

\[
 \begin{pmatrix}0&2&2\\2&0&2\\2&2&0\end{pmatrix}.                 \tag{5.1}
\]

The three double-cover deck involutions are actual automorphisms. The
double-cover divisor relation
\(H_i+\iota_i^*H_i=2H_j+2H_k\), together with
\(\iota_i^*H_j=H_j\) and \(\iota_i^*H_k=H_k\), shows that the involution
fixing \(H_j,H_k\) acts as the reflection in

\[
 r_i=H_j+H_k-H_i,\quad r_i^2=-4,\quad (r_i,r_j)=4\ (i\ne j).
                                                                    \tag{5.2}
\]

Indeed the three ideal vertices are the isotropic rays \(\mathbf R_{>0}H_i\),
and \((r_i,r_j)^2=r_i^2r_j^2=16\) for \(i\ne j\). The mirrors are therefore
pairwise asymptotic and bound a finite-area ideal triangle in
\(\mathbb H^2\). Also \(x^2=4(ab+ac+bc)\) for
\(x=aH_1+bH_2+cH_3\), so the lattice represents no \((-2)\)-vector and the
ample cone is the entire positive cone. Hashimoto--Lee record (5.1), the
three deck involutions, and
\(\operatorname {Aut}(X)=C_2*C_2*C_2\cong\Pi(2)\), with \(\Pi(2)\) of index
six in \(PGL_2(\mathbf Z)\), in Sections 1 and 5 of
[*Free automorphism groups of K3 surfaces with Picard number
3*](https://arxiv.org/abs/2306.16147). Thus this is an arithmetic
finite-covolume reflection group, and the intrinsic actual-automorphism lane
is genuinely nonempty.

## 6. Hostile classification audit

| lane | verdict |
|---|---|
| minimal rational/ruled | impossible for \(n\ge2\), since \(\rho\le2\) |
| K3 | possible; \(d(n+1)\le20\); Wehler supplies an actual-automorphism \(H^2\) carrier |
| Enriques | possible only within rank \(10\) |
| abelian | \(d(n+1)\le4\) over \(\mathbf C\) |
| bielliptic | rank \(2\), hence no \(n\ge2\) |
| properly elliptic/general type | impossible for canonical carriers by (1.3)--(1.4) |
| smooth pseudo-automorphism | same as automorphism |
| Cremona/arbitrary marking isometry | fails the actual-automorphism lane; enters the canonical lane only if every condition in (1.2a) holds |

Large arithmetic automorphism group is not enough. A generic unnodal
Enriques surface has finite-index group \(W_{E_{10}}(2)\), but its
bielliptic involutions act as
\({\rm id}_U\oplus(-{\rm id}_{E_8})\), so their action on the **whole**
numerical lattice is not a rational hyperplane reflection. This does not by
itself exclude a smaller common invariant trace block on which their
restrictions are reflections: the definition constrains \(u_i|_T\), not the
action on \(T^\perp\). No finite-covolume common trace block for these
generators is established here;
see Martin, [*Automorphisms of unnodal Enriques
surfaces*](https://arxiv.org/abs/1908.00049). Totaro's Hesse rational
surface admits \(GL_2(\mathbf Z[\zeta_3])\), projectively of finite covolume
on \(\mathbb H^3\), but its action on that full rank-four Lorentzian factor
is the orientation-preserving Bianchi action, not a hyperplane-reflection
group. Again, this full-factor observation does not rule out a smaller
invariant trace block; none is certified here. See
[*Algebraic surfaces and hyperbolic
geometry*](https://www.math.ucla.edu/~totaro/papers/public_html/algebraic.pdf),
Section 10.

For nonminimal rational surfaces Picard rank is unbounded and standard Weyl
groups act by Cremona transformations on markings. Surface classification
alone gives no rank bound after minimality and biregularity are dropped. An
unbounded family in which all simple hyperbolic reflections lift to actual
automorphisms of one rational surface is not established here.

## 7. Maximality and the exact universality gap

Reflection maximality lives inside \(O(T)\) and does not imply surface
minimality or control an exceptional negative complement. If a realization
has a finite invariant zero-dimensional centre, equivariant blowup preserves
the \(T\)-action and its reflection maximality.

On a fixed K3, maximality means every further integral reflection of \(T\)
extending across discriminant glue, fixing the period, and preserving the
common ample cone already lies in \(\Gamma\). This is a Torelli/glue
condition, not nef-wall maximality.

The next falsifiable lemma for the actual-automorphism lane has an
indispensable bounded-complement gate:

> **Minimal automorphism lift lemma.** Prove from raw reflection-maximal
> Coxeter data either that \(T_{\mathbf Q}=\operatorname {NS}(X)_{\mathbf Q}\),
> or that the orthogonal Neron--Severi complement has uniformly bounded rank
> and discriminant. Then enumerate the resulting primitive K3 embeddings,
> complement actions, and discriminant-compatible extensions to \(H^2\);
> decide whether the simple reflections fix a common period and ample
> chamber. Prove that every member of a stated raw Coxeter class has a
> survivor, or exhibit the first exact no-survivor Gram.

Once a particular negative-definite integral complement is fixed, its
orthogonal group and the relevant discriminant glue are finite. But fixed
\(T\) alone does **not** give finitely many Picard enhancements:
Noether--Lefschetz complements can have unbounded discriminant. The missing
statements are bounded complement and **essential surjectivity from a raw
Coxeter class**. Until both are proved, (2.1) is a nonempty conditional
pushback, not an ABSW-type theorem for all arithmetic polytopes.

Finally, (2.1) bounds field degree, not field discriminant or a finite field
list. That still requires the independent integral Gram-envelope/interval
packet.
