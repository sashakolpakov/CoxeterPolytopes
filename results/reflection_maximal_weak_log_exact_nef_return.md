# A nef weak-log class removes the hidden-wall debt

**Status (2026-08-10): proved project-native exact-return theorem, after
independent hostile audit.**
The theorem generalizes the strong-Coble return argument from
\(K+B/2=0\) to the full weak-log lane

\[
                       D=-\left(K_S+{B\over2}\right)\quad\text{nef}.
\]

It does not infer nefness from reflection maximality.  A compiler must
construct the surface, the boundary and the proposed facet curves, and must
give an independent nef or semiample certificate for \(D\).  Once those
geometric gates hold, relative reflection maximality proves that the proposed
chamber is the entire nef cone; no enumeration of all remaining negative
curves is needed.

## 1. Setup

Let \(S\) be a smooth projective rational surface over an algebraically
closed field.  Let

\[
 B=B_1+\cdots+B_b                                      \tag{1.1}
\]

be a possibly empty disjoint union of smooth rational curves with
\(B_i^2=-4\).  Put \(L=\operatorname {NS}(S)\) and use the adapted host

\[
 M=\left\langle L,{B_1\over2},\ldots,{B_b\over2}\right\rangle,
 \qquad (\ ,\ )_M=2(\ ,\ )_L.                         \tag{1.2}
\]

Assume that

\[
                  D=-\left(K_S+{B\over2}\right)       \tag{1.3}
\]

is nef.  This is a genuine geometric premise.  It may be certified by a
base-point-free linear system, a morphism, or an independently complete
Mori calculation; membership in a merely proposed wall cone is not enough.

Let \(\mathcal C_S^+\) be the component of the positive cone containing the
ample classes.  Let \(P\) be the ideal closure of a finite-sided
finite-volume Coxeter chamber in \(\mathbb P(\mathcal C_S^+)\).  Suppose its oriented
facets are represented by distinct effective irreducible negative curves
\(C_1,\ldots,C_f\) on \(S\), so

\[
 P=\{[x]\in\mathbb P(\overline{\mathcal C_S^+}):
                  (x,C_j)\geq0\text{ for every }j\}.   \tag{1.4}
\]

Assume every displayed facet reflection preserves \(M\), and let \(W(P)\)
be their reflection group.  The group hypothesis is

> \(W(P)\) has no proper discrete reflection overgroup inside
> \(O^+(M)\).

Two stronger hypotheses also suffice:

1. absolute maximality among discrete hyperbolic reflection groups; or
2. maximality among arithmetic **reflection overgroups in the real
   hyperbolic isometry group**, allowing a commensurable change of integral
   host.

Maximality only inside \(O^+(L)\) need not imply maximality inside
\(O^+(M)\) when \(M\ne L\).  Likewise, maximality of a displayed facet list,
root span, or lattice overclosure is not the group hypothesis used here.

### 1.1 What counts as a nef certificate for \(D\)

Because \(D\) can be half-integral, a finite certificate may clear its
denominator.  Any of the following is sufficient:

1. an integer \(m>0\) for which \(mD\) is Cartier and
   \(|mD|\) is base-point-free;
2. a morphism \(g:S\to T\) and a nef divisor \(H\) on \(T\) with
   \(mD=g^*H\); or
3. an independently proved complete curve-cone calculation on which
   \(D\) is nonnegative.

Checking only \(D.C_j\ge0\) on the proposed facet list says merely that
\(D\) belongs to the proposed dual cone.  Using the desired equality
\(P=\mathbb P(\operatorname {Nef}S)\) to promote that membership to nefness
would be circular.  Neither reflection maximality nor Lemma 2.1 removes
this obligation, because the proof of the lemma itself uses \(D.C\ge0\)
for every hidden curve.

## 2. Nefness classifies every negative curve

### Lemma 2.1 (weak-log negative alphabet)

Every irreducible curve \(C\) with \(C^2<0\) is exactly one of:

1. a boundary component \(B_i\), with square \(-4\);
2. a smooth rational \((-2)\)-curve satisfying
   \(D.C=B.C=0\); or
3. a smooth rational \((-1)\)-curve satisfying

   \[
                     D.C+{B.C\over2}=1.               \tag{2.1}
   \]

The reflection in its wall preserves \(M\) in every case.

### Proof

Suppose first that \(C\) is not a component of \(B\).  Distinct effective
irreducible curves have nonnegative intersection, so \(B.C\geq0\).  Since
\(D\) is nef, \(D.C\geq0\).  Equation (1.3) and adjunction give

\[
\begin{aligned}
 2p_a(C)-2
   &=C^2+K_S.C\\
   &=C^2-D.C-{B.C\over2},
\end{aligned}
\]

or equivalently

\[
             C^2=2p_a(C)-2+D.C+{B.C\over2}.           \tag{2.2}
\]

If \(p_a(C)\geq1\), the right side is nonnegative.  Thus a negative curve
has \(p_a(C)=0\).  Since an integral projective curve has
\(p_a=g(\widetilde C)+\sum_p\delta_p\), this means that \(C\) is smooth and
rational, and (2.2) shows
\(C^2\geq-2\).  The square \(-2\) case forces both nonnegative correction
terms to vanish.  The square \(-1\) case gives (2.1).

For a \((-2)\)-curve, every \(B_i.C=0\), so its integral reflection on
\(L\) fixes every half-boundary.  For a \((-1)\)-curve, the doubled form on
\(M\) gives

\[
 s_C(x)=x+2(x.C)C\quad(x\in L),
 \qquad
 s_C(B_i/2)=B_i/2+(B_i.C)C,                            \tag{2.3}
\]

which preserves \(M\).  Finally, \(B_i/2\) has square \(-2\) in the
doubled form, and disjointness gives

\[
 s_{B_i/2}(x)=x+(B_i.x)B_i/2                          \tag{2.4}
\]

while fixing every other half-boundary.  \(\square\)

The lemma uses only actual nefness of \(D\), adjunction and boundary
disjointness.  It does not assume that the displayed curves generate the
Mori cone.

## 3. Reflection saturation gives exact nef return

### Theorem 3.1 (maximal weak-log return)

Under the hypotheses of Section 1,

\[
                         \boxed{P=\mathbf P(\operatorname {Nef}S)}.
\tag{3.1}
\]

### Proof

Let \(R(M)\) be generated by every integral hyperbolic reflection of \(M\).
The group \(O^+(M)\) is discrete in the real hyperbolic isometry group.
For any such reflection \(s\), the group
\(\langle W(P),s\rangle\) therefore remains discrete inside \(O^+(M)\) and
is generated by reflections.  It contains the finite-covolume group
\(W(P)\), so it is a finite-index arithmetic overgroup.  Relative
maximality therefore gives

\[
                              W(P)=R(M).               \tag{3.2}
\]

Choose an ample divisor \(A\) on \(S\).  Every displayed facet curve is
effective, so \(A.C_j>0\).  The explicit choice of
\(\mathcal C_S^+\) in (1.4) therefore puts \([A]\) in the interior of \(P\);
this positive-component orientation is part of the carrier return, not an
automatic property of an unoriented Gram matrix.

Let \(C\) be any irreducible curve.  If \(C^2<0\), Lemma 2.1 says that its
reflection belongs to \(R(M)=W(P)\).  A reflection mirror of \(W(P)\) cannot
cross the interior of its fundamental chamber.  Since \(A.C>0\), the chosen
side is

\[
                         x.C\geq0\quad(x\in P).        \tag{3.3}
\]

If \(C^2\geq0\), then \(A.C>0\) places \(C\) in the closed future causal
cone.  Every class in the closed cone over \(P\) lies in the same future
component.  The Lorentzian reverse Cauchy--Schwarz inequality again gives
\(x.C\geq0\).

Thus every class in the cone over \(P\) is nef.  Conversely a nonzero nef
class \(N\) has \(N^2\ge0\), \(A.N>0\), and nonnegative intersection with
the displayed effective facet curves, so it lies in the right side of
(1.4).  If \(N^2=0\), use \(N+\epsilon A\) and take
\(\epsilon\to0\).  The zero class is restored before projectivization.
This proves equality of the ideal closures
(3.1).  \(\square\)

## 4. Compiler consequence

For a reflection-maximal raw input, a rational-surface compiler in this lane
must still prove:

1. nonemptiness of the cluster/incidence locus;
2. disjoint smooth rational \((-4)\) boundary components;
3. effectivity and irreducibility of every proposed facet curve;
4. exact labelled Gram and orientation recovery; and
5. actual nefness of \(D=-(K+B/2)\), preferably through a finite semiample
   or morphism certificate.

After those five gates, Theorem 3.1 supplies automatically:

- the complete negative-curve alphabet;
- exclusion of every hidden nef wall;
- literal same-\(P\) return; and
- the complete ideal closure of the nef chamber.

This is stronger than the strong-Coble special case because \(D\) may be
nonzero.  If \(D=0\), equation (2.1) specializes to \(B.C=2\) on every
\((-1)\)-curve and recovers the earlier Coble theorem.  If \(B=0\), the
theorem says that an independently certified weak-del-Pezzo surface with a
relative-reflection-maximal proposed chamber has no hidden wall.

Together with the weak-incidence orbifold calculation, a complex carrier in
this lane satisfies

\[
                              \rho(S)\le26,
 \qquad                         \dim P\le25.           \tag{4.1}
\]

The dimension theorem still requires the boundary contraction and orbifold
\(c_2\) hypotheses recorded there; Theorem 3.1 is the exact-return step.

### 4.1 Exact non-Coble calibration

The theorem is nonvacuous outside the strong-Coble locus.  Blow up the torus
fixed point on the positive section of \(\mathbb F_4\), away from its
negative section \(B\).  With \(F\) the pulled-back fibre, \(E\) the
exceptional curve and \(G=F-E\), the complete toric Mori cone is

\[
                       \overline {NE}(S)=\operatorname {Cone}(B,E,G).
\]

In the basis \((B,F,E)\),

\[
 K=(-2,-6,1),\qquad
 D=-\left(K+\frac B2\right)=\frac32B+6F-E.             \tag{4.2}
\]

The exact pairings are

\[
 (D.B,D.E,D.G)=(0,1,1/2),                               \tag{4.3}
\]

so the independent toric cone proves that \(D\) is nef.  On the adapted
host the three roots \(B/2,E,G\) have Coxeter angles
\((\pi/2,\pi/3,0)\), giving the reflection-maximal
\((2,3,\infty)\) triangle.  The theorem correctly returns its already known
full nef cone even though \(B\ne-2K\).  This calibration is replayed in the
[hostile audit](reflection_maximal_weak_log_exact_nef_return_hostile_audit.md).

## 5. Exact limits

This theorem does **not** assert that:

- maximality produces the surface, the boundary, or the weak-log divisor;
- membership of \(D\) in a proposed cone proves that \(D\) is nef;
- an arbitrary numerical Stage-0 packet is effective;
- every arithmetic chamber has a weak-log rational carrier; or
- the adapted host is canonical before the boundary has been constructed.

If \(D\) is not independently nef, equation (2.2) has no sign control and
hidden curves of lower square or positive genus may occur.  Their
reflections need not preserve \(M\), so maximality cannot repair the gap.
The nef certificate is therefore exactly the additional geometric property
used by the proof, not cosmetic input.

Relative maximality must also be checked on the **adapted host** \(M\).
Maximality on the raw Neron--Severi lattice alone is insufficient whenever
some half-boundary reflection exists only after passing to \(M\).

The strong-Coble special case and its independent audit are
[here](reflection_maximal_coble_exact_nef_return.md) and
[here](reflection_maximal_coble_exact_nef_return_hostile_audit.md).  The
rank consequence is the
[weak-incidence orbifold theorem](rational_wrong_log_orbifold_rank_gate.md)
with its
[hostile audit](rational_wrong_log_orbifold_rank_gate_hostile_audit.md).
