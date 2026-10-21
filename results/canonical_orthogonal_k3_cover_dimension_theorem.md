# Canonical-orthogonal and K3-cover dimension theorem

**Status (2026-08-10): exact project-native theorem, a class-wide strong-Coble
corollary, and thirteen frozen carrier calibrations.**  No
Nikulin--Vinberg dimension theorem, ABSW, Prasad, or volume formula is used.

The actual-automorphism hypothesis is stronger than the rank argument needs.
The exact datum consumed by that argument is

\[
 \boxed{\text{one full trace block }T\hookrightarrow K_X^\perp
        \text{ on one minimal surface}.}                         \tag{0.1}
\]

Over the complex numbers this gives

\[
                     [k:\mathbf Q](n+1)\leq20,                  \tag{0.2}
\]

and in positive characteristic it gives the sharp unconditional cohomological
ceiling

\[
                     [k:\mathbf Q](n+1)\leq22.                  \tag{0.3}
\]

Pairing isometries, correspondences, monodromy, and automorphisms are
different ways to certify (0.1); they are not part of its minimal statement.
This distinction turns the four frozen Coble chambers into genuine minimal
K3-cover carriers and gives a direct class-wide corollary:

\[
 \boxed{\dim \mathbf P(\operatorname {Nef}S)\leq19
        \text{ for every complex strong Coble surface},}         \tag{0.4}
\]

The same K3-cover argument by itself gives the weaker cohomological ceiling
\(21\) in characteristic \(p\ne2\).  Independently, the general Coble bound
\(b\leq10\) gives the sharper class-wide ceiling \(19\) in every
characteristic different from two.  None of these statements uses the
classification of finite-automorphism Coble surfaces.

## 1. The exact weakest hypothesis used by the rank proof

Let \(k\) be totally real of degree \(d\), and let \((V,q)\) have dimension
\(n+1\) over \(k\), signature \((1,n)\) at one real embedding, and negative-
definite conjugates at the other embeddings.  Put

\[
 T=\operatorname {Res}_{k/\mathbf Q}V,
 \qquad Q_T=\operatorname {Tr}_{k/\mathbf Q}q.
\]

Thus

\[
 \dim_{\mathbf Q}T=d(n+1),\qquad
 \operatorname {sig}(T)=(1,d(n+1)-1).                           \tag{1.1}
\]

### Definition 1.1

A **canonical-orthogonal full trace carrier** is a minimal smooth projective
surface \(X\) and an injective positive similarity

\[
 j:T\longrightarrow\operatorname {Num}(X)_{\mathbf Q},\qquad
 (j(x),j(y))=cQ_T(x,y),\quad c\in\mathbf Q_{>0},                 \tag{1.2}
\]

such that

\[
                              j(T)\perp K_X.                    \tag{1.3}
\]

Allowing a positive common scale is necessary for finite covers and loses no
Coxeter information.  No facet effectivity, operator, group relation,
faithfulness of an ambient action, or preservation of an ample chamber is
part of Definition 1.1.  The injection must contain the **whole** restriction
of scalars.  A distinguished real eigenspace alone would erase the factor
\(d\).

### Theorem 1.2

For \(n\geq2\), a complex canonical-orthogonal full trace carrier has

\[
             d(n+1)\leq\rho(X)\leq20.                           \tag{1.4}
\]

Over an algebraically closed field of characteristic \(p>0\), the same
hypotheses give

\[
             d(n+1)\leq\rho(X)\leq22.                           \tag{1.5}
\]

**Proof.**  Hodge index says that \(\operatorname {Num}(X)_{\mathbf Q}\)
has one positive direction.  The Lorentzian block \(j(T)\) consumes it, so
\(j(T)^\perp\) is negative definite.  If \(\kappa(X)=-\infty\), minimality
makes \(X\) rational or ruled and \(\rho(X)\leq2\), impossible because
\(d(n+1)\geq3\).  If \(\kappa(X)\geq0\), then \(K_X\) is nef and
\(K_X^2\geq0\).  Condition (1.3) places \(K_X\) in a negative-definite
space, so \(K_X\equiv0\).  Surface abundance/classification puts \(X\) in
Kodaira dimension zero.  Over \(\mathbf C\), the numerical rank ceilings are

```text
K3 20,   Enriques 10,   abelian 4,   bielliptic 2.
```

In positive characteristic, K3 surfaces have \(\rho\leq b_2=22\), Enriques
surfaces have numerical rank \(10\), abelian surfaces have \(\rho\leq6\),
and the remaining Kodaira-zero types are smaller.  This proves (1.4)--(1.5).
\(\square\)

This is the minimal hypothesis for the **canonical-class proof**, not a
claim of logical necessity for every isolated surface.  Dropping minimality
permits arbitrarily many blowups; dropping (1.3) permits the canonical class
to occupy the negative complement; and replacing \(T\) by one real summand
destroys the field-degree factor.

## 2. How operators certify canonical orthogonality

Suppose the facet root lines \(kr_i\) span \(V\).  It is enough to choose a
spanning subset of them and, for each chosen root, a pairing isometry

\[
 u_i\in O(\operatorname {Num}(X)_{\mathbf Q}),\qquad
 u_i(K_X)=K_X,\qquad u_i|_{j(T)}=\operatorname {Res}(s_{r_i}).  \tag{2.1}
\]

The \((-1)\)-eigenspace of \(\operatorname {Res}(s_{r_i})\) is the rational
space underlying \(kr_i\), of dimension \(d\).  Hence

\[
 (K_X,j(ar_i))=(u_i(K_X),u_i(j(ar_i)))=-(K_X,j(ar_i))=0
\]

for every \(a\in k\).  The root lines span, so \(K_X\perp j(T)\).  Notice
what was not used: the \(u_i\) need not satisfy Coxeter relations on
\(j(T)^\perp\), generate a faithful ambient group, preserve an effective
cone, or come from automorphisms.

The same argument applies to a codimension-two correspondence or a
cohomological operator only after checking all four properties in (2.1):
preservation of the numerical lattice, preservation of the pairing,
fixing \(K_X\), and the required restriction.  An arbitrary algebraic
correspondence does not automatically have them.

For a smooth proper family of surfaces, monodromy preserves the intersection
form and the relative canonical class.  To apply (2.1), the full trace block
must additionally be a monodromy-stable algebraic sub-local-system.  Merely
having a local system of the right rank does not put it in
\(\operatorname {Num}(X)\).

Finally, if \(K_X\equiv0\), condition (1.3) is automatic.  A K3, Enriques,
abelian, or bielliptic carrier therefore needs only the full trace
similarity injection for the rank theorem; no reflection lift is required.

## 3. The canonical K3 cover of every strong Coble surface

Let \(S\) be a strong Coble surface over an algebraically closed field of
characteristic different from two.  Here **strong Coble** includes
rationality:

\[
 |-K_S|=\varnothing,\qquad
 |-2K_S|=\{B\},\qquad
 B=B_1+\cdots+B_b
\]

with the \(B_i\) disjoint smooth rational curves.  Adjunction gives
\(B_i^2=-4\), hence

\[
 K_S^2=-b,\qquad \rho(S)=10+b.                                 \tag{3.1}
\]

The second equality in (3.1) uses rationality and Noether's formula
\(\rho(S)=10-K_S^2\); it would not follow from the displayed linear systems
alone.

Take the double cover \(\pi:X\to S\) defined by the line bundle
\(L=-K_S\) and branched over \(B\in|2L|\).  The branch is smooth, so \(X\)
is smooth, and

\[
 K_X=\pi^*(K_S+L)=0.                                           \tag{3.2}
\]

Moreover

\[
 \pi_*\mathcal O_X=\mathcal O_S\oplus\mathcal O_S(K_S),
 \qquad \chi(\mathcal O_X)=2.
\]

Since \(S\) is rational, Riemann--Roch gives
\(\chi(\mathcal O_S)=\chi(\mathcal O_S(K_S))=1\), which proves the displayed
Euler characteristic.  The cover is connected because its branch is
nonempty.  Serre duality and \(K_X=0\) then give \(q(X)=0\), so \(X\) is a
K3 surface.

Put \(N=\operatorname {Pic}(S)\) and

\[
 \widetilde N=\langle N,B_1/2,\ldots,B_b/2\rangle\subset N_{\mathbf Q}.
\]

If \(R_i\subset X\) is the ramification curve over \(B_i\), define

\[
 \jmath:\widetilde N\to\operatorname {NS}(X),\quad
 \jmath|_N=\pi^*,
 \qquad \jmath(B_i/2)=R_i.                                    \tag{3.3}
\]

Equivalently, \(\jmath\) is the restriction of rational pullback to
\(\widetilde N\).  The relation \(\pi^*B_i=2R_i\) shows that its image is
integral, while \(\pi_*\pi^*=2\) proves injectivity.  Projection formula
and \(R_i^2=-2\) give

\[
 (\jmath u,\jmath v)_X=2(u,v)_S.                               \tag{3.4}
\]

Consequently

\[
                 \widetilde N(2)\hookrightarrow\operatorname {NS}(X)
                                                                    \tag{3.5}
\]

is an integral isometric injection of rank \(10+b\).  This is exactly the
scaled extended-Picard lattice on which the boundary reflections become
integral in the frozen Coble certificates.

The K3 Picard bounds now give, without a finite-automorphism classification,

\[
 \begin{array}{ll}
 \mathbf C:&10+b\le20,\quad \dim\mathbf P(\operatorname {Nef}S)=9+b\le19,\\
 \operatorname {char}p\ne2:&10+b\le22,\quad
              \dim\mathbf P(\operatorname {Nef}S)=9+b\le21.
 \end{array}                                                    \tag{3.6}
\]

Thus the cover argument independently recovers \(b\leq10\) over
\(\mathbf C\), and gives only \(b\leq12\) in positive characteristic.
For comparison, the general algebraic-geometric Coble theorem gives

\[
                         b\leq10                              \tag{3.7}
\]

in every characteristic different from two; see Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
Section 3, citing Dolgachev--Kond\=o, Corollary 9.1.5.  Hence the best known
class-wide projectivized-nef ceiling is \(19\) also in positive
characteristic.  The finite-automorphism classification further improves
that ceiling to \(17\) on its locus.  Equation (3.6), rather than these
stronger imported facts, is the direct cohomological output of the
polygon/polytope-to-K3-cover route.

## 4. Facet lifts: correspondence is not the same as monodromy

The frozen Coble rows have only three physical facet types.  Their cover
lifts are exact:

| downstairs facet | class upstairs | reflection lift |
|---|---|---|
| \(B_i/2\) | ramification root \(R_i\), \(R_i^2=-2\) | \(s_{R_i}\) |
| \(E^2=-1\) | \(e=\pi^*E\), \(e^2=-2\) | \(s_e\) |
| \(C^2=-2,\ C\mathbin{\cdot}B=0\) | \(\pi^*C=c+c'\), \(c^2=c'^2=-2,\ c\mathbin{\cdot}c'=0\) | \(s_cs_{c'}\) |

The last row follows because the unramified double cover of
\(C\simeq\mathbf P^1\) splits.  If \(x\) is invariant under the deck
involution, then \(x\mathbin{\cdot}c=x\mathbin{\cdot}c'\), and therefore

\[
             s_cs_{c'}(x)=s_{c+c'}(x).                          \tag{4.1}
\]

For any divisor root \(\delta^2=-2\), the algebraic correspondence

\[
                    \Delta_X+\delta\times\delta                \tag{4.2}
\]

acts on divisor cohomology by
\(x\mapsto x+(x\mathbin{\cdot}\delta)\delta\), the root
reflection.  Compositions of (4.2) give the paired lifts.  Thus every facet
in all four frozen Coble rows has a genuine algebraic correspondence lift;
these are not merely matrices in an abstract orthogonal group.

Over \(\mathbf C\), a \((-2)\)-reflection is a global K3 monodromy
operator: the K3 monodromy group is \(O^+(H^2(X,\mathbf Z))\).  See
Proposition 5.5
of Daniel Huybrechts, [*Lectures on K3 Surfaces*](https://doi.org/10.1017/CBO9781316594193).
Therefore the complex Coble and Enriques lifts below are **generatorwise
actual monodromy elements**.  This does not prove that the chosen lifts of
all simple reflections satisfy every Coxeter relation on the anti-invariant
orthogonal complement.  No single family or coherent Coxeter-monodromy
representation is archived for those rows.

In characteristic three, (4.2) remains an algebraic correspondence and its
action on the numerical lattice remains the required isometry.  The frozen
artifacts do not construct crystalline or étale
Picard--Lefschetz families.
Calling those actions monodromy would therefore overstate the result.

## 5. Honest frozen calibration

The following table records the strongest proved action level, separately
from the rank hypothesis.  `gen.` means generatorwise complex-K3 monodromy;
it does not mean a coherent representation.

| carrier | \(d(n+1)\) | minimal \(K\equiv0\) carrier | correspondence | complex K3 monodromy | coherent Coxeter lift | automorphisms |
|---|---:|---|---|---|---|---|
| Wehler `H2/3` | \(3\) | the K3 itself | graphs | yes | yes | yes |
| Bugaenko `H8/11` | \(2\cdot9=18\) | rank-18 K3 | integral paired roots | yes | yes | no |
| Coble II `H10/30` | \(11\) | canonical K3 cover | root/paired-root | gen. | not certified | no |
| Coble I `H11/26` | \(12\) | canonical K3 cover | root/paired-root | gen. | not certified | no |
| Coble VI char. 3 `H14/60` | \(15\) | supersingular K3 cover | root/paired-root | not applicable | not certified | no |
| Coble MII char. 3 `H17/48` | \(18\) | supersingular K3 cover | root/paired-root | not applicable | not certified | no |
| seven complex Enriques full-nef `H9` types | \(10\) each | Enriques itself, or its K3 cover | paired roots upstairs | gen. | not certified | no |

The Bugaenko row is stronger than a bare abstract Neron--Severi action.  Its
rank-18 trace lattice is the full \(\operatorname {NS}\) of a complex K3,
and every field reflection is the product of two orthogonal K3 root
reflections.  Those products extend by the identity on
\(\operatorname {NS}^{\perp}\); hence their Coxeter relations and
faithfulness survive on full \(H^2\).  They move the nef chamber and are not
surface automorphisms.  See the frozen
[K3 eigen-nef theorem](bugaenko_k3_eigen_nef.md).

For the seven Enriques rows, \(K_Y\equiv0\) and
\(\operatorname {Num}(Y)=U+E_8(-1)\) has rank ten, so the sharper intrinsic
bound \(d(n+1)\leq10\) applies and all seven `H9` chambers attain equality.
A nodal curve on \(Y\) lifts under the étale K3 cover to two disjoint
exchanged roots; their paired reflection restricts to its wall reflection.
The seven exact chambers are stored in
[`enriques_finite_automorphism_nef_chambers.json`](../data/known_examples/enriques_finite_automorphism_nef_chambers.json).

The four Coble ranks are the invariant K3-cover lattice ranks

```text
b=1 -> 11,   b=2 -> 12,   b=5 -> 15,   b=8 -> 18.
```

In particular, their original rational blowups are not being smuggled into
a minimal-surface theorem.  The carrier used here is the smooth canonical
K3 cover.

## 6. What this proves, and the remaining global obstruction

This gives a nonempty multi-row theorem with three genuinely different
geometric mechanisms:

```text
actual automorphisms        : Wehler;
coherent K3 monodromy       : Bugaenko;
K3-cover correspondences    : four Coble rows and seven Enriques rows.
```

It also proves the class-wide strong-Coble ceiling (3.6).  What it does not
prove is that a general compact or finite-volume arithmetic Coxeter
polytope admits any canonical-orthogonal carrier.  Reflection maximality
alone does not construct a K3/Enriques cover, a strong-Coble boundary, or a
full trace injection.  The global theorem target is now sharper:

> From a declared raw maximal arithmetic Coxeter class, construct a minimal
> \(K\)-trivial surface (or prove canonical orthogonality on another minimal
> surface) carrying the **entire** restriction-of-scalars trace block by a
> bounded algebraic packet.

Once that compiler exists, (0.2)--(0.3) give the dimension bound immediately.
Field discriminant and a finite field list still require the independent
integral Gram/interval packet; Picard rank alone does not supply them.

## 7. Machine replay

The compact ledger is
[`canonical_orthogonal_k3_cover_calibration.json`](../data/enumerations/canonical_orthogonal_k3_cover_calibration.json).
The certifier and tests hash every direct frozen input, verify all trace and
cover ranks, replay the scaled pullback pairing, and check the single-root
and paired-root reflection identities:

```text
python3 scripts/certify_canonical_orthogonal_k3_cover_bound.py
python3 -m pytest -q tests/test_canonical_orthogonal_k3_cover_bound.py
```

The source reports for the four Coble rows are the
[Type-II `H10`](coble_type_ii_b3_plane_inverse.md),
[Type-I `H11`](coble_type_i_two_boundary_full_ns_h11.md),
[Type-VI `H14`](coble_type_vi_char3_five_boundary_full_ns_h14.md), and
[Type-MII `H17`](coble_type_mii_full_ns_h17.md) certificates.  Their common
cover construction is developed in the
[K3-cover bridge](coble_k3_cover_polytope_to_surface_bridge.md).
