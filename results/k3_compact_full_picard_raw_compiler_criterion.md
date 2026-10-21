# A finite iff compiler for compact full-Picard K3 eigen-nef returns

**Status (2026-08-11): proved as a finite decision theorem for a pinned exact
raw arithmetic realization.** The compiler can return `EMPTY`. Therefore this
is not the missing essential-surjectivity assertion that every raw maximal
arithmetic chamber has a K3 carrier, and it is not Property 4 by itself.

The gain is an exact separation of the remaining compact nonrational K3
problem. Once the raw chamber is given by exact algebraic wall lines, the
existence of a primitive full-trace, full-Picard complex-K3 carrier whose
complete eigen-nef section is that chamber is decidable by a finite list:

1. a finite simultaneous physical trace-root normalization census;
2. a finite even-overlattice and primitive K3-embedding census; and
3. a finite exhaustive test for every ambient square-minus-two root whose
   restriction can meet the compact chamber.

If one candidate passes, periods, the K3 Weyl theorem and the K3 cone theorem
construct the surface. Every returned facet is then supported by a smooth
rational square-minus-two curve represented by the same vector in the trace
module. Thus the physical in-trace star used in
[the compact host theorem](k3_compact_physical_star_trace_host_finiteness.md)
is a conclusion, not an extra premise, on every accepted full-Picard branch.

Raw reflection maximality does **not** generally delete the third test. An
ordinary K3 root reflection has a one-dimensional rational minus eigenspace
and need not preserve the distinguished eigenspace; a native field reflection
has rational minus-eigenspace dimension \([k:\mathbf Q]\) and need not
preserve the candidate lattice. Maximality is a valid shortcut only after a
separate same-host native-reflection condition has been checked.

## 1. Exact raw input and desired output

Let \(k\) be a totally real field of degree \(d>1\), with distinguished real
embedding \(\tau _0\), and let \(V\) be an \(m\)-dimensional \(k\)-vector
space, \(m=n+1\ge3\). A **pinned raw compact arithmetic packet** is

\[
 \mathcal P=(k,\tau _0,V,[g_0],\{\ell_a^+\}_{a\in A}),       \tag{1.1}
\]

where:

- \([g_0]\) is a \(k\)-similarity class of symmetric forms with signature
  \((1,m-1)\) at \(\tau _0\) and negative-definite conjugates at every other
  real place;
- the \(\ell_a\) are oriented \(k\)-rational negative lines whose
  distinguished hyperplanes are the irredundant facets of an exact compact
  Coxeter polytope \(P\subset\mathbf H^{m-1}\); and
- the algebraic Gram entries, chosen halfspaces, Coxeter/Vinberg predicates,
  and, when desired, raw reflection maximality are part of the exact
  presentation.

This is still surface-free and lattice-free. The pinning is only the exact
finite presentation of an arithmetic polytope; it does not choose an integral
normalization, a trace module, a Picard lattice, a K3 embedding or a period.
An abstract diagram without exact dotted distances does not specify this
existence problem.

A **successful full-Picard K3 return** consists of

\[
 (S,\mathcal O,M,g,B,j,E,p),                             \tag{1.2}
\]

where \(\mathcal O\) is an order in \(k\), \(M\subset V\) is a finite
projective rank-\(m\) order module, \(g\in[g_0]\), and

\[
 B(x,y)=\operatorname {Tr}_{k/\mathbf Q}g(x,y)           \tag{1.3}
\]

makes \(N=M_{\mathbf Z}\) even integral. Every raw facet reflection preserves
\(M\). The embedding

\[
             j:N\hookrightarrow\Lambda_{K3}
\]

is primitive, the projective K3 surface has
\(\operatorname {NS}(S)=N\), and, after a marked Weyl framing, the selected
point satisfies

\[
       p\in E\cap\operatorname {Amp}(S),
       \qquad [p]\in P^\circ,                            \tag{1.4a}
\]

while

\[
 \boxed{
 P=\mathbf P\bigl(\operatorname {Nef}(S)\cap E\bigr)
  =\overline{\mathbf P\bigl(\operatorname {Amp}(S)\cap E\bigr)}.}
                                                               \tag{1.4}
\]

Here \(E=V^{\tau _0}\) is the distinguished real summand for the conjugated
self-adjoint order action. Equality in (1.4) is the exact nonrational
replacement for literal equality with the full nef cone. Condition (1.4a)
is part of the return contract, not an inferred decoration: it rules out an
eigenspace lying inside a root wall.

## 2. The simultaneous physical normalization census is finite

Choose a facet of \(P\) and its full facet star. Compactness implies that the
star normals span \(V\). Enumerate the finitely many \(m\)-subsets that are
\(k\)-bases. Fix one such subset and representatives

\[
             r_1,\ldots,r_m,
 \qquad h_i=g_0(r_i,r_i).                               \tag{2.1}
\]

The nonzero-pairing graph of this basis has one Lorentz block. Its
scaling-invariant pair and chordless Cartan cycles generate \(k\). In every
integral trace carrier these cycles are traces of products of the rank-one
operators \(1-s_{r_i}\), hence are algebraic integers in its coefficient
order. Let

\[
       \mathcal O_0=\mathbf Z[\text{those cycles}].      \tag{2.2}
\]

If the cycles are not integral or do not generate \(k\), this basis branch is
rejected. Otherwise \(\mathcal O_0\) is a full order in \(k\). It is enough
to bound the census: any successful coefficient order lies in the finite
interval

\[
             \mathcal O_0\subseteq\mathcal O
                    \subseteq\mathcal O_k.              \tag{2.2a}
\]

Enumerate all multiplicatively closed additive intermediate groups in the
finite quotient \(\mathcal O_k/\mathcal O_0\). Thus every possible order is
retained, while arbitrary unbounded conductors need not be guessed.

### 2.1 The finite diagonal energy shell

Suppose a successful carrier exists and choose the actual physical roots
\(\delta_i\) on these \(m\) wall lines. Write

\[
          g=\lambda g_0,\qquad \delta_i=a_i r_i,
          \qquad \lambda,a_i\in k^\times.              \tag{2.3}
\]

The physical square condition says

\[
 u_i:=g(\delta_i,\delta_i)=\lambda a_i^2h_i,
 \qquad \operatorname {Tr}_{k/\mathbf Q}(u_i)=-2.      \tag{2.4}
\]

Every \(u_i\) is totally negative. Moreover, for
\(a\in\mathcal O_0\),

\[
 \operatorname {Tr}(a u_i)=B(a\delta_i,\delta_i)\in\mathbf Z,
\]

so \(u_i\in\mathcal O_0^\vee\), the order codifferent. Consequently every
diagonal belongs to the finite set

\[
 \mathscr E_{-2}(\mathcal O_0)=
 \left\{u\in\mathcal O_0^\vee:
   \tau(u)<0\ \forall\tau,
   \operatorname {Tr}_{k/\mathbf Q}u=-2\right\}.       \tag{2.5}
\]

Indeed the positive place energies \(t_\tau=-\tau(u)\) satisfy

\[
                         \sum_\tau t_\tau=2,            \tag{2.6}
\]

so the Minkowski image of \(u\) lies in a compact box. A codifferent is a
full Euclidean lattice, and its intersection with that box is finite and
effectively enumerable.

### 2.2 Common square classes reconstruct the field Gram

For a tuple
\((u_1,\ldots,u_m)\in\mathscr E_{-2}^m\), equations
(2.3)--(2.4) have a solution precisely when the totally positive elements

\[
                         u_i/h_i                        \tag{2.7}
\]

belong to one common coset of \(k^{\times2}\). Equivalently,

\[
        \frac{u_i/h_i}{u_1/h_1}\in k^{\times2}
        \quad(2\le i\le m).                            \tag{2.8}
\]

This is an exact number-field square test. Gauge \(a_1=1\), put
\(\lambda=u_1/h_1\), and choose the unique signs of the square roots
\(a_i\) that respect the oriented distinguished normals. The complete field
Gram is forced:

\[
                     G_{ij}=\lambda a_i a_jg_0(r_i,r_j). \tag{2.9}
\]

The gauge loses no marked host. In an ungauged solution, multiplication of
\(V\) and \(M\) by \(a_1^{-1}\), together with multiplication of the form by
\(a_1^2\), is an oriented isometry of trace packages that fixes every raw
wall line. Thus it is exactly the declared marked raw-similarity equivalence.

Thus a fixed raw packet leaves only finitely many physical star Grams. Test
each one for exact raw oriented similarity, admissible signatures, and the
codifferent/integrality conditions. The broader narrow-star estimate

\[
 |\tau_0(G_{ij})|<14,\qquad |\tau(G_{ij})|\le2
 \quad(\tau\ne\tau_0)                                  \tag{2.10}
\]

is an independent finite-box check useful for a class-wide host bound, but
(2.5)--(2.9) already make the normalization census finite for one fixed raw
\(P\).

### 2.3 Only finitely many integral trace modules remain

For one surviving Gram put

\[
                  F_0=\bigoplus_{i=1}^m
                         \mathcal O_0\delta_i.          \tag{2.11}
\]

Its complete integral Gram in a chosen \(\mathbf Z\)-basis of
\(\mathcal O_0\) is explicit. Every successful trace lattice satisfies

\[
                         F_0\subseteq M\subseteq F_0^\#. \tag{2.12}
\]

The quotient \(F_0^\#/F_0\) is finite. Enumerate its subgroups and pair them
with the finitely many orders in (2.2a). For each \(M\) compute its multiplier
order

\[
       \mathcal O_M=\{a\in k:aM\subseteq M\}.           \tag{2.12a}
\]

This is an effectively computed overorder of \(\mathcal O_0\) contained in
\(\mathcal O_k\). Retain a pair \((\mathcal O,M)\) exactly when
\(\mathcal O_0\subseteq\mathcal O\subseteq\mathcal O_M\), the lattice is
even and integral, \(M\) is projective over that \(\mathcal O\), and every
raw facet field reflection preserves \(M\). This is a finite computation in a
discriminant module. A successful module need not be projective over
\(\mathcal O_0\), nor does projectivity have to ascend to
\(\mathcal O_M\); enumerating every intermediate order is what makes the
filter exhaustive.

For each nonbasis facet line \(\ell_a\), enumerate

\[
       \mathscr R_a(M)=
       \{\delta\in M\cap\ell_a:B(\delta,\delta)=-2\}/\{\pm1\}. \tag{2.13}
\]

This set is finite: the trace form restricted to the \(\mathbf Q\)-space
underlying a totally negative \(k\)-line is negative definite. Reject the
candidate unless every \(\mathscr R_a(M)\) is nonempty. This is simultaneous
in-trace root normalization of **all** facets, not equality of projected
lines with unrelated Picard roots.

> **Finite-host lemma.** From a pinned compact raw packet one can compute a
> finite list \(\mathfrak H(P)\) containing, up to marked raw similarity,
> every even integral full-trace host in which every facet line contains a
> square-minus-two vector. The list can be empty.

This is the exact step unavailable in the projected-root scaling family:
there the physical roots lie outside the trace block, so equations
(2.4)--(2.6) do not constrain the escaping trace scale.

## 3. Primitive K3 embedding is a finite exact gate

For \(N=M_{\mathbf Z}\) first check

\[
                \operatorname {rk}N=dm\le20,
 \qquad \operatorname {sig}N=(1,dm-1).                 \tag{3.1}
\]

Compute the finite discriminant quadratic form

\[
                         (A_N,q_N)=(N^\#/N,q_N).        \tag{3.2}
\]

There is a primitive embedding
\(N\hookrightarrow\Lambda_{K3}\) if and only if there is an even lattice
\(T\) of signature

\[
                         (2,20-dm)                      \tag{3.3}
\]

with discriminant form

\[
                         q_T\simeq-q_N.                 \tag{3.4}
\]

Necessity is discriminant anti-glue in the unimodular K3 lattice. Conversely
the graph of an anti-isometry \(A_N\simeq A_T\) glues \(N\oplus T\) to the
unique even unimodular lattice of signature \((3,19)\). Exact primary Jordan
decomposition and the local even-lattice existence conditions decide
(3.3)--(3.4), including all nonstable-range cases. If marked embedding
components are wanted, enumerate the finitely many complementary genera and
the finite anti-isometry double cosets.

Failure of a stable-range length inequality is not rejection; only failure
of the exact discriminant-form gate is.

## 4. Compactness makes the complete hidden-root test finite

Fix a surviving \(N\) and its distinguished summand \(E\). Let
\(p\in P^\circ\) be an exact algebraic point. A K3 root is any

\[
                         \delta\in N,
 \qquad B(\delta,\delta)=-2.                            \tag{4.1}
\]

There can be infinitely many roots in \(N\), but only finitely many can have
a restricted mirror meeting compact \(P\).

### 4.1 All-place boundedness

Use the orthogonal real decomposition

\[
 N_{\mathbf R}=E_{\tau_0}\oplus
       \bigoplus_{\tau\ne\tau_0}E_\tau.                \tag{4.2}
\]

If the proper restricted mirror of \(\delta\) contains \(x\in P\), then
\(\delta_{\tau_0}\perp x\). Hence the distinguished component lies in the
negative-definite space \(x^\perp\), while every nonidentity summand is
already negative definite. The identity

\[
 \sum_\tau B_\tau(\delta_\tau,\delta_\tau)=-2           \tag{4.3}
\]

expresses \(2\) as a sum of nonnegative place energies. Each component energy
is at most \(2\). As \(x\) ranges over compact \(P\), the negative forms on
\(x^\perp\) remain uniformly definite. Thus all such Minkowski components
lie in one bounded ellipsoid. The lattice \(N\) has finite intersection with
it.

Roots with zero distinguished projection must also be included. Their
restricted hyperplane is all of \(E\), so no point of \(E\) can be ample;
they are a terminal rejection, not an invisible harmless root. They lie in
the same finite search because their nonidentity energies sum to \(2\).

### 4.2 An explicit one-ellipsoid bound

For an executable bound, let

\[
 R=\max_{x\in P}d_{\mathbf H}([p],[x]),                 \tag{4.4}
\]

which is computable from the finitely many vertices of compact \(P\). If a
proper root mirror meets \(P\), the hyperbolic point-to-hyperplane formula
gives

\[
 \frac{|B(p,\delta)|^2}{B(p,p)}\le2\sinh^2R.            \tag{4.5}
\]

The real quadratic form

\[
 Q_p(z)=\frac{2B(p,z)^2}{B(p,p)}-B(z,z)                \tag{4.6}
\]

is positive definite: writing \(z=ap+z_\perp\) gives
\(Q_p(z)=a^2B(p,p)-B(z_\perp,z_\perp)\). Every relevant root therefore
satisfies

\[
                         Q_p(\delta)\le2\cosh(2R).      \tag{4.7}
\]

Equation (4.7) is a terminating lattice-point enumeration. Exact algebraic
linear programming on a Klein slice then decides whether
\(\delta^\perp\cap P^\circ\) is nonempty and whether the hyperplane equals a
prescribed facet.

Define \(\Delta_P(N)\) to be the resulting finite signed-root list. The
complete-root gate is:

1. no root has zero distinguished projection;
2. no member of \(\Delta_P(N)\) has a proper mirror meeting \(P^\circ\); and
3. every facet hyperplane of \(P\) is the restriction of at least one member
   of \(\Delta_P(N)\).

These conditions are necessary and sufficient for \(P^\circ\) to be a
connected component of the complete restricted K3 root arrangement. The
second condition puts \(P^\circ\) in one component; the third prevents that
component from extending through any facet. Local finiteness supplies the
converse without a hidden limiting wall.

## 5. Periods, Weyl framing and physical curves add no further gate

Suppose the primitive-embedding and complete-root gates pass. Put
\(T=j(N)^\perp\). Choose an oriented K3 period in \(T\) outside every extra
Noether--Lefschetz hyperplane. Period surjectivity gives a marked K3 surface
with

\[
                         \operatorname {NS}(S)=N.       \tag{5.1}
\]

For \(dm<20\) this is the complement of a countable hyperplane arrangement in
the type-IV domain. For \(dm=20\), the positive rank-two complement is already
the period plane, and exact Picard equality is automatic.

Let \(D\) be the ambient root chamber containing \(p\). Section 4 gives

\[
                              D\cap E=P^\circ.           \tag{5.2}
\]

The K3 Weyl group is transitive on positive root chambers. Re-mark, or
conjugate the self-adjoint \(k\)-action, by a Weyl element taking \(D\) to the
ample chamber. Under the resulting marking the selected \(p\) is ample, so
(1.4a) holds, and the K3 cone theorem gives (1.4). This is a Weyl framing, not
a claim that the Weyl word is induced by an automorphism.

At a generic point of each restricted facet, at least one ambient simple
root wall supports \(D\), even if the restricted facet comes from an ambient
face of codimension greater than one. Indeed, local finiteness leaves only
finitely many ambient simple inequalities there. Their restrictions are
nonnegative on \(D\cap E\); any nonzero restriction vanishing at a generic
facet point vanishes on that whole codimension-one facet. At least one is
nonzero, because the prescribed supporting root has nonzero restriction and
is a nonnegative combination of the ambient simple roots. After Weyl framing
the effective sign of such a simple root is the class of a smooth rational
square-minus-two curve; this is the usual
Riemann--Roch/root-effectivity and K3 nef-cone theorem. Because
\(\operatorname {NS}(S)=M_{\mathbf Z}\), that very curve vector belongs to
the trace module. Equality of its distinguished normal line with the raw
\(k\)-line descends over \(k\), so it is a rescaling of the native raw normal.
Consequently every facet in any chosen star is represented by the same
in-trace vector as its physical curve, and a compact star supplies a spanning
subset.

The canonical requirement has no extra content on this branch:
\(K_S=0\). The self-adjoint order action can be recorded by a rational
divisor-product correspondence; no surface automorphism or action on the
transcendental lattice is asserted.

## 6. The exact compiler theorem

### Theorem 6.1 (finite iff/EMPTY compiler)

For every pinned exact compact admissible arithmetic packet (1.1), the
procedure in Sections 2--4 terminates and returns a finite set
\(\operatorname {K3Comp}(P)\) of marked lattice, embedding and root-chamber
packages. It satisfies

\[
 \boxed{
 \operatorname {K3Comp}(P)\ne\varnothing
 \Longleftrightarrow
 P\text{ has a primitive full-trace, full-Picard complex-K3 return (1.4)}.}
                                                               \tag{6.1}
\]

Every accepted arithmetic package constructs a K3 return by Section 5.
Every such return appears in the list, up to the declared marked raw
similarity and primitive-embedding equivalence.

#### Necessity

In an actual full-Picard return, each facet of the exact eigen-nef section is
supported by an ambient simple square-minus-two root. A compact facet star
spans \(V\), so choose a root basis from it. Reflection-word traces put the
field-generating cycles in the coefficient order, (2.4)--(2.6) put the root
norms in the finite energy shell, and trace duality gives (2.12). Thus its
host occurs in the finite census. The Picard embedding supplies the
discriminant complement, and exact eigen-nef equality forbids every hidden
root in Section 4.

#### Sufficiency

A surviving module is an even hyperbolic trace lattice with physical roots
on every raw facet line. The discriminant complement and period theorem give
an exact-Picard projective K3. The complete-root gate identifies \(P\) with
the restricted ambient root chamber, and Weyl framing identifies that chamber
with the eigen-nef section. Root effectivity supplies the physical curves.

The theorem is an exact existence **decision** for one supplied raw packet.
It does not assert that the finite set is nonempty.

## 7. What maximality can and cannot replace

For a K3 root \(\delta\in N\), the ordinary Picard reflection is

\[
       r_\delta^B(x)=x+B(x,\delta)\delta.               \tag{7.1}
\]

It preserves \(N\), but its rational minus eigenspace has dimension one and
it generally mixes the real field summands. The reflection of the
distinguished eigenspace is instead the native field reflection

\[
 s_\delta^g(x)=x-
       \frac{2g(x,\delta)}{g(\delta,\delta)}\delta,      \tag{7.2}
\]

whose rational minus eigenspace is \(k\delta\), of dimension \(d\). Equation
(7.2) need not preserve \(M\). Therefore neither ordinary K3-root maximality
nor abstract maximality of the displayed facet group excludes a hidden
restricted K3 wall.

There is one safe shortcut. Suppose, for every nonzero-projection root in
the finite window \(\Delta_P(N)\), (7.2) preserves the same marked \(M\), and
suppose \(P\) is reflection-maximal relative to that full native host. A
hidden cut would then adjoin a new reflection inside the discrete arithmetic
group \(O(M,g)\), contradicting relative maximality. Under these two
hypotheses the interior-cut part of Section 4 is automatic. Zero-projection
roots and physical facet witnesses must still be tested.

For comparison, an exact finite-volume **noncompact** full-trace K3 return in
dimension at least two has \(k=\mathbf Q\). Then \(E=N_{\mathbf R}\), (7.1)
is the relevant host reflection, and full-host reflection maximality excludes
a hidden square-minus-two cut. The proved
[bounded closed-star theorem](noncompact_k3_bounded_spanning_root_lemma.md)
supplies a full physical root basis with off-diagonal Gram entries in
\(\{0,\ldots,13\}\). Hence \(R\subseteq N\subseteq R^*\) gives a finite
rational Picard-host census before the same discriminant/period construction.
This makes the pinned, full-host-maximal K3 decision compiler finite on the
cusped lane as well. It remains separate from the compact nonrational
argument: cusp-local data alone do not bound the host, and neither theorem
proves raw Property-4 nonemptiness.

## 8. Exact remaining Property-4 obstruction

The compiler can stop with a rigorous `EMPTY` certificate at any of three
finite gates:

\[
 \begin{array}{c|l}
 \text{gate}&\text{exact obstruction}\\ \hline
 \mathrm{PHYS}&\text{no simultaneous codifferent trace-}(-2)
   \text{ normalization and even full-trace module},\\
 \mathrm{EMB}&\text{no K3 discriminant-form complement},\\
 \mathrm{ROOT}&\text{a zero-projection or chamber-cutting hidden }(-2)
   \text{ root}.
 \end{array}                                             \tag{8.1}
\]

The first nonautomatic datum is `PHYS`: a common square-class-compatible
normalization of a spanning compact star by literal in-trace roots, extended
simultaneously to all facets in one even integral module. Projected physical
walls do not supply it. Once `PHYS`, `EMB`, and `ROOT` pass, there is no
remaining period, Torelli, root-effectivity or scheme-theoretic nonemptiness
condition on the K3 branch.

Universal Property 4 would require the additional theorem

\[
 P\text{ raw maximal in the declared compact nonrational class}
 \Longrightarrow \operatorname {K3Comp}(P)\ne\varnothing
 \quad\text{or a realized alternative AG branch}.      \tag{8.2}
\]

Nothing above proves (8.2). The immediate rank rejection
\(d(n+1)>20\), an empty physical energy/square-class census, a failed
discriminant complement, or an enumerated hidden root can all make the K3
list empty. Thus the result closes the finite **compiler/decision** debt for
this lane while leaving the essential-surjectivity/nonemptiness debt exactly
visible.

## 9. Reproduction

- [certificate generator](../scripts/certify_k3_compact_full_picard_raw_compiler_criterion.py)
- [exact JSON certificate](../data/enumerations/k3_compact_full_picard_raw_compiler_criterion.json)
- [regression tests](../tests/test_k3_compact_full_picard_raw_compiler_criterion.py)

Run:

```text
python3 -m scripts.certify_k3_compact_full_picard_raw_compiler_criterion --check-output
pytest -q tests/test_k3_compact_full_picard_raw_compiler_criterion.py
```

The executable replay checks the finite golden codifferent energy shell, the
explicit positive-definite hidden-root ellipsoid on a Lorentz lattice, the
rank ceiling, the exact gate order, and the maximality/nonemptiness scope.
Number-field lattice enumeration, discriminant-form existence, the K3 period
map and the K3 cone theorem remain cited mathematical algorithms/theorems;
they are not represented as having been proved by the small replay.
