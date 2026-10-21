# A fixed Jacobian has projective K3 torsors of every cusp index

**Status (2026-08-11): proved algebraic-geometric blindness theorem, with an
exact lattice replay.**  Even the strongest data attached to one K3 cusp do
not bound its divisibility.  There is one fixed Jacobian elliptic K3 surface

\[
 S\longrightarrow\mathbf P^1,
 \qquad \operatorname{NS}(S)\simeq
 U\oplus E_8(-1)^2,
\]

and, for every integer \(d\ge2\), a projective K3 torsor
\(X_d\to\mathbf P^1\) with relative Jacobian \(S\) such that

\[
 \operatorname{NS}(X_d)\simeq U(d)\oplus E_8(-1)^2,
 \qquad \operatorname{div}(F_d)=d.
\]

All members have the same two \(II^*\) fibres and four \(I_1\) fibres, the
same physical
\(\widetilde E_8\sqcup\widetilde E_8\) incident packet, Euler number 24,
Picard rank 18, Jacobian Mordell--Weil rank zero, and discriminant length two.
The incident curves have full local parabolic rank: their common orthogonal
is exactly \(\mathbf QF_d\).  Nevertheless the numerical multisection index
is the arbitrary integer \(d\), and every transverse divisor has degree in
\(d\mathbf Z\).

Thus

\[
\boxed{
 \begin{gathered}
 \text{Kodaira/ADE data, Euler 24, Shioda--Tate, MW rank zero,}\\
 \text{bounded discriminant length, and primitive K3 embedding}\\
 \text{do not bound cusp divisibility or a transverse physical-root degree.}
 \end{gathered}}
\]

This is a **cusp-local and selected-fibration no-go**.  The two affine
\(E_8\) packets are not asserted to be the complete walls of a finite-volume
chamber.  No \(X_d\) is asserted to have finite automorphism group or a
finite-volume complete nef chamber.  Hence this is not a counterfamily to
the global noncompact host theorem.

The exact noncompact field-collapse theorem has already reduced the returned
coefficient field to \(\mathbf Q\).  Every lattice here is integral and hence
already lies over that collapsed field.  The family shows that the field
collapse itself supplies no further control of the integral cusp index.

## 1. The fixed Jacobian

Choose a very general Jacobian elliptic K3 surface \(S\to\mathbf P^1\) with
singular fibres

\[
                       2II^*+4I_1.                 \tag{1.1}
\]

This is the generic member of the two-parameter family with two marked
\(II^*\) fibres.  For example, the Weierstrass shape

\[
 y^2=x^3+a t^4x+t^5(b_0+b_1t+b_2t^2)
\]

has \(II^*\) fibres at zero and infinity and a residual quartic
discriminant; for a general choice it has four simple residual zeros.  The
trivial lattice is \(U\oplus E_8(-1)^2\), and very generality removes every
extra Noether--Lefschetz class.  Therefore

\[
 \operatorname{NS}(S)=U\oplus E_8(-1)^2,
 \qquad T(S)=U^2.                                  \tag{1.2}
\]

Shioda--Tate gives

\[
 18=2+8+8+\operatorname{rank}MW(S),
\]

so \(MW(S)\) has rank zero.  It actually has trivial torsion: a \(II^*\)
fibre has trivial component group, and the height formula leaves no
nonzero torsion section.  The Euler calculation is

\[
                    10+10+4\cdot1=24.             \tag{1.3}
\]

Nothing below varies \(S\).

## 2. Ogg--Shafarevich produces every index projectively

For a complex K3 surface there is a canonical identification

\[
 \operatorname{Br}(S)simeq
 \operatorname{Hom}(T(S),\mathbf Q/\mathbf Z).
 \tag{2.1}
\]

In particular,

\[
 \operatorname{Br}(S)[d]\simeq(\mathbf Z/d)^4.     \tag{2.2}
\]

Write the first copy of \(U\subset T(S)=U^2\) as
\(\mathbf Ze_1+\mathbf Zf_1\).  Choose

\[
 \alpha_d(ae_1+bf_1+ce_2+ef_2)=\frac a d\pmod{\mathbf Z}.
 \tag{2.3}
\]

It has exact order \(d\).  Here one must use the *algebraic*, not the analytic,
Tate--Shafarevich group.  Huybrechts defines the algebraic group
\(\Sha_{\rm alg}(S)\) as the set of projective K3 surfaces with a genus-one
fibration and a specified isomorphism of their relative Jacobian with \(S\),
and proves

\[
                     \Sha_{\rm alg}(S)\simeq\operatorname{Br}(S). \tag{2.4}
\]

(The analytic group is larger and typically parametrizes nonprojective K3
surfaces.)  Thus the algebraic Brauer class \(\alpha_d\) gives a **projective
K3 surface** \(X_d\), a genus-one fibration
\(X_d\to\mathbf P^1\), and an isomorphism of its relative Jacobian with the
fixed \(S\).  Equivalently, Meinsma--Shinder's equation (4.1), Definition 4.2,
and the surjectivity in Theorem 4.13 produce precisely this projective
elliptic-K3 torsor.  This is an existence theorem for every Brauer class, not
just for classes already known to come from a surface.

The exact period--index statements needed here are recorded in
Meinsma--Shinder, Lemma 4.4:

\[
 \operatorname{ind}(X_d)=\operatorname{ord}(\alpha_d)=d,
 \qquad
 T(X_d)\simeq\ker(\alpha_d:T(S)\to\mathbf Z/d).
 \tag{2.5}
\]

The kernel in (2.5) has basis

\[
                  de_1, f_1, e_2, f_2,
\]

and hence

\[
                         T(X_d)\simeq U(d)\oplus U. \tag{2.6}
\]

This proves both projectivity and the exact order/index statement; no
informal appeal to an abstract genus-one torsor is being made.

## 3. Exact Neron--Severi complement

The lattice in (2.6) has signature \((2,2)\) and discriminant group
\((\mathbf Z/d)^2\), of length two.  The stable embedding theorem used here
is Huybrechts, Theorem 14.1.12 (Nikulin, Theorem 1.14.4): an even lattice
\(T\) of signature \((t_+,t_-)\) has a unique primitive-embedding orbit in an
even unimodular lattice \(\Lambda\) of signature \((l_+,l_-)\) if

\[
 t_+<l_+,
 \qquad t_-<l_-,
 \qquad \ell(A_T)+2\le \operatorname{rk}\Lambda-\operatorname{rk}T.
 \tag{3.1}
\]

For \(T=U(d)\oplus U\) and \(\Lambda=\Lambda_{K3}\), both signature
inequalities are strict and

\[
       \ell(A_T)+2=4\le22-4=18.                    \tag{3.2}
\]

There is also an explicit representative, so the complement is not inferred
from determinant or discriminant order.  In

\[
 \Lambda_{K3}=U_1\oplus U_2\oplus U_3\oplus E_8(-1)^2
\]

with bases \(e_i,f_i\), embed \(T=U(d)\oplus U\) using the four vectors

\[
        x=-e_1+df_2,\qquad y=e_2,\qquad e_3,\qquad f_3. \tag{3.3}
\]

The Gram matrix of (3.3) is \(U(d)\oplus U\).  Its coordinate matrix has a
maximal minor \(-1\), using the rows \(e_1,e_2,e_3,f_3\), so the embedding is
primitive.  Direct solution of the four orthogonality equations gives its
orthogonal complement:

\[
  \mathbf Ze_1+\mathbf Z(df_1+e_2)+E_8(-1)^2
       \simeq U(d)\oplus E_8(-1)^2.                 \tag{3.4}
\]

For completeness, the displayed complement is primitive too: the rows
\(e_1,e_2\) and the sixteen \(E_8\) rows give an identity maximal minor.

For every complex K3 surface, \(T(X_d)=\operatorname{NS}(X_d)^\perp\) is a
primitive sublattice of the unimodular lattice
\(H^2(X_d,\mathbf Z)\simeq\Lambda_{K3}\).  Apply the unique-orbit statement
(3.1)--(3.2) to this actual embedding and to (3.3).  An ambient isometry
carries one to the other and therefore carries their orthogonal complements
to each other.  Hence

\[
 \boxed{\operatorname{NS}(X_d)\simeq
        U(d)\oplus E_8(-1)^2}.                    \tag{3.5}
\]

In particular,

\[
 |\operatorname{disc}\operatorname{NS}(X_d)|=d^2,
 \qquad A_{\operatorname{NS}(X_d)}\simeq(\mathbf Z/d)^2. \tag{3.6}
\]

The hosts are pairwise nonisometric, yet their discriminant groups always
need only two generators.  This pinpoints the discriminant-form issue:
**bounded length** is blind to \(d\), while fixing the entire discriminant
form is not blind because its order already contains \(d^2\).

The decomposition can be aligned with the *actual* fibre packet, not merely
with an abstract copy of the Picard lattice.  Choose eight finite \(E_8\)
components in each \(II^*\) fibre and let

\[
                         R=E_8(-1)^2\subset\operatorname{NS}(X_d).
\]

Because \(R\) is unimodular, orthogonal projection to \(R^*=R\) gives an
integral splitting

\[
             \operatorname{NS}(X_d)=R\oplus K,\qquad \operatorname{rk}K=2.
                                                               \tag{3.7}
\]

The primitive fibre \(F_d\) lies in \(K\), and its divisibility is the
multisection index \(d\).  Complete \(F_d\) to a basis \(F_d,D\) of \(K\).
For some integer \(a\), its Gram matrix is

\[
                         \begin{pmatrix}0&d\\d&2a\end{pmatrix}. \tag{3.8}
\]

Here is a short check that (3.8) is really \(U(d)\) with \(F_d\) as its
first standard vector.  The K3 gluing identity gives
\(q_K=q_{\operatorname{NS}}=-q_T\simeq q_{U(d)}\).  The Smith group of
(3.8) is \((\mathbf Z/d)^2\), so \(d\mid2a\).  Replacing \(D\) by
\(D+kF_d\) changes \(a\) by \(kd\); hence the only possible nonzero residue
is \(a=d/2\) when \(d\) is even.  In that exceptional residue, with dual
generators \(P=(D-F_d)/d\) and \(Q=F_d/d\),

\[
                         q(pP+qQ)=\frac{-p^2+2pq}{d}\pmod {2\mathbf Z}.
                                                               \tag{3.9}
\]

Every isotropic element then has even \(P\)-coordinate, so isotropic
elements cannot generate \((\mathbf Z/d)^2\).  This contradicts the two
isotropic standard generators of \(q_{U(d)}\).  Thus \(a\equiv0\pmod d\),
and (3.7) is the actual marked splitting

\[
       (\operatorname{NS}(X_d),F_d,R)
       \simeq (U(d)\oplus E_8(-1)^2,u,E_8(-1)^2).  \tag{3.10}
\]

## 4. The physical cusp packet is already full locally

This fibre assertion does not follow from equality of Euler numbers alone.
Huybrechts, *Lectures on K3 Surfaces*, Chapter 11, Section 4.1, proves for a
genus-one K3 fibration \(X\to\mathbf P^1\) that

\[
                       J^k(X)_t\simeq X_t            \tag{4.1}
\]

for every integer \(k\) and every closed base point \(t\).  His proof uses
the absence of multiple fibres on a K3 and constructs a local analytic or
etale section through a smooth point; consequently \(X\) and \(J^k(X)\) are
locally isomorphic over a cover of the base.  Taking \(k=0\) and
\(J^0(X_d)=S\) shows that the Kodaira fibre schemes in (1.1), including
their component multiplicities, are unchanged.  Thus the two \(II^*\)
fibres on every \(X_d\) supply eighteen actual irreducible smooth rational
\((-2)\)-curves.

For either fibre, let \(C_0,C_1,\ldots,C_8\) be its components and use the
highest-root marks

\[
                    (2,3,4,6,5,4,3,2).
\]

Then

\[
 F_d=C_0+2C_1+3C_2+4C_3+6C_4+5C_5+4C_6+3C_7+2C_8.
 \tag{4.2}
\]

Each nine-curve Gram is affine \(E_8\), of rank eight.  The two vector
spans meet in the common fibre line, so their combined vector span has rank
17, their Gram has rank 16, and

\[
 \left\langle C_i\text{ in both }II^*\text{ fibres}\right\rangle_{\mathbf Q}^{\perp}
                         =\mathbf QF_d.             \tag{4.3}
\]

Equivalently, their total finite ADE rank is

\[
                         8+8=16=\rho(X_d)-2.        \tag{4.4}
\]

This is exactly the full-rank **local parabolic** or MW-rank-zero condition.
It is one dimension short of spanning \(\operatorname{NS}(X_d)_{\mathbf Q}\).
Equations (4.3)--(4.4) do not certify that the complete root chamber has
finite volume: a transverse wall and all other edge closures remain global
data.

## 5. Divisibility and every transverse physical curve escape

In the basis of (3.5), \(F_d=u\).  Its pairing row is

\[
                         (0,d,0,\ldots,0),
\]

so

\[
             \operatorname{div}_{\operatorname{NS}(X_d)}(F_d)=d. \tag{5.1}
\]

This agrees with the geometric index in (2.5).  Consequently every divisor
\(D\) satisfies

\[
                              F_d.D\in d\mathbf Z.  \tag{5.2}
\]

In particular, if \(C\) is any horizontal irreducible curve, and especially
if it is a transverse physical \((-2)\)-curve, then

\[
                              F_d.C\ge d.           \tag{5.3}
\]

Thus the frozen invariant tuple in Sections 1 and 4 cannot imply either the
existence of a transverse physical root of uniformly bounded degree or a
uniform bound on all such roots: for \(d>B\), (5.2) rules out every degree in
\(\{1,\ldots,B\}\).  This conclusion does not require claiming that one
particular abstract root is irreducible.

The lattice nevertheless makes the missing direction visible.  If \(\theta\)
is the highest root in the first \(E_8(-1)\), then

\[
                             c=v-\theta
\]

has

\[
                     c^2=-2,
 \qquad F_d.c=d.                                   \tag{5.4}
\]

It has nonnegative pairing with every component in the displayed affine
packet.  Formula (5.4) is an exact transverse **lattice-root** witness.  It is
not labelled a physical facet without a separate global chamber proof.

Even the automorphisms visible from the chosen fibration do not repair the
escape.  Its translation subgroup is \(MW(S)=0\).  An automorphism preserving
\(F_d\) acts on the base preserving the colored six-point discriminant locus,
and its generic-fibre kernel is finite.  Hence

\[
                         \operatorname{Aut}(X_d,F_d)
                         \quad\text{is finite}.     \tag{5.5}
\]

This is only the stabilizer of the selected fibration, not the full
automorphism group of \(X_d\).

## 6. Exact separation of the global hypotheses

For Picard rank at least three, a complete K3 root chamber has finite volume
if and only if the full K3 automorphism group is finite.  This is a global
condition on **all** root walls and all genus-one fibrations.  It is not a
second local number that can be inserted into Shioda--Tate.

The published Nikulin--Vinberg classification implies that only finitely many
rank-at-least-three Picard lattices have this property.  Therefore full finite
volume certainly bounds all cusp divisibilities after that classification is
imported.  But that finiteness is the reflection-lattice theorem whose
independent replacement is at issue; the later algebraic interpretation as
finite automorphism group does not reprove it.

Same-host reflection completeness is also insufficient by itself.  Every K3
nef cone is the chamber of the complete \((-2)\)-mirror arrangement of its
own Picard lattice.  Saturating reflections after the host is known does not
compare the different lattices in (3.5), and it does not create a bounded
transverse curve.

The exact algebro-geometric audit is therefore:

| input | what it controls | why it misses \(d\) |
|---|---|---|
| Kodaira fibres and Euler 24 | local singular-fibre topology | fixed under every torsor twist |
| Shioda--Tate, MW rank zero | rank of the vertical frame | torsor order is transverse glue |
| multisection/torsor theory | identifies the missing parameter | \(\operatorname{Br}(S)\) has elements of every order |
| Picard rank \(18\) | ambient dimension | constant throughout the family |
| primitive K3 embedding | period realizability | explicit for every \(d\) |
| discriminant length \(2\) | number of glue generators | exponent and order remain unbounded |
| finite \(\operatorname{Aut}(X_d,F_d)\) | selected-fibration symmetries | does not control other fibrations/walls |
| complete finite-volume nef chamber | all global walls | not satisfied or claimed by this family |

## 7. What a successful global theorem adds

A cusp-adapted transverse-wall theorem would immediately kill the escape.  If
one actual transverse physical curve \(C\) satisfies

\[
                         C.C_i\le B
\]

for all nine components of either \(II^*\) fibre, then (4.2) gives

\[
                 d\le F_d.C\le30B.                 \tag{7.1}
\]

Thus the missing input is exactly a bounded physical wall outside
\(F_d^\perp\), or more generally a bounded full-spanning physical-root
packet.  Equation (7.1) also shows why such a theorem is genuinely global:
none of Euler, Shioda--Tate, or the affine diagram supplies \(C\).

There is an important quantifier distinction.  A theorem producing some
full spanning packet with pairwise intersections at most \(B\) need not put
the nine preselected \(II^*\) components in that packet.  Therefore one may
not substitute \(B=13\) into (7.1) unless the global theorem is additionally
cusp-adapted.

The frozen bounded-spanning-root theorem still bounds this family
indirectly if \(X_d\) is assumed to have a finite-volume complete chamber.
For \(\rho=18\), its selected root basis has Gram matrix \(G\) with diagonal
\(-2\) and off-diagonal entries at most 13.  If \(R\) is the lattice it
spans, then

\[
 d^2=|\operatorname{disc}\operatorname{NS}(X_d)|
      \le |\operatorname{disc}R|=|\det G|
      \le (4+17\cdot13^2)^9=2877^9.                \tag{7.2}
\]

The last inequality is Hadamard's row-norm bound.  Thus the global theorem
does close the escape without classifying elliptic fibres, but its automatic
corollary is (7.2), not the sharper conditional value \(d\le390\).

## 8. Sources and replay

The Ogg--Shafarevich, period--index, and transcendental-kernel inputs are
stated explicitly in:

- Reinder Meinsma and Evgeny Shinder,
  [*Derived Equivalence for Elliptic K3 Surfaces and Jacobians*](https://doi.org/10.1093/imrn/rnae061),
  IMRN 2024, equations (2.6)--(2.7), equation (4.1), and Lemma 4.4;
- Daniel Huybrechts,
  [*Brilliant families of K3 surfaces: Twistor spaces, Brauer groups, and
  Noether--Lefschetz loci*](https://doi.org/10.5802/afst.1741), Section 1.2
  and Corollary 3.3, for the distinction between algebraic/projective and
  analytic/nonprojective Tate--Shafarevich families;
- Daniel Huybrechts,
  [*Lectures on K3 Surfaces*](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf),
  Chapter 11, Section 4.1, equation (4.1), for fibrewise and local-etale
  equality with the relative Jacobians, and Chapter 14, Theorem 1.12, for
  Nikulin's unique primitive-embedding criterion;
- V. V. Nikulin,
  [*Integral symmetric bilinear forms and some of their geometric
  applications*](https://www.mathnet.ru/eng/im1677), Theorem 1.14.4, for the
  original stable primitive-embedding theorem.

The exact replay is
[`certify_noncompact_k3_fixed_jacobian_index_blindness.py`](../scripts/certify_noncompact_k3_fixed_jacobian_index_blindness.py),
with frozen output
[`noncompact_k3_fixed_jacobian_index_blindness.json`](../data/enumerations/noncompact_k3_fixed_jacobian_index_blindness.json)
and regressions in
[`test_noncompact_k3_fixed_jacobian_index_blindness.py`](../tests/test_noncompact_k3_fixed_jacobian_index_blindness.py).
It verifies the primitive unit minor, the exact orthogonal complement, both
affine \(E_8\) relations, the vector/Gram/common-orthogonal ranks, the
discriminant and divisibility formulas, and the transverse lattice root.  It
does not pretend to machine-prove Ogg--Shafarevich theory or a finite-volume
claim.

Run

```bash
python3 scripts/certify_noncompact_k3_fixed_jacobian_index_blindness.py \
  --check data/enumerations/noncompact_k3_fixed_jacobian_index_blindness.json
pytest -q tests/test_noncompact_k3_fixed_jacobian_index_blindness.py
```
