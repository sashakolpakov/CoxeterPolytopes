# Independent hostile audit: fixed-Jacobian K3 torsors do not bound cusp index

**Verdict (2026-08-11): PASS, exactly in the asserted cusp-local and
selected-fibration scope.** One fixed Jacobian elliptic K3 surface has
projective genus-one K3 torsors of every index \(d\geq2\), with

\[
 T(X_d)\simeq U(d)\oplus U,
 \qquad
 \operatorname {NS}(X_d)\simeq U(d)\oplus E_8(-1)^2,
 \qquad
 \operatorname {div}(F_d)=d.
\]

The two \(II^*\) fibres and four \(I_1\) fibres remain unchanged. Their two
physical affine \(E_8\) packets have finite-root rank
\(16=\rho(X_d)-2\) and common orthogonal line \(\mathbf QF_d\). Thus all of
the selected-cusp data in the owner theorem are genuinely fixed while the
integral transverse glue \(d\) escapes.

This is not a finite-volume family and not a family of maximal reflection
chambers. The separate \(B=13\) spanning theorem is fully compatible with
it. Under an additional complete finite-volume hypothesis its automatic
consequence is

\[
                         d^2\leq2877^9,
\]

not \(d\leq390\). The latter needs a stronger cusp-adapted assertion putting
one bounded transverse wall together with the nine preselected \(II^*\)
components.

## 1. Frozen bytes audited

The independent audit applies to these exact owner artifacts:

```text
14520f6cbbe6dc95f559ebb36d1232dd572adaf54e55457ed434f6108f39852b  results/noncompact_k3_fixed_jacobian_index_blindness.md
50ca142bd2f3d3d8ad7fb469ab31a56e9d5fcfdfde60999c1191f78fb75e76ef  results/noncompact_k3_fixed_jacobian_index_blindness_hostile_scope_audit.md
bb72c96a3aed2b02289856167938a643e06938f1fc6345aa9a3674bb6c2bca6b  scripts/certify_noncompact_k3_fixed_jacobian_index_blindness.py
1f7bf7c64ef44043e6df6700558c8a610a30b01a8f37a3e1c3f31c49bdcfcc47  data/enumerations/noncompact_k3_fixed_jacobian_index_blindness.json
e736f24793b87e1865aa3c8638159209a3a3082a3fe7a0c97cb196fb00ad3408  tests/test_noncompact_k3_fixed_jacobian_index_blindness.py
```

The corrected global input is frozen separately:

```text
de3b78a9070dc7035a8ecf51c5aa697b0736ad1559f0db8867b6629eb3cc4446  results/noncompact_k3_bounded_spanning_root_lemma.md
70ae2c1e09ddd0fe8cd1999087172ed14499dbd47709ebbc2236aa77ad985bde  results/noncompact_k3_bounded_spanning_root_lemma_independent_hostile_audit.md
```

The proof below reconstructs the mathematical argument. The owner certifier
and archive are regression evidence, not evidence for the imported
algebraic-geometric theorems.

## 2. The fixed Jacobian is a genuine two-\(II^*\) K3

Consider

\[
 y^2=x^3+a t^4x+t^5(b_0+b_1t+b_2t^2).                 \tag{2.1}
\]

At \(t=0\), the orders of \(A,B,\Delta\) are \((4,5,10)\). After setting
\(s=1/t\) and applying the standard K3 Weierstrass rescaling, one gets

\[
 A_\infty=a s^4,
 \qquad
 B_\infty=s^5(b_2+b_1s+b_0s^2),
\]

so the same orders occur at infinity. The remaining discriminant factor is

\[
             4a^3t^2+27(b_0+b_1t+b_2t^2)^2,           \tag{2.2}
\]

a quartic with four simple nonzero roots for a general choice. Hence the
fibre list is exactly

\[
                             2II^*+4I_1,
\]

and its Euler sum is \(10+10+4=24\).

The section, fibre, and two finite \(E_8\) component lattices give
\(U\oplus E_8(-1)^2\) in the Neron--Severi group. The corresponding
lattice-polarized K3 moduli space has dimension two, and \(2.1\) is the usual
two-parameter family after Weierstrass and base rescaling. Away from its
Noether--Lefschetz loci,

\[
 \operatorname {NS}(S)=U\oplus E_8(-1)^2,
 \qquad T(S)=U^2.                                    \tag{2.3}
\]

This is also the standard \(U\oplus E_8\oplus E_8\)-polarized family studied
by [Clingher--Doran](https://arxiv.org/abs/math/0602146). Shioda--Tate gives
Mordell--Weil rank zero. There is no torsion: a \(II^*\) fibre has trivial
component group, and its identity smooth part is additive and has no
nonzero complex torsion.

## 3. Every order occurs projectively

For a complex projective K3 surface, the canonical description is

\[
 \operatorname {Br}(S)
 \simeq\operatorname {Hom}(T(S),\mathbf Q/\mathbf Z).
                                                               \tag{3.1}
\]

Write \(T(S)=U_1\oplus U_2\), with first basis \(e_1,f_1\), and define

\[
 \alpha_d(ae_1+bf_1+ce_2+ef_2)=a/d\pmod {\mathbf Z}.            \tag{3.2}
\]

This homomorphism has exact order \(d\). The important existence statement
is not an informal appeal to the analytic Tate--Shafarevich family.
[Meinsma--Shinder, Section 4](https://eprints.whiterose.ac.uk/id/eprint/233188/1/rnae061.pdf)
works in the complex projective elliptic-K3 category, identifies
\(\operatorname {Br}(S)\) with the Ogg--Shafarevich torsor group, and proves
surjectivity in Theorem 4.13. Thus every \(\alpha_d\) is represented by a
projective genus-one K3 surface \(X_d\to\mathbf P^1\) with relative
Jacobian \(S\).

Lemma 4.4 of the same paper states precisely

\[
 \operatorname {ind}(X_d)=\operatorname {ord}(\alpha_d)=d,
 \qquad
 T(X_d)\simeq\ker(\alpha_d:T(S)\to\mathbf Z/d).        \tag{3.3}
\]

The kernel has basis \(de_1,f_1,e_2,f_2\), so its Gram matrix is

\[
                           U(d)\oplus U.                \tag{3.4}
\]

This closes both possible loopholes: the torsor exists as a projective K3,
and its index is exactly \(d\), rather than merely divisible by the period.
Huybrechts's discussion of the algebraic and analytic families gives the
same distinction; see
[*Brilliant families of K3 surfaces*](https://doi.org/10.5802/afst.1741),
Corollary 3.3.

## 4. The primitive K3 complement is exact

The lattice \(T_d=U(d)\oplus U\) has signature \((2,2)\) and discriminant
group \((\mathbf Z/d)^2\), so \(\ell(A_{T_d})=2\). Nikulin's stable
primitive-embedding theorem applies because

\[
 2<3,
 \qquad 2<19,
 \qquad \ell(A_{T_d})+2=4\leq22-4=18.                 \tag{4.1}
\]

It gives one primitive-embedding orbit in the K3 lattice. This is
[Theorem 14.1.12 in Huybrechts's lectures](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf),
equivalently
[Nikulin, Theorem 1.14.4](https://www.mathnet.ru/eng/im1677).

There is also a direct representative. In

\[
 \Lambda_{K3}=U_1\oplus U_2\oplus U_3\oplus E_8(-1)^2,
\]

embed \(T_d\) by

\[
                    -e_1+df_2,\quad e_2,\quad e_3,\quad f_3.    \tag{4.2}
\]

The rows \(e_1,e_2,e_3,f_3\) give a maximal minor of absolute value one,
so the embedding is primitive. Direct orthogonality gives the complement

\[
 \mathbf Ze_1+\mathbf Z(df_1+e_2)+E_8(-1)^2
 \simeq U(d)\oplus E_8(-1)^2.                           \tag{4.3}
\]

This displayed complement is primitive as well. Uniqueness of the orbit
therefore identifies the actual complement of \(T(X_d)\), proving

\[
 \operatorname {NS}(X_d)\simeq U(d)\oplus E_8(-1)^2.   \tag{4.4}
\]

In particular,

\[
 |\operatorname {disc}\operatorname {NS}(X_d)|=d^2,
 \qquad
 A_{\operatorname {NS}(X_d)}\simeq(\mathbf Z/d)^2.      \tag{4.5}
\]

Thus bounded discriminant length is genuinely blind to \(d\), while the
full discriminant group and form are not.

## 5. The actual fibre marking survives the complement calculation

An abstract Picard isometry would not by itself identify the displayed
\(U(d)\) vector with the geometric fibre. Let \(R=E_8(-1)^2\) be generated
by eight finite components from each actual \(II^*\) fibre. Since \(R\) is
unimodular, its embedding splits integrally:

\[
             \operatorname {NS}(X_d)=R\oplus K.         \tag{5.1}
\]

The primitive fibre \(F_d\in K\) has divisibility \(d\). Completing it to a
basis \(F_d,D\) gives

\[
 \operatorname {Gram}(K)=
 \begin{pmatrix}0&d\\ d&2a\end{pmatrix}.               \tag{5.2}
\]

The Smith group is \((\mathbf Z/d)^2\), so \(d\mid2a\). Replacing \(D\) by
\(D+kF_d\) changes \(a\) by \(kd\). For odd \(d\), this already gives
\(a\equiv0\pmod d\). For even \(d\), the sole competitor is
\(a\equiv d/2\pmod d\), with Gram

\[
                         \begin{pmatrix}0&d\\d&d\end{pmatrix}.
\]

In its dual basis \(P=(D-F_d)/d\), \(Q=F_d/d\), the discriminant form is

\[
             q(pP+qQ)=\frac{-p^2+2pq}{d}\pmod {2\mathbf Z}.     \tag{5.3}
\]

Isotropy forces \(p\) even. Hence all isotropic elements lie in a proper
index-two subgroup and cannot generate \((\mathbf Z/d)^2\). By contrast,
\(q_{U(d)}\) has the two standard isotropic generators. The exceptional
binary lattice therefore has the wrong discriminant form, and

\[
 (\operatorname {NS}(X_d),F_d,R)
 \simeq(U(d)\oplus E_8(-1)^2,u,E_8(-1)^2).              \tag{5.4}
\]

The independent executable checks this parity invariant for every even
\(2\leq d\leq50\), in addition to the explicit primitive complements. No
determinant-only inference remains.

## 6. The fibre packet is physical but only locally complete

Huybrechts proves in
[*Lectures on K3 Surfaces*, Chapter 11, Section 4.1](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf)
that

\[
                         J^k(X)_t\simeq X_t             \tag{6.1}
\]

for every integer \(k\) and every closed base point. For singular fibres,
the proof constructs a local analytic or etale section through a smooth
point and uses the absence of multiple fibres on a K3 surface. Taking
\(k=0\) proves that each \(X_d\) has the same fibre schemes as \(S\), not
merely the same Euler number.

Each \(II^*\) fibre consequently supplies nine actual irreducible smooth
rational square-\((-2)\) curves. With the affine coefficient followed by
the highest-root marks,

\[
                 (1,2,3,4,6,5,4,3,2),                 \tag{6.2}
\]

their weighted sum is \(F_d\). The two packet vector spaces have dimension
nine and intersect in \(\mathbf QF_d\), so their union has vector-span rank
17 and Gram rank 16. In the rank-18 Picard space,

\[
 \left\langle\text{the two affine }E_8\text{ packets}\right\rangle_{\mathbf Q}
 =F_d^\perp,
 \qquad
 \left(\cdots\right)^\perp=\mathbf QF_d.               \tag{6.3}
\]

This is the strongest possible local parabolic rank, but it is still one
vector short of spanning the Picard space. It is not a complete
finite-volume wall system.

## 7. Divisibility is the escaping transverse parameter

In the marked \(U(d)\) plane, the pairing row of \(F_d=u\) is
\((0,d)\). Therefore

\[
 \operatorname {div}(F_d)=d,
 \qquad
 F_d.D\in d\mathbf Z
 \quad\text{for every }D\in\operatorname {NS}(X_d).     \tag{7.1}
\]

This agrees with the multisection index from (3.3). Any horizontal physical
curve, if present, has degree at least \(d\). Thus no bound depending only
on the fixed Jacobian, fibre packet, Euler number, Shioda--Tate rank,
Picard rank, primitive K3 embedding, or discriminant length can guarantee a
bounded-degree transverse physical root.

There is a formal lattice root of degree \(d\): if \(v\) is the second
\(U(d)\) vector and \(\theta\) is the first \(E_8\) highest root, then

\[
                         c=v-\theta,
 \qquad c^2=-2,
 \qquad F_d.c=d.                                      \tag{7.2}
\]

It pairs nonnegatively with the displayed affine packet. The owner theorem
correctly does not call it an irreducible curve or a nef-cone facet.

The selected-fibration stabilizer is finite, but this also remains local.
Its translation part is the rank-zero, torsion-free \(MW(S)\); the remaining
generic-fibre kernel is finite; and its base image preserves a six-point
discriminant set. This gives no information about the full automorphism
group, whose elements may move \(F_d\).

## 8. Exact compatibility with the \(B=13\) global theorem

Suppose, only conditionally, that one \(X_d\) had a finite-sided,
finite-volume, complete full-Picard root chamber. The global theorem then
selects 18 physical roots with Gram matrix \(G\), diagonal entries \(-2\),
and off-diagonal entries in \(\{0,1,\ldots,13\}\). Every Euclidean row norm
satisfies

\[
 \|G_i\|_2^2\leq4+17\cdot13^2=2877.
\]

If \(R\) is their full-rank lattice, then

\[
 |\det G|
 =[\operatorname {NS}(X_d):R]^2
  |\operatorname {disc}\operatorname {NS}(X_d)|
 =[\operatorname {NS}(X_d):R]^2d^2.                   \tag{8.1}
\]

Hadamard's inequality gives

\[
             d^2\leq|\det G|\leq2877^{18/2}=2877^9.   \tag{8.2}
\]

The exponent and constant are therefore correct.

There is a sharper statement under a stronger, cusp-adapted hypothesis. If
one actual transverse curve \(C\) has \(C.C_i\leq13\) for all nine
components of one selected \(II^*\) fibre, then the marks in (6.2) sum to
30 and

\[
                    d\leq F_d.C\leq30\cdot13=390.       \tag{8.3}
\]

But the \(B=13\) theorem selects an arbitrary spanning packet in a narrow
closed star. It does not place the nine preselected fibre components in
that packet. Hence (8.3) is not an automatic consequence. Numerically the
distinction is visible already because \(391^2<2877^9\).

The fixed-Jacobian family therefore calibrates, rather than contradicts,
the global theorem: local parabolic completeness misses \(d\); one global
spanning packet bounds \(d\) coarsely; a packet adapted to the selected cusp
would give the sharper value 390.

## 9. Final scope table

| assertion | independent verdict |
|---|---|
| projective torsor for every \(d\geq2\) | PASS; Ogg--Shafarevich surjectivity in the projective elliptic-K3 category |
| index and fibre divisibility equal \(d\) | PASS; period--index theorem plus the marked \(U(d)\) plane |
| \(T(X_d)=U(d)\oplus U\) | PASS; exact Brauer kernel |
| \(\operatorname {NS}(X_d)=U(d)\oplus E_8(-1)^2\) | PASS; unique primitive embedding and explicit complement |
| actual fibre-aligned \(U(d)\oplus E_8(-1)^2\) splitting | PASS; unimodularity plus the even-\(d\) discriminant-form exclusion |
| unchanged physical \(2\widetilde E_8\) packet | PASS; fibrewise relative-Jacobian theorem |
| local invariants bound \(d\) | false; this family disproves the implication |
| the family has finite-volume complete chambers | not claimed or proved |
| the family is reflection-maximal | not claimed or proved |
| \(B=13\) automatically gives \(d\leq390\) | false; cusp adaptation is missing |
| finite volume conditionally gives \(d^2\leq2877^9\) | PASS; exact index identity and Hadamard |
| universal Property 4 | not claimed |

## 10. Independent reproduction

The independent replay is
[audit_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py](../scripts/audit_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py),
with regressions in
[test_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py](../tests/test_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py).
It imports no owner code. It checks the two \(II^*\) Weierstrass orders, both
primitive unit minors, exact orthogonal complement, discriminants, fibre
divisibility, both affine \(E_8\) relations and ranks, the even-\(d\)
discriminant-form obstruction, the owner archive self-hash and negative
scope flags, and the corrected \(2877^9\) Hadamard constant.

Run

```bash
python3 scripts/audit_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py --check-pins
pytest -q tests/test_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py
```
