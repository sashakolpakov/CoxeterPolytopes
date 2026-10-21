# Independent hostile audit: compact physical stars do finite-ize the trace host

**Verdict (2026-08-11): PASS, with the stated literal same-vector and
same-host scopes essential.**  I audited the frozen owner theorem without
using its certifier as mathematical evidence.  A spanning compact
minimal-aperture star represented by primitive square-\((-2)\) vectors in
the trace module does eliminate all five proposed integral escapes:

\[
 \text{common trace scale},\quad
 \text{order conductor},\quad
 \text{root-line ideals},\quad
 \text{same-span overlattice index},\quad
 \text{trace-host isometry type}.
\]

It does **not** bound an orthogonal Picard extension or a K3 complement, and
a physical curve whose wall merely projects to the same eigenspace line is
not enough.  The independent
[projected-root scaling counterfamily](k3_compact_projected_root_scaling_no_go.md)
shows that this distinction is sharp.

## 1. Frozen owner bytes audited

The audit applies to exactly these files:

```text
5cc86e7225c11697c83af0ec34fa5c773e8e969c47107eabdd2a327c23c1ea3a  results/k3_compact_physical_star_trace_host_finiteness.md
c4a52c76f28192eada049fe66031fd62bb35f6366332f92abc54cec995d58d55  scripts/certify_k3_compact_physical_star_trace_host_finiteness.py
cc2d48b699aeebd1aaa10ae6b4b8249addde2eac681c9bb475aa19fcd3e14e90  data/enumerations/k3_compact_physical_star_trace_host_finiteness.json
7f38e30ed6fb20af8c54a266c9b8ec4cad59ddccd90f1fd451d40f9565dbd3ec  tests/test_k3_compact_physical_star_trace_host_finiteness.py
```

The proof below rederives the implication directly.  The owner certificate
is useful regression evidence, not a substitute for the argument.

## 2. Energy audit: square \(-2\) locks the trace kernel

Let \(g\) be the \(k\)-valued trace kernel and let
\(\delta_1,\ldots,\delta_m\) be the selected physical star basis.  The key
hypothesis is literal:

\[
             \delta_i\in M,\qquad
             B(\delta_i,\delta_i)=-2,                 \tag{2.1}
\]

and the same \(\delta_i\) is the native field-normal vector for its facet.
At the distinguished place it is a negative wall normal; all other places
are negative definite.  Therefore

\[
 t_{i,\tau}=-\tau(g(\delta_i,\delta_i))>0,
 \qquad \sum_\tau t_{i,\tau}=2.                       \tag{2.2}
\]

This is the missing scale lock.  It bounds every diagonal contribution by
two.  For a pair of star roots, the narrow-facet estimate gives

\[
 |\tau_0(g_{ij})|<7\sqrt{t_{i,\tau_0}t_{j,\tau_0}}<14, \tag{2.3}
\]

while definite-place Cauchy--Schwarz gives

\[
 |\tau(g_{ij})|\le\sqrt{t_{i,\tau}t_{j,\tau}}\le2.
                                                               \tag{2.4}
\]

The independent trace check is also correct:

\[
 |B(\delta_i,\delta_j)|
 <14+
 \sqrt{(2-t_{i,\tau_0})(2-t_{j,\tau_0})}
 \le16.                                                     \tag{2.5}
\]

This last integer bound alone would not control the field form; (2.2)--
(2.4), combined with the codifferent lattice, do.

## 3. Conductor audit: the cycles bound the order from below

For \(T_i=1-s_{\delta_i}\), the projective-module trace of a product of the
rank-one operators is the corresponding directed Cartan cycle.  It lies in
the declared coefficient order \(\mathcal O\).  The compact Lorentz-block
argument supplies a finite field-generating packet \(\mathcal C\) of such
integral cycles, with degree bounded by the K3 rank and every conjugate
bounded by \(14^m\).

There are consequently only finitely many possible packets and fields.  For
each one,

\[
          \mathcal O_0=\mathbf Z[\mathcal C]
          \subseteq\mathcal O\subseteq\mathcal O_k.   \tag{3.1}
\]

Because \(\mathcal C\) generates \(k\), \(\mathcal O_0\) is a full order.
The finite additive group \(\mathcal O_k/\mathcal O_0\) has only finitely
many subgroups, so it has only finitely many intermediate subrings.  This
rules out an unbounded order conductor.  Notice that this part uses exact
integral reflection-module traces; a field bound by itself would not suffice.

## 4. Codifferent audit: only finitely many field Grams

Self-adjointness and integrality give, for all \(a\in\mathcal O\),

\[
 B(a\delta_i,\delta_j)=\operatorname {Tr}_{k/\mathbf Q}
       (a g_{ij})\in\mathbf Z.
\]

Hence

\[
                          g_{ij}\in\mathcal O^\vee.   \tag{4.1}
\]

For each of the finitely many orders, \(\mathcal O^\vee\) is a discrete
full Minkowski lattice.  Equations (2.2)--(2.4) place \(g_{ij}\) in a fixed
bounded box, whose intersection with that lattice is finite.  Thus the
entire \(k\)-Gram matrix on the star basis belongs to a finite set.  In
particular, rewriting the form as \(g=-\alpha h\) cannot reintroduce a
common scale: the reflection host depends on \(g\), and \(g\) has already
been bounded.

## 5. Root-line and overlattice audit

Set

\[
                         F=\bigoplus_i\mathcal O\delta_i. \tag{5.1}
\]

The physical star spans \(V\), so the chosen \(m\) vectors are
\(k\)-independent and \(F\) has full additive rank.  The possible order,
field Gram and an integral basis of \(\mathcal O\) determine the full
\(\mathbf Z\)-Gram of \(F\), including all root-line ideal information.

Every integral host containing this frame satisfies

\[
                          F\subseteq M\subseteq F^\#. \tag{5.2}
\]

The discriminant group \(F^\#/F\) is finite, so only finitely many
intermediate additive lattices occur.  Moreover

\[
                  [M:F]^2\mid |\det F|.               \tag{5.3}
\]

Projectivity over \(\mathcal O\), evenness, reflection stability and K3
embeddability are filters on this finite list.  This proves that neither a
fractional root-line ideal nor an integral same-span saturation can escape.

## 6. Physical-root location and the full-Picard corollary

The owner theorem does not infer \(\delta_i\in M\) from equality of a real
wall line; it assumes it.  That is logically correct.  If
\(\operatorname {NS}(X)=M_{\mathbf Z}\), however, the condition becomes
automatic.  An ambient K3 wall root \(\delta\) then lies in \(M\).  If its
distinguished projection is proportional to a native normal \(r\), choose
a nonzero coordinate \(\delta_\ell\).  The ratio
\(a=r_\ell/\delta_\ell\in k\) has distinguished value equal to the real
proportionality scalar, and injectivity of the real field embedding gives

\[
                              r=a\delta.               \tag{6.1}
\]

Native reflections are unchanged by this rescaling.  Thus exact full-Picard
eigen-nef return with actual K3 facet walls supplies the physical-star
hypothesis.

If the trace block is only primitive in a larger Picard lattice, the root
may instead have the form

\[
                  \delta=t+u,qquad
                  0\ne u\in T_{\mathbf Q}^{\perp},     \tag{6.2}
\]

possibly through discriminant gluing.  Its restricted wall sees only \(t\).
There is then no identity
\(\operatorname {Tr}g(t,t)=-2\), and the energy proof cannot start.

The projected-root counterfamily makes (6.2) explicit for the fixed compact
golden \((2,5,5)\) triangle.  It has trace form \(sQ\), primitive trace
embedding, physical square-\((-2)\) facets and exact eigen-nef return for
every \(s\ge2\), while

\[
 [S_s:T_s\oplus K_s]=s^3,
 \qquad |\det S_s|=64s^3(s-1)^3.                       \tag{6.3}
\]

So physical projection alone genuinely fails; this is not merely a gap in
the proof.

## 7. Chamber-finiteness composition

The owner corollary is valid only for maximality in the same complete native
\(k\)-reflection host.  Once the finite host list is known, any additional
native integral reflection lies in the same discrete arithmetic group.
Same-host reflection maximality therefore forces the displayed facet group
to equal the group generated by all native host reflections.  Local
finiteness and reflection adjacency give one orbit of chambers for that
complete arrangement, so each oriented host contributes at most one
unmarked normalized chamber class.

For \([k:\mathbf Q]>1\), an ordinary K3 root reflection has rational minus
eigenspace of dimension one, whereas a native \(k\)-reflection has rational
minus eigenspace of dimension \([k:\mathbf Q]\).  Maximality for the first
species cannot silently replace maximality for the second.

## 8. Final scope table

| proposed conclusion | audit verdict |
|---|---|
| finite coefficient fields | PASS; inherited from compact narrow-star cycles |
| finite coefficient orders | PASS; \(\mathcal O_0\subseteq\mathcal O\subseteq\mathcal O_k\) |
| finite trace kernels/common scales | PASS; square-\((-2)\) energy plus codifferent box |
| finite projective root modules and same-span overlattices | PASS; \(F\subseteq M\subseteq F^\#\) |
| finite native chambers | PASS only with exact return and same-host \(k\)-reflection maximality |
| finite ambient NS extensions or K3 complements | **not claimed and false without another bound** |
| same result from merely projected physical roots | **false**, by the golden scaling/glue family |
| universal Property 4 | **not claimed** |

The theorem is therefore sound and sharp in the direction relevant to the
compact finite-host implication: the nontrivial AG property needed is not
merely a physical curve above each wall, but a spanning minimal star of
literal physical \((-2)\)-vectors inside the full trace module.

## 9. Audit reproduction

The independent regression test pins the owner bytes, replays the energy,
order, codifferent and discriminant-index mechanisms, and cross-checks the
projected-root no-go:

```text
pytest -q tests/test_k3_compact_physical_star_trace_host_finiteness_independent_hostile_audit.py
```
