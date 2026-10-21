# A finite P-first semiampleness certificate for weak-log exact return

**Status (2026-08-10): proved finite local-algebra certificate theorem and
nonzero ruled calibration; the PB-Atlas B1/B2a handoff below is an exact
adapter contract, not an implemented end-to-end adapter.**  Let

\[
                  D=-\left(K_S+\frac12B\right)
\]

on a smooth rational cluster surface.  This note gives a finite certificate
that a specified multiple \(mD\) is base-point-free, hence that \(D\) is
nef.  It supplies exactly the missing nef premise in the
[reflection-maximal weak-log exact-return theorem](reflection_maximal_weak_log_exact_nef_return.md).

The certificate does **not** infer nefness from the proposed Coxeter
halfspaces.  It uses explicit affine surface charts, local trivializations of
\(\mathcal O_S(mD)\), returned global sections, and exact localized
Bezout/saturation identities.  Incidence nonemptiness, boundary geometry,
facet effectivity, reflection maximality, and raw surface-free exhaustivity
remain separate hypotheses.

The machine artifacts are:

- [schema](../schemas/pfirst-weak-log-semiampleness.schema.json);
- [localized-unit/F4 checker](../scripts/certify_pfirst_weak_log_semiampleness.py);
- [frozen archive](../data/enumerations/pfirst_weak_log_semiampleness.json); and
- [tests](../tests/test_pfirst_weak_log_semiampleness.py).

## 1. The local algebra certificate

Let \(k\) be an exact ground field and let \(L\) be a line bundle on a
projective surface \(S\).  Suppose a certified finite affine cover has rings

\[
 A_\alpha=
 \left(k[z_1,\ldots,z_N]/(I_{\alpha,1},\ldots,I_{\alpha,r_\alpha})\right)
 [g_{\alpha,1}^{-1},\ldots,g_{\alpha,c_\alpha}^{-1}].       \tag{1.1}
\]

This presentation also covers a parameterized PB-Atlas incidence component:
the parameter equations are included among the \(I_{\alpha,\ell}\), and all
declared open factors are among the \(g_{\alpha,h}\).  Put

\[
                         g_\alpha=\prod_h g_{\alpha,h}.
\]

Fix a local frame \(e_\alpha\) of \(L\) on every chart.  Let
\(s_1,\ldots,s_t\) be returned global sections and write

\[
                  s_j=f_{\alpha j}e_\alpha,
       \qquad f_{\alpha j}={n_{\alpha j}\over
                  \prod_hg_{\alpha,h}^{a_{\alpha jh}}}
                  \in A_\alpha .                            \tag{1.2}
\]

All denominators in (1.2) are units.  The finite certificate on chart
\(U_\alpha\) is an integer \(N_\alpha\ge0\) and polynomials
\(q_{\alpha j},r_{\alpha\ell}\) satisfying the literal identity

\[
 g_\alpha^{N_\alpha}
  =\sum_jq_{\alpha j}n_{\alpha j}
   +\sum_\ell r_{\alpha\ell}I_{\alpha,\ell}.                \tag{1.3}
\]

Equivalently,

\[
 1\in(n_{\alpha1},\ldots,n_{\alpha t})A_\alpha,
 \quad\text{or}\quad
 ((I_\alpha,n_{\alpha1},\ldots,n_{\alpha t}):g_\alpha^\infty)
       =k[z_1,\ldots,z_N].                                  \tag{1.4}
\]

The archive includes the nonconstant control

\[
                      x=y+(x-y)                              \tag{1.5}
\]

in \(\mathbf Q[x,y]_x\).  Thus the checker exercises genuine localization:
the two returned numerators are \(y\) and \(x-y\), neither a constant unit
before localization.

### Theorem 1.1 (finite semiampleness certificate)

Assume:

1. the charts (1.1) are certified to cover \(S\);
2. the frames \(e_\alpha\) and transition units certify one line bundle
   \(L\);
3. the local expressions (1.2) agree under those transitions and hence
   define global sections; and
4. identity (1.3) holds on every chart.

Then the returned sections generate \(L\), so \(L\) is base-point-free.  If
\(L=\mathcal O_S(mD)\) for a positive integer \(m\), then \(D\) is nef.

For fixed returned sections and a fixed finite cover, the existence of the
identities (1.3) is also necessary for base-point-freeness.

### Proof

After inverting \(g_\alpha\), equation (1.3) gives \(1\) in the ideal of the
local section functions.  The evaluation map

\[
       \langle s_1,\ldots,s_t\rangle_k\otimes\mathcal O_{U_\alpha}
                    \longrightarrow L|_{U_\alpha}           \tag{1.6}
\]

is therefore surjective.  The affine charts cover \(S\), so (1.6) is
surjective globally.  The resulting morphism
\(\phi:S\to\mathbf P^{t-1}\) satisfies
\(L\simeq\phi^*\mathcal O(1)\).  Hence \(L.C\ge0\) for every irreducible
curve \(C\), and

\[
                            D.C={L.C\over m}\ge0.             \tag{1.7}
\]

Conversely, if these fixed sections generate on \(U_\alpha\), their
numerators generate the unit ideal after localization.  Clear the finitely
many denominators in a unit expression to obtain (1.3) for some
\(N_\alpha\).  This is the elementary localization lemma; the
Nullstellensatz may be used to discover the empty base locus, but replay
requires the stronger exact identity.  \(\square\)

No algebraic closure is needed to check (1.3).  If emptiness is first proved
geometrically after extending the field, faithful flatness descends the unit
ideal; the serialized certificate is still an identity over the declared
exact coefficient field.

## 2. Certifying the affine cover and the global sections

Merely listing affine charts is not a cover proof.  For a PB-Atlas blowup
packet the required downstream adapter contract is inductive.  The frozen
checker in this note does not yet construct that cover from a live B1/B2a
packet.

### 2.1 Cover contract

Start with the three standard affine charts of \(\mathbf P^2\), or the
standard toric charts of \(\mathbb F_e\).  At a chronological B1 centre with
local ideal \((u,v)\), replace the center chart by the two Rees charts

\[
 A[t]/(v-ut),\qquad A[s]/(u-vs).                             \tag{2.1}
\]

They are the two standard opens of
\(\operatorname {Proj}A[U,V]/(vU-uV)\), and cover its exceptional
\(\mathbf P^1\).  Retain the complement of the centre.  Iterating (2.1)
along the B1 chronological centre records gives a finite affine cover of the
final blowup, including free and satellite centres.  The selected B2a chart
determines the local substitutions used by the section engine, but it may not
stand in for the other blowup chart: both opens, or an exact refinement of
them, belong to the cover certificate.

### 2.2 Section contract

Let

\[
         mD=\pi^*L_0-\sum_i\mu_iE_i^{\rm total}.             \tag{2.2}
\]

A future B1/B2a adapter must perform five exact checks.

1. The integer \(m\) clears every coefficient denominator of \(D\), and the
   B1 class row agrees coefficientwise with (2.2).
2. The endpoint coefficient vector of every returned section uses the B1
   standard section columns: plane monomials for \(\mathbf P^2\), or the
   summands

   \[
       H^0(\mathbb F_e,aS+bF)
          =\bigoplus_{i=0}^aH^0(\mathbf P^1,\mathcal O(b-ie))
                                                                  \tag{2.3}
   \]

   on a ruled endpoint.
3. The exact B2a jet matrix, or a ruled analogue with the same row contract,
   kills each coefficient vector.  Expected dimension is never substituted
   for this kernel calculation.
4. At centre \(p_i\), the strict transform is divided by the exceptional
   coordinate to exactly the B1 multiplicity \(\mu_i\).  The B2a low-jet
   deletion guarantees that the quotient is regular, not Laurent.
5. The resulting local functions satisfy the transition identities on all
   overlaps.  These identities either are reconstructed from the single
   endpoint section, or are returned as further cleared-denominator
   polynomial identities.

The current Stage-B2a engine implements this section/jet transport only for
plane packets.  The schema records that limitation rather than feeding a
ruled packet through plane formulas.  A ruled chart provider must implement
(2.3).  The calibration below supplies such a provider torically for
\(\operatorname {Bl}_p\mathbb F_4\).

The present checker implements the sparse localized identity verifier and
the complete toric \(\mathbb F_4\) provider.  It pins the B1/B2a interfaces,
but it does **not** read an arbitrary B1 packet, call B2a, solve its jet
kernels, construct every Rees chart, or verify all transition identities.
Those operations are the specified adapter debt; metadata saying that they
were performed is not accepted as their proof.

## 3. Exact serialized schema and checker

The schema identifier is

```text
coxeter.pb_atlas.pfirst_weak_log_semiampleness.v1
```

The reusable per-carrier contract specified here must pin its B1/B2a packet
or alternate provider and contain:

```text
surface_cover
  endpoint affine charts
  chronological blowup/refinement records
  chart rings Q[parameters,coordinates]/I localized at g

line_bundle
  denominator m
  class mD in B1 chronological coordinates
  local frames and transition units

global_sections
  endpoint coefficient vectors
  exact jet-kernel witnesses
  local numerator and denominator records

unit_certificates[chart]
  sparse chart relations I_l
  sparse localizers g_h
  sparse section numerators n_j
  saturation exponent N
  relation multipliers r_l
  section multipliers q_j
```

Sparse coefficients are rational strings and exponent arrays; booleans,
floating approximations, malformed denominator powers, and exponent-length
mismatches are rejected.  The checker multiplies sparse polynomials over
\(\mathbf Q\) with `Fraction`, verifies (1.3) coefficientwise, validates the
JSON schema, and pins the B1/B2a and weak-log exact-return dependencies.
The frozen archive instantiates the full surface-cover part of this contract
only in the toric calibration of Section 4.  Its `pb_atlas_adapter` block is
a versioned contract specification, not a claim that arbitrary B1/B2a
packets were consumed.

For a universal incidence component, its ideal belongs among the chart
relations and every collision/nonvanishing factor belongs to the localizer.
Thus one identity certifies the whole declared constructible open, not a
finite sample of points.  Nonemptiness of that open is still a separate
incidence certificate.

## 4. Nontrivial calibration on \(\operatorname {Bl}_p\mathbb F_4\)

Use the smooth complete fan of \(\mathbb F_4\) with cyclic rays

\[
 (0,1),\quad(-1,0),\quad(0,-1),\quad(1,4).                  \tag{4.1}
\]

The first ray is the negative section \(B\).  Blow up the torus fixed point
where the positive section meets the second fibre.  Star subdivision inserts

\[
                       (-1,-1)=(-1,0)+(0,-1).                \tag{4.2}
\]

In cyclic order label the new rays

\[
\begin{array}{c|ccccc}
\text{ray}&B&G&E&P&F\\ \hline
v_\rho&(0,1)&(-1,0)&(-1,-1)&(0,-1)&(1,4)\\
D_\rho^2&-4&-1&-1&3&0.
\end{array}                                                  \tag{4.3}
\]

Here \(E\) is exceptional and \(G=F-E\).  In the total-transform basis
\((B,F,E)\),

\[
 (\ .\ )=
 \begin{pmatrix}-4&1&0\\1&0&0\\0&0&-1\end{pmatrix},
 \qquad K=(-2,-6,1).
\tag{4.4}
\]

For boundary \(B\),

\[
 D=-\left(K+{B\over2}\right)
   ={3\over2}B+6F-E,
 \qquad L=2D=3B+12F-2E.                                    \tag{4.5}
\]

The invariant-divisor coefficients of \(L\), in the order (4.3), are

\[
                           (3,0,-2,0,12).                    \tag{4.6}
\]

Thus its character polytope consists of \((x,y)\in\mathbf Z^2\) satisfying

\[
 y\ge-3,\quad -x\ge0,\quad -x-y\ge2,\quad -y\ge0,\quad
 x+4y\ge-12.                                                \tag{4.7}
\]

The Cartier data on the five affine cone charts are

\[
\begin{array}{c|ccccc}
U_\sigma&U_{B,G}&U_{G,E}&U_{E,P}&U_{P,F}&U_{F,B}\\ \hline
m_\sigma&(0,-3)&(0,-2)&(-2,0)&(-12,0)&(0,-3).
\end{array}                                                  \tag{4.8}
\]

Every vector in (4.8) satisfies all five inequalities (4.7).  Return the
four distinct character sections

\[
             \chi^{(0,-3)},\quad\chi^{(0,-2)},\quad
             \chi^{(-2,0)},\quad\chi^{(-12,0)}.              \tag{4.9}
\]

On \(U_\sigma\), use \(\chi^{m_\sigma}\) as the local frame of
\(\mathcal O(L)\).  The section \(\chi^{m_\sigma}\) then has local
representative \(1\).  Consequently every one of the five charts has the
literal unit identity \(1=1\), while the checker independently verifies
that all other local character differences lie in the corresponding dual
semigroup.  This proves from the explicit affine cover—not from a proposed
Mori cone—that \(L=2D\) is base-point-free.

As numerical cross-checks,

\[
 D^2=8,\qquad (D.B,D.E,D.G)=\left(0,1,{1\over2}\right),
 \qquad L^2=32.                                           \tag{4.10}
\]

Thus this is a genuinely nonzero weak-log calibration.  It is not the Coble
case \(D=0\).

## 5. Handoff to weak-log exact return

Suppose a P-first packet has independently certified:

1. a nonempty smooth rational cluster surface;
2. a disjoint smooth rational \((-4)\) boundary;
3. every proposed facet as an effective irreducible curve with the exact
   labelled Gram and orientation;
4. preservation of the adapted host by its facet reflections; and
5. relative, absolute, or arithmetic-reflection maximality in the precise
   sense required by the weak-log theorem.

If it also passes the finite certificate of Theorem 1.1 for
\(D=-(K+B/2)\), then \(D\) is nef.  The independently hostile-audited
[weak-log theorem](reflection_maximal_weak_log_exact_nef_return.md) now
classifies every hidden negative curve as an adapted integral root and uses
reflection saturation to prove

\[
                         P=\mathbf P(\operatorname {Nef}S).  \tag{5.1}
\]

This closes the nef/exact-return gate **for that supplied carrier packet**.
It does not prove that every raw arithmetic Coxeter chamber produces such a
packet, that every PB-Atlas incidence component is nonempty, or that
maximality manufactures effectivity.  In particular, inequalities
\(D.C_j\ge0\) against only the proposed facets are never accepted as a
semiampleness certificate.

## 6. Replay and dependency pins

Run

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_pfirst_weak_log_semiampleness.py --check
/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_pfirst_weak_log_semiampleness.py
```

The archive pins the current B1 and B2a schemas/cores.  It also pins the
weak-log exact-return theorem at
`5336e6f58377f0f5ab8220b4760d4472a578a05afd0888f3c27cde1281b221e6`
and its independent hostile audit at
`3c948a26b1841dc580a6fd0430c5415c47cf7a7f749a09e73340aa97f5f548c3`.
