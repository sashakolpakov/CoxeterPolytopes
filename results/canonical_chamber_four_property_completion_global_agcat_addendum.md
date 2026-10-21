# Global AGCat addendum to the canonical four-property completion audit

**Status (2026-08-12): terminal scope-preserving integration addendum.**
Low-dimensional descent and the assembled global composition have
independent hostile PASS; the full relevant regression suite reports
`31 passed`.  This file does not modify or supersede the
mathematical verdict of the frozen
[surface four-property audit](canonical_chamber_four_property_completion_audit.md)
(SHA-256
`f1c56a0ade1b37a498b759262b3309c97f96a54ae00900d11d4035a490a05e36`)
or its
[independent hostile audit](canonical_chamber_four_property_completion_independent_hostile_audit.md)
(SHA-256
`8ee64106a10655c9126990365e78e1546f7d50806a31ebab7b5936ee6baf80a0`).
It records the logically different global theorem proved after those bytes
were frozen.

## 1. Two carrier questions that must not be conflated

The canonical audit tests the strong surface question:

\[
 P\longmapsto
 (S,\operatorname {NS}(S),\text{effective divisors},\operatorname {Nef}(S))
 \longmapsto P,                                                   \tag{1.1}
\]

with exact effectivity, MMP/canonical fidelity, hidden-wall exclusion, and an
exhaustive surface carrier for every raw maximal input.  Its negative verdict
on Property 4 remains unchanged.  The three-branch surface disjunction is
false, universal projected-root K3 repair is false over the all-finite-volume
class, and a mixed divisor/action or different surface category is still
required.

The new theorem answers the broader AG-category question with a different
carrier:

\[
 P\longmapsto
 (\mathbf P(V),Q_q,H_1,\ldots,H_f,\sigma_0,C_P)
 \longmapsto P.                                                   \tag{1.2}
\]

Here all algebraic data are defined over a displayed totally real field
$k$; $Q_q$ is a smooth projective quadric, the $H_i$ are ordered polar
hyperplanes, and $(\sigma_0,C_P)$ is finite real framing data.  Morphisms are
field-isomorphism-compatible projective isomorphisms preserving the quadric,
each ordered hyperplane, and the real framing.  Quadric polarity and the
oriented component return the identical chamber.  Thus (1.2) is a genuine
exact algebraic-arrangement functor, but it is not a claim that $P$ is the
full nef cone of a surface.  The structural reduction in ABSW Lemma 5.2
citing Vinberg Lemma 7 places arithmetic reflection inputs in the admissible
quadratic/lattice model; it supplies no finiteness conclusion.

## 2. What is now globally proved

Integral rescaling-invariant Cartan-cycle functions on (1.2), together with
real-place signature constraints, give a finite global field set through the
following exhaustive branches:

- the dimension-`>=6` minimality-14 arrangement funnel;
- exact dimension-`4/5` descent to `H3` or `H2`;
- the compact-`H3` general-fork determinant coordinate; and
- the `H2` genus-zero orbifold-curve theorem, with noncompact cusps forcing
  the rational field.

The geometric dimension ceiling `n<=995` leaves finitely many field/rank
pairs.  For each pair, the connected bounded wall frame with `A=62`, finite
integral cycle packets, diagonal-gauge reconstruction, Burnside, and the
full-order conductor sandwich give finitely many full projective native
hosts.  Reflection maximality gives $W(P)=R(L)$, and a finite generic
gallery gives one orbit of complete mirror chambers.  Consequently there
are finitely many unmarked normalized reflection-maximal finite-volume
arithmetic Coxeter chambers and side-reflection groups.

This replaces the old conditional downstream use of Nikulin's
fixed-field/fixed-rank reflective-lattice finiteness conclusion.  It does not
replace or repair surface Property 4.

## 3. Attribution and surviving open claims

| statement | verdict after the global theorem |
|---|---|
| exact functor to a genuine algebraic category and exact return | **proved globally** by the real-framed quadric arrangement |
| global finite set of minimal Vinberg fields | **proved** from integral arrangement coordinates plus Fekete/orbifold arithmetic |
| finiteness of maximal chambers and reflection groups | **proved** after `n<=995`, `A=62`, reconstructed order theory, and maximality |
| universal exact surface/nef-cone Property 4 | **still false/open in the precise sense of the frozen audit** |
| direct `AG-DEG` implication from maximality to a bounded Hilbert/log degree | **open** |
| proof using algebraic geometry alone | **not claimed** |

The correct description is therefore **AG-mediated hybrid finiteness**, not
AG-only finiteness and not universal surface realization.  The exact evidence
chain is the
[global composite theorem](agcat_global_maximal_arithmetic_coxeter_finiteness.md),
its [independent composition audit](agcat_global_maximal_arithmetic_coxeter_finiteness_independent_hostile_audit.md),
[terminal synthesis](exact_agcat_functor_and_ag_finiteness_pushback_synthesis.md),
[low-dimensional descent](low_dimensional_minimality14_field_descent_to_h3.md),
its [independent hostile audit](low_dimensional_minimality14_field_descent_to_h3_independent_hostile_audit.md),
[compact-`H3` theorem](compact_h3_general_fork_arrangement_field_finiteness.md),
and [native-host reconstruction](ag_field_funnel_bounded_frame_native_host_closure.md).

The low-dimensional owner/audit SHA-256 pins are respectively
`8157c0779705f5f265e75126d10361441bf2e2b504e522d917cbda848b75f019`
and
`74c32b41b593419951ecb6db0c43c3364134cc033e77038175156bffb385b5da`.
The global composite owner/audit SHA-256 pins are respectively
`441deb5a379fc59c6b4c1e1e86212f5f84db4e603a4b37259200db43a4508dfd`
and
`6fbf3f5954d66390e82916b3fd68720ba58ed106c681e7ff7196f8e30667a44e`.
