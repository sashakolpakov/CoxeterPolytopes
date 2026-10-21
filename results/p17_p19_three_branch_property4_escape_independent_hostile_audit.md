# Independent hostile audit of the `P17`--`P19` three-branch escape claim

**Verdict (2026-08-11): PASS in the exact standard-host/relative-maximal
scope.**

The direct rational and same-line K3 obstructions survive their stated
scope.  An earlier producer draft reused canonical MMP orbits computed under
the primitive physical scaling and therefore failed to quantify over
boundary doublings.  The frozen producer now removes that circularity and
enumerates every adjunction colouring before MMP.  This independent
wall-relation replay confirms its result more strongly: the primitive
square-\((-2)\) facets alone span the full host in every row, so weak-log
adjunction forces \(K=0\), contradicting every remaining square-\((-1)\)
physical ray.

Two outer qualifications remain.  Complete reflection return on the declared
host proves relative host maximality, not absolute or maximal-arithmetic
reflection maximality.  And a single pinned integral host does not by itself
exhaust every commensurable normalization allowed by a broader raw-diagram
compiler.

## 1. What was audited

The producer is
[`p17_p19_three_branch_property4_escape.md`](p17_p19_three_branch_property4_escape.md),
with its script, certificate, and tests.  Its three claims are:

1. no direct rational full-nef carrier;
2. no certified rational weak-log carrier; and
3. no literal full-rank K3 carrier with rational same-line physical
   \((-2)\)-roots.

The input is explicitly the standard wall module

\[
                         L_n=I_{1,n},\qquad n=17,18,19,     \tag{1.1}
\]

not every integral realization of the unlabelled diagram.

The audited producer quartet is pinned by SHA-256:

```text
6bc82285d0ca42c4e7eeaed40d2b6302182083557e88db9c57e20e842d67b4b0  report
0755edff0bf042970ba19c87b01a79b2a44f5e21fd36662e570b2282953ad5e7  certificate generator
005b24f6ab15ee4114f25194fc38f605e456ecc693faaca327ad1119a03810af  archive
3c57d95d4ec1ac858277716c8d3d44000d4d582afd62cb588619980ac829d746  tests
```

## 2. Direct rational branch: PASS in inherited scope

The producer pins the existing marking-independent MMP/Abel certificates.
Those enumerate respectively \(8,24,120\) canonical markings for a carrier
whose displayed physical curves have the frozen standard scaling, and every
branch reaches the same degree-one Abel contradiction on an integral cubic.

This is valid inherited evidence for emptiness of that direct physical
species.  It does not independently upgrade the standard chamber to absolute
reflection maximality or enumerate other commensurable wall hosts.

## 3. Weak-log proof: historical defect removed

The historical defect was precise.  The old draft invoked the same
\(8,24,120\) canonical MMP orbits and then checked that no primitive
norm-one facet \(r\) had \(K.r=1\).  But those MMP orbits were obtained by
contracting physical \((-1)\)-curves.  In a putative weak-log escape, a
norm-one ray may instead be represented by a reduced boundary curve

\[
                              B=2r,\qquad B^2=-4.            \tag{3.1}
\]

Then \(r\) need not be effective and cannot be contracted as a \((-1)\)-curve.
Consequently the old canonical orbit did not exhaust the new physical
scalings.  The frozen producer no longer uses that argument: its Lemmas
3.1--3.2 enumerate all \(4,32768,32\) physical adjunction colourings and
reject them against the complete wall relations.  The independent
calculation below replays the same quantifier from the frozen wall source.

## 4. Independent exact check: the norm-two walls already span

Assume a certified weak-log carrier on the exact standard host.  Its divisor

\[
                            D=-(K+B/2)                       \tag{4.1}
\]

is nef.  Let \(C\) be any displayed nonboundary curve on a primitive
norm-two ray.  Nefness, nonnegative intersection of distinct curves, and
adjunction give

\[
 0\le D.C=-K.C-B.C/2,
 \qquad K.C=2p_a(C)\ge0.                                  \tag{4.2}
\]

Thus \(p_a(C)=0\), \(B.C=0\), and \(K.C=0\).  The rational reflection
integrality/adjunction gate rules out a nontrivial physical multiple on a
primitive norm-two ray, so these are the actual primitive \((-2)\)-classes.

The exact raw wall matrices give

\[
\begin{array}{c|c|c|c}
 n&\rho&\#\{r_i^2=-2\}&
 \operatorname {rank}\langle r_i:r_i^2=-2\rangle\\ \hline
17&18&20&18\\
18&19&22&19\\
19&20&45&20.
\end{array}                                                \tag{4.3}
\]

Hence (4.2) forces \(K=0\) in the whole Néron--Severi space.  But every
primitive norm-one facet must be either a physical \((-1)\)-curve, with
\(K.r=-1\), or a doubled boundary \(B=2r\), with \(K.r=1\).  Each row has at
least one such facet.  Contradiction.

For redundancy, the independent executable audit enumerates all sign
assignments

\[
 K.r_i\in\{-1,+1\}\quad(r_i^2=-1),\qquad
 K.r_i=0\quad(r_i^2=-2)                                   \tag{4.4}
\]

against every exact wall relation.  It rejects

\[
                             4,\quad32768,\quad32          \tag{4.5}
\]

assignments, with zero survivors.  This confirmation is pre-MMP and does not
use the old canonical orbits.

## 5. K3/minimal trace branch: PASS

For a rational same-line physical root on every facet, one common similarity
\(c\) and rational line scalars \(\lambda_1,\lambda_2\) on norm-one and
norm-two rays would satisfy

\[
 c\lambda_1^2=2,qquad 2c\lambda_2^2=2,qquad
 \left(\frac{\lambda_1}{\lambda_2}\right)^2=2.             \tag{5.1}
\]

This is impossible over \(\mathbb Q\).  Thus the literal physical K3
`PHYS` gate is empty.  This matches the declared same-line/common-similarity
contract and does not apply to uncontrolled projections from a larger Picard
lattice.

The third CCR branch is worded more broadly than “K3.”  The frozen producer
now closes that wording explicitly: exact return gives full Picard ranks
\(18,19,20\), finite-orbit descent cannot lower rank, and canonical
orthogonality forces \(K\equiv0\).  The complex minimal-surface
classification then leaves K3 as the only possible Néron--Severi host:
Enriques has numerical rank ten, and abelian or bielliptic surfaces have
Picard rank at most four.  Thus the same-line K3 calculation closes the exact
third branch in the stated scope.

## 6. Maximality and normalization: unresolved outer quantifiers

The producer says the rows are “genuine raw maximal arithmetic chambers.”
Its cited exact evidence proves that each is the complete reflection chamber
of the declared host \(I_{1,n}\), hence maximal **relative to that host**.
It does not prove absence of a reflection overgroup on another commensurable
host.  The distinction is material: the separate `Gamma16` maximality proof
required a genuine all-discrete-overgroup argument.

The current CCR theorem explicitly permits a branch-specific relative-host
maximality hypothesis, so the frozen standard-host theorem is still a
valid counterexample to automatic nonemptiness in that relative formulation.
It is not yet a strict counterexample to “raw maximal” if that phrase means
absolute or maximal-arithmetic reflection maximality.

Likewise, the standard-host coloring repair does not enumerate all rational
similarities, saturated wall lattices, or commensurable adapted hosts which a
diagram-only outer compiler might allow.  Its exact conclusion is therefore

> for each pinned \(I_{1,n}\) realization, the direct, weak-log, and literal
> same-line K3 species are empty; relative host maximality alone cannot make
> one of them nonempty.

The compact \(\Delta(2,4,6)\) theorem avoids both qualifications: its
maximality is absolute by a self-contained area argument, and rational-host
return already fails by anisotropy/compactness.

## 7. Replay

```bash
python3 -m py_compile \
  scripts/certify_p17_p19_weak_log_coloring_independent_audit.py

PYTHONPATH=src:. python3 \
  scripts/certify_p17_p19_weak_log_coloring_independent_audit.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_p17_p19_weak_log_coloring_independent_audit.py
```

Certificate:
[`p17_p19_weak_log_coloring_independent_audit.json`](../data/enumerations/p17_p19_weak_log_coloring_independent_audit.json).
