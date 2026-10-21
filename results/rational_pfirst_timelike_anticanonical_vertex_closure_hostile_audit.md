# Hostile audit: timelike anticanonical-vertex closure

**Audit status (2026-08-11): PASS for \(I_{1,n}\), \(3\le n\le8\), with
ordinary \((-1)/(-2)\) physical facets.  The report does not cover
twice-\((-4)\) walls or \(-K\notin P\).**

This audit checks
results/rational_pfirst_timelike_anticanonical_vertex_closure.md without
using a reflection or weak-del-Pezzo classification.

## 1. Standard simplex and finite volume

For each \(3\le n\le8\), the roots

\[
E_n,\quad E_i-E_{i+1}\ (1\le i<n),\quad
H-E_1-E_2-E_3
\]

form a determinant-\(-1\) basis of \(I_{1,n}\).  Their squares and canonical
pairings are

\[
(-1,-2,\ldots,-2),\qquad(-1,0,\ldots,0).
\]

Exact dualization returns one ideal vertex and otherwise finite vertices.
The distinguished finite vertex is

\[
F=-K=(3,-1,\ldots,-1),\qquad F^2=9-n>0.
\]

Thus the chamber is finite-volume and \(-K\) is genuinely timelike, not a
disguised cusp.

## 2. Reflection-maximality descent

The integrality reduction is valid.  In a unimodular lattice, a primitive
reflecting vector has divisibility one, so its square divides \(2\).  Only
\(-1\) and \(-2\) occur.

Signed exceptional permutations and overall root sign put such a vector in
the form

\[
(a,-m_1,\ldots,-m_n),\qquad
a\ge0,\quad m_1\ge\cdots\ge m_n\ge0.
\]

For \(s=m_1+m_2+m_3\), negative norm rules out \(s\le a\).  The critical
inequality

\[
\sum m_i^2\le m_1^2+m_2^2+(n-2)m_3^2\le s^2
\]

has exact gap

\[
(9-n)k^2+4k(A+B)+2AB
\]

after writing \(m_1=k+A,m_2=k+B,m_3=k\).  It is nonnegative throughout
the audited range.

If \(s=3a\), the norm equation gives the Cremona simple root.  Otherwise
reflection in that root strictly lowers the absolute time coefficient.
At time zero, only a coordinate root or a two-coordinate difference
remains.  Therefore every integral host reflection is generated.

This proves relative maximality in the full integral isometry group, not
only maximality of a root list or bounded search.

## 3. Canonical rigidity

Ordinary physical colors independently imply

\[
F.C=1\text{ on \((-1)\) facets},\qquad
F.C=0\text{ on \((-2)\) facets}.
\]

These signs place \(F\) in the homogeneous wall cone but do not, before
reduction, choose its positive sheet.  The proof correctly avoids that
orientation shortcut.  Relative maximality first identifies the side group
with the full reflection group.  Choose an oriented reflection isometry
\(w:P\to P_n\).  Then \(G=wF\) satisfies every standard oriented
inequality, is characteristic, and has square \(9-n\).

For a characteristic vector of square \(9-n\) in the homogeneous standard
inequality cone,
write

\[
a=2\ell+3,\qquad m_i=2k_i+1.
\]

The norm and chamber conditions are

\[
\ell(\ell+3)=\sum k_i(k_i+1),\qquad
\ell\ge k_1+k_2+k_3.
\]

The exact comparison gap is

\[
(9-n)k(k+1)+4k(A+B)+2AB+2A+2B.
\]

For \(n\le8\), equality forces all variables to vanish.  Thus \(F\) is the
unique characteristic vector of the required square in that cone.
Consequently \(G=F\): the translating reflection fixes \(F\) and \(K=-F\).
Only after this conclusion is the standard finite vertex transported back
to the chosen positive component.  Thus transporting a marked carrier
preserves canonical data without an unproved initial sheet choice.

## 4. Independent nef certificate on the carrier

The flex pencil has local ideal

\[
(z^3,z-x^3-xz^2).
\]

After the first three basepoint blow-ups it is \((x^6,w)\), up to local
units and a regular coordinate change.  Truncating after \(n\) blow-ups,
\(3\le n\le8\), leaves residual ideal

\[
(x^{9-n},w).
\]

This is a point base locus, not a fixed curve.  The transformed pencil has
class

\[
3H-\sum_{i=1}^nE_i=-K_{S_n}.
\]

A pencil with no fixed curve certifies nefness directly: for an irreducible
curve \(C\), select a member not containing \(C\), and compute a
nonnegative intersection.  This avoids the circular inference
\(-K\in P\Rightarrow-K\) nef.

The chronological centers also directly realize every facet:

- the last exceptional \(E_n\);
- the strict exceptional transforms \(E_i-E_{i+1}\); and
- the strict flex line \(H-E_1-E_2-E_3\).

## 5. Full nef return

With \(F=-K\) independently nef, adjunction for any irreducible negative
curve gives

\[
C^2=2p_a(C)-2+F.C.
\]

It forces a smooth rational \((-2)\)-curve with \(F.C=0\), or a smooth
rational \((-1)\)-curve with \(F.C=1\).  Both wall reflections are integral.
Since the side group contains every integral reflection, no hidden
negative mirror crosses the fundamental chamber interior.  Ample
orientation selects the correct side.  Nonnegative-square curves follow
from the Lorentz future-cone inequality.

The reverse containment uses only effectiveness of the displayed facets.
The exact-nef equality is therefore sound.

## 6. Audit boundary

The proved implication is:

\[
\begin{gathered}
3\le n\le8,\quad L=I_{1,n},\quad
\text{ordinary physical facets},\\
\text{relative reflection maximality}
\Longrightarrow
\text{canonical rigidity and an exact rational carrier}.
\end{gathered}
\]

The following are not proved:

- a chamber with a twice-primitive \((-4)\) wall;
- canonical escape with \(-K\) outside the chamber;
- the higher-rank case inside this timelike argument (the combined ordinary
  theorem now excludes it by the spacelike wall-cone lemma);
- a nonunimodular adapted host; or
- arbitrary surface existence from bare lattice data.

The first, second, fourth, and fifth items remain OPEN; the ordinary
higher-rank item is closed by the combined theorem.
