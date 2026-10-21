# The norm-six OG10 compiler adds no projected mirror

## Global theorem

For the rank-preserving wall-first compiler

\[
H_6(k)=\begin{pmatrix}F&3k\\3k^t&6\end{pmatrix},
\qquad p_0=F^{-1}k,
\]

the exact Schur gate is

\[
0<q=2p_0^2<\frac43,
\qquad 0<q'<\frac43.
\]

The open Minkowski rectangle contains only

\[
\boxed{q=1}.
\]

Indeed, writing \(q=a+b\phi\), the two embedding bounds force
\(|b|<8/(3\sqrt5)<2\) and \(|a|<4/3+\phi<3\); exact enumeration of the
resulting fifteen pairs leaves only `1`.

For the normalized covectors \(k=(A,E_8w,C)\),

\[
y=2p_0=(-A/\phi,2w,C)\in F.
\]

Consequently

\[
y^2=4p_0^2=2q=2,
\qquad p_0^\perp=y^\perp.
\]

Thus every admissible norm-six/divisibility-three wall in this compiler
projects to a norm-two mirror already present in the complete arrangement of
\(F\).  This is a global projected-arrangement collapse, not a finite-box
observation.

## Which old roots occur?

Not every norm-two root occurs.  The exact converse is

\[
y\text{ occurs}\iff Fy/2\in F^\vee
\iff y\in F\cap2F^\vee.
\]

The `x` and `u` diagonal blocks already contain a factor two.  Since `E8` is
unimodular modulo two, the remaining condition is precisely

\[
y_{E_8}\in2E_8.
\]

In the frozen 71-root prefix this parity subset consists of root indices

\[
8,\ 9,\ 61,\ 70.
\]

Only indices `8,9` occur in the original 61-wall prefix.  Direct exact
pairing with the archived ray gives, respectively,

\[
-\frac{2\phi}{13},\quad
\frac{2-24\phi}{13},\quad
\frac{-22-32\phi}{13},\quad
\frac{-18-12\phi}{13}.
\]

All four are strictly negative, so the two newly exposed parity roots at the
next completed level do not truncate \(v_{61}\).

## Geometric consequence and scope

The compiler may still provide a different **ambient OG10 divisor lift** of
an old mirror: a selected trace class can have Beauville square `-6` and
ambient divisibility `3`.  What it cannot do is add a new hyperplane to the
fixed \(H^9\) arrangement.  The projected chamber is therefore the same old
complete norm-two chamber \(P_9\), already proved to have infinite volume.

No ambient \(H^{10}\) volume conclusion follows merely from this
codimension-one equality—finite-volume ideal polytopes can have
infinite-length sections.  Any ambient conclusion must explicitly invoke the
previously proved facet/arithmetic implication for the relevant carrier, or
supply a new Vinberg/stabilizer argument.

The source-independent certificate is
`scripts/independent_og10_golden_norm_six_projection_collapse.py`, with frozen
data in `data/candidates/og10_golden_norm_six_projection_collapse.json`.
