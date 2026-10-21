# Exact Vinberg Transcript for Bugaenko's Compact H8 Chamber

Status: independently recomputable exact wall enumeration.

## Lattice, normalization, and control vector

Let

\[
L=\langle-2\phi\rangle\oplus E_8
\quad\text{over}\quad \mathcal O=\mathbb Z[\phi],
\qquad \phi^2=\phi+1,
\]

and take the control vector (p=(1,0,\ldots,0)), of norm (-2\phi).
This choice is intrinsic to the displayed orthogonal splitting and lies on the
distinguished hyperbolic sheet.

The Smith invariants of (L) are eight units and (2), up to the unit
\(-\phi\).  Since (2) is inert in \(\mathbb Z[\phi]\), possible primitive
root norms modulo squared units are represented by (1,2,4).  Evenness
excludes (1).  If a primitive root (r=(x,y)) had norm (4), reflection
divisibility would give (E_8y\in2\mathcal O^8), hence (y\in2E_8).  Reducing

\[
-2\phi x^2+E_8(y)=4
\]

modulo (4) then gives (x^2\in2\mathcal O).  Inertness of (2) forces
(x\in2\mathcal O), contradicting primitivity.  Thus every root hyperplane
has a unit multiple of norm (2), and the search loses no hyperplanes by using
norm-two representatives.

For (r=(x,y)) of norm (2), oriented so that (p) is on its negative side,
(x>0) and Vinberg's exact distance parameter is

\[
\delta(r)=-\frac{(p,r)^2}{(p,p)(r,r)}=\phi x^2.
\]

## The height-zero cone

At (x=0), write (y=a+b\phi), with (a,b\in E_8).  The equation
(E_8(y)=2) is equivalent to

\[
E_8(a)+E_8(b)=2,
\qquad 2(a,b)+E_8(b)=0.
\]

Because (E_8) is positive definite and even, this forces (b=0) and
(E_8(a)=2).  Hence the height-zero roots are exactly the 240 ordinary
(E_8) roots.  Exact closure under the eight simple reflections recomputes
all 240 and verifies that the oriented roots (-e_1,\ldots,-e_8) delimit a
fundamental cone.

## Positive-height enumeration

Put (z=E_8y).  The cone inequalities are precisely (z_i\ge0), and the norm
equation becomes

\[
z^T E_8^{-1}z=2+2\phi x^2=:N_x.
\]

At the conjugate embedding, positivity gives (x'^2\le\phi).  Exhaustive
coefficient enumeration under (0<x\le7) and this conjugate bound gives,
in increasing height order,

\[
1,\ \phi,\ 1+\phi,\ 2\phi,\ 1+2\phi,\ 2+2\phi,
\ 1+3\phi,\ 2+3\phi.
\]

For each (x), a coordinate of (z) satisfies (z_i^2\le2N_x) at both real
embeddings, because the minimum of (z^TE_8^{-1}z) with (z_i) fixed is
(z_i^2/2).  This gives a finite, explicitly archived coordinate alphabet.
The recursive search uses exact rational Schur complements to reject a prefix
only when its minimum possible remaining energy exceeds (N_x) at one of the
two embeddings.  Linear inequalities from earlier walls provide the other
exact pruning rule.

The eight exact solution counts are

\[
1,\ 1,\ 0,\ 0,\ 0,\ 0,\ 0,\ 1.
\]

Together with the eight cone roots, the three generated positive-height roots
are exactly the eleven walls in the independently verified chamber
certificate.  In particular, the six intervening empty levels are proved
empty rather than skipped.

## Termination and independent check

The search transcript and polyhedral verifier are separate.  The former
exhausts every admissible root through the last generated wall in Vinberg
height order.  The latter recomputes all projective extreme rays from the
eleven inequalities and proves that they are 41 finite vertices, with no cusp,
hyperideal ray, or opposite-sheet ray.  Vinberg's stopping criterion therefore
applies: once the current acute chamber is compact, no later root can be a new
simple wall.

The archived transcript is
`data/known_examples/bugaenko_vinberg_transcript.json`.  Reproduce it with

```sh
PYTHONPATH=src python3 scripts/build_bugaenko_vinberg_transcript.py
```

and verify both enumeration and compact termination with
`tests/test_bugaenko_vinberg.py`.
