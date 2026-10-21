# Exact Norm-Two Vinberg Prefix for the H11 Field Lattice

## Outcome

For

\[
H_{11}=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}\oplus E_8
\quad\text{over }\mathbb Z[\phi],
\]

the exact Vinberg search has exhausted every norm-two root height whose
positive control coordinate satisfies \(0<x\leq20\).  There are 29 admissible
positive heights and 95 walls through this bound.  The resulting cone is
**not** yet compact or finite-volume: it has an exactly certified spacelike
extreme ray.

This is a bounded-prefix result.  It neither proves reflectivity nor proves
non-reflectivity of \(H_{11}\).

The full machine-readable transcript is
`data/candidates/h11_vinberg_prefix_bound20.json`.  It contains every root,
every coordinate box, and every recursion counter.  A smaller bound-five
regression transcript remains at
`data/candidates/h11_vinberg_prefix.json`.

## Control vector and algebraic reduction

Take

\[
p=(2,-1,0^8),\qquad D=1+4\phi.
\]

Then

\[
H_{11}p=(-D,0^9),\qquad (p,p)=-2D.
\]

For \(r=(x,y,z)\), put \(t=x+2y\) and \(w=E_8z\).  Direct completion of
the coupled binary block gives

\[
2(r,r)=t^2-Dx^2+2w^TE_8^{-1}w.
\]

Thus the norm-two equation is the positive-definite problem

\[
t^2+2w^TE_8^{-1}w=4+Dx^2,
\qquad t\equiv x\pmod{2\mathbb Z[\phi]}.
\]

This is the useful generalization of the Bugaenko enumerator: the coupled
binary block becomes one additional \(A_1\) coordinate and a parity gate,
while the same exact Schur-complement recursion handles the \(E_8\) weights.
The Vinberg height is

\[
\delta(r)=-\frac{(p,r)^2}{(p,p)(r,r)}=\frac{Dx^2}{4}.
\]

At the conjugate embedding, \(D'=5-4\phi<0\).  The norm equation therefore
forces

\[
(x')^2\leq\frac4{4\phi-5}=\frac{4D}{11}.
\]

Together with \(0<x\leq B\), this gives a finite, explicitly proved
coefficient box for every control-coordinate scan.

## Height zero

When \(x=0\), the equation is

\[
2y^2+E_8(z)=2.
\]

The two terms divided by two are totally nonnegative algebraic integers
summing to one.  If either is nonzero, its two embeddings lie in \((0,1]\),
so its nonzero integral norm forces both embeddings to equal one.  Hence one
term is one and the other zero.  Comparing the two coefficients of
\(E_8(a+b\phi)=2\) then forces \(b=0\).  The height-zero roots are consequently
the two ordinary \(A_1\) roots and the 240 ordinary \(E_8\) roots.

The nine oriented simple walls are

\[
(0,-1,0^8),\qquad (0,0,-e_i)\quad(1\leq i\leq8).
\]

Their chamber conditions are exactly \(t\geq0\) and \(w_i\geq0\).

## Exact enumeration

For fixed \(x\), every coordinate \(u\in\{t,w_1,\ldots,w_8\}\) obeys
\(u^2\leq4+Dx^2\) at both real embeddings.  Prefix Schur complements give
the exact minimum remaining quadratic energy.  Earlier walls give the
monotone pruning inequalities

\[
tt_j+2w\cdot z_j\leq Dxx_j.
\]

All hot-loop arithmetic is performed as integer pairs \((a,b)\) representing
\(a+b\phi\); the sign test at either real embedding is exact.  No
floating-point value makes an enumeration or pruning decision.

The 29 positive levels have solution counts

```text
1,1,1,0,1,0,3,5,1,5,0,1,5,2,1,1,9,2,4,8,0,0,9,0,3,11,0,6,6
```

Together with the nine height-zero simple roots, these give 95 walls.  The
recursion visited 65,621,346 prefix nodes; all per-depth visited, parity,
Schur-pruned, wall-pruned, and leaf-mismatch counters are archived.

## Independent hyperideal obstruction

The independent exact cone verifier uses incremental double description.  It
was cross-checked against exhaustive \(9\)-wall-subset enumeration on the two
triangle fixtures, Bugaenko's compact \(\mathbb H^8\) chamber, and the H11
bound-five cone.  In particular, both algorithms give 54 finite and 172
hyperideal rays at bound five.

At bound twenty, a dual floating-point hull was used only to *propose* an
active set.  The accepted witness was then recomputed from scratch over
\(\mathbb Q(\sqrt5)\).  Facets

\[
(0,3,4,5,6,7,8,11,85)
\]

have full rank nine and define the ray

\[
\begin{aligned}
v={}&(1,-\tfrac12,\tfrac12+2\phi,
\tfrac{11}4+\tfrac{11}4\phi,
\tfrac{55}{14}+\tfrac{55}{14}\phi,
\tfrac{22}7+\tfrac{22}7\phi,\\
&\tfrac{33}{14}+\tfrac{33}{14}\phi,
\tfrac{11}7+\tfrac{11}7\phi,
\tfrac{11}{14}+\tfrac{11}{14}\phi,
\tfrac{55}{28}+\tfrac{55}{28}\phi).
\end{aligned}
\]

Its pairings with all 95 walls are nonpositive, with exactly the nine listed
facets active, and

\[
(v,v)=\frac{81}{28}-\frac{25}{14}\phi
=\frac{81-50\phi}{28}>0.
\]

The last inequality is exact: \(81/50>\phi\), since
\((81/50)^2-81/50-1=11/2500>0\).  Therefore this is a genuine hyperideal
extreme ray, so Vinberg's finite-volume stopping criterion does not apply at
the archived bound.

## Precise remaining gap

Every admissible root with \(0<x\leq20\) has been processed in increasing
Vinberg height.  Roots with \(x>20\) have not.  A later wall may cut the
certified hyperideal direction, as happened to earlier witnesses during this
run.  Consequently the correct status is

`norm_two_vinberg_prefix_x_le_20_nonfinite_current_cone_reflectivity_undecided`.

For the K3 construction, evenness identifies the complete restricted
trace-root arrangement with these norm-two mirrors up to unit rescaling.  It
does not turn this bounded prefix into a complete chamber.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/build_h11_vinberg_transcript.py --identity-bound 20
PYTHONPATH=src python3 scripts/diagnose_h11_vinberg_prefix.py
python3 -m pytest tests/test_h11_vinberg.py tests/test_incremental_chamber.py
```
