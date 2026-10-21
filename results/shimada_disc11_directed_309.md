# A 309-wall K3-generated Coxeter subsystem in dimension nine

**Status (2026-08-05).** Exact theorem and reproducible checkpoint.  This
strengthens the certified wall lower bound for the discriminant-11 K3
eigen-nef chamber from 181 to 309.  The recorded 309-wall outer polyhedron is
still infinite-volume: an exact hyperideal projective ray survives.  Thus the
result is not a finite-volume or compactness theorem for the complete chamber.

## The algebro-geometric object

Let (X_{11}) be the marked singular K3 surface used in the coupled
discriminant-11 construction, with

\[
 \operatorname{NS}(X_{11})\cong
 U+T_{11}(-1)+E_8(-1)^2.
\]

The self-adjoint golden action cuts out a rank-ten real-multiplication
summand whose field form is

\[
 H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}+E_8,
 \qquad \operatorname{sign}(H)=(9,1).
\]

Projectivizing its negative cone gives (mathbb H^9).  Intersecting this
summand with the K3 nef cone produces the geometric chamber under study.
Every wall accepted below is obtained from an integral physical class
(C\in\operatorname{NS}(X_{11})) with (C^2=-2), located in Shimada's
complete induced-nef atlas.  Its restriction and exact normalization gives a
root (r\in\mathbf Z[\phi]^{10}) with (H(r,r)=2).  In this precise sense,
the divisor lattice of one algebraic surface is the wall generator; the
Lorentzian signature is supplied by Hodge theory rather than imposed after
the fact.

The enhanced automorphisms

\[
 V_4=\{1,A,B_+,B_-\}\subset \operatorname{Aut}(X_{11},\tau)
\]

preserve the complete eigen-nef chamber.  Consequently one physical curve
wall automatically supplies its entire (V_4)-orbit.

## Exact theorem

There is an explicitly listed (V_4)-stable set

\[
 R_{309}\subset \mathbf Z[\phi]^{10},\qquad |R_{309}|=309,
\]

with the following properties.

1. Every (r\in R_{309}) has (H(r,r)=2), has the ample orientation, and is
   induced by an integral physical Néron--Severi class of square (-2).
2. Distinct roots have nonpositive inner product.  All
   (inom{309}{2}=47,586) pairs have recognized Coxeter type:

   \[
   2680\ (m=2),\qquad 776\ (m=3),\qquad 460\ (m=5),
   \qquad 43,670\ \text{ultraparallel}.
   \]

   There are no unrecognized finite angles.
3. The half-space system is (V_4)-invariant and is a certified rational-curve
   facet lower bound for the complete K3 eigen-nef chamber.  Its canonical
   ordered-root hash is

   ```text
   d5b9030d49c16fd55f2186ff115ad23eb6afe99f489a62682049691cb069d197
   ```

4. Its finite outer cone is not compact and not finite-volume.  An exact
   feasible hyperideal ray is

   \[
   \begin{aligned}
   v={}&\bigl(-1,
   (-164+415\phi)/1021,
   (-941+196\phi)/1021,
   (-840-16\phi)/1021,\\
   &(-1781+180\phi)/1021,
   (-1781+180\phi)/1021,
   (-1781+180\phi)/1021,\\
   &(-1781+180\phi)/1021,
   (-1781+180\phi)/1021,
   (-1680-32\phi)/1021\bigr),
   \end{aligned}
   \]

   with

   \[
   H(v,v)=\frac{5514638-3408230\phi}{1042441}>0.
   \]

   Its active facets are
   ({4,5,6,8,9,11,13,249,253}), every other wall pairing is
   nonpositive, and the complete pairing hash is

   ```text
   2588b7add212ce29cc7df2021c255568d7524a754330d0176868d9a89dbd642a
   ```

Hence (R_{309}) is a new explicit, surface-generated Coxeter subsystem in
hyperbolic dimension nine, but it is not yet the desired finite-volume
example.

## Construction and proof

The proof begins with the 97 exact Vinberg walls already transported to
physical rational-curve walls on (X_{11}).  Closing them under (V_4)
gives 181 distinct roots.  Two directed curve searches give 189.  Starting
from that frozen base, 30 more directed searches each produce a new orbit of
four roots:

\[
189\longrightarrow193\longrightarrow\cdots\longrightarrow309.
\]

At each step a floating linear program is used only to propose an active set
and an escape direction.  It is not part of the proof.  The accepted step is
then checked exactly in (mathbf Q(\phi)):

- the proposed ray satisfies every current wall inequality and is
  hyperideal;
- the dyadic point immediately before the light cone is strict on all current
  walls;
- Shimada's exact gallery ends at a `_minustwo` boundary;
- the boundary functional uniquely normalizes to an integral norm-two field
  root;
- the corresponding physical Néron--Severi vector is integral of square
  (-2);
- the crossing point is timelike, lies on the new mirror, and is strict on
  every prior wall;
- the full (V_4)-orbit adds four new roots, remains pairwise acute, and cuts
  the chosen escape-ray orbit.

The complete transcript contains 2,217 Shimada auxiliary transitions.  A
standalone replay reconstructed all 30 galleries, physical roots, field
normalizations, crossing faces, (V_4)-orbits, and terminal pairings.  Its
sealed payload hash is

```text
133759a4077a3c4fa543fb259228c8b9d4f7b20512ea82adb0610f6de521a346
```

The transcript and checkpoint are in
[`data/candidates/shimada_disc11_directed_309.json`](../data/candidates/shimada_disc11_directed_309.json).
The verifier and resume implementation are in
[`scripts/probe_shimada_disc11_hyperideal_escape.py`](../scripts/probe_shimada_disc11_hyperideal_escape.py).

An independent replay is run by

```bash
PYTHONPATH=src:. python3 scripts/probe_shimada_disc11_hyperideal_escape.py \
  /path/to/disc11.txt \
  --verify-checkpoint data/candidates/shimada_disc11_directed_309.json
```

Only after this verification may the search be extended:

```bash
PYTHONPATH=src:. python3 scripts/probe_shimada_disc11_hyperideal_escape.py \
  /path/to/disc11.txt \
  --resume-checkpoint data/candidates/shimada_disc11_directed_309.json \
  --additional-orbits N
```

## What this says algebraically, and what remains open

The useful algebraic object is not merely the abstract lattice
(operatorname{NS}(X)).  It is the package

\[
(X_{11},\operatorname{NS}(X_{11}),\operatorname{Nef}(X_{11}),
\tau,V_4,h),
\]

where effectivity distinguishes the actual curve walls, the ample class
chooses their orientation, and the automorphisms propagate them.  The
directed search is therefore an algebraic curve-discovery procedure, not a
blind enumeration of Lorentzian vectors.

This does not settle the global finiteness question.  The complete chamber
may contain further curve walls that cut the displayed hyperideal ray.  The
separate fixed-curve theorem still gives the geometric decision criterion:
finite enhanced covolume is equivalent to finiteness of the
(operatorname{Aut}(X_{11},\tau))-orbit of the unique genus-one curve fixed
by (A).  The 309-wall theorem advances the constructive boundary while
leaving that orbit, and hence full finite volume, open.
