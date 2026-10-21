# Hostile audit of the Type-I downstairs anti-Gram census

**Verdict:** `PASS/FREEZE`, with a strict strengthening: the rank-ten residue
is eliminated.

The hostile archive is
[`enriques_type_i_downstairs_anti_gram_census_hostile_audit.json`](../data/enumerations/enriques_type_i_downstairs_anti_gram_census_hostile_audit.json),
rebuilt by
[`audit_enriques_type_i_downstairs_anti_gram_census_hostile.py`](../scripts/audit_enriques_type_i_downstairs_anti_gram_census_hostile.py).

## Independent switching and inertia replay

The production certifier gauges sorted breadth-first-tree edges negative.
The hostile replay instead gauges sorted depth-first-tree edges positive.
Starting again from Shimada's ten intersection-one and three
intersection-two edges, it independently generates all

\[
2^{10}3^3=27,648
\]

anti-Grams.  The different gauge again gives exactly eighteen classes.  For
every class, the complete diagonal-sign orbit is generated directly and its
size agrees with `2^(12-c)`, where `c` is the number of support components.
The orbit-size census is

```text
 512:  2 orbits
1024:  6 orbits
2048: 10 orbits.
```

The audit does not call the project's congruence-inertia routine.  It
factorizes each exact characteristic polynomial and uses exact Sturm root
counts.  All eighteen inertias and all weighted multiplicities agree with
the production archive.  There is one negative-semidefinite rank-nine orbit
and one negative-semidefinite rank-ten orbit before the rootlessness gate.
The rank-nine orbit agrees exactly, after sheet switching, with the explicit
geometric package.

## The rank-ten killer

In the unique rank-ten orbit, the eleventh and twelfth anti-generators are
orthogonal:

\[
\delta_{11}^2=\delta_{12}^2=-4,
\qquad \delta_{11}\cdot\delta_{12}=0.
\]

Their downstairs roots satisfy

\[
r_{11}+r_{12}=(8,4,8,14,12,10,8,6,4,2)\in2L_{10}.
\]

For any independent choices
`eps_11,eps_12 in {+1,-1}`, half-sum integrality gives

\[
e=\frac{\epsilon_{11}\delta_{11}+\epsilon_{12}\delta_{12}}2
  =b_{11}^{\epsilon_{11}}+b_{12}^{\epsilon_{12}}
   -\frac{u_{11}+u_{12}}2
\]

as an integral anti-invariant vector in the saturation of the generated
span.  All four sign choices have

\[
e^2=\frac{-4-4}{4}=-2.
\]

Therefore the saturated anti-lattice is not rootless, and `L2` eliminates
the entire rank-ten orbit.  This is a glue obstruction, not a numerical
signature heuristic.

Before saturation, the audit also identifies the generated quotient lattice
exactly as

\[
E_8(-2)\oplus\langle-4\rangle^2.
\]

It is itself rootless; the root appears precisely because the half-sum glue
forces an index-two saturation vector.  This distinction is why checking
only the displayed Gram would have left a false candidate alive.

## Strengthened theorem and scope

The frozen statement is now:

> For the Type-I `12A` wall system, the 27,648 two-sheet anti-Grams reduce,
> after negative-semidefiniteness, half-sum glue, and saturated rootlessness,
> to the unique rank-nine geometric orbit up to independent sheet swaps.

Thus there is no surviving rank-ten period branch, and a Picard-rank-19
assumption is no longer needed to force the anti-Gram.  What remains is a
different compiler layer: recover the surviving rank-nine primitive
embedding and glue orbit, period orientation, and `L4` certificate from the
downstairs datum without importing OSCAR's known carrier.  No new surface or
abstract Coxeter polytope is claimed.
