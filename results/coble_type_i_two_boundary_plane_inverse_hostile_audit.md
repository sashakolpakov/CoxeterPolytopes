# Hostile audit: the `H11/26` plane inverse passes after a specific RR repair

**Verdict (2026-08-09):
`PASS_AFTER_H11_SPECIFIC_RIEMANN_ROCH_REPAIR`.**  The line--quintic
construction, the strong-Coble conclusion, and the literal same-marked-`P`
conclusion are correct.  The producer's original root-seed-only handoff did,
however, skip one direction: complete enumeration of integral cusp lifts and
exclusion of hidden curves do not by themselves prove that each of the ten
listed but unconstructed lifts is effective.  On every one of the four
`P2` branches, a short Riemann--Roch argument proves effectivity; the acute
wall Gram then proves irreducibility and smoothness.

This qualification does **not** invalidate the full `C_CobVin` incidence
theorem.  A geometric point of that compiler is defined to include a smooth
irreducible curve representing every raw wall.  The gap concerned only the
weaker use of that theorem after constructing a root seed plus four of the
fourteen `(-1)` walls.

The helper-free executable is
[`audit_coble_type_i_two_boundary_plane_inverse_hostile.py`](../scripts/audit_coble_type_i_two_boundary_plane_inverse_hostile.py),
with frozen
[`JSON`](../data/enumerations/coble_type_i_two_boundary_plane_inverse_hostile_audit.json)
and
[`tests`](../tests/test_coble_type_i_two_boundary_plane_inverse_hostile_audit.py).
It imports no producer function and performs no `cdd` enumeration.

## 1. Independent jet calculation

The audit rebuilds the four proximity blocks directly in blow-up charts:

```text
A: (1,1,1)             3 rows
B: (1,1)               2 rows
C: (2,2,2)             9 rows
D: (2,2,2)             9 rows
                         -------
                         23 rows on 21 quintic coefficients.
```

For the independently generated matrix `M(a,b)`, deleting zero-based rows
`11,20` gives

\[
             -64(a-b)(a+b).
\]

Symbolic ranks are

```text
rank M(a,a)   = 20,
rank M(a,-a)  = 20,
rank M(0,0)   = 20,
rank M(1,2)   = 21.
```

To exclude exceptional parameter values hidden by a function-field rank,
the audit also computes covering families of `20 x 20` minors.  Their gcds
are the nonzero constants `32` on `a=b=t` and `16` on `a=-b=t`.
Consequently both determinantal components have constant rank twenty.

The first two equalities say that every maximal minor vanishes on both
components; the displayed minor is nonzero off their union.  Therefore the
rank-drop locus over `Q` is exactly `a^2-b^2=0`, not merely a locus suggested
by one determinant.

The rejected component was also strengthened from a sample calculation to
a coefficientwise identity.  On `a=t,b=-t`, a kernel vector is

\[
 y\big((t-2)x^2+(2t-1)xy-2tx+ty^2-2ty+t\big)^2.
\]

The conic is never the zero polynomial, so every point of this component,
including `t=0`, is nonreduced.  At `a=b=1`, the other component gives
exactly the quintic printed in the producer.

## 2. Quintic and boundary geometry

An independent `Singular` run gives:

```text
absolute nonunit components:       1
absolute exponent:                 1
projective singular radical degree: 2
minimal singular primes:           2
singular-scheme length:            10
support:                            C=(0:1:1), D=(-1:2:1).
```

Thus the quintic is absolutely irreducible and has no unreported projective
singularity, including at infinity.  Direct local blow-ups at both points
give multiplicity sequence `(2,2,2)`.  The terminal quadratic cones have
discriminants `32` and `-32`; one final blow-up separates their two branches.
Hence both singularities are `A5`, each with delta invariant three.  Since a
plane quintic has arithmetic genus six, the strict transform is a smooth
geometrically integral rational curve.

The boundary-line restriction is

\[
                         Q(x,0,1)=16x^3(x-1)^2.
\]

The gradients at `A=(0,0)` and `B=(1,0)` are respectively `(0,-1)` and
`(0,-4)`, so these are smooth contact points.  Their contact orders `3+2`
exhaust Bezout intersection with the line.  After the declared five blow-ups
the two strict transforms are therefore disjoint.  The root-line restrictions

```text
AC:  -y(y-1)^4
AD:   2x(x+1)^4
BCD:  x^2(x-1)(x+1)^2
```

independently confirm the required proper and infinitely-near centers.

In the canonical `I_(1,11)` plane basis,

\[
 B_1=H-E_1-\cdots-E_5,
 \qquad
 B_5=5H-E_1-\cdots-E_5-2E_6-\cdots-2E_{11}.
\]

The audit recovers

\[
 B_1^2=B_5^2=-4,
 \quad B_1B_5=0,
 \quad K B_1=K B_5=2,
 \quad B_1+B_5=-2K.
\]

If `D~-2K` is effective, then `D.B_i=-4`, so each reduced boundary is a
component of `D`; their sum already has class `D`.  Thus
`|-2K|={B_1+B_5}`.  If `A~-K` were effective, equality
`2A=B_1+B_5` would contradict the odd component multiplicities of the
reduced right-hand side.  Therefore `|-K|` is empty and the surface is strong
Coble.

## 3. The precise transfer gap

The valid hidden-curve argument is a necessity statement.

1. The ten actual `(-2)` curves and two crossing `(-1)` roots give the
   finite-index twelve-root Type-I Coble--Mukai diagram.
2. Kondō's Propositions 3.5 and 3.7 make this the complete effective
   irreducible root set.  This excludes hidden `(-2)` curves and hidden
   `(-1)` curves meeting the two boundaries once each.
3. Any remaining actual `(-1)` curve meeting one boundary twice produces a
   nef primitive isotropic class `f=E+B_i/2`.  It must be an ideal ray of the
   root chamber, so the complete raw cusp computation lists its class.

This proves

```text
actual negative curve  ==>  one of the 26 raw wall classes.
```

It does not prove the reverse implication for an integral cusp-lift class.
Kondō's genus-one/Halphen results produce the fibration associated to `f`,
but bare lattice integrality of `f-B_i/2` is not, without a further local or
effectivity argument, the assertion that this class is represented by a
smooth irreducible exceptional curve.  The producer explicitly constructed
only the two crossings and two same-boundary lifts, leaving ten such classes.

## 4. H11-specific Riemann--Roch repair

The hostile audit reconstructs all four determinant-one plane bases directly
from the raw Gram, walls, plane class, and contraction order.  In each branch
the ten missing raw cusp lifts have the same plane-degree multiset

```text
2,2,2,2,3,3,4,4,5,5.
```

Let `E` be any one of them and let `H` be the pullback of a general line.
The raw arithmetic gives

\[
 E^2=-1,
 \qquad K.E=-1,
 \qquad \chi(\mathcal O(E))
       =1+\frac{E(E-K)}2=1.
\]

Moreover `d=H.E` lies between two and five, while

\[
                    H.(K-E)=-3-d<0.
\]

Since `H` is nef on an iterated blow-up of `P2`, `K-E` cannot be effective.
Serre duality and Riemann--Roch now give `h0(E)>=1`.  This proves effectivity
for all ten classes on all four branches without guessing their plane
equations.

Effectivity is still not irreducibility.  Write an effective representative
as `E=sum m_C C`.  Since `E^2=-1`, some component satisfies `E.C<0`; hence
`C^2<0`.  The valid hidden-curve direction above puts the class of `C` among
the raw walls.  But the independently recomputed off-diagonal raw wall Gram
lies in `[0,5]`.  Therefore no raw wall distinct from `E` can have negative
pairing with `E`, and `[C]=[E]`.  Positive plane degree forces multiplicity
one and leaves no residual effective divisor.  Thus `E` is irreducible.
Adjunction gives arithmetic genus zero, so it is a smooth rational
`(-1)`-curve.

We have now proved the two logically distinct halves:

```text
all 26 listed wall classes are represented by the required curves;
every irreducible negative curve has one of those 26 classes.
```

Consequently the full nef cone is literally the marked raw chamber.  The
independent cone certificate remains `508` finite, `42` ideal, and zero
hyperideal rays, hence finite-volume and noncompact.

## 5. Lineage verdict

The producer's final qualifications are accurate.  Kondō already constructed
and classified the two-boundary Type-I surface, and already proved the
twelve-root Coble--Mukai chamber.  The abstract `H11/26` full-nef chamber was
also already frozen in this project before this inverse.  The defensible new
claim is only:

> a project-new raw-wall-derived `P2` line--quintic model and a repaired exact
> same-marked-`P` inverse for that already known project chamber.

Worldwide priority for this plane presentation remains open.  It is not a
new surface, a new abstract polytope, a compact example, or a class-wide L4
theorem.

Primary AG source: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially Propositions 3.4, 3.5, 3.7, 3.8, Lemma 3.11, Example 4.1,
Proposition 4.2, and Theorem 1.1.

## 6. Bounded replay

```bash
PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 scripts/audit_coble_type_i_two_boundary_plane_inverse_hostile.py \
  --check-output data/enumerations/coble_type_i_two_boundary_plane_inverse_hostile_audit.json

PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 -m pytest -q tests/test_coble_type_i_two_boundary_plane_inverse_hostile_audit.py
```

The standalone certificate peaked at about `0.082 GiB` in the hostile run.
