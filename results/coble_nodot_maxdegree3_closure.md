# The complete subcubic nodot theorem, and why it is still not strict L4

**Status (2026-08-09): exact class-wide carrier theorem on a natural raw
denominator, with two positive same-polytope returns.  The result closes
essential surjectivity in this denominator but deliberately does not claim
strict L4: the terminal algebraic geometry does not sharpen the P-side
classification.**

## Result

Consider twelve proposed norm-minus-two Coble--Mukai roots.  Require only
that

1. their nonorthogonality support is connected on twelve vertices and has
   maximum degree at most three;
2. every nonzero pairing is `1` or `2`, with at most three pairing-two
   pairs;
3. the root Gram has exact inertia `(1,9,2)`;
4. its saturated quotient embeds in the even unimodular lattice `E10`;
5. the root cone has finite volume; and
6. adjoining the characteristic boundary and every primitive integral cusp
   lift gives the complete one-boundary raw chamber of `D_1_split`.

This is a surface-free condition.  It mentions no Coble surface, K3 period,
Kondō type, known Gram, incidence solution, or compiler success.

> **Theorem.** Every member of this denominator is integrally isometric to
> exactly one of the Type-I or Type-II root chambers.  Its complete cusp lift
> is respectively the `H10/22` or `H10/30` chamber.  The exhaustive split-K3
> compiler constructs a complex one-boundary strong Coble surface and returns
> the same marked full nef chamber.  Thus every raw member has an exact
> carrier; there are exactly two positive polytope isometry classes.

The two polytopes and surfaces were known.  What is new in this project is the
closed raw-denominator theorem and the exact direction

```text
arbitrary P in the declared class
  -> exact finite Coxeter census
  -> exhaustive K3/Coble compiler
  -> complex strong Coble surface
  -> the same marked P.
```

## Exact census

A connected twelve-vertex graph of maximum degree three has between eleven
and eighteen edges, so the edge schedule is exhaustive.  Across all four
pairing-two counts the exact denominator is

```text
10,456,309 weighted graph trials
     1,038 Lorentzian corank-two marked presentations
       837 weighted root-diagram isometry classes
         2 finite-volume root classes
         2 complete cusp-lifted full chambers
         2 exact same-P Coble returns.
```

The two new tail computations are:

| support edges | marked Lorentzian presentations | weighted classes | finite-volume classes |
|---:|---:|---:|---:|
| `11--13` | 156 | 130 | 1, Type I |
| `16--18` | 408 | 302 | 0 |

In the low tail, pairing-two counts one and two contribute no cone without a
hyperideal extreme ray.  Pairing-two count three contributes exactly one:
its root cone has `14 finite + 9 ideal + 0 hyperideal` vertices and passes the
literal Vinberg test.  Its weighted diagram is Type I.  Every one of the 302
high-tail classes has an exact hyperideal extreme ray before an `E10`
saturation is imposed.

The frozen fourteen-edge census contributes 194 classes and no survivor.
The four frozen fifteen-edge shards contribute 211 classes and exactly one
survivor, Type II.  Their two integral `E10` markings give the same full
`H10/30` chamber.  These middle archives also supply independent exact
double-description and literal Vinberg checks; “noncompact” is nowhere
confused with “finite volume.”

## Exact carriers and all-carrier exhaustiveness

The Type-I and Type-II full cusp lifts are not merely matched to surface
names.  The one-boundary K3 compiler exhausts every sheet allocation of the
twelve roots, every saturated anti-invariant lattice, every even overclosure,
every primitive K3 graph glue, and the complete root-regular irreducibility
gate.  Its two passing branches construct one-boundary strong Coble surfaces.

Kondō's effective-root theorem and the primitive-isotropic/Halphen argument
then exclude every hidden root and every hidden `(-1)` wall.  Consequently

\[
                 \mathbf P(\operatorname{Nef}S)=P
\]

with the original integral orientation.  Conversely, any exact carrier of
the declared species appears in the finite compiler.  Hence this is genuine
class-wide essential surjectivity, not a classification of only the
carrier-positive subset.

## Why this does not earn strict L4

The five logical gates have the following scores.

| gate | verdict |
|---|---|
| independent raw class | pass |
| a carrier for every raw member | pass |
| exhaustive carrier compiler | pass |
| exact same-`P` return | pass |
| independent AG input | pass: Kondō supplies existence and classifies the one-boundary finite-Aut surfaces |
| new nonthin theorem on `P` caused by AG | **fail** |

The failure in the last row is substantive.  Exact P-side cone enumeration
has already isolated Type I and Type II before Kondō is invoked.  Algebraic
geometry supplies the carriers, their interpretation, hidden-wall
exhaustion, automorphism groups, and moduli lineage, but it does not remove a
third raw polytope, lower a dimension bound, or otherwise sharpen the
polytope classification.  The honest rating is therefore

```text
L3 + class-wide scoped essential surjectivity;
not strict L4 under the nonthin-payoff requirement.
```

This is still the strongest current positive benchmark for the requested
direction: both positive inputs are chosen by a natural raw Coxeter grammar,
not by surface success, and every input returns exactly.

## A sharp no-go: automaticity on all of `D_CobVin` is false

There is also an explicit reason the theorem cannot simply delete the
subcubic/twelve-root hypothesis.  The classical `P10` chamber
`R12-E8A1x2-C2-H10-12` passes all of the broad raw one-boundary arithmetic
tests:

```text
unique characteristic norm -1 boundary;
boundary orthogonal complement = E10;
ten roots spanning the complete finite-volume E10 chamber;
one primitive ideal cusp;
the unique exceptional wall is exactly f-b;
the complete twelve-wall package regenerates from P.
```

Nevertheless it cannot be the full nef cone of a complex one-boundary strong
Coble surface.  Finite root-reflection index would force finite automorphism
group, while Kondō's exhaustive complex one-boundary list has Types I and II,
both with twelve effective irreducible roots; `P10` has ten.  Therefore

\[
                 R_{\rm CobVin}(\mathbf C)\subsetneq D_{\rm CobVin}.
\]

This is a concrete counterexample to the maximal incidence-automaticity
lemma, not merely a missing proof.  Any future strict-L4 subdenominator must
exclude `P10` for a structural P-side reason while remaining broad enough
that the final AG theorem has something nontrivial left to say.  Adding
“twelve roots” does exclude it, but the present finite census then isolates
the two answers before AG, which explains the exact bottleneck.

The displayed strict containment is intentionally a **complex** statement.
It uses Kondō's characteristic-zero one-boundary classification exactly as
replayed by the hostile audit.  It does not claim an obstruction in
characteristic two, in every positive characteristic, or for a different
surface/carrier species.

## Remaining theorem target

The smallest honest next target is now:

> Find a natural raw condition broader than the answer-isolating subcubic
> census, not defined by incidence-stack nonemptiness, such that every input
> has a geometric Coble incidence branch, but several P-side possibilities
> remain before Kondō is applied.  Kondō must then eliminate or bound some of
> those possibilities.

Equivalently, one needs a structural incidence-realization theorem which is
stronger than local Vinberg/cusp arithmetic (because of `P10`) and softer
than an exhaustive graph grammar that has already computed the answer.
Mnëv--Sturmfels universality rules out a theorem based only on local expected
dimension or bounded incidence neighbourhoods; a global realizable design,
period theorem, or a different carrier species is required.

## Reproducibility

The complete new-tail ledger is
[`coble_nodot_maxdegree3_closure.json`](../data/enumerations/coble_nodot_maxdegree3_closure.json),
generated by
[`certify_coble_nodot_maxdegree3_closure.py`](../scripts/certify_coble_nodot_maxdegree3_closure.py)
and tested by
[`test_coble_nodot_maxdegree3_closure.py`](../tests/test_coble_nodot_maxdegree3_closure.py).
Every new root cone has a frozen exact ray digest and a replayed strict
timelike interior certificate.  The middle lanes, full cusp lifts, K3
compiler, uniform hidden-wall theorem, and `P10` hostile audit are pinned by
SHA-256.

Replay with

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_maxdegree3_closure.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_maxdegree3_closure.py
```

The independent AG source is Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the main classification theorem and the effective-root and
genus-one-fibration results used by the exact-return compiler.
