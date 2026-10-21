# Compact Coxeter dimension at most23: closing both dimension24 frontiers

20 October 2026. Complete global composition by root, held for
separate mathematical review. Paper only: no scientific, Lean or
TeX execution. Earlier24/26/27 sources are preserved.

**Theorem.** Every compact hyperbolic Coxeter polytope has dimension
at most23, once the exact component proofs and this composition
have passed their complete source-bound audits.

## 1. Exact input boundary

All paths below are under results/. Use the complete written
dimension24 proof and its necessary equality-dimension frontier:

`compact_dimension24_farthest_forest_composition_20261020.md`,
SHA256 `7dd4495aacc4344b2433b55a2b11c0fa3f613785483c31441c3a6aed179b986a`.

Its independently reviewed component geometry and small-defect
accounting exclude dimensions25--29; the explicit external published
ceiling29 supplies the higher-dimensional tail. At a farthest vertex
of a hypothetical dimension24 polytope, it leaves only:

* A20 and four spectator ranks, with one short A1 or one isolated
  rank-two7 defect. Fully ordinary spectators are already excluded.
* D23 and one isolated short A1 of type I.

These are necessary configurations of the SAME actual farthest
vector x, with ||x||^2<1 and additive orthogonal component energies.
No other frontier or metric-face inheritance is assumed.

Use two additional written theorems:

| Theorem | Source and SHA256 |
| --- | --- |
| Actual A20 projection>11/14, arbitrary finite spectators | `compact_a20_actual_projection_eleven_fourteenths_20261020.md`, `1720c1d864ff624fad0dd11db9312afa6034249e6b36c9caa0710a44f049c673` |
| D23 cannot coexist with isolated short A1 I at a farthest vertex | `compact_d23_isolated_a1_coupled_prefix_obstruction_20261020.md`, `530d16ad49e603da0efd844dc2fc47af5066a6a35b39e73403d601448a6fbba6` |

Both are actual-component statements with full wall feasibility
and arbitrary spectator scope, not edge-only relaxations. Their
proofs and source-bound reviews are internal mathematical inputs,
not attributed to published literature or silently external axioms.

## 2. The A20 configuration is impossible

The new A20 theorem gives actual energy>11/14. If there is one
short A1, discard that single rank. The other three spectator
ranks retain their strict ordinary rate, so

    ||x||^2>11/14+3/14=1.

If instead there is an isolated rank-two defect, its cost is
>9/100 and the remaining two ranks are ordinary. Therefore

    ||x||^2>11/14+9/100+2/14=713/700>1.

These are exactly the two permitted exceptional alternatives,
independent of the large component. Both contradict the common
actual norm. If all four ranks were ordinary the contradiction
would be stronger still.

## 3. The D23 configuration is impossible

Its sole spectator is an isolated short A1 I, exactly the hypothesis
of the new coupled D23 theorem. That theorem handles all squared
short-wall heights0<h<1/2. It uses the full D4 inventory at index19,
the same-vector prefix lower bound and the actual-wall height gap;
all three height intervals force the total norm above1. There is
no remaining alternative in this configuration.

Thus neither necessary farthest diagram can occur. Dimension24
is impossible, and the already proved global24 tail gives dimension
at most23. This applies to arbitrary compact Coxeter polytopes,
without arithmeticity. The internal arithmetic compact ceiling20
is already stronger and is unchanged; no compact example above8
or sharpness claim is made.

## 4. Verification and manuscript scope

This source is the short final implication. The actual geometry,
finite classifications, exact inequalities and small exception
isolation remain proved in the explicitly bound companion sources.
Separate source-bound reviews are required for both new inputs and
this composition. Byte/hash replay, manuscript compilation and Lean
verification are different checks and do not replace those proofs.

Keep one canonical manuscript, updating its bound and adding the
two frontier arguments. The preserved24/26/27 notes are supporting
stages rather than additional papers. The new logic has not yet
been compiled in Lean; compilation must respect the serial resource
gate. Lowering23 further, improving the other nonsharp classes and
constructing new examples remain active research objectives.
