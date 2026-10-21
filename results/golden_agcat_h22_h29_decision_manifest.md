# Golden `H22--H29`: what AGCat proves, and what must still be enumerated

## Verdict

The proposed nonexistence conclusion does **not** currently follow from the
AGCat functor.

What is proved is

\[
 \boxed{\text{only finitely many reflection-maximal arithmetic chambers occur}}
\]

after combining the exact real-framed arrangement carrier with classical
hyperbolic dimension/frame bounds and the cycle/order/conductor argument.
For the fixed slice

\[
 K=\mathbf Q(\sqrt5),\qquad 22\le n\le29,
\]

this gives a finite decision problem.  It does not say that the finite set is
empty.  No complete cycle-packet, reflection-order, conductor-lattice and
Vinberg-host enumeration has yet been run for any of the eight native ranks
`23,...,30`.

The new machine manifest

[`golden_agcat_h22_h29_decision_manifest.json`](../data/enumerations/golden_agcat_h22_h29_decision_manifest.json)

freezes that distinction and materializes the first universal finite stage.
Its verdict is deliberately `OPEN`, not a bounded-search synonym for
nonexistence.

## Two AG statements which cannot be conflated

### 1. Universal real-framed arrangement AGCat

Every finite-volume arithmetic Coxeter chamber in the Vinberg model maps to

\[
 (\mathbf P(V),Q_q,H_1,\ldots,H_f,\sigma_0,C_P),
\]

and polarity plus the real framing returns the identical chamber.  This is a
universal exact functor, but the carrier is the original quadratic form and
its polar hyperplanes written as an algebraic arrangement.  Its exact return
is therefore an encoding theorem, not an independent dimension obstruction.

Integral cycle functions on this carrier are useful: with the external
dimension ceiling, the connected `A=62` frame theorem, Burnside and the
conductor sandwich, they prove finiteness of maximal native hosts.  The
implication is

```text
AG arrangement + hyperbolic geometry + arithmetic of orders
       => a finite host set,
```

not

```text
AG arrangement => the host set is empty.
```

### 2. Coherent full-`H^2` golden K3 carrier

If a chamber has a coherent same-chamber complex-K3 carrier with a unital
self-adjoint quadratic action on all of `H^2`, then the period plane consumes
two additional field slots and

\[
       2(n+3)\le b_2(\mathrm{K3})=22,
       \qquad n\le8.
\]

That is a genuine Hodge-theoretic obstruction.  However, essential
surjectivity from arbitrary golden Coxeter chambers to this K3 category is
not proved.  In fact, for dimensions `22--29` the proposed K3 carrier is
impossible on rank grounds; this does not imply that the original chambers
are impossible.  Using it that way would assume the missing universal arrow
whose consequence one is trying to prove.

Thus the universal AGCat carrier has universality but no independent
dimension leverage, while the coherent K3 carrier has strong leverage but no
universal existence theorem.  The machine manifest records both flags.

## The newly materialized universal golden edge box

Let `m=n+1` and select the connected `m`-wall frame supplied by the bound
`A=62`.  For every unordered pair put

\[
 b_{ij}=a_{ij}a_{ji}\in\mathbf Z[\varphi].
\]

At the distinguished place and the positive-definite conjugate place,

\[
 0\le \sigma_0(b_{ij})<62^2,
 \qquad 0\le \sigma_1(b_{ij})<4.                  \tag{1}
\]

The second inequality is strict for a frame pair.  Equality in
Cauchy--Schwarz for the positive-definite conjugate form would make the two
embedded roots proportional, contradicting their independence.  The weak
closed Cauchy box has one extra element, `4`; compactness also independently
excludes the corresponding parallel selected facets.

The certifier enumerates every golden integer in (1), with exact radical
signs and no floating-point interval test:

| finite stage | exact count |
|---|---:|
| weak closed Cauchy superset | 6,878 |
| exact independent-frame box | 6,877 |
| Cauchy-equality value removed | 1 (`4`) |
| invalid subparallel non-Coxeter values | 2 |
| valid Coxeter pair alphabet | **6,875** |
| dotted values in that alphabet | 6,869 |

The two rejected values are

\[
             2-\varphi,\qquad 3-\varphi.
\]

For a finite Coxeter label `r`,
`b=4 cos^2(pi/r)` and
`Q(b)=Q(cos(2pi/r))`.  Degree at most two first restricts to
`r in {2,3,4,5,6,8,10,12}`; labels `8` and `12` generate
`Q(sqrt(2))` and `Q(sqrt(3))`.  Hence the exact golden finite-angle alphabet
is

\[
\begin{array}{c|cccccc}
r&2&3&4&5&6&10\\ \hline
b&0&1&2&1+\varphi&3&2+\varphi.
\end{array}
\]

The parallel value `4` is absent.  Every retained value with distinguished
value greater than four is a possible dotted product.  The JSON stores all
6,877 exact independent-frame box elements and all 6,875 filtered elements,
together with independent SHA-256 digests; it separately records the 6,878
element weak closed superset and its unique equality-boundary value.

### Exact local cycle-packet tranches

The implementation now goes one step beyond an edge alphabet.  It exhausts
every complete triangle whose three edges are finite, nonorthogonal Coxeter
relations over the golden field.  There are 35 unordered label multisets from
`{3,4,5,6,10}`.  If the edge products are `c12,c23,c31`, cyclic integrality
and symmetrizability require an integral golden integer

\[
 b_{123}^2=c_{12}c_{23}c_{31};
\]

the acute orientation chooses the negative root.  The exact determinant is

\[
 \det G=8-2(c_{12}+c_{23}+c_{31})+2b_{123}.           \tag{3}
\]

Only ten of the 35 multisets have an integral golden three-cycle.  Requiring
negative determinant at the selected place and positive determinant at the
conjugate place leaves exactly five:

\[
 (3,3,5),\quad(3,5,5),\quad(4,4,5),\quad(5,5,5),
 \quad(5,10,10).                                      \tag{4}
\]

This is a complete local cycle-packet census, modulo all triangle vertex
permutations, and therefore a real constraint for the higher-rank packet
CSP.

The local triangle tranche is now complete even after restoring every dotted
edge.  The 6,874 nonzero edge values give 54,158,527,500 unordered raw edge
multisets.  An exact element-squareclass XOR join reduces these to 4,768,141
integral golden square-product packets; 2,633,885 have positive conjugate
determinant, and every one of those has Lorentzian selected restriction.
There is one raw selected-degenerate packet, but its conjugate determinant is
also zero, so no selected-degenerate packet survives the conjugate
positive-definite gate.  The complete certificate and stream digests are in
[`golden_agcat_full_triangle_packets.md`](golden_agcat_full_triangle_packets.md).

The queue now also has its first exact four-node join.  Modulo \(S_4\), all
connected finite/orthogonal assignments give 2,380 orbits, of which 327 are
simultaneously cycle-compatible and 119 pass both real places.  With exactly
one dotted edge, moving that unique edge to \(01\) and quotienting its
order-four stabilizer gives 14,287,520 canonical assignments after ranging
over all 6,869 dotted values.  The simultaneous triangle/four-cycle
squareclass gate leaves 532,734 and the two-place Gram gate leaves 180,402.
This includes 3,277 retained diamonds and 1,168 retained complete \(K_4\)
packets, so overlapping cycles are genuinely checked.  The exact theorem and
stream digests are in
[`golden_agcat_four_node_at_most_one_dotted.md`](golden_agcat_four_node_at_most_one_dotted.md).
Its clean-room hostile replay independently reconstructs the edge alphabet,
squareclasses, orbit quotients, determinant expansion, all 48 sign strata,
four stream hashes, and the full-rank path counts; it is terminal `PASS` in
[`golden_agcat_four_node_at_most_one_dotted_hostile_audit.md`](golden_agcat_four_node_at_most_one_dotted_hostile_audit.md).

There is also an exact global path family in every native rank 23--30: one
terminal dotted product \(c\), all remaining path edges labelled 3, and
leading determinants \(D_j=2j-(j-1)c\).  The number satisfying conjugate
positive definiteness ranges from 3,590 in rank 23 to 3,552 in rank 30;
\(c=2+2\varphi\) works uniformly.  Positive definiteness at the conjugate
place also makes each ambient arithmetic orthogonal lattice cocompact.  These
are therefore cocompact arithmetic host/frame seeds, not finite-volume
simplices or completed reflective hosts.  In particular they prove that
frame compatibility by itself does not empty the golden queue.

The next bounded local tranche is closed as well: every connected four-node
tree with exactly two dotted edges.  From 471,865,955 canonical assignments,
97,172,716 pass both real places: 25,023,281 stars, 81 paths with disjoint
dotted ends, and 72,149,354 paths with adjacent dotted edges.  A clean-room
audit caught the thin 81-packet disjoint locus, correcting an earlier blanket
inertia rejection, and independently reproduced the final exact count.  See
[`golden_agcat_two_dotted_four_node_trees.md`](golden_agcat_two_dotted_four_node_trees.md)
and its
[`independent hostile audit`](golden_agcat_two_dotted_four_node_trees_independent_hostile_audit.md).

This still is not a high-rank packet census.  Four-node subsets with two or
more dotted edges on cycles, subsets with at least three dotted edges, and
their joins into connected ranks 23--30 remain open, as do the gauge, order,
conductor, and Vinberg-host stages.

## Exact rank ledger and finite debt

For every cycle length `r=2,...,30`, the manifest records the finite AGCat
Minkowski box

\[
 \mathcal S_{K,r}=\{x\in\mathbf Z[\varphi]:
 |\sigma_0(x)|<62^r,\ |\sigma_1(x)|\le2^r\}.
\]

It also gives a proved integral coefficient rectangle containing the box.  If
`x=a+b phi`, then one may take

\[
 |b|\le {62^r+2^r\over2}-1,
 \qquad
 |a|\le2^r+{62^r+2^r\over2}-1.                       \tag{2}
\]

The eight native ranks require 204 rank/cycle-box applications in total:

| hyperbolic dimensions | native ranks | required cycle lengths |
|---|---|---|
| `22,...,29` | `23,...,30` | `2,...,m` in rank `m` |

The exact remaining queue is:

1. enumerate connected frame graphs and complete compatible packets beyond
   the closed triangle, at-most-one-dotted four-node, and exactly-two-dotted
   four-node-tree local tranches;
2. implement the golden prime-ideal gauge and reconstruct each reflection
   representation;
3. compute the full reflection order, its conductor, and every stable lattice
   between the conductor bounds;
4. run a complete replayable Vinberg enumeration for every retained host;
5. canonicalize the hosts and certify the final list empty or exhibit a
   survivor.

There are presently zero universal complete packet manifests, zero universal
order/conductor manifests and zero complete host transcripts in these eight
ranks.  Those zeroes, rather than the existence of the finite bounds, are why
nonexistence remains open.

An empty maximal-completion list would suffice for all compact chambers in
this slice, but this bridge uses a specific standard argument rather than
mere terminology.  Fix a compact arithmetic side group `W` in one of the
eight dimensions, inside the one fixed ambient Lie group
`Isom(H^n)`.  No common integral lattice or preselected commensurator for its
overgroups is assumed.  For every discrete reflection group `R` containing
`W`, the covolume formula makes `R` a lattice of finite index and gives

\[
 [R:W]={\operatorname{vol}(W\backslash\mathbf H^n)
        \over\operatorname{vol}(R\backslash\mathbf H^n)}.           \tag{5}
\]

The Kazhdan--Margulis positive lower covolume bound in the fixed dimension
makes the set of these indices a bounded nonempty subset of the positive
integers.  Choose a member `R_max` of maximum index.  If it had a strict
reflection overgroup `R'`, then

\[
 [R':W]=[R':R_{\max}][R_{\max}:W]>[R_{\max}:W],                    \tag{6}
\]

a contradiction.  The other candidate overgroups need not be nested; the
maximum-index argument compares all of them through the fixed subgroup `W`.
The finite extension preserves dimension, cocompactness and arithmeticity.
The minimal Vinberg field also persists under finite index: Vinberg's least
adjoint ring/field and cyclic Cartan-product criterion identify it as a
commensurability invariant in this setting.  These two ingredients, including
the same max-index argument and its independent hostile audit, are written
out in
[`higher_dimensional_coxeter_orbifold_ag_quotient_no_go.md`](higher_dimensional_coxeter_orbifold_ag_quotient_no_go.md)
[and its hostile audit](higher_dimensional_coxeter_orbifold_ag_quotient_no_go_hostile_audit.md),
and
[`unconditional_maximal_arithmetic_reflection_chamber_finiteness.md`](unconditional_maximal_arithmetic_reflection_chamber_finiteness.md).

Consequently a hypothetical nonmaximal compact golden chamber in `H22--H29`
would force a maximal compact golden chamber in the same dimension and
field.

There are three sets here and they should not be confused.  The computational
candidate queue reconstructs **all full native hosts** `R(L)` arising from
retained frame/order/lattice packages, and every candidate must receive a
complete Vinberg classification.  The direct terminal set is the subset of
**compact reflective full-host survivors**.  That terminal set must be empty;
the other full hosts need not themselves be absent—they may be nonreflective,
noncompact, or have a non-finite-volume reflection chamber.

This direct route does not actually require maximalizing first.  For any
hypothetical compact `W`, the invariant-lattice construction already gives

\[
                  W\subset R(L),\qquad [R(L):W]<\infty.
\]

Thus `R(L)` would itself be a compact reflective full-host survivor.  Emptying
that survivor set excludes `W`.  Reflection-maximal completions provide a
useful equivalent cross-check: for a hypothetical `R_max`, the same
construction gives `R_max subset R(L)`, and maximality forces equality.  No
direct survivor-empty certificate or maximal-completion empty certificate has
yet been produced.  The missing step is computational emptiness, not a second
carrier theorem.

## Existing golden host lanes

Throwing the golden field back into the direct search gives four certified
but incomplete families.

* `'<-2phi> + D_n'`, `n=22,...,29`: norm-two controls through `phi^2` and
  exact current cones.  Every prefix has many hyperideal rays, but the
  nonunimodular mirror alphabet and later levels are incomplete.
* `'<-2phi> + N(A1^24)'`: the mirror norm class is complete.  There are 24
  walls at `x=1`, 1,035 at `x=phi`, two accepted `x=phi^2` orbits with
  18,216 walls, and an accepted 6,072-root cutter family at
  `x=phi^3=1+2phi`.  A replacement hyperideal extreme ray survives; the full
  `phi^3` packet and later levels remain open.
* `'<-a> + Leech'`: all 47 squarefree coefficient-ideal directions of norm at
  most 128 have complete local mirror classes and strict witnesses against
  their first necessary packets.  Realization, acceptance, nonunit square
  multiples, later levels and coefficients beyond the norm bound remain.
* `'<-a> + S23'`, where `S23` is the determinant-three odd-Leech section:
  for `a=phi` and `a=2phi`, the complete first packets have respectively
  553,104 and 53,682 walls.  Both exact prefixes retain strict hyperideal
  points; later packets remain open.

These are real golden arithmetic search examples, but none is a compact or
finite-volume survivor and none is an exhaustive field/rank decision.

## Reproduction

```console
PYTHONPATH=src python3 scripts/certify_golden_agcat_h22_h29_decision_manifest.py
PYTHONPATH=src python3 scripts/certify_golden_agcat_h22_h29_decision_manifest.py --verify
PYTHONPATH=src:. pytest -q tests/test_golden_agcat_h22_h29_decision_manifest.py
```

The focused suite reports `6 passed` and the frozen manifest rebuild reports
`PASS`.
