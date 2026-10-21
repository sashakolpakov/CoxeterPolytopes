# Unconditional compact dimension29: exact savings targets

The current unconditional compact rule has maximum87 at rank29, while a valid universal vertex bound84 would exclude dimension29. This packet makes the required reductions explicit: **3/4 per component of ranks4,7,8,9 and1/8 per rank10 component** would suffice. The largest proposed component reduction,3/4, is the smallest possible for a nonnegative rank-only scalar target. These reductions are **not a proved face weighting**. The compact dimension bound remains29.

The [producer](../scripts/fresh_compact29_savings_frontier_20261018.py), [saved result](fresh_compact29_savings_frontier_20261018.json), and [standalone scalar checker](../scripts/check_compact29_savings_saved_20261018.py) use no arithmetic-field hypothesis. They preserve the [frozen original dimension29 reduction](fresh_compact_dimension_bound_20261015.md) and its [independent review](fresh_compact_dimension_bound_20261015_review.md). The original2,126-forest classification and79 local completions are inherited results, not new discoveries here.

## Complete scalar domain

For the published unconditional Vinberg rule, the weight of an irreducible finite component is

\[
w(r)=\begin{cases}
r(r-1)&1\le r\le4,\\
\binom r2&5\le r\le7,\\
24&r=8,\\27&r=9,\\r+19&r\ge10.
\end{cases}
\]

The original packet establishes these weights and their face inequalities; the new diagnostic does not reprove either. Every one of the4,565 rank partitions of29 is checked. The only forests with weight above84 use ranks4,7,8,9,10. Their complete finite-type lists are A/B/D/F/H at rank4, A/B/D/E at ranks7 and8, and A/B/D at ranks9 and10. Repeated equal-rank factors are unordered. In particular no restriction on dihedral labels has been inserted: every finite rank2 factor has deficit4 from3r and so cannot occur in a forest above84.

Write \(x_r\ge0\) for a hypothetical saving on each whole rank\(r\) component. The exact necessary and sufficient scalar conditions for lowering all these forests to84 are:

| Component ranks | Old weight | Typed forests | Required saving |
|---|---:|---:|---|
|4+4+4+4+4+9|87|378|\(5x_4+x_9\ge3\)|
|4+4+4+7+10|86|420|\(3x_4+x_7+x_{10}\ge2\)|
|4+4+4+8+9|87|420|\(3x_4+x_8+x_9\ge3\)|
|4+4+7+7+7|87|300|\(2x_4+3x_7\ge3\)|
|4+7+8+10|86|240|\(x_4+x_7+x_8+x_{10}\ge2\)|
|4+7+9+9|87|120|\(x_4+x_7+2x_9\ge3\)|
|4+8+8+9|87|150|\(x_4+2x_8+x_9\ge3\)|
|7+7+7+8|87|80|\(3x_7+x_8\ge3\)|
|9+10+10|85|18|\(x_9+2x_{10}\ge1\)|

There are1,448 typed forests of weight87,660 of weight86, and18 of weight85, totaling2,126. The maximum87 occurs in six rank patterns. Every saved typed forest is compared literally with the original complete list.

Setting \(x_4=x_7=x_8=x_9=3/4\), \(x_{10}=1/8\), and all other savings zero satisfies every row. Forests initially at most84 remain at most84 because the proposed savings are nonnegative. Thus the conditional rank-only maximum is exactly84. If all proposed savings were at most\(t\), the eighth row would force \(4t\ge3\). Hence \(t\ge3/4\), attained by the stated vector. At that cap, the last row forces \(x_{10}\ge(1-3/4)/2=1/8\), also attained. This is a complete elementary minimax calculation, not an optimizer output or a geometric existence argument for the rule.

## Why this is the relevant threshold

For a hypothetical compact29-polytope, let \(\sigma(v)\) be the original weight and \(s(F)=\sigma(F)-(5-k_F)\ge0\) its face slack. The frozen derivation gives

\[
\sum_v(\sigma(v)-84)
=5\sum_{i=0}^{14}(14-i)(15-i)h_i+\sum_Fs(F)\ge1050.
\]

Here \(h_0=1\), all\(h_i\ge0\), and Dehn–Sommerville symmetry has already been used. Since each positive vertex excess is at most3, at least350 vertices must have weight above84. Setting weights on faces with at least five sides to zero retains the face inequalities, so the same necessary count applies to small-face weight above84. A new valid face rule with vertex bound84 would contradict the displayed positive balance. This packet supplies no such rule.

## Concrete next step and the existing local obstruction

The [frozen small-face packet](fresh_compact29_small_faces_20261015.md), with its [independent review](fresh_compact29_small_faces_20261015_review.md), gives exact separate local completions for all21 E7 pairs,28 E8 pairs, and30 rank4 pairs. In particular every positive-weight angle in \(E_7+E_7+E_7+E_8\), of total weight87, has an individual triangle or quadrilateral completion with the required finite corners, Lorentz inertia, and no parabolic subdiagram. Those79 realizations are not simultaneous completions.

Therefore the next useful unconditional compact calculation must couple multiple face stars at one actual vertex, or use a global incidence/charge argument. A concrete first target is the \(E_7^3+E_8\) forest: select two positive-angle stars sharing original support, retain all original rank29 roots, and test the mutual products of the newly introduced normals in the same Lorentz space. Independent one-star success cannot establish their joint compatibility; a finite collection of failed selected joins cannot establish universal incompatibility. Any proposed exclusion must cover all valid stars, including continuous dotted contacts, or be stated only for its declared finite subdomain. This packet stages that target without claiming the existing79 witnesses are exhaustive.

Vinberg's published obstruction to uniform nonnegative diagram-local weights with constant below3 is inherited through the original report. It does not exclude rank29-specific full-forest charges or coupled incidences. Arithmetic field restrictions used by the separate compact-arithmetic21 proof are not available in this general lane.

The compactH9 construction lane is separate: its current exact66-facet seed remains infinite-volume, and the four-finite-duty single-cap family is excluded by its own pending-review packet. A split cap with changed incidence is the next construction direction; no construction result is used in this dimension diagnostic.

## Reproduction and verification scope

```sh
python3 scripts/fresh_compact29_savings_frontier_20261018.py --verify
python3 scripts/check_compact29_savings_saved_20261018.py
```

The producer checks all rank partitions, typed forests, source bindings, and rational scalar inequalities. The standalone checker independently enumerates rank multiplicities, reconstructs the complete typed domain, checks every saved row and both minimax witnesses, and verifies the explicit no-new-theorem flags. Neither program reconstructs the published face weights, geometric79-star witnesses, or the original h-vector proof. Generation374e87, separate fresh replay aec311, and saved-checker execution15e4ea all terminated with PASS. Metadata check61ca34 passed9 source bindings,7 local links and2 Python ASTs.

Root completed a full mathematical producer/report/standalone-checker read: the complete4,565 rank partitions and2,126 typed forests, nine inequalities, minimax3/4 and rank10 1/8 targets, strict balance and350-heavy-vertex interpretation, and need for coupled-star information are coherent. The reviewer did not execute a program or reprove the inherited79 geometric witnesses or h-vector argument. The packet is final/frozen; its unproved scalar target is explicitly separate from any new compact dimension theorem.
