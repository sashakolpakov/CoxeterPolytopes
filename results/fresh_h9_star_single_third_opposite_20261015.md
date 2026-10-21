# A single-contact third leaf with opposite single-contact anchors

The [checker](../scripts/fresh_h9_star_single_third_opposite_20261015.py) and [certificate](fresh_h9_star_single_third_opposite_20261015.json) exclude the star assignment saved at [26,951 visits](fresh_h9_star_single_third_20261014_search_batch_013.json) as part of a complete unbounded label family. The owner block has any finite label at least 4; the other three blocks have any finite labels at least 7. The third leaf has one right own contact and one own label **at least 5**. The two anchors have opposite single own contacts, both label 3. All noncore labels and dotted distances remain free.

This is an obstruction to a **compact Coxeter realization of that family**, not an exclusion of the entire star incidence. It requires no arithmetic-field assumption and no upper label cutoff. The proof uses 54 positive-kernel face guards, three projection budgets, and 16 strict-sign polynomial certificates with **165 positive and three negative coefficients**.

## Contact domains from actual faces

Normalize the roots to squared norm 1. Set
\[
S=(3,4;6,7;9,10;12,13).
\]
The third leaf is 14, with right contact to 3 and cosine \(\beta\ge\cos(\pi/5)\) at 4. Anchor 11 contacts only root 13 in the last block, while anchor 17 contacts only 12, both with cosine \(1/2\).

The 24 actual affine 4–3–4 path guards force the four blocks mutually right. A projection budget retaining root 4 forces the third leaf right to every foreign block: its retained own cost exceeds \(1/2\), while any foreign nonright contact costs at least \(1/2\). Each anchor has retained own cost \(1/4\), so it permits at most one foreign nonright contact, necessarily label 3. Actual affine 6–3 paths prohibit contact with either high middle block. If an anchor contacts the owner block, the same path forces that block's internal label to be 4 or 5.

The remaining restrictions are derived afresh for opposite anchors:

- The mutual anchor label is at most 5, by an actual 3–6 path at trial label 6.
- Each third-anchor label is at most 3, by an actual 4–4 path.
- Two owner foreign contacts are impossible. If they use the same owner root, an actual affine star with bonds 4,3,3 occurs. If they use different roots, an actual affine path with labels 3,3,4,3 occurs.
- Two third-anchor contacts are impossible, by the actual 4,3,3 star centered at the third leaf.
- A nonright mutual anchor contact cannot coexist with a third-anchor contact or an owner foreign contact. In each case an actual five-root branched affine diagram provides a positive null vector.
- If the unique foreign contact hits root 4, both third-anchor contacts are right. One is forced right by the actual all-3 triangle; the other is forced right by an actual 3,4,3,3 path.

All 54 kernels are strictly positive and their supports are checked against saved actual vertices. Raising labels or adding further acute contacts decreases the kernel quadratic form, so the relaxed affine diagrams contradict positive definiteness. The guards are conditional only on explicitly stated contacts; the old shared-root cycle arguments are not reused.

There are exactly **22 remaining discrete contact patterns**, reduced to **14 distinct residual formulas** by naming the unique foreign anchor first. Both anchor roles and their retained-root budgets are checked literally. With no foreign contact, the mutual label may be 2,3,4,5 when both third-anchor contacts are right; otherwise the mutual contact is right and exactly one third-anchor contact can have label 3. A foreign contact forces the mutual contact right. At root 3, it permits either no third-anchor contact or one of the two label-3 contacts; at root 4 it permits neither. Foreign cases use owner label exactly 4 or 5.

## The residual Gram matrix

Write \(c\) for the owner cosine, \(t\) for the last-block cosine,
\(q=1-c^2>0\), and \(D=1/(1-t^2)\). If there is no foreign contact, set \(F=T=0\). Otherwise
\[
F=\frac1{4q},\qquad
T=\begin{cases}\beta c/(2q)&\text{at root 3},\\
\beta/(2q)&\text{at root 4}.
\end{cases}
\]
Let \(\mu\) be the mutual anchor cosine and \(\delta_x,\delta_y\) the third-anchor cosines. Projection off the positive eight-root support gives
\[
R=\begin{pmatrix}
1-D/4-F&-\mu-tD/4&-\delta_x-T\\
-\mu-tD/4&1-D/4&-\delta_y\\
-\delta_x-T&-\delta_y&1-\beta^2/q
\end{pmatrix}.
\]
The factor \(t\) in the off-diagonal entry is essential: the anchors use opposite last-block roots. The shared-anchor matrix would be incorrect here.

The checker forms this matrix for each domain and clears the positive denominator:
\[
P=q(1-t^2)\det R.
\]
The resulting expression is polynomial and has degree at most two in \(t\). For the ordinary regions it uses nonnegative shifts
\[
\beta=(1+\sqrt5)/4+x,\qquad
c=1/\sqrt2+y
\]
when the owner label is unrestricted; a foreign contact instead fixes \(c\) to the exact label-4 or label-5 cosine. Each coefficient of a monomial in \(x,y\) is expressed in the degree-two Bernstein basis on the stated interval for \(t\).

## The exceptional integer-label split

Thirteen of the fourteen formulas have strictly positive coefficients on \(9/10\le t\le1\). One formula changes sign within that continuous interval: owner label 4, foreign contact at root 3, \(\delta_x=0\), \(\delta_y=1/2\), and \(\mu=0\). At third label 5 it is
\[
P(t)=\frac{3(1+\sqrt5)}{32}t^2
-\frac{\sqrt{10}+\sqrt2}{64}t
-\frac{5+3\sqrt5}{64}.
\]
A blanket continuous positivity claim would be false. The finite integer labels instead give a complete three-way split:

| Third label | Last label | Exact parameter region | Sign |
|---|---|---|---|
| at least 6 | at least 7 | \(\beta\ge6/7,\ 9/10\le t\le1\) | positive |
| 5 | at least 8 | \(\beta=\cos(\pi/5),\ 23/25\le t\le1\) | positive |
| 5 | 7 | \(\beta=\cos(\pi/5),\ 9/10\le t\le901/1000\) | negative |

These regions are rigorous. With \(u=\cos(2\pi/7)>1/2\), the identity
\[
T_7(u)-1=(u-1)(8u^3+4u^2-4u-1)^2
\]
and strict monotonicity of the cubic on \([1/2,1]\) bound its root between \(31/50\) and \(311801/500000\). Thus
\(9/10<\cos(\pi/7)<901/1000\).
The other bounds follow from
\(\cos^2(\pi/8)=(2+\sqrt2)/4>(23/25)^2\) and
\(\cos^2(\pi/6)=3/4>(6/7)^2\).

The split gives **16 certificates in total**. Exact rational-square comparisons in \(\mathbb Q(\sqrt2,\sqrt5)\) verify 165 strictly positive Bernstein coefficients and three strictly negative ones. Each region has a nonzero constant-in-\(x,y\) coefficient at every Bernstein index, so its determinant sign is strict throughout. No case is dropped at a closed endpoint.

In every possible integer-label case, \(\det R\ne0\). But three projected normals lie in the two-dimensional orthogonal complement of the actual PD8 support in a ten-dimensional Lorentz space. Their Gram determinant must vanish. This contradiction proves the stated compact-family obstruction.

## Verification and continuation

The frozen predecessor archive and batch 013 were freshly replayed before this work; both passed. Primary generation and saved replay pass. The [independent review](fresh_h9_star_single_third_opposite_20261015_review.md) also passes generation and saved replay. Its separate `Fraction` implementation reconstructs all 54 conditional kernels, three budgets, 22-to-14 pattern reduction, 16 polynomial regions, 165 positive and three negative coefficients, and the exact cosine bounds, without producer or SymPy arithmetic imports.

```sh
python scripts/fresh_h9_star_single_third_opposite_20261015.py --verify
```

The owner and last-root reversals are checked against the complete vertex and minimal-nonface lists; they yield four source predicates. This proof supplies no new polytope and does not claim the whole star incidence is closed.

## Literal transfers and resumed traversal

The separate [transfer checker](../scripts/fresh_h9_star_single_third_opposite_20261015_gates.py) and [certificate](fresh_h9_star_single_third_opposite_20261015_gates.json) require the independent theorem review and the frozen [twelve-core transfer review](fresh_h9_star_core_transfers_20261014_review.md). They recheck each complete induced core and all **684 new path, star, fork and budget support occurrences**. The independent reviewer separately reconstructs these twelve maps, their support occurrences and the resulting predicate set. No whole-nerve symmetry is assumed for the local maps.

Four source orientations through twelve placements add 48 predicates, bringing the total from 349 to **397**. Binder generation and saved replay pass. All original files and the batch 013 cursor remain unchanged.

The [resumable search](../scripts/fresh_h9_star_single_third_opposite_20261015_search.py) resumed batch 014 from 26,951 visits with all nine pending frames and their later siblings. It stopped at the first outside necessary-label survivor after **33,064 new visits in 138.455 seconds**, using less than its 100,000-visit allowance. The [saved batch](fresh_h9_star_single_third_opposite_20261015_search_batch_014.json) records **60,015 cumulative visits**, nine pending frames, and 161 applications of the new opposite-single-third predicates. Visits are traversal visits, not a count of distinct complete assignments. The saved replay passes with exactly the same traversal and stopping state.

The next state has block labels `(5, >=7, >=7, >=7)`. Leaf 2 has own contacts `(3,3)`; leaf 14 has `(2,4)`. The remaining three blocks each have opposite single contacts: leaves 5/15, 8/16, and 11/17 have `(2,3)` and `(3,2)`, respectively. Category `>=7` retains every finite integer label in that range. This lies outside the new theorem because its third nonright contact is 4, not at least 5.

The state has not been proved metrically feasible. Its unassigned meeting contacts and continuous dotted distances still need a Lorentzian rank/signature completion. A concrete next task is to analyze the same residual matrix at owner label 5 and third label 4, deriving fresh literal contact domains before using any polynomial sign test. The whole star incidence remains open; no new search or gate is included beyond this saved frontier.
