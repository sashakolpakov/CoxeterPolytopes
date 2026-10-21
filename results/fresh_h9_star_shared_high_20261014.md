# Shared-root high-label family in star case0

**The shared single-contact family has no compact H9 Coxeter realization for any finite last-block label.** The first three blocks retain both leaves with own contacts `(3,3)`. Leaves11 and17 are both right to12 and both have π/3 contact to13. Their own missing triples `(2,3,m)` force the last internal label `m>=7`. The literal swap of12 and13 covers the other orientation. Unlisted finite labels and dotted distances remain free.

The [exact checker](../scripts/fresh_h9_star_shared_high_20261014.py) and [certificate](fresh_h9_star_shared_high_20261014.json) cover7,840 residual cases:7,392 have a strictly positive determinant throughout the continuous high-label range, and448 have exact nonpositive actual-vertex witnesses. This excludes the [saved103-visit necessary assignment](fresh_h9_star_lastblock_20261014_search_batch_004.json), not the whole incidence.

The reviewed two-dense gate excludes any first-three internal label at least6. Their own `(3,3)` missing triples force labels at least4, so these three labels reduce to4 or5. Actual four-root supports force all cross-block contacts right by the relaxed affine path4–3–4. The eight block roots therefore have a positive direct-sum Gram.

An actual finite support retains root13 and all six foreign roots for each anchor. Its own projection cost is1/4. Every foreign nonright contact contributes at least1/2, because each foreign inverse diagonal is at least2. Thus each anchor can have at most one foreign contact, necessarily π/3: seven possibilities. The actual triangle `{11,13,17}` forces the anchors mutually right.

Leaf14 has the same retained-own-root budget, but the actual triangle `{12,13,14}` contains the high internal bond. A second nonright bond would give a three-node path with determinant strictly below `1−3/4−1/4=0`. Consequently14 is right to both last roots and has only five foreign patterns: none, or one π/3 contact to one of the four middle-block roots. Two fresh exact actual-vertex witnesses exclude a label4 between14 and either anchor. All other unspecified contacts in these witnesses are right. The last internal label is relaxed to4, so the failures persist for every larger label.

Put
\[
D=\frac{1}{\sin^2(\pi/m)}>4,
\]
since `m>=7` implies `cos(pi/m)>cos(pi/6)=sqrt3/2`. No estimate for `cos(pi/7)` or upper label bound is needed.

After projection off the first six block roots, write the residual Gram of11,17,14 as
\[
R(D)=\begin{pmatrix}
a-D/4&b-D/4&e\\
b-D/4&d-D/4&f\\
e&f&q
\end{pmatrix}.
\]
The last-block contribution is rank one because both anchors contact the same root and14 contacts neither. Hence
\[
\det R(D)=B+AD,\qquad
A=\frac{(e-f)^2-(a+d-2b)q}{4},\quad
B=\det\begin{pmatrix}a&b&e\\b&d&f\\e&f&q\end{pmatrix}.
\]
All coefficients lie in `Q(sqrt2,sqrt5)`; the variable D is allowed to be any real number greater than4 and need not belong to that field.

The complete domain has `8*7*7*5*4=7,840` cases: eight internal choices, seven patterns for each anchor, five for14 and four remaining contact choices. Every slope A is exactly positive. In7,392 cases the exact boundary value `B+4A` is positive, so the determinant stays positive for all D>4. Each remaining case has a saved nonpositive actual-vertex LDL pivot with the last internal label relaxed to4. These448 witnesses are verified directly in the exact field; the finite-diagram classifier only helps locate them.

A nonzero three-by-three residual determinant, together with the positive eight-root block, gives rank11 in a ten-dimensional Lorentzian space. The alternative actual-vertex failures contradict compactness directly. No zero-polynomial branch or high-label tail remains in this specified family. This is a continuous rank/vertex argument, with no cyclotomic field-membership restriction.

Generation and saved replay are run with:

```sh
python scripts/fresh_h9_star_shared_high_20261014.py --verify
```

The primary uses a frozen exact field helper; independent arithmetic scope is stated separately in its review.

## Saved continuation

The [gate binder](../scripts/fresh_h9_star_shared_high_20261014_gates.py) appends the two literal root orientations, giving158 proven predicates. The [search](../scripts/fresh_h9_star_shared_high_20261014_search.py) resumes the103-visit cursor without replaying its explored prefix. [Batch005](fresh_h9_star_shared_high_20261014_search_batch_005.json) stops after two additional visits,105 cumulative, with nine pending frames.

The next necessary assignment has opposite single own π/3 contacts: leaf11 has profile `(2,3)` and leaf17 has `(3,2)`, with last internal label `m>=7`. The first three blocks remain B2 with all own contacts3. This state lies outside the present shared-root theorem and outside the earlier opposite-single theorem, whose nonright own labels must be at least4. No metric realization or whole-incidence exclusion is claimed.

The binder and resumed batch have separate saved replays:

```sh
python scripts/fresh_h9_star_shared_high_20261014_gates.py --verify
python scripts/fresh_h9_star_shared_high_20261014_search.py --batch 5 --verify
```

The [independent review](fresh_h9_star_shared_high_20261014_review.md) passes generation and saved replay. It separately reconstructs the exact residual arithmetic, the affine determinant identity and continuous interval bounds, and the literal domain/actual-vertex witnesses. Its report distinguishes this arithmetic audit from the primary's reuse of the frozen field helper.
