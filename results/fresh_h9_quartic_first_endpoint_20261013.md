# Three certified nearest-mirror steps in the quartic host

**Subsequent result:** [the next-frontier report](fresh_h9_quartic_next_frontier_20261013.md) uses the first two complete steps below to connect exact chamber symmetries and prove that this full H5/H9 arithmetic host is nonreflective. Its full reflection chamber has infinite volume; further searches in these same lattices cannot produce a finite-volume polytope.

The complete interval census certifies the first two proposed edge steps and finds a **previously missed norm-$(2+T)$ mirror at the third edge**. Its four symmetry images retain all 182 preceding walls, producing a coherent alternative cone with **186 walls in H9 and 76 in H5**. Both cones still have exact strict spacelike escapes. No finite-volume or compact example is claimed.

Here $K=\mathbb Q(T)$, $T=2\cos(\pi/8)$, $T^4-4T^2+2=0$, and the unchanged lattice is

\[
L=\langle-(T-1)\rangle+F_4+F_4+\langle1\rangle.
\]

The [earlier report](fresh_h9_quartic_mixed_20261013.md) and [structural certificate](fresh_h9_glued_host_20261013_structure.json) establish the ring, signatures, unimodular Gram matrix and complete necessary primitive norm pool

\[
1,\quad 2+T,\quad T^2,\quad 2T^2-T^3,\quad 2.
\]

There is no new-lattice or new-group claim. The improvement is a complete first-mirror calculation, including the nonsquare residual classes omitted by the previous unit grids.

## Complete finite intervals

For a mirror of norm $q$ meeting a segment between time-normalized timelike points $k_0,k_1$, write its time coordinate as $x$ and put $\alpha=T-1$. Positive-spatial Cauchy–Schwarz and convexity give the exact bound

\[
\frac{\alpha x^2}{q}\le H,
\qquad H=\max_{j=0,1}\left(\frac{\alpha}{-k_j^2}-1\right).
\]

At the other three embeddings the entire lattice is positive definite, giving $x_\sigma^2\le-q_\sigma/\alpha_\sigma$. Integral coordinates in these four bounds form a finite alphabet. Enumeration uses the exact trace identity

\[
\operatorname{Tr}\big((a+bT+cT^2+dT^3)^2\big)
=4(a+2c)^2+8(b+3d)^2+8c^2+8d^2.
\]

The finite edge reflection group fixes the segment pointwise, so any crossing root can be folded into its dominant sector. Its four integral pairing magnitudes $p_i$ satisfy $p_i^2\le q e_i^2$ at every embedding. **Every** such bounded pairing is enumerated, not just Coxeter label values. The edge-span component is fixed by the inverse edge Gram. The remaining binary component has positive norm at every embedding; with each bounded time coordinate fixed, an exact quadratic equation enumerates all its solutions in K. No assumption that the residual norm is a square, or that one relative unit generates all orbits, enters this census. Norm-zero residuals are handled separately as mirrors containing the entire line. Every retained normal is primitive and integral, and every reflection coefficient is integral.

The producer also exhausts the omitted positive F4 block. For every norm in the pool its bounded dominant vectors with $S\le q$ at every embedding have only $S=0$ or $S=q$. A crossing root with a nonzero H5 component requires $S<q$ at every embedding, so it has zero omitted block. The equality case gives only the 24 pure-F4 mirrors containing the edge. This lifts the complete H5 census to **all H9 mirrors meeting the segment**.

Finite reflection-orbit comparisons verify the complete starting, edge and endpoint stabilizers. The incident roots form a positive definite Coxeter chamber for each complete stabilizer. Thus these are genuine local corners and adjacent edge steps of the full lattice reflection chamber, not merely finite vertices of a selected subgroup cone. This does not prove that every archived wall is a simple wall of that full chamber, nor does it complete its edge frontier.

| Step | H5 edge | Integral profiles | Nonsquare residual profiles | Quadratics | H9 mirrors: edge / start / next | Result |
|---|---|---:|---:|---:|---|---|
| 1 | $B_2+A_1+A_1$ | 23 | 3 | 2,229 | 30 / 33 / 41 | Old norm-2 proposal is first |
| 2 | $B_3+A_1$ | 8 | 2 | 1,206 | 34 / 41 / 49 | Old norm-2 proposal is first |
| 3 | $F_4$ | 5 | 2 | 7,105 | 48 / 49 / 49 | New norm-$(2+T)$ mirror is first |

The first two open segments contain no mirrors. The third contains exactly one, at affine parameter $T/2$ between the old starting corner and the old norm-one proposal. The three source-linked complete certificates are [first](fresh_h9_quartic_first_endpoint_20261013.json), [second](fresh_h9_quartic_first_endpoint_20261013_second.json), and [third](fresh_h9_quartic_first_endpoint_20261013_third.json). Each records exact bounds, complete alphabets, profiles, roots, source hashes, stage identity and stabilizer counts.

## The new mirror and corrected cone

In the H5 coordinates $(x;F_4;z)$, the new primitive normal is

\[
\begin{aligned}
x&=-77-42T+131T^2+71T^3,\\
F_4&=(-90-48T+152T^2+82T^3,\;-158-84T+267T^2+144T^3,\\
&\qquad -156-86T+269T^2+146T^3,\;-94-52T+162T^2+88T^3),\\
z&=-1-2T+T^2+T^3.
\end{aligned}
\]

Its norm is $q=2+T$ and its original-basepoint height is

\[
\alpha x^2/q=-12513-6772T+21361T^2+\frac{23121}{2}T^3
\approx120836.0712.
\]

The earlier norm-one proposal had height approximately 3,174,864.1373. It remains a valid alternative subgroup construction, but is not the first wall on this edge. The corrected time-one endpoint has norm

\[
482-643T-111T^2+172T^3<0.
\]

The [corrected 186-wall construction](fresh_h9_quartic_first_endpoint_20261013_completion.json) adds the four images under F4 reversal and block interchange to the preceding 182 walls. Its [separate Fraction audit](fresh_h9_quartic_first_endpoint_20261013_completion_independent.json) checks all 17,205 H9 pairs and 2,850 H5 pairs, the full norm and reflection-integrality gates, strict common interiors, actual finite vertices and exact strict escapes. A changed-coordinate fake root is rejected.

The displayed H5 escape has norm $-114+69T+203T^2-112T^3>0$. The original H9 escape survives, with norm

\[
-\frac{56876}{289}+\frac{34263}{289}T+
\frac{201861}{578}T^2-\frac{55639}{289}T^3>0.
\]

The strict perturbations and their future orientation are checked independently. These certify infinite volume of the displayed cones; they do not exclude completion by further lattice mirrors.

## Corrected frontier and replay

The corrected H5 corner has incident wall indices **[1, 2, 70, 72, 74]** in the corrected 76-wall artifact. Exact existing-wall lookahead is archived in [the frontier certificate](fresh_h9_quartic_first_endpoint_20261013_frontier.json). Releasing wall 74 returns along the just-certified edge to the preceding finite corner. The other four directions reach the light cone before any of the existing 76 walls:

- support [2, 70, 72, 74], releasing wall 1;
- support [1, 70, 72, 74], releasing wall 2;
- support [1, 2, 72, 74], releasing wall 70;
- support [1, 2, 70, 74], releasing wall 72.

These four edges are the next unresolved root-census tasks. Their starting point, directions and exact existing-wall intersections are saved; no new root search on them was performed this round.

[Allcock's edgewalking framework](https://arxiv.org/html/2110.03784v1#S9) motivates the first-mirror requirement. The new calculation proves that requirement for these three finite intervals over this quartic ring. It does not assert a generalization of the paper's integer-lattice algorithm or an exhaustive walk of the full chamber.

The independent implementation uses separate Fraction field arithmetic, a different affine free coordinate, an independently implemented complete quadratic-tower square test, and a separate F4 block census. All three [first](fresh_h9_quartic_first_endpoint_20261013_independent.json), [second](fresh_h9_quartic_first_endpoint_20261013_second_independent.json), and [third](fresh_h9_quartic_first_endpoint_20261013_third_independent.json) replays pass. Another agent reviewed the finite-interval bounds, folding and H9 lift without finding an issue.

```sh
python scripts/fresh_h9_quartic_first_endpoint_20261013_independent.py --stage 1 --verify
python scripts/fresh_h9_quartic_first_endpoint_20261013_independent.py --stage 2 --verify
python scripts/fresh_h9_quartic_first_endpoint_20261013_independent.py --stage 3 --verify
python scripts/fresh_h9_quartic_first_endpoint_20261013_completion_independent.py --verify
```

The producer is [first_endpoint](../scripts/fresh_h9_quartic_first_endpoint_20261013.py), with explicit `--stage 1`, `--stage 2`, or `--stage 3`. The corrected orbit [producer](../scripts/fresh_h9_quartic_first_endpoint_20261013_completion.py) and existing-wall [frontier calculation](../scripts/fresh_h9_quartic_first_endpoint_20261013_frontier.py) run with `PYTHONPATH=src:scripts`. There are no unfinished searches in this batch.
