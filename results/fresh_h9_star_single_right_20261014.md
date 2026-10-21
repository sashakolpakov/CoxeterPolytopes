The star 0 rank obstruction remains valid with **both cut leaves 15 and 16 unrestricted**. Their four own contacts may be right or any other finite Coxeter labels. In particular, this excludes the saved 415-visit state with facet 16's profile (2,5), without treating a decrease from label 3 to label 2 as monotone.

The hypotheses retained are:

- the positive leaves 2, 5 and 8 have both own labels at least 3;
- facet 14 has own labels exactly 3/3;
- all four internal block labels are at least 4.

The own profiles of the last anchors 11 and 17, all other meeting labels, and all dotted distances remain unrestricted. This is a label-family exclusion on the fixed star incidence, not a whole-incidence exclusion.

The [checker](../scripts/fresh_h9_star_single_right_20261014.py) and [certificate](fresh_h9_star_single_right_20261014.json) append one predicate to the frozen 160, giving 161 predicates. They bind the [previous dependency widening](fresh_h9_star_dense_frontier_20261014.md), all five primary arithmetic proofs and their separate independent reviews. No new determinant bank is required. The separate [independent dependency and symmetry review](fresh_h9_star_single_right_20261014_review.md) passed generation and saved replay, reproducing all 4,956 needed witness occurrences, 23 bound reductions and 161 predicate images.

The common rank core is

\[
S\cup\{11,14,17\},\qquad S=\{3,4,6,7,9,10,12,13\}.
\]

None of the removed pairs \((6,15),(7,15),(9,16),(10,16)\) belongs to this eleven-normal core. Every projection entry, determinant and polynomial in the 155,680 inherited rank cases is therefore unchanged. The argument does not relax the rank of the full Gram matrix.

The essential additional check concerns domain and vertex failures. All 4,956 needed stored actual-vertex witnesses have a literal support or a specified leading principal support containing none of those four pairs. Consequently their already verified nonpositive principal Grams are also unchanged. The checker reconstructs each support from the complete literal vertex list and tests this absence individually.

For the all-own-3 proof, only the budgets of the three core leaves 11, 14 and 17 enter the rank calculation. Its needed right/upper-contact guards are the 24 cross-block guards, the mutual 11–17 guard and the two third-contact upper guards. The three guards on unused leaf pairs and five unused leaf budgets can be discarded. That bank has no zero-determinant cases requiring further witnesses. For the other four proofs every stored witness is retained and checked: 96, 708, 479 and 3,643 occurrences respectively.

All 17 arithmetic high-block reductions and all six reductions in the last-block composition use lower-bound premises independent of the four removed pairs. Hence the complete last-internal-at-least 4 partition is unchanged. The remaining six noncore contacts belong to positive leaves 2, 5 and 8 and retain their lower bound 3; the earlier reviewed Perron-monotonicity argument still applies to increases of those labels. The internal lower bounds remain explicit rather than being inferred from a changed missing triple.

This proves a family with two entire unused cut-leaf profiles free. A right contact on a positive leaf, a changed facet 14 profile, or an internal label 2/3 remains outside this new predicate. The saved 415-visit cursor will be resumed only after the dependency certificate is reviewed, preserving all pending siblings.

A separate [literal symmetry binder](../scripts/fresh_h9_star_single_right_20261014_symmetry.py) checks the involution

\[
(0\ 1)(2\ 17)(3\ 12)(4\ 13)(5\ 16)(6\ 9)(7\ 10)(8\ 15)(11\ 14)
\]

against all 378 actual vertices and the complete minimal-nonface clutter. It reverses the four blocks and exchanges their positive and cut leaves. Every pair in every fixed-label or lower-bound predicate is mapped explicitly. Of the 161 images, twelve are additional distinct predicates, giving 173 in total. This transfer is an incidence relabeling, so it preserves rank, signature and the continuous dotted distances as well as finite labels; it introduces no label cutoff. No automorphism of another catalogue case is asserted.

The [first resumed batch](fresh_h9_star_single_right_20261014_search_batch_008.json) passed generation and saved replay: 786 new visits, 1,201 cumulative, nine pending frames. The unrestricted-cut predicate pruned 122 states. It stopped at internal labels (4,4,5,4), with leaves 8/16 carrying opposite single profiles (2,4)/(4,2), while every other own profile remains3/3. This state is outside the two-cut-leaf theorem because positive leaf8 now has a right contact.

A separate [finite-core domain refinement](fresh_h9_star_single_right_20261014_finite_core.md) addresses that entire new family. For all 16 internal4/5 choices it removes every noncore own restriction, using two new literal affineF4 witnesses inside the rank core and the existing nonzero rank bank. Its source and data are separate; it is not a retroactive change to the preceding frozen theorem or batch.

The final [composition certificate](fresh_h9_star_single_right_20261014_finite_core_gates.json) contains **205 predicates**: the 160 inherited predicates, the unrestricted two-cut-leaf theorem, twelve new reversal images, and sixteen finite-core families with sixteen reversal images. Both new proof scopes have independent source/domain reviews and saved replays; the symmetry compositions also pass generation and saved replay.

| Saved continuation | Predicates | New visits | Cumulative visits | Pending frames |
|---|---:|---:|---:|---:|
|[Batch008](fresh_h9_star_single_right_20261014_search_batch_008.json)|173|786|1,201|9|
|[Batch009](fresh_h9_star_single_right_20261014_finite_core_search_batch_009.json)|205|68|1,269|9|

Batch009 uses the [new continuation wrapper](../scripts/fresh_h9_star_single_right_20261014_finite_core_search.py) and resumes the exact 1,201-visit cursor. Twelve states were pruned by the finite-core predicates or their reversals. It preserves the next outside assignment:

- internal block labels are (4,4,5,>=7), where the last category means **every finite integer label at least 7**, with no finite cutoff;
- leaves 8 and 16 have opposite single profiles (2,4) and(4,2);
- leaves 11 and 17 have the shared single profile (2,3), right to root 12 and at angle pi/3 to root 13;
- all remaining own profiles are 3/3.

The positive leaf8 right contact violates the first new theorem's lower bound. The high last-internal category and its single own contacts violate the finite-core theorem's hypotheses. The exact saved assignment is therefore retained, not classified as excluded or as a metric realization. Its omitted meeting contacts, continuous dotted distances and Lorentzian rank/signature conditions remain unresolved. All pending DFS siblings are preserved; no prior traversal was restarted.

Both saved batch replays passed, reproducing their counts, rejection records, survivors and all pending cursors. Primary proofs, symmetry compositions, independent reviews, syntax and report links pass their checks. This phase is frozen; no jobs remain.
