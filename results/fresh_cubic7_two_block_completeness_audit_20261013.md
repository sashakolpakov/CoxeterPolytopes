# Completeness review of the thin-alphabet two-block search

The reviewed enumeration is complete for its stated support restriction and height bound. No omission was found in the norm alphabet, coefficient bounds, zero/equal-block cases, or block placements. A separate height-30 calibration agrees with the full four-block `Host.candidates` enumerator on **all 20 scalar-zero roots supported on at most two blocks**: eight have one nonzero block and twelve have two.

This is a mathematical code review plus an enumeration calibration. It reuses the project cubic arithmetic, independently brute-forces coefficient boxes without the thin enumerator's interval pruning, and compares two different block-joining algorithms. It does not claim an independent height-1000 census or geometry replay. Those saved geometry certificates are checked by the separate Fraction verifier.

## Exact coefficient bounds

The three `DUAL` elements in [the thin enumerator](../scripts/fresh_cubic7_target_alphabet_20261013.py) are exactly trace-dual to `1,τ,τ²`; the calibration checks the full identity matrix `Tr(d_i τ^j)=δ_ij`. Therefore an integral element `x=a+bτ+cτ²` satisfies

\[
 |a_i|=|\operatorname{Tr}(d_i x)|
 \leq\sum_e |d_i^{(e)}|\,|x^{(e)}|.
\]

The saved rational bounds strictly enclose all nine conjugate absolute values of the dual basis. The power-of-two embedding bounds enclose the square roots of the requested caps, so the resulting integer coefficient box is complete. The subsequent intervals for `bτ+cτ²` can overestimate its range but cannot remove a valid integer constant coefficient. Every retained element is checked against the original caps at all three embeddings.

The calibration verifies the three exact isolating intervals, all dual inequalities, and six thin-alphabet comparisons: zero caps, a negative cap, and both time/weight alphabets at height 30 for norms `2` and `2β`. Its independent oracle exhausts every coefficient of a certified box and performs no interval pruning.

## Why the scalar-zero norm alphabet is exactly 2 and 2β

The [existing full-host audit](fresh_cubic7_dihedral_independent_audit_20261013.md) proves that primitive reflection norms can be normalized to the safe pool `1,2,4,β,2β,4β` for `α=2τ`. The ring is `Z[τ]`, the field has class number one, and every totally positive unit is a square. Thus the normalization loses no primitive mirror.

With the scalar coordinate zero, the quadratic form is even, so norms `1` and `β` cannot occur. For `q=4` or `4β`, reflection integrality forces every block coordinate to be divisible by 2: the block matrix is invertible modulo 2 because its determinant is `β`, a unit modulo the inert prime 2. The norm equation then forces the time coordinate to be divisible by 2. Dividing the vector by 2 would give a scalar-zero vector of norm `1` or `β`, again contradicting evenness. The calibration directly checks that squaring has no nonzero kernel in the eight-element residue field.

The remaining ideals `(2)` and `(2β)` are squarefree, so an integral vector of either norm is automatically primitive. A common nonunit divisor of its coordinates would force the squared divisor ideal to divide the norm ideal.

## The chamber, energy bounds, and support coverage

Write `u,w` for the covector coordinates in one block. Its two original zero-height walls impose `-u≤0` and `D(u-w)≤0`, with `D>0`, hence precisely `0≤u≤w`. The complete zero-height chamber was independently established in the existing audit. Equality and the zero block must be included; the reviewed [two-block producer](../scripts/fresh_cubic7_two_block_search_20261013.py) includes both.

The inverse block formula

\[
 y=((2u+Cw)/\beta,(Cu+2w)/\beta)
\]

is exact because `4-C²=β`. Every block is positive definite at every real embedding. Thus each block energy is nonnegative, and the norm equation implies energy bounds `q(H+1)` at the distinguished embedding and `q` at the other embeddings. Cauchy's inequality gives `u²,w²≤2E`, exactly the producer's weight caps. The same norm equation yields the time caps `Hq/α` and `-q/α` at the other two embeddings. All coordinate and reflection-integrality conditions are applied after these complete bounds.

The energy join enumerates every pair of allowed blocks. Twelve ordered placements in four block positions cover every support of size two, while the included zero block covers supports of size one and zero. Duplicate placements are harmless because saved roots use exact ray keys. No exchange within a block is omitted: the chamber inequality has already selected the allowed orientation; the coordinate exchange is one of the original chamber reflections.

At height 30, the calibration independently joins a Cartesian product of two block lists, then compares it with `Host.candidates`, which joins four sorted block multisets and the scalar coordinate. It runs all six full-host norm classes and then restricts to zero scalar and at most two nonzero blocks. The sets coincide exactly, including all one-block cases. No extra roots occur in the discarded norm classes.

The height-1000 search's subsequent cross-compatible selection is greedy. Completeness of the input dictionary does not mean that every possible compatible subset was tested, nor that roots with three/four blocks or nonzero scalar were included. The separately appended norm-one Pell cap is explicitly outside this bounded dictionary.

## Reproduction

The [calibration script](../scripts/fresh_cubic7_two_block_completeness_audit_20261013.py) writes the [PASS record](fresh_cubic7_two_block_completeness_audit_20261013.json), including exact alphabet caps and counts, per-norm control counts, and the matched roots.

```sh
PYTHONPATH=src:. python scripts/fresh_cubic7_two_block_completeness_audit_20261013.py
```

The completed run took about three seconds. The [two-block search report](fresh_cubic7_two_block_report_20261013.md) supplies the larger saved dictionaries and their exact infinite-volume conclusions.

## Follow-up source review: equal/free blocks with the scalar coordinate

The same review was extended, without rerunning the searches, to [the equal-block/scalar producer](../scripts/fresh_cubic7_scalar_search_20261013.py) and [the mixed-block/scalar producer](../scripts/fresh_cubic7_mixed_scalar_search_20261013.py). Their enumeration logic is complete for their declared slices and height bounds. This statement concerns input enumeration; it does not assert that greedy selection covers every compatible subset or that any resulting cone has finite volume.

An equal-coordinate block has vector `(a/D,a/D)`, covector `(Da,Da)`, and norm `2a²`. Since `D` is an integral unit, integral `a` parameterizes the entire integral equal-coordinate line, and this line is saturated. The slice consisting of that block, time and scalar has norm `2a²+z²−2τx²`. For normalized norms `β` or `2β`, ambient reflection integrality forces each of `x,a,z` to be divisible by `β`: the factors `2`, `τ` and `D` are units at the prime `(β)`. This contradicts valuation one of the proposed norm at that prime. Together with the existing exclusion of primitive norms `4,4β`, it leaves precisely the norm classes `1,2`. For these classes the omitted separate reflection-divisibility filters are automatic for integral `x,a,z`; the producer also asserts full reflection integrality on every reconstructed root.

The equal-block/scalar bounds are `a²≤q(H+1)/2` and `z²≤q(H+1)` at the distinguished place, and `a²≤q/2`, `z²≤q` elsewhere. The time bounds are unchanged. Every term in the spatial energy is nonnegative at every embedding, so these inequalities exclude no solution. The square lookup for `z` loses no solutions because `z≥0` at the distinguished place makes `z↦z²` injective. Zero equal blocks and zero scalar coordinates are included. All four placements of the equal block are tested, including the pure-scalar case.

The mixed slice adds one free dihedral block and correctly restores all four primitive norm classes `1,2,β,2β`: divisibility of a free block's covector by `β` need not force the block coordinates themselves to be divisible by `β`. Its separate divisibility filters are exactly those for time, equal block, scalar and free-block covector coordinates. The free block again includes every `0≤u≤w`, including equality and zero. The energy dictionary maps each value `2a²+z²` to a **list of all** `(a,z)` pairs, so collisions do not discard distinct solutions. The exact join with the free-block energy solves the full norm equation. Twelve ordered placements cover every declared support, including degenerations with either block zero or equal. The bounds are the same positive-energy and Cauchy bounds proved above.

The replacement-dictionary provenance is also correct. The scalar producer imports `compatible_roots` from the earlier original-25 two-block replacement artifact, not that artifact's greedily selected wall list. It joins this full dictionary with all its own original-25-compatible scalar roots. All **283** saved union rows were checked directly against their `origin` and `origin_index`, including exact vector, norm and height equality. The mixed producer in turn imports that entire `combined_replacement_dictionary` and adds every newly enumerated original-25-compatible mixed root; it does not condition the union on compatibility with the intervening larger cones. Ray deduplication is exact, and sorting by exact height precedes each greedy rebuild from the original 25 walls.

Conditional completions of previous cones are separately filtered by their saved compatibility flags. Their `selected_indices` index the filtered candidate list, whose mapping back to the main dictionary is saved as `candidate_indices`. The original-25 rebuild indexes its complete union directly. These distinct index conventions are consistent with the saved provenance. The Pell cap is a separate, explicitly tested addition; a failed cross-angle gate does not discard the earlier scalar candidates or alter the completeness of the underlying dictionaries.

This follow-up is a source-level proof and saved scalar-provenance check, not a second run of either height-1000 census. The mixed artifact was still being produced when the source review finished; independent saved-root, angle and geometry replay is handled separately.
