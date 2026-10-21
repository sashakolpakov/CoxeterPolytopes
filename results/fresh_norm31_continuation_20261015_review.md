# Independent review of the discriminant-31 nonreflectivity witness

**PASS for this literal lattice.** The [separate checker](../scripts/fresh_norm31_continuation_20261015_review.py)
reconstructs the exact field arithmetic, all 28 wall and facet witnesses,
the complete 11 finite corners, and the infinite-order integral chamber map
in the [primary report](fresh_norm31_continuation_20261015.md). Its
[certificate](fresh_norm31_continuation_20261015_review.json) binds the frozen
inputs and their transitive manifests. It imports no producer code.

The conclusion is nonreflectivity of

    E8 ⊥ [[−2√2,−1], [−1,2]] over Z[√2].

It excludes a finite-volume Coxeter chamber whose reflections all preserve
this lattice. It is not a dimension bound or a claim about every lattice
in its commensurability class. No new polytope was obtained.

The independent calculation uses rational pairs for Q(√2), exact signs,
Gaussian determinants, and positive-pivot elimination. It verifies the
physical determinant −1−4√2 and positive definiteness at the other real
embedding. Every saved root has norm 2 or 18+10√2 and satisfies integral
reflection coefficients. All pairwise angles, strict center inequalities,
and explicit relative interior facet points are checked directly.

For completeness of the finite-corner census, this checker builds the
finite-pair graph and enumerates its maximal cliques, then every 9-subset
of those cliques. This differs from the producer's fixed-size recursion.
Exactly 35 candidate supports occur; full positive definiteness leaves
the same 11 finite vertices. Their saved rays satisfy all wall inequalities
with precisely nine equalities. The ten finite edges are reproduced.
An ideal or hyperideal census is not needed for the obstruction.

The matrix and inverse satisfy every entry of TᵀGT=G and TT⁻¹=T⁻¹T=I,
and all entries lie in Z[√2]. The determinant is −1 and the trace is
9+5√2>10. The checker verifies all nine outward root images, the positive
apex scale 18+14√2, both corners' complete sign lists, and both exact
basepoint decompositions with positive apex and nonpositive wall coefficients.
Thus the claim uses the actual embedded sectors, not merely isomorphic diagrams.

The mathematical review of the enumeration finds its completeness argument
sound. Evenness gives q=2b for a primitive reflector. Reflection integrality
and the adjugate show b divides the prime discriminant ideal. The PID and
totally positive unit-square arguments leave exactly the two norm classes.
The zero sector is completely determined by the rational even E8+A1 shell.
Both-embedding bounds give finite coefficient boxes for all controls and
all nonnegative sector coordinates. The equation, parity, and long-root
divisibility conditions are exact. Prefix Schur minima are necessary lower
bounds at both embeddings, and prior-wall pruning has nonnegative remaining
physical coefficients. Hence these prunes cannot discard a simple wall.
The full Vinberg search is replayed separately; this review does not claim
an independently reimplemented root enumeration.

For a finite corner P of the genuine prefix, the verified decomposition
H=aP+Σbᵢnᵢ has a>0 and bᵢ≤0. A later simple wall r has (H,r)<0 and
(nᵢ,r)≤0, which force (P,r)<0. Thus both complete corner sectors persist
in the full reflection chamber. The integral future isometry permutes
all lattice mirrors and maps one complete outward sector to the other,
so it stabilizes that full chamber. Its trace proves infinite order.

The admissible integral isometry group is discrete; the positive definite
conjugate form supplies compactness at the second embedding. Distinct powers
of the infinite-order chamber stabilizer give infinitely many disjoint
translates of a sufficiently small interior ball. The chamber therefore
has infinite volume. Equivalently, the full reflection subgroup has infinite
index in the integral isometry group. This proves the claimed lattice
obstruction without inferring nonreflectivity from an unfinished search.

Replay: `python3 scripts/fresh_norm31_continuation_20261015_review.py --verify`.
