# Arithmetic20 internal-field translation: independent mathematical review

**Verdict: PASS at the stated finite-data and technical-supplement boundaries.**
Reviewed 21 October 2026, without importing or executing a mathematical checker,
Lean, or TeX. No further mathematical or exposition correction is required in
the reviewed insertion. A bounded build of this successor remains a separate
obligation; the old 17-page PDF is not evidence of its typesetting.

## Bound source and scope

The reviewed manuscript is
`manuscripts/compact_arithmetic_dimension_bound.tex`, SHA-256
`f8200c5dfbc76f8c9c7fa8969ecd33e879ff0f3d7909c00f9c4751110c4a6ce4`
(72,879 bytes, 1,441 lines). Its map is
`docs/compact_arithmetic20_proof_map.md`, SHA-256
`443023db750f72a145e40cfdd8de1e8dd2033a3fd84691f90e2988380be912c7`.
I read the complete new section, lines 1013–1259, and the field composition,
lines 1261–1335, against both original mathematical arguments:

* `results/fresh_arithmetic_high_label_20261016.md`, SHA-256
  `2ea56d166c4a29d8e0da397be99a927fe9849053838f7a67409b964c3dbe7157`;
* `results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.md`, SHA-256
  `549a1152c68a8a70f55f404b9b6c3df64116e172092e196b460763613d85d9c1`.

I also read the needed current face/projection and strict counting definitions
(lines 79–139), the exact triangle implementation's field, sign, square-root,
alphabet and triangle routines, and the retained TeX supplement's two-field
section. I did not author these field proofs or their TeX translation. Earlier
shared work on counting and other conditional kernels does not supply the new
field exclusions as premises in this review.

## Mathematical checks

**Actual field, rather than a defining overfield (1018–1033).** The determinant
expansion uses cyclic products, so it does not incorrectly place individual
unit-normal Gram entries in the ground field. Nonphysical embeddings extend
to the totally real entry field. A nonzero conjugate of the Lannér determinant
cannot vanish, and positive semidefiniteness therefore gives strict positivity.
A proper determinant subfield would admit another real embedding of the actual
field fixing its negative determinant. This proves the required equality.

**Arbitrary-white geometric reduction (1035–1081).** Both minimal pieces contain
their opposite black pair. Absence of a length-one or length-two marked corner
forbids every cross black contact; a shared white also gives such a corner.
An order-two piece would then be orthogonal to the other indefinite piece.
Thus only two triangles on distinct whites, joined by their white bond, remain.
Their opposite pairs are positive definite proper subdiagrams. After eliminating
the four black roots, costs greater than two leave two diagonals below minus
one and an unchanged off-diagonal of square below one. The resulting negative
definite plane contradicts Lorentz index one. Extra white roots cannot repair
that principal obstruction. The six-normal figure correctly suppresses the
opposite-pair internal products and allows only the indicated white contacts;
it is not a claim that every depicted support is realizable.

**Eight nongolden fields (1083–1134).** The three-node test excludes a label at
least six from a Lannér diagram of order at least four, including parallel or
dotted proper pairs. The remaining cosine entries lie in the biquadratic field,
contradicting the actual-field lemma. The label-seven two-cycle and the right
triangle determinant give exactly the remaining high-bond alternatives, including
the necessary dense 6/8 case. A distinguished bond in a finite corner is an
entire dihedral component. Both continuous projection-cost estimates are strict:
the cubic calculation gives cos²(pi/9)>7/8, while two contacts on a 6/8 bond
give cost at least 1/[2(1-c)]>2. The triangle-star argument supplies the negative
direction and forces the three black nodes themselves to be the minimal block.
At least two nonright black bonds pay the required two units. The tree degree
formula and the isolated-dihedral extra unit give at most 2r−2 per component,
hence at most 40 for every rank-21 forest. The strict 40f0+550 counting lemma
then contradicts this upper bound.

**Exact biquadratic field and finite proposition (1136–1223).** The published
arithmetic Lannér classification is explicitly retained for orders four/five;
it is not inferred from an arbitrary displayed entry field. Order two remains
the full unsampled algebraic-integer family and is excluded by the geometric
reduction before a finite triangle list is used. The cyclotomic degree and
group argument exhaust all finite labels and gives precisely the seven-letter
alphabet. The 84 unordered triples with repetition are therefore a complete
domain. The written successive quadratic-root algorithm includes zero mixed
coefficients, both norm-root signs, nonzero denominator conditions and exact
verification, with a fixed-depth rational-square base. This proves termination
and square-root-test completeness. Galois stabilizers test exact ground-field
generation; the displayed square-two determinant and the three conjugate signs
are the arithmetic Gram tests, with integral cyclic products and positive proper
minors. These are the same tests as the standalone rational implementation
`c2c2067d4269ba7fdb80029a5d2c193ea5d5ba9823866c87865449a6c7097024`.
The counts 13/37/28/4/2 and the two surviving triples remain an explicit finite
certificate output, not a newly executed 84-case computation in this review.

**Five cost positions and payments (1225–1259).** I independently expanded the
five costs from (d_p+d_q+P)/(4−d_m); all displayed expressions agree. Four exceed
two. In the remaining opposite-right position, the white's label-eight contact
and the compulsory nonright white bridge form a nonpositive three-node finite
corner, so the position is impossible. Thus the same principal Schur argument
and distance-two weights work without an ambient white-rank cutoff. Both listed
black triangles pay at least two and the same rank-21 upper bound is valid.

**Thirteen branches and higher dimensions (1261–1335).** The table covers the
rational field, two retained all-rank field theorems, eight distinguished-bond
fields, the exact biquadratic field and the golden field. The degree/Galois
checks establish both hypotheses for all eight fields; the coprime-conductor
argument handles K9 separately. The biquadratic theorem is never applied to its
golden proper subfield. The dimension-14 field list is used at 21; the shorter
restriction is used only at 22–29. The published ceiling 29 and the retained
dimension-specific technical supplement close the higher tail. The proof still
has result 20 and imports neither later bound-18 research nor parked refinements.

## Preservation, limits and readiness

The preservation receipt
`results/compact_arithmetic20_internal_field_translation_preservation_20261021.json`
has SHA-256 `66a7de5a9d8f306e9caed7baee53884d9f532a7df7164bf32b61638a2c4523d2`.
Its four source/PDF/build/map archives round-trip to the stated predecessors.
The intermediate 5a7de355 source is separately preserved in
`results/source_archive/compact_arithmetic20_before_positive_cost_arithmetic_correction_20261021.tex.gz`
(compressed SHA-256
`b573a8f27cb9fb5679697bbe58ccb813e99ad5aa61cdd7fb7220ac02159f115b`).
Its exact diff to the reviewed successor changes only the false printed sum
`1+.7+.4+.3=2` to the correct `1+7/10+2/5+3/10=12/5>2`.

This is a full mathematical review of the insertion and its composition,
not a new audit or regeneration of every inherited golden queue, field catalogue
or forest table. The Q(sqrt2), K7 and 22–29 arguments remain explicit inherited
technical-supplement theorems; published structural results and named exact finite
outputs remain inputs at their declared boundaries. The current map accurately
states these limits and the pending successor build. A targeted check of the
older E8 triangle-domain and consumer discussion (492–518, 746–786) found no
additional untranslated proof-critical argument required for this insertion;
its literal finite tables and auxiliary report are still visibly identified.
No end-to-end Lean, fresh checker execution, external refereeing or successor
PDF validation is claimed.
