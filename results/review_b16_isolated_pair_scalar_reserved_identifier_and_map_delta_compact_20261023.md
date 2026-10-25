# B16 scalar kernel: different-author repair-delta review

23 October2026. Reviewer: finite57_composition_audit. **STATIC DELTA
PASS, no remaining correction.** No compiler, import, test, worker,
enumeration or mathematical job was launched for this review.

Root authored the implementation and this repair. I did not author
either. I contributed/reviewed earlier written B16 and actual-point
inputs; those mathematical interfaces remain external and are not
represented as independently formalized geometry by this report.

## Exact successor and preserved predecessor

|Object|Bytes/lines|SHA256|
|---|---|---|
|formal/CoxeterProofs/B16IsolatedPairScalarKernel.lean|12112/282|aa1f20bed7ae1dffc50423c72d5828eaa96cd9742dbc20a132020edca6207138|
|docs/b16_isolated_pair_scalar_kernel_logic_20261023.md|10717/174|3dbcda28df4f5362a955b0e4e781bc83a554fc7438e5f474be4c65870a415983|
|results/lean_b16_isolated_pair_scalar_verification_20261023/failed_compile_source_preservation.json|2252|bab0e81fc275938ced4ea0be674c20af357baafa39cbae075562d48ba8253881|

I read the ENTIRE current282-line source and current map, the ENTIRE
original different-author precompile report and root clearance, and
the complete preserved compiler diagnostics/stdout. Read-only gzip
decoding verified all FOUR manifest records' compressed hashes/bytes
and decoded hashes/bytes, including source, map, diagnostics and stdout.

The decoded predecessor source is12102 bytes, SHA256
`b92f33dcbe030fb2489e8e162b8b8dfdeac568fbc3dd66339eeb651927efcdd6`;
its gzip is2989 bytes,
`c2c95a48d0e422ed592f874fe3446afa991cb911f4af2f1afec1b26a37e5cb56`.
The decoded predecessor map is9904 bytes,
`0fc165f58f097b26ac015a75f82d2f83e1a3350f0de7a583ab5bd5dd812fa733`;
its gzip is4542 bytes,
`dd07cb9710f30d0e16ca9dfdce4813b0e51240f86130515a19fef2b617641bc8`.

## Exact two-line source delta and unchanged proof boundary

The COMPLETE decoded-predecessor/current diff is precisely TWO local
identifier replacements, at the cubic and quadratic prefix proofs:

    have prefix := OrderedRing.mul_pos
    have middle_gain := OrderedRing.mul_pos

The current bytes equal the predecessor with exactly those two
replacements; the10-byte increase is exactly two five-byte changes.
No other character, statement, theorem signature, hypothesis, type,
coefficient, expression, tactic, import or footer request changed.
The renamed local propositions are not referenced by identifier in
later source; their positive-product facts remain in the `grind` context.
`middle_gain` is an ordinary identifier, with no conflict in either
local proof. This addresses the observed reserved-token parser errors
without changing the intended proof or adding assumptions.

The original FULL static report
`results/review_b16_isolated_pair_scalar_kernel_precompile_full_source_and_map_replay_20261023.md`
freshly hashes to
`98c956a5fa3c2bb6be7f95559ddd457d9e4bfad53da21c8d0eaf966c7bea6812`.
The original root clearance
`results/review_b16_isolated_pair_scalar_kernel_root_precompile_acceptance_20261023.md`
hashes to
`e54e699fe2dc87523696f163bca6c5ee6e43f404889fb803b223ddd85acb9fa3`.
Both bind the exact preserved b92/0fc predecessor; their execution
credit remains zero and they are not silently rebound to the successor.

The twelve selected declarations, all ordinary external scalar/data
premises, generic ordered-ring metric types and integer valuation
premises are unchanged. The source still imports ONLY Std. There is
no new axiom, sorry, instance, native_decide, mathematical fallback or
assumed final gain. I freshly matched ALL nine unchanged mathematical/
data map pins and the unchanged v4.32.2 toolchain pin
`2bdc48adfa58d0017e538a0ad117c5d73d35deec879978f909406a80c8037273`.

## Failed chronology and the complete map delta

The preserved full compiler diagnostics are4543 bytes,
`55e6e36ecf86e03bffd52375f7c52610a2ee85b91731825ccdfdc0b10917b210`.
They report unexpected reserved token `prefix` at lines221 and254,
with consequent unsolved goals. Two unused-section-variable warnings
are retained. The corresponding two failed prefix declarations print
`sorryAx`; those diagnostic prints are FAILED output, not a successful
object or an actual separate axiom audit. The816-byte runner stdout,
`9319ca8e586befb502265c93039bae13b10d019dc646423d9c29316bdebf65f5`,
records child exit1. No later compilation result is asserted here.

Relative to the preserved0fc map, the final map changes only the
source byte/hash binding and honest failed-attempt/repair chronology.
Its complete mathematical input, scalar algebra, external-assumption
and twelve-declaration tables are unchanged. The intermediate97baa
map retained a stale final sentence saying no Lean job had run. I
flagged that chronology defect; root corrected ONLY that sentence.

The intermediate map was byte-preserved in
`results/lean_b16_isolated_pair_scalar_verification_20261023/map_before_chronology_correction.gz`,
4825 bytes,SHA256
`aa75af0cc3a346acb6c3f4b4fddab902a5e08ecf76b89eca3f9932b8932511c0`.
Its independently decoded10670 bytes hash to
`97baa05b659dab48b94b95cca0497af9b79192bb86bbeb5855159bc7e77b2330`.
The sole final delta now states that no SUCCESSFUL compilation or
separate actual audit is credited. Thus source-only preparation,
observed failure, prior written/runtime mathematics and eventual Lean
success are distinguished accurately in the bound successor.

**Verdict:** source aa1f20/map3dbc are suitable for root's next named
bounded compilation if independently admitted, followed ONLY after
success by a SEPARATE actual twelve-declaration axiom audit. This
static report supplies zero compilation/audit PASS credit. Unchanged
resource gates, one shared-lock process,2048MiB/one thread and the
prohibition on rerunning successful older dependencies remain in force.
