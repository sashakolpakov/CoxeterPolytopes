# Independent source review of the arithmetic15 layout proposal

25 October 2026. Reviewer `/root/ra_release_completion_audit`, different
from the proposal author. **Source-semantic clearance: PASS.** Root
acceptance and any bounded compilation/rendered inspection are separate.

Reviewed exact files:

- Authoritative preimage: `manuscripts/compact_arithmetic_dimension_bound_consolidated.tex`,
  SHA256 `821fc73f18021dee86d411148c7278c2f672346d539ce25952db3d108ec9df88`.
- Proposal patch: `/private/tmp/arithmetic15_remaining_layout_frac_proposal_20261025.patch`,
  70,505 bytes, SHA256 `1b9b6bfbae984e6d57d76484ef8a1fff29860b193e671bc0933497d4ee84f00e`.
- Candidate: `/private/tmp/arithmetic15_remaining_layout_frac_proposal_20261025.tex`,
  695,927 bytes, SHA256 `f3484a664a9e6f93ad1a8ac8d02f217c781c6b5aa3a9578cdf83dbf17683f05c`.
- Metadata: `/private/tmp/arithmetic15_remaining_layout_frac_proposal_20261025.json`,
  43,504 bytes, SHA256 `7ef883ed7712978543cd105eda28028b128af4729a3d5b45a7f8ea5b3b94a2e6`.

The patch is exactly the unified difference between the authoritative
preimage and candidate. The authoritative source remained at its preimage
hash throughout this review.

## All 361 grouped infix fractions

I independently scanned the original source, respecting TeX control
sequences, escaped braces and comments. It contains exactly 361 relevant
`\over` tokens. Each belongs to one explicit balanced brace group, with
exactly one top-level infix fraction in that group. Every ordered metadata
entry agrees byte for byte with that original group and its numerator and
denominator. None contains a cross-fraction `\left`/`\right` pair,
explicit group command or environment, and none changes a top-level font
or math-style switch.

For every entry, the actual replacement is exactly

    {\frac{numerator.strip()}{denominator.strip()}}

where `strip()` removes only the outer whitespace of the two recorded
payload strings. The original outer group remains. The numerator and
denominator tokens, their order and all internal grouping are preserved.
Thus this is the ordinary LaTeX spelling of the same grouped fraction,
with no scope or operator-precedence change. No infix `\over` remains
in the candidate. Reconstructing all 361 replacements independently
leaves precisely the eight local edits below. The bound patch and the
metadata's original/component strings together record every exact pair.
The independent control-word scan includes 237 fractions whose denominator
starts immediately with a digit; it does not use a regex word boundary
after `\over`. The candidate has zero actual `\over` control words.

## Eight local edits

1. The single end-hole energy `12569249/15000000` moves from inline math
   to a displayed fraction, with the same numerator and denominator.
2. The three squared radical margins become a two-line gathered display;
   all three values remain in their original order.
3. The three adjacent intervals receive separate math delimiters and an
   ordinary conjunction. Their endpoints and order are unchanged.
4. The four cubic-polynomial signed fractions become a table of the four
   previously stated arguments and their numerators/denominators. I read
   the full original polynomial/sign paragraph: the table preserves the
   argument order and each exact signed fraction. No value was recalculated.
5. The two separately paid actual-left-H4 totals become two aligned rows;
   all summands, signs and strict inequalities are unchanged.
6. The three continuous lower expressions become separate gathered rows.
   The following table retains the same three square-root bounds and their
   respective gain lower bounds. Calling the radical bounds upper bounds
   exposes the direction used by the existing circle-inequality lower seal;
   it supplies no new estimate or different pairing. The full-interval,
   no-sampling and no-I-spending statements remain unchanged.
7. The printed-proof-home table receives wrapping, ragged-right columns.
   Every row, label, mathematical statement and proof reference is retained.
8. The source-interface longtable receives a width tied to the text block
   and ragged-right wrapping. Its contents and row order are unchanged.

The required `array` package was already loaded. The complete ordered lists
of labels, references and theorem/proposition/lemma/corollary/proof
environment boundaries agree between preimage and candidate. No hypothesis,
case split, mathematical token in the fraction payloads, theorem conclusion,
computation premise, verification scope or source dependency was changed.

## Conclusion and limits

The exact patch and candidate above are cleared for root review and a
necessary bounded layout build. No source-semantic repair was identified.
This clearance does not establish the visual result or a new mathematical
claim. I made no authoritative edit, imported no scientific worker, and ran
no tests, mathematical calculations, enumeration, Lean/TeX compilation or
PDF process. Only source/byte comparisons and this private review note
were produced.
