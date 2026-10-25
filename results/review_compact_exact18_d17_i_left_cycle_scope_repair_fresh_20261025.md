# D17 actual-I: independent review of the one-hunk F1 repair

25 October 2026. Reviewer: `audit_finite57_algorithms`. Source only.

I compared the complete immutable 464-line predecessor
`results/source_archive/compact_exact18_d17_actual_i_transcription_20261025/compact_exact18_d17_actual_i_working.tex`
(21953 bytes, SHA
`b777eebcc8c912b9829d30eaf8f52057ad4b42c39c9b1ca5c2e4f474bd963172`)
with the complete current 468-line fragment
`manuscripts/revisions/parts/compact_exact18_d17_actual_i_working.tex`
(22177 bytes, SHA
`6a38015b2519e3a58b19262fb57dbb9124064e7fa67dffd4d2ca9e1e8587934d`).
The actual source difference is exactly the single hunk recorded in
`results/compact_exact18_d17_i_left_cycle_scope_repair_20261025.patch`,
SHA `fffc7aab41b4f43b63b70902d5330a6d288d6bb7254336b6336dd42bf4444075`.

The successor restricts the own-three all-three cycle to a simple-three
left contact. It keeps the far-end-three star/path/null cases separate
and sends the other positive near/next/end-four cases to their direct
ratio bounds. For end-four at either end, the smallest positive own
is g=1/2, Delta<=1/15 and the least numerator is
g+1/(15 sqrt(2)), so the printed lower ratio is

    sqrt(15)/2 + 1/sqrt(30) > 2.11 > 1.554.

A foreign cost only increases this ratio. Thus the repaired end-four
own-three branch is strictly strong without an affine-cycle claim.
This is precisely F1 of the complete outside report
`results/review_compact_exact18_d17_actual_i_transcription_fresh_20261025.md`,
SHA `40a9414fd244193dcaf39be9c62d56f08022658bff6a1b2bdf5ede22921e67e3`.
No other argument, numerical bound, domain or dependency changed.

The complete source-review assessment in that 321-line report therefore
transfers to this corrected successor with F1 closed. The shared-kernel
authorship disclosure, named common geometric inputs and limits of that
assessment remain in force. Root acceptance and same-family integration
are separate actions. This is not a new global dimension bound, build,
runtime result or Lean proof. No worker/import, test, enumeration,
scientific process or compilation was run; only source/diff/hash I/O
was used to inspect the actual repair.
