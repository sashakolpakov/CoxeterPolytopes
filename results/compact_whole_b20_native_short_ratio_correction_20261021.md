# Whole B20: correction of the zero-own short-contact ratio

21 October 2026. Author correction following root's full mathematical
review. The predecessor's proposed THREE-wall reduction was unsupported:
it dropped a factor1/sqrt2 in one geometric ratio. This is a substantive
coverage correction, not a wording-only clarification.

The exact last-long ratio, already correctly displayed in that source,
is

    R=(gamma+delta/sqrt2)/sqrt(1-delta^2).

At zero own contact gamma0 and short label m, delta=cos(pi/m), hence

    R=cot(pi/m)/sqrt2.

In particular m6 gives sqrt(3/2), not sqrt3; m7 is also below the
claimed uniform strong threshold. Neither is excluded by the supplied
proof. At m8 the ratio is exactly1+1/sqrt2>33/20, and it increases
with m. Thus the proved direct-strong tail is m>=8. Zero-own short
labels5,6,7 must all remain, together with own3/short3 and own3/short4.
The corrected source claims only this FIVE-type necessary reduction,
still pending complete review; no whole-B20 or global exclusion follows.

Before editing, I preserved the complete predecessor at
`results/source_archive/compact_whole_b20_before_zero_own_short_ratio_correction_20261021.md.gz`:

* Original16308 bytes, SHA256
  `acd331054227f6b7ae6f0c860eb4fae705df333b696428fe8b0d81b9f70707cb`.
* Gzip7490 bytes, SHA256
  `d33ffdf275d914ac062b4f5556c946c697d000d874cc0b51f0dcfb65eb4cf1eb`.
* Decompression was checked byte-for-byte against that original.

Current successor:
`results/compact_whole_b20_continuation_20261021.md`,16495 bytes/337 lines,
SHA256 `ea689009a12b4c2d8d107ad397e5a6718f27450e194157e09695eabdd9a1628a`.
The inspected exact delta changes the header/theorem/outro to five types,
corrects the zero-own threshold paragraph, and lists the additional actual
short-edge cosines C20=cos(pi/m)*sqrt20 for m5,6,7. The other actual
point, coefficient-payment and reference arguments are unchanged. No job,
worker import, Lean/TeX run, frozen input or guide edit occurred.
