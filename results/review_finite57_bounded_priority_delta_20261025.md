# Finite57: bounded-tail priority wording delta

25 October 2026. Source-only review by the bounded-tail audit author.
No worker, import, enumeration, test, compilation or certificate arithmetic
was executed for this delta.

The live same-family revision has SHA256
8069e413481fa03333bd81eb20c65de9ed42a7ba208427179cbe711bdde312a0.
Its preserved predecessor is
results/source_archive/finite57_before_bounded_priority_wording_20261025/finite_volume57_consolidated.tex.gz.
A literal source comparison has exactly one hunk, at the concluding
paragraph of the dimension59--65 section.

The successor says that **every newly reduced affine corner**, including
a nonzero fractional reduction, receives its whole-face theorem first.
Only when the actual face has no such corner are its ideal weights said
to agree pointwise with the preceding layer. Finite sparse changes are
handled afterwards. This is precisely the priority distinction required
in section6 of
results/audit_finite57_bounded_tail_layers_fresh_20261025.md.
It avoids the predecessor's overly narrow wording “a new affine zero”.

I accept this source delta at the scope of that recommendation. The delta
does not change the marked arrays, their selectors, the finite outputs,
or any computational or geometric premise. In particular, all five
binding gates in section9 of the bounded-tail report remain open exactly
as recorded there. This is a wording correction and source review, not a
new execution receipt or a complete runtime verification of the tail.
