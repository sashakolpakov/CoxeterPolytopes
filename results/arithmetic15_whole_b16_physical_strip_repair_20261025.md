# Whole B16: exact physical-strip wording repair

25 October2026. SOURCE ONLY. No worker/import, enumeration, arithmetic
checker, test, Lean, TeX or build ran. The six canonical pairs and
root-owned mains remain untouched. This is a proposed precision repair
awaiting independent delta review and root acceptance.

The original whole-B16 freeze72979324 remains byte-exact. Its only new-body
issue F1, identified in the complete outside report
review_arithmetic15_whole_b16_fresh_20261025.md
(SHA cd4cb9a3df7da08d01ba7d9b2e50d66a612239e2ced3c4fa984ba6a8abab3b16),
is that the eleven-row table's actual physical cutoff needs to be printed
immediately before the table, alongside the coarse bound used for degree.

The one inserted paragraph now states both exact native/edge strips:
far sqrt2<=S<(16sqrt17-1)/8<65/8 and
near2+sqrt2<=S<(16sqrt17+1)/8<67/8.
It derives them from the two already established C_S formulas and C_S<16,
and explains explicitly why far4+3sqrt2 is outside the actual domain.
No row, numerical value, elimination or algorithm predicate changed.
Primary proof home remains4fb20275 Section5; historical sources unchanged.

The growing A12-universal preimage4ce85c46 was preserved at
results/source_archive/arithmetic15_before_whole_b16_physical_strip_repair_20261025/arithmetic15_missing_foundations_working.tex.
The complete live-preimage diff consists of precisely this one insertion.
In particular the earlier separately corrected actual-I far-next3
wording remains unchanged, and all later A-chain bodies remain unchanged.

Fixed successor archive:
results/source_archive/arithmetic15_whole_b16_physical_strip_repair_20261025/.

* whole_b16_section.tex:635 lines/20583 bytes,
  SHA2562e26857715fc9293c5009244ccffb07f829cc6e73c8fc9963e72dedfa50071c7.
  This is the reviewed original whole-B16 section plus ONLY F1.
* whole_b16_section_predecessor.tex: exact original7297 section6000–6617.
* arithmetic15_missing_foundations_working.tex:9726 lines/324843 bytes,
  SHA256913fe38c7c06a5d368ff253c7efaeae81598bde1219b0a8fa1ee0b835527a00a.
  This whole-part archive includes later bodies outside the whole-B16
  review; no whole-file successor acceptance is requested.
* results/arithmetic15_whole_b16_physical_strip_repair_20261025.patch:
  complete section-only minimal diff,
  SHA256db2efc19b4ebf13737e1d0443929206008027bc90e811bb650bfaa4535223013.

Integration must take ONLY the repaired whole-B16 section onto the
already corrected accepted main. It must not replace a main with an old
whole-part prefix. No independent review, root acceptance, execution or
publication-ready status is claimed merely by preparing this successor.
