# Independent paper review: repeated B2 corners

20 October 2026. Reviewer: batch_guard. Source/math reading only; no
import, computation, proof job or build was performed.

Reviewed in full:
`results/compact_b2_matching_corner_drift_exclusion_20261020.md`,
SHA256 `4170c8c1e4e1641148b2287d8b177d1e09fd9befeacd7cde9ae629b25183fdd2`.

**PASS within the stated scope.** No mathematical flaw was found.

I independently checked the B2 inverse Gram and the actual finite-edge
Schur formula. Positive definiteness permits a contact to another whole
B2 only at one node, with label3 and cost1/2, and at most one such other
block. The zero-contact, two-contact and off-block branches in Section1
are complete. Native contact1 is correctly excluded using compactness.

The positive-definite projection argument is valid: its coefficients
are positive, all omitted wall inequalities remain strict by acuteness,
and the resulting point is future timelike. For a finite one-contact
wall it proves both proposed endpoints actual. It therefore justifies
using one wall for both edge directions, rather than independently
assigning hypothetical endpoints.

I rederived both lengths in (3), the k>=6 half-tanh bound, both strict
k5 rational estimates, and the positive-height obstruction to two
k5 exception walls in different blocks. The complete block sums and
the final arithmetic are correct: at r6 the squared lower bound is
1369/32, and its ratio to r increases thereafter.

Every released white diagram is the same A1+B2^(r-1), so the argument
applies to ANY one globally shared positive white weight. It excludes
B2^r vertices from the undirected drift witness for r>=6. It neither
excludes all finite vertex types nor proves a new uniform dimension
ceiling. No arithmeticity premise is used.
