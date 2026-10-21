# Bounded cavity scouts and the selected replacement

These numerical scouts use midpoint coordinates from the frozen exact
release34 residual4 model. They are discovery diagnostics, not exact
nonexistence certificates. Every candidate is tested against all native
old contacts, including projected walls redundant for points of S6.

The [right-S6 scout](../scripts/fresh_nonarithmetic_release34_cavity_20261017.py)
tests the two-endpoint frames `(7,12,47)`, `(8,12,52)`, `(30,47,52)` with
labels2..32,48,64,128. Its [result](fresh_nonarithmetic_release34_cavity_20261017.json)
has2,148 candidates and149,472 compatible dotted pairs between groups,
but no finite new-new pair. Every first/second-group pair is obtuse.
The [continuous gate](fresh_nonarithmetic_release34_cavity_gate_20261017.md)
certifies the large-label portion without a cutoff.

The [seven-white-profile scout](../scripts/fresh_nonarithmetic_release34_cavity_a2_20261017.py)
adds the six possible single A2 attachments of labels3,4,5, right to H4.
It tests framing labels2..16,24,32,64; its
[result](fresh_nonarithmetic_release34_cavity_a2_20261017.json) has588
candidates and again no finite compatible pair. Nonzero A2 profiles survive
only in the third group. The exact
[white-contact input](fresh_nonarithmetic_release34_cavity_20261017_white_input.json)
is rebuilt by the source's `--rebuild-input` option. Native white offsets
are included in old-contact and new-new-contact tests.

The [replacement-event scout](../scripts/fresh_nonarithmetic_release34_replace12_47_20261017.py)
preserves12's contacts with7 and8 and imposes a finite event with47 or52,
using labels2..128,192,256,512. Retaining all65 walls gives
[no candidate](fresh_nonarithmetic_release34_replace12_47_20261017.json).
After releasing54, the
[separate result](fresh_nonarithmetic_release34_replace12_47_20261017_drop54.json)
has126 replacements, starting with label6 to52; each has16 intrinsic open
ends after deletion costs. The subsequent bounded replacement47 events
produce no accepted pair. No all-label exclusion is inferred.

The simplest label6 replacement has been promoted to the
[exact64-wall construction and complete boundary](fresh_nonarithmetic_release54_replace12_20261017.md).
It closes a compact heptagon and improves17 open ends to16. That exact
certificate, rather than a numerical tolerance test, supports the geometric
claim.

Generation and saved numerical replays passed. These bounded negative
searches leave other labels, contact patterns, replacements and larger
clusters open.
