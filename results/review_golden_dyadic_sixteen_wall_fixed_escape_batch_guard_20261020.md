# Independent pre-execution review: fixed sixteen-wall escape

20 October 2026. Source and paper review only. I did not import or execute
the worker, enumerate rays, or run a proof process.

**PASS within the stated finite-prefix scope.** This review binds:

- `scripts/check_golden_dyadic_sixteen_wall_cone_20261020.py`, SHA256
  `5539759638a792ac56bd4f762fec2f3ae2c1bd4c1c0751cb297fd272208f28c9`;
- preparation `golden_dyadic_sixteen_wall_cone_preparation_20261020.md`,
  SHA256 `6cdeff97990a342b3f5e7afc5026c89cf17bc0f747b853a563fd80e618d9233d`;
- paper `golden_dyadic_sixteen_wall_caps_20261020.md`, SHA256
  `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`;
- its coordinate JSON, SHA256
  `80f18b60003375876f5bae9c122c9466ff1eb1c69fed582d87c7facf50fd5745`.

I read the complete worker and its definitions-only quadratic-field and
matrix-inverse helper. Input pairs in `1,phi` are independently converted
to the helper's `1,sqrt(5)` representation. Exact signs and the ring test
use the correct ring of integers. The ten displayed embedding columns
are compared entrywise with the input and their full Gram is checked.
The time/five-coordinate/four-D4-root basis is a complete original
integral basis. Checking each reflection on these ten vectors proves
lattice preservation; no field-similarity or primitivity shortcut is
used. All sixteen roots, 120 products, native angle alternatives and
strict center are reconstructed independently of the supplied ledger.

The old Gram positive-pivot calculation and inverse justify the corner
coordinate change. Both the strict point and the fixed escape round-trip
through it. The escape has the asserted square `47-29phi>0`, positive
original time, and all sixteen nonpositive products. The nine active
full rows are actually reduced to rank nine; the future row is strict.
This establishes the stated extreme ray without an omitted DD stage.

The written infinite-volume implication is sound: perturbing the escape
toward the strict timelike point gives a strict spacelike feasible point
on the same time side. Its segment to the strict timelike point crosses
the future light cone strictly inside every halfspace, yielding an open
ideal-boundary patch. Rank nine is useful extra evidence, not a premise
needed for that implication.

The worker checks source pins before and after, enforces the bounded
parent-owner interface, refuses replacement when writing, and recomputes
the entire certificate when verifying. It has no cone/height search or
fallback expansion. This is a pre-execution gate; it asserts no actual
generation or replay status. Its conclusion concerns the sixteen chosen
halfspaces only. Later mirrors and the full integral R5 chamber remain
unresolved; neither a compact example nor full-host nonreflectivity is
claimed.
