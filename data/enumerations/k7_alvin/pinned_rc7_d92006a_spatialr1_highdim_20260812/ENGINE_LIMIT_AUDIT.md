# Pinned RC7 high-dimensional engine-limit audit

Image: `coxeter-alvin:rc7-spatial-r1-d92006a`  
Immutable ID: `sha256:b995880fcd16c6940af225e0c6e5d8734be495f434d8e66640774a4714252896`  
AlVin commit: `d92006a842cb4668e7c89e8c8e2d30337e4b0f26`

This is a source-level launch audit, not a termination or reflectivity result.

## Dimension and index storage

- `/opt/AlVin/alvin.cpp:39` sets the dimension to the runtime diagonal length
  minus one.  Root vectors and Coxeter matrices are dynamically sized from
  that value (`alvin.cpp:68-117`, `rcyclotomic7integer_alvin.cpp:68-73`).
- The dimension and loop indices are `unsigned int`; no fixed H22/H29 array or
  bit-mask cap was found in AlVin or its bundled CoxIter.  Ambient ranks 23 and
  30 are therefore far below the representational limit.

## RC7 norm-factor gate

- `/opt/AlVin/rcyclotomic7integer_alvin.cpp:454-517` collects distinct prime
  factors of the diagonal, adds the rational prime 2 and two fundamental
  units, and throws if the resulting distinct-factor count exceeds 8.
- It then scans at most `2^f - 1` raw factor subsets.  Repeating the coefficient
  1 twenty-one or twenty-eight times introduces no new factor, so dimension
  does not enlarge this subset scan.  The p7 and inert2 forms remain safely
  below the eight-factor gate (at most five collected entries before
  deduplication effects), hence at most 31 raw subsets.

## Practical limits

- The recursive coordinate enumeration and CoxIter finite-covolume checks are
  dimension-dependent and combinatorial.  They are the practical runtime and
  memory risks; the source exposes no reliable a priori completion-time bound.
- A complete AlVin termination still requires independent exact replay.  No
  finite-volume, compactness, candidate, or non-reflectivity claim is made for
  a live or externally interrupted process.

## Form arity gate

`tests/test_run_cubic7_alvin_live.py` proves that each H22 form has exactly 23
diagonal entries with 21 unit slots, and each H29 form has exactly 30 entries
with 28 unit slots.  The exact RC7 semantic conversion is checked as well.
