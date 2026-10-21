# Independent golden `D24` / `D24+` prefix audit

An implementation-independent exact replay confirms the norm-two Vinberg
prefix of

\[
\langle-2\varphi\rangle\oplus D_{24}
\quad\hbox{and}\quad
\langle-2\varphi\rangle\oplus D_{24}^{+}
\]

through controls `x=1, phi, phi^2`.  The audit uses standalone rational-pair
arithmetic for `Q(phi)`, enumerates all four residues of `O/2O`, imposes the
two real-place bounds exactly, checks lattice membership, and independently
replays every pairing with the `D24` simple walls and earlier accepted walls.
It imports neither the production search nor its exact-arithmetic class.

The complete dominant candidate counts at the three controls are `1,3,4`;
the accepted counts are `1,2,1`, for both `D24` and `D24+`.  All three
nonzero half-spinor residues have zero candidates through this prefix, so the
two lattice results genuinely coincide here.  The four accepted spatial
supports are

- `x=1`: `(phi,phi)`;
- `x=phi`: `(phi,phi,phi,phi)` and `(2phi)`;
- `x=phi^2`: `(1+phi,1+phi,phi,phi,phi,phi)`.

Each accepted vector belongs to the `D24` root coset and has Lorentzian norm
two.  In the displayed order their exact pairing matrix is

```text
 2          0          0          0
 0          2       -2phi         0
 0       -2phi         2       -2-2phi
 0          0       -2-2phi       2
```

The other four candidates have archived positive-pairing witnesses: the
mixed candidate at `x=phi` pairs by `phi` with the `x=1` wall; at `x=phi^2`,
two candidates pair by `1` with that wall and the remaining candidate pairs
by `phi` with the accepted `(2phi)` wall.  Thus the accepted count is not an
artifact of ordering inside an equal-height packet.

The coefficient box is exhaustive.  For doubled spatial coordinates
`z=A+B phi`, the largest target gives `|z|<12` and `|z'|<4`; hence
`|B| sqrt(5)<16`, so `|B|<=7`, and `|2A+B|<16`, so `|A|<=11`.  The audit
enumerates the slightly larger box `-12<=A<=12`, `-7<=B<=7`.

This certificate is deliberately limited to norm two and controls through
`phi^2`; it makes no claim about later shells or other mirror norms.

- Audit: `scripts/audit_golden_d24_prefix_independent.py`
- Tests: `tests/test_golden_d24_prefix_independent_audit.py`
- Archive: `data/enumerations/golden_d24_prefix_independent_audit.json`
