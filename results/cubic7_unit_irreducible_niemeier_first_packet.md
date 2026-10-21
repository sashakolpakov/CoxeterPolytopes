# Cubic-seven unit coefficient: irreducible Niemeier first packets

The two irreducible rank-24 Niemeier root systems were screened for

\[
\langle-(\theta-1)\rangle\oplus N,
\qquad N=D_{24}^{+}\text{ or }N(A_{24}).
\]

Both have complete canonical primitive mirror norms `1,2`.  Their first
positive wall is the norm-one affine highest-root wall at control
`2+5theta+2theta^2` and height `5+14theta+6theta^2`, with scale
`2theta+theta^2`.

The actual glue bases and root systems are used, not coordinate-axis proxy
walls.  For `D24+`, the spatial vector is the D24 root
`(2theta+theta^2)(e1+e2)` and the affine marks
`(1,2^21,1,1)` reconstruct the highest root.  The 25 exact alcove rays are
2 finite, 0 ideal, and 23 hyperideal.

For the A24 Niemeier lattice, the spatial vector is
`(2theta+theta^2)(e1-e25)`.  The 25 exact alcove rays are 1 finite, 0 ideal,
and 24 hyperideal.

Thus neither irreducible Niemeier unit seed closes at its first packet.  The
statement is deliberately scoped to this exact packet; later walls can cut
hyperideal rays.

- Backends: `src/coxeter_polytopes/search/cubic7_d24_plus_unit.py`,
  `src/coxeter_polytopes/search/cubic7_a24_plus_unit.py`
- Audit: `scripts/audit_cubic7_unit_irreducible_niemeier.py`
- Tests: `tests/test_cubic7_d24_plus_unit.py`,
  `tests/test_cubic7_a24_plus_unit.py`
- Archive:
  `data/enumerations/cubic7_unit_irreducible_niemeier_first_packet.json`
