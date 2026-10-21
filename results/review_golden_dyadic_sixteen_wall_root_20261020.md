# Root review of the fixed sixteen-wall escape

20 October 2026. Complete source and mathematical review before execution.
No new compact example or full-lattice nonreflectivity follows.

Reviewed scientific note `golden_dyadic_sixteen_wall_caps_20261020.md`
(SHA256 `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`),
coordinate input (SHA256
`80f18b60003375876f5bae9c122c9466ff1eb1c69fed582d87c7facf50fd5745`),
preparation (SHA256
`6cdeff97990a342b3f5e7afc5026c89cf17bc0f747b853a563fd80e618d9233d`),
and complete fixed checker `check_golden_dyadic_sixteen_wall_cone_20261020.py`
(SHA256 `5539759638a792ac56bd4f762fec2f3ae2c1bd4c1c0751cb297fd272208f28c9`).

The ten embedding columns recover the literal original lattice and
the diagonal M6 plus positive H3 metric. The checker uses a complete
original integral basis for its 160 reflection tests. Its 120 pair
tests recompute products in both systems, rather than merely comparing
the supplied ledger. The old nine-wall Gram has positive elimination
pivots, and its inverse and duals are checked independently. All sixteen
strict-point products are negative, and its square is negative.

For the fixed escape, its M6 square is `21-13phi=-phi^-7`, and its
positive H3 square is `26-16phi=2phi^-6`. Their sum
`47-29phi` is strictly positive. All sixteen wall products are
nonpositive, and the nine active rows have rank nine. The checker
does not enumerate extreme rays or make a completeness assumption.
The spacelike witness alone suffices for noncompactness of the
finite prefix. Perturbing it towards the strict future timelike point
preserves positive square and makes every wall product strictly
negative. The segment between these points crosses the future light
cone in the strict interior of every halfspace. This supplies an open
ideal-boundary patch and hence infinite hyperbolic volume.

The source's bounded affine section and fixed midpoint section are
additional hand calculations, not prerequisites for this implication.
The checker does not claim to verify those unused section enumerations.
In particular the full integral chamber can have additional walls
cutting the escape. Its compactness remains open.

Execution contract: exclusive certificate creation; fresh separate
recomputation for replay; 128 KiB output cap; exact source hashes
before and after; required parent watchdog and shared lock; no search,
cone enumeration, automatic retry, or heavy dependency. Independent
pre-execution review `review_golden_dyadic_sixteen_wall_fixed_escape_batch_guard_20261020.md`
(SHA256 `771e08e469c6483c03a701f5c6619acad6ba54aef572e8ab997ad0a55bdb84f1`)
also passes. Root authorizes one bounded generation and one separate
replay, subject to the existing resource gate. Execution is recorded
separately; it is not asserted by this pre-execution review.
