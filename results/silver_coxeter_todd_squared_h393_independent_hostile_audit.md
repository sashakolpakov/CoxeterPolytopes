# Independent hostile audit: silver K12^2 through H393

## Verdict

**PASS.** A clean-room implementation, importing no project arithmetic or
search backend, reproduces the complete 16-class chronology through
`H=393+278sqrt(2)`: 128 controls at
87 heights and exactly six nonempty
necessary rows. The explicit norm-two cutter is primitive, reflective, and
acute to the complete earlier prefix. The `7/4` dual-direction replacement
is strictly spacelike and strict against every necessary row through the
bound. This remains a bounded infinite-volume certificate, not Vinberg
termination or a global nonreflectivity result.

## Independent checks

- The literal K12 Gram has determinant 729; `G h=e1`, `h^2=4/3`, and `3h`
  is integral of norm 12. Independent PARI enumeration recovers coordinate
  minima `0,4,4,10` and cutter shell maxima `36+26sqrt(2)` and
  `40+30sqrt(2)`.
- The exact control census has 50 controls
  through the first packet, 76 strictly between
  it and the cutter, and 2 at the cutter.
- Direct coordinate multiplication gives `Gram(D,E)=(28,84,46)`,
  `Gram(U,V)=(380,204,278)`, `r^2=2`, and `(h.U,h.V)=(22,16)`.
  All reflection pairings are integral. An odd rational coefficient of `U`
  excludes divisibility by the unique prime over two, proving primitivity.
- Complete shell maximisation makes the cutter strictly acute to both first
  packet components. Necessary four-vector Schur censuses give intervening
  maxima `0,0,-7-4sqrt(2)`; all four simultaneous norm-one maxima are strict.
- For a distinct simultaneous norm-two root, the two-embedding Cauchy box
  leaves only positive pairings `1,sqrt(2),2`. The first two force forbidden
  coefficient-lattice squares below the K12^2 minimum four; `2` forces equal
  spatial vectors. Thus no distinct positive same-height pairing exists.
- The replacement norm is `205/12 - 12*sqrt(2)` and its tight
  wall margin is `1/2 - 1*sqrt(2)`.

Owner claims and the independently reconstructed summary agree exactly. The
control/Gram stream hash is
`d2638caf5af400fa4e586bd2b27478f38f6a0498dccec0b760e5cf0a2c9a85c4`.

## Reproduction

```bash
python3 scripts/independent_hostile_audit_silver_coxeter_todd_squared_h393_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_h393_independent_hostile_audit.py
```
