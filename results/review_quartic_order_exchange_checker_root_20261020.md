# Root source audit: fixed quartic order-exchange checker

20 October2026. Source review and syntax check only; execution pending.

Worker: `scripts/check_quartic_order_exchange_period_20261020.py`,
SHA256 `a30ea0e6b2b8e5e71f94ab3fd9dc5acab58a320a0a624c42929cceb9b08a9340`.
Arithmetic helper: `scripts/fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.py`,
SHA256 `122cfeea5190517d9203a61d933fbe964485652be04c5fdca157e492eca2c811`.

**PASS for the fixed arithmetic scope.** The worker uses exact rational
coefficients in the biquadratic field. The helper's catalogue does not
execute on import. Integral-basis conversion and its integrality test
match the full ring from the written preflight.

All loops are fixed five-dimensional matrix operations. The explicit
block Gram inverse avoids indefinite Gaussian pivots. Both frame
inverses, both period inverses, the metric, every oriented simple-root
image, future points and the exact trace are checked. Initial and final
source digests bind the same inputs; certificate creation refuses
replacement and fresh replay recomputes all identities.

Full arrangement, actual corners, finite classification, free chamber
action and the arithmetic lattice theorem remain explicit written
inputs. This worker is not a geometric completeness proof or a search.
Run generation and replay separately through the parent-lock bounded
wrapper, with its59s/1GiB/one-thread caps and40GiB/6GiB start gates.
No new Lean module or TeX build is needed for this research-family
obstruction; the previously checked conditional chamber-period kernel
does not itself formalize this geometric instance.
