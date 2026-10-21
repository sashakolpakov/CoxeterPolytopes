# Root source review: exact first-edge worker

20 October 2026. Full source and mathematical-interface review, before
execution. No certificate counts or successful execution are asserted here.

* Worker `scripts/check_compact_h9_sqrt3_diagonal_first_edge_20261020.py`,
  SHA256 `56149b65859e61b980eabc60248d7fb4541e40b7a488a949d423bc53995abcd8`.
* Current route `results/compact_h9_sqrt3_diagonal_construction_route_20261020.md`,
  SHA256 `16732ee8454f20051c524cbaeb8a56ec8b1b9dbb916f7644d6b08adfd1f60297`.

PASS for a bounded, complete first-edge census within the explicit written
geometric interface. The route successor corrects only the unit interval's
closed lower endpoint and clarifies unit-square normalization. The earlier
root mathematical review of route1ac6604b remains applicable to these
unchanged derivations.

The K-arithmetic is exact. The sign function handles opposite-sign
coefficients by the nonzero integer a^2-3b^2; field quotient integrality
uses conjugation and an integer denominator. No floating point decision
occurs. Reflection checks inspect all entries of the actual ten-by-ten
matrix, so allowing extra nonprimitive candidates does not introduce
nonintegral formal mirrors.

The coefficient boxes cover the whole reduced strip. From
0<z<48 and |z'|<5, its rational coefficient is an integer in[-2,26].
Using sqrt3>3/2 bounds the other coefficient within[-1,17]. Likewise
0<=X<35 and |X'|<4 give the displayed[-1,19] and[-1,12] boxes. The
ten R values satisfy0<R<160 and0<R'<=18. The square dictionary is an
exhaustive two-dimensional enumeration, not a sampled optimization.
At X>=0,z>0, the strict reduced-strip condition is equivalent to
z^2<epsilon^2*R, or X^2<(1+sqrt3)^2*R; imposing both loses no point.

All ten equations, both placements of the nonisolated profiles, every
power0 through4, and the exact certified upper-mirror cutoff are tested.
Coordinate congruences reconstruct a and b integrally; b>=0 enforces
the released zero wall as well. The other eight zero-wall products are
checked against the exact intended gaps. The chosen norm and complete
reflection matrix are checked independently of the profile constructor.

For simultaneous candidates the prefix dominance sign is correct.
The vector a_left*r_right-a_right*r_left has zero time coordinate and
zero spatial sum. Writing it as sum d_i*(-e_i+e_(i+1)), its coefficients
are d_i=a_right*prefix(r_left)-a_left*prefix(r_right). Nonnegative d_i
make the left inequality imply the right inequality inside the A8
sector. The worker requires a common dominant simple inequality and
checks projective equivalence of all such choices.

The first-collision conclusion uses standard discrete reflection-chamber
theory: the actual first wall has the listed finite rank-nine corner,
and is included by the complete primitive norm inventory. Every extra
accepted candidate is itself an integral reflection mirror oriented
away from the initial point, so none can cross this ray before its
actual first mirror. This is an explicitly written input, not a global
chamber-completion theorem supplied by the program.

Generation refuses to overwrite evidence. Replay rebuilds the full
enumeration and compares the entire JSON object, including the pinned
sources and census digests. This is a separate execution of the same
implementation, not an independent algorithm. The external wrapper must
retain the common lock and all resource caps. No compactness or finite
volume result follows from a single completed edge.
