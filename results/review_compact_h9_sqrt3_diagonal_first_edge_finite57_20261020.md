# Independent source review: the square-root-of-three diagonal edge

20 October 2026. **PASS for source and mathematical coverage, before
execution.** No import, enumeration, checker run, compilation or typesetting
job was performed in this review. No outcome or candidate count is asserted.

Sources read in full:

* `results/compact_h9_sqrt3_diagonal_construction_route_20261020.md`, SHA256
  `16732ee8454f20051c524cbaeb8a56ec8b1b9dbb916f7644d6b08adfd1f60297`.
* `scripts/check_compact_h9_sqrt3_diagonal_first_edge_20261020.py`, SHA256
  `56149b65859e61b980eabc60248d7fb4541e40b7a488a949d423bc53995abcd8`.

## Inventory and finite reduction

For O=Z[sqrt(3)], the primitive-reflector norm divisibility and the unit
argument give exactly the six necessary representatives used by the worker,
modulo squares of O-units. These are not six classes modulo arbitrary
squares in the field. The finite rank-nine extensions of the retained A8
give the six isolated rows and four contact rows, with the two placements
specified in each contact row. Recomputing
`R=9q-k(9-k)g^2` gives every stored R and no missing placement.

The upper mirror construction uses the actual ninth power of
`U=(2+s)+(1+s)*sqrt(s)`. Its integral coefficients have relative norm one;
the divisibilities for a=X9/3 and b=(Z9+1)/9 give an integral reflector of
norm 2b. The worker verifies its full reflection matrix, rather than
assuming primitivity or that its square belongs to the primitive pool.

The strip reduction needs no fundamental-unit assertion: U and its inverse
preserve O-pairs. For z,X positive, reducing w=z+sqrt(s)X into
`sqrt(R)<=w<U*sqrt(R)` uses a nonnegative power to return to the original
solution. The upper mirror has ratio U^(9/2), so powers 0 through4 cover
every solution at or before its collision. The exact cross multiplication
in `census()` has the correct orientation and includes equality.

The enlarged boxes are complete. From `0<z<48, |z'|<5`, integer coordinates
of z satisfy `-2<=A<=26, -1<=B<=17`; from
`0<=X<35, |X'|<4`, those of X satisfy
`-1<=C<=19, -1<=D<=12`. These follow from 3/2<s<7/4 and the ten R bounds.
The square index is an exact two-coordinate lookup. The further strip and
conjugate tests remove no reduced solution; the strict upper strip boundary
and the possible X=0 lower boundary are treated correctly.

## Arithmetic, roots and minimum

The pair arithmetic, sign test for opposite-sign coefficients, conjugation
and integral quotient are exact. In `representatives()`, division by s
correctly maps A+B*s to B+(A/3)*s. `build_root()` enforces X in3O and
z-kg in9O, a>0, b>=0, the exact norm, and nonpositive products with all
nine original walls, including the released one. Its reflection test checks
every entry `2*r_i*(G*r)_j/q` for O-integrality. Thus a retained vector is
an actual integral reflection, even when it is not primitive.

Keeping such nonprimitive extras is safe: they cannot create a spurious
earlier mirror. Every primitive norm representative remains covered, and
every extra accepted vector still supplies a genuine mirror of the same
integral reflection arrangement. The certificate explicitly records this
handling; it does not certify that all retained vectors are primitive.

The collision comparison uses a/z with positive denominators. For tied
candidates L and R, the prefix tested in `simple_dominates(L,R)` is

    d_i = a_R sum_(j<=i)y_(L,j) - a_L sum_(j<=i)y_(R,j).

These are exactly the coefficients of the old simple roots in
`a_L*r_R-a_R*r_L`. Nonnegative prefixes therefore mean L's inequality
implies R's throughout the old A8 sector. The final total coordinate is zero
because collisions tie. Requiring one candidate to dominate all ties, and
mutual domination among all such choices, correctly identifies a single
projective simple wall if these assertions pass.

The first-endpoint conclusion still has explicit written geometric inputs:
the integral reflection arrangement is locally finite, its initial sector
is the full B9 chamber, and a finite first corner extending A8 has one of
the inventoried finite diagrams. Hence its genuine new simple wall is in
the enumerated set. Conversely every accepted candidate is an actual mirror,
so none can precede the arrangement's first collision. This establishes the
minimum, rather than merely finding a local finite root. The simultaneous
prefix test resolves the chosen sector. The retained A8 walls plus the
selected wall are positive definite by their finite profile; equivalently
the displayed endpoint has square
`-s*(z^2-9s*a^2)/z^2<0` and is future directed.

No independent native-angle census with every other integral mirror is
claimed by the source. The actual simple-sector interpretation uses the
preceding reflection-chamber and finite-diagram facts. It is not inferred
merely from the candidate's norm equation.

## Execution and scope

The source has an internal45-second alarm, bounded JSON, exclusive creation,
source hashes before and after the computation, and canonical JSON replay.
It acquires no inner global lock: the root wrapper must provide the shared
lock, disk/RAM gates and process-tree memory cap. Replay recomputes the same
finite census; it is not an independently implemented arithmetic checker.

Subject to an actual successful execution, the result certifies one first
finite edge in this host. It neither enumerates the remaining chamber nor
certifies finite volume or a compact H9 example.
