# Independent review of the D5 fork quadrilateral packet

18 October 2026. **PASS: written local proof, exact witness source, and
fresh source-bound replay.** This review concerns only the right-corner
pattern in
[the packet](fresh_compact21_d5_fork_quad_local_witness_20261018.md), not
every D5 fork quadrilateral or a global dimension improvement.

## Geometry and exact witness

I read the complete report and producer. The seven-root Gram matrix has
the claimed field-valued scaling of the last root; replacing its square
two normalization by square `2Y` does not change its hyperplane or the
opposite Cartan square. The positive original D5 corner and negative
six-root determinant give physical inertia `(5,1,1)` at full corank one.
The conjugate matrix is positive semidefinite with a positive definite
six-root principal minor, giving inertia `(6,0,1)`. The source checks the
generalized Cartan entries in the full ring `Z[phi]`, rather than only
checking rational denominators.

The Schur computation removes the three white roots. Four actual
timelike corner rays satisfy every inactive black inequality strictly
and have negative mutual products. Their convex hull is therefore the
claimed compact quadrilateral. A small negative white-dual perturbation
of the ray sum gives a strict timelike point of the full acute chamber.
This establishes an actual local face, not merely positive principal
submatrices. The 127-subset census excludes affine subsets but does not
establish finite volume. Its 32 skipped subsets all contain the strictly
indefinite two-root pair, which justifies those skips.

## Complete right-corner arithmetic branch

The necessary corank equation is

`(2X-5)(2Y-5)=9`.

Both opposite white-augmented pairs must be indefinite. A positive
definite pair would project a strict center to an actual face: the
inverse of its positive definite Gram M-matrix is nonnegative, so all
other wall inequalities remain strict. A positive semidefinite singular
pair would give an ideal face, contrary to compactness. Consequently
`X,Y>5/2`; zero contacts and the alternate sign branch are not omitted.

Cartan-cycle integrality gives `X,Y in Z[phi]`, without presuming that
every square-two root is field-valued. Conjugate positive semidefiniteness
gives `0<=X',Y'<=8/5`. For `U=2X-5`, `V=2Y-5`, these imply
`UV=9`, positive physical values and `-5<=U',V'<=-9/5`.

The elementary unit argument is valid. The polynomial `t^2-t-1` is
irreducible modulo 3, so a norm divisible by 3 forces divisibility by 3
in the integer ring. Taking norms in `UV=9` then reduces `U` to
`3^a phi^r`, with `a=0,1,2`. The short fundamental-interval argument
indeed proves that every positive unit is a power of phi. Negative norm
forces odd `r`; integrality of `(5+U)/2` forces `r=0 mod 3` because phi
has order three in the residue field modulo 2. The conjugate interval
leaves exactly `(a,r)=(0,-3),(2,3)`. Thus the unordered pair is exactly

`{X,Y}={(3+sqrt(5))/2,(23+9sqrt(5))/2}`.

This covers all exponents, including negative exponents. It is not a
bounded-height experimental classification. One of the two new sides
is the earlier L3/L4 triangle wall. Its products with either E8 triangle
type have squared values `(23+9sqrt(5))/4` and `3+sqrt(5)/2`, respectively.
Both fail integrality in `Z[phi]`. At the specified D5+2E8 vertex this
proves `t_E+8z<=8` for the designated triangle count and this specific
quadrilateral indicator.

## Verification and boundary

A fresh root process ran the producer with `--verify` successfully in
terminal chunk `40ce6c`. The whole gated wrapper took 0.136 seconds;
the producer reported peak RSS 27,983,872 bytes. It held the shared
exclusive computation lock, used its 59-second alarm, and reproduced
the saved certificate exactly. Before starting, the resource check
reported 200.03 GiB free disk and 10.47 GiB available memory.

Producer SHA-256:
`8cd609e79cf06f1bbea59faf57f9338fe35f06407784ce612efa4bc1e3eb59e0`.

This is independent mathematical/source review and a fresh replay of
the author's checker. It is not a second independently implemented
checker, external refereeing, or Lean verification. The unit
classification is a reviewed written proof; the script checks the
explicit witness and local arithmetic consequences. Other corner
patterns, universal replacement payments, and a complete compact
extension remain unproved. No existing vertex allowance or dimension
ceiling is reduced by this packet alone.
