# A coupled affine-A20 repair at the exposed H22 cusp

No finite-volume example results. This calculation tests a new local repair,
outside the eleven transported wall families excluded previously. An exact
dual certificate shows that every wall in this proposed repair at the fixed
cusp leaves the current global escape. Repairs at other cusps, other incident
diagrams, and nonincident cutting walls remain outside the conclusion.

## Construction being tested

At the lower placement's null ray `delta4`, the finite components have types
`A5 + A3 + A3 + A4 + A1`; a separate affine `A1` is already present. There
are four unused positive tangent directions. One can try to join the five
finite chains into a single affine `A20` cycle by inserting five norm-two
walls. Together with the old affine `A1`, that would give a bounded local
link of rank 21 if the full compatible geometric realization is achieved.

An individual bridge pairs to `-1` with an endpoint of each of two distinct
chains and to zero with every other incident old wall. There are exactly
32 such endpoint pairs. Its residual vector varies on a sphere in the four
unused directions; adding an arbitrary real multiple of `delta4` does not
change its norm or incident pairings. The calculation covers this entire
continuous family, not a sampled list of directions or shifts.

## Exact obstruction

Let `X` be the already certified strict spacelike escape of the 5,397,548-wall
lower placement, and let `B` be the positive definite Gram matrix of the
16 finite-component normals and old wall 26. For a bridge, its prescribed
pairing vector `p` has two entries `-1` and the rest zero. Put

\[
\rho=2-p^{\mathsf T}B^{-1}p>0.
\]

For each endpoint pair, the certificate specifies one or two old walls and
positive coefficients `lambda_j`, with

\[
z=X-\sum_j\lambda_j w_j,\qquad (z,\delta_4)=0.
\]

Write `z_E` for the column of its pairings with the 17 basis normals. Define

\[
\kappa=p^{\mathsf T}B^{-1}z_E,\qquad
\zeta=z^2-z_E^{\mathsf T}B^{-1}z_E.
\]

All 32 certificates verify exactly in `Q(sqrt(11))` that

\[
\kappa<0,\qquad \zeta\ge0,\qquad
\kappa^2-\rho\zeta>0.
\]

The Lorentz form has signature `(22,1)`. Cauchy–Schwarz in the positive
quotient of `delta4`-perpendicular space therefore gives, for every such
bridge `r` having nonpositive pairings with the old walls,

\[
(X,r)\le (z,r)\le\kappa+\sqrt{\rho\zeta}<0.
\]

Thus this escape survives any finite collection of these bridges at this
fixed cusp, including every coherent five-bridge affine-A20 completion.
The proof only needs acute inequalities against the original 210 walls.
There is no restriction on the field, finite angle labels to nonincident
walls, real shifts, or tangent directions of a proposed bridge.

The 32 numerical optimizations were used to find certificates. They are
not part of the proof: the verifier reconstructs the rational Gram matrix,
checks its signature and the exhaustive endpoint list, then verifies every
coefficient and sign by rational arithmetic and squaring.

## Reproduction and limits

Run `python3 scripts/fresh_h22_coupled_cusp_20261013_verify.py --verify`.
The [certificate](fresh_h22_coupled_cusp_20261013.json) records all 32 dual
inequalities and hashes its inputs. The [verifier](../scripts/fresh_h22_coupled_cusp_20261013_verify.py)
does not call the numerical optimizer.

This is an obstruction to a specified local completion cutting a specified
global escape. It is not nonexistence of H22 Coxeter polytopes, not exclusion
of other cusp-completion diagrams, and not an obstruction to adding further
walls away from this cusp or transporting repairs to other cusps.
