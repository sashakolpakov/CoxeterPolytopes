# Root review: wall-height propagation and the mixed N5 case

20 October 2026. Independent complete written review of
`compact_corner_wall_height_propagation_20261020.md`, SHA256
`1c2b4a1e56509112aaf9dd680db4181b0217075d94261c9f4ceefd7ecf3c528c`.
PASS. No execution, uniform dimension ceiling or Lean verification.

For a nonincident actual wall the contact vector is nonnegative. The
inverse finite Gram is strictly positive within each connected component.
Any distinct nonincident wall either contacts one of those components,
giving spatial product at least half the minimum displayed inverse
coordinate, or has zero spatial product, giving a negative actual-wall
product of at most-1/2. Acuteness proves the stated height floor in both
cases. Incident old walls are correctly excluded from its scope. The
Schur formula then transfers the height floor to every other edge.

For N5 I independently checked h=(4phi-5)/12 and mu=1/4. Each of the
three replacements is positive definite; acute projection makes its
endpoint actual. The A1 and A2-mate cosh values are correct. Their
half-tangent squares exceed1/6 and1/5 by exactly the two displayed
rational phi bounds. The spectator height estimate is strictly above
40/81, giving every untouched direction a half-tangent square above1/10.
The rank-two energy lower bound yields more than1/15 per untouched
rank. The total exceeds one already at n=13.

The final source corrects an illustrative example identified during
review: the genuine N5 local model is the compact three-dimensional
simplex[3,5,3]. The four-dimensional[5,3,3,3] example used elsewhere has
a different mixed profile. This correction changes no inequality.

The result excludes only this particular profile at a farthest A1/A2
vertex, rather than every local occurrence or every mixed forest. Its
constant may decay on long finite components. Those limitations are
retained, and no compact dimension-12 theorem follows here.
