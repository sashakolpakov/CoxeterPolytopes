# Root review of actual spherical-face projection pressure

20 October 2026. Full paper review of
`compact_farthest_spherical_face_projection_pressure_20261020.md`, SHA256
`abfcb96b6eb74aefdc6498cda362d81b1d4a1f9c32ad17fe3b464cff26fba565`.
No computation, build, end-to-end formalization or dimension change.

For c=a*p-v, the old selected projection w has c+w orthogonal to
the old J, norm delta, and time pairing-a. Projection of p onto
J+c gives q0=(1+a^2/delta)*p-(a/delta)*(v-w). Normalization therefore
gives exactly the stated v_J,R_J and the right side H(a,delta).
Inverse positivity makes w a nonnegative old-normal combination.
Since the center z is strictly interior, x.oldnormal>0, so x.w>0
whenever its contact vector is nonzero. The old actual mixed-projection
lemma supplies feasibility; an arbitrary abstract Gram is not enough.

For the isolated two-wall specialization, a^2+delta=beta^2(D-1).
Dividing H by beta*sqrt(D) gives the stated strict floor S. With
t=tan(theta), its square is exactly
(D-1)/D*(sqrt(D-1)-t)/(sqrt(D-1)+t).
Comparison with the LONG half-tangent S_long requires
S_long^2=(sqrt(D)-1)/(sqrt(D)+1), with the square retained.
After cancellation, S>S_long is equivalent to the displayed bound
t<F(sqrt(D)). Differentiating F^2 yields a positive numerator
2q(4q^3+18q^2+18q+5). At sqrt6 its lower comparison is
331>116sqrt6; the comparison to tan^2(pi/5) uses
sqrt5>67/30. Both elementary exact inequalities hold.

For each selected one-sided A3, its inverse Gram column is
(3/2,1,1/2), giving precisely the H4 residual and correction. Far
contacts do remain in the entire vector v and in a; this result
does not replace their spatial contribution by beta*omega_i.
For two sides, the correction norms and cross product are respectively
(3/2)beta^2 and(3/4)beta^2. Expansion of the average therefore gives
||y_bar||^2=a^2+1-(15/8)beta^2, as stated. The strict Cauchy bound
uses the same farthest-point x of norm less than one.

The A14 barrier comparison is also exact: its minimum central formal
tau^2 is149769/280000>1/2, while D<=15/2 gives the corrected long
floor below1/2. Thus the informal hope that the corrected scalar floor
alone would exclude A14 is explicitly withdrawn. The actual spherical
face vector constraints contain genuinely additional information, but
the source does not claim they already exclude a whole vertex type.

Verdict: PASS for the stated identities and necessary inequalities.
Root co-derived the general projection mechanism; this review is not
external refereeing. Complete actual-profile coverage and a simultaneous
energy contradiction remain required for a numerical dimension bound.
