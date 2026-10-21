# Independent review: central-ball completion volume cost

20 October 2026. **PASS within the stated conditional geometry.**
This was an independent full written mathematical/source review; no
scientific program, build, Lean job or geometric enumeration was run.

Reviewed source:
`results/compact_ra4_central_ball_completion_volume_cost_20261020.md`,
SHA256 `df5e75b3cfd8466b9b35f36b44f6c6b37f94b25afcdb276e7a443039b1d6da87`.
The two inherited input hashes were also checked against the actual files:
the compact-RA4 manuscript `afe80767…ad976f` and collar interpretation
`c1a3f94c…42669` agree with the source's pins.

I found no mathematical correction to request. In particular:

* The twenty specified outgoing walls really have unit normals
  `m_v=sqrt(phi)*v-phi*n`. The central wall, its intrinsic point O and
  both displayed projected-normal formulas have the stated norms and
  signs. The three ultraparallel checks are strict, so every actual
  neighbor of the central wall is an additional wall.
* Summing the twenty projected vertex walls proves `-<ell,O>>0` for
  every actual neighbor ell: equality would make ell orthogonal to a
  spanning Lorentz four-space. The face-average decomposition has the
  correct coefficients. Its positive O term proves `<ell,b_j><0` for
  all twelve base face walls, without guessing the unknown neighbors.
  The compact right-angled normal gap then gives `-<ell,b_j>>=1`.
* The dodecahedral support bound gives the stated E. For `t<=7/4`, an
  actual zero product with some `m_v` is forced by the normal gap. In
  the resulting orthogonal frame I independently obtain
  `a=4sqrt(2)S/(X-9)`, `t=h sqrt(17)S/(X-9)`,
  `y_j<=lambda*S`, and `c^2=phi*lambda^2`. The exact loose bounds in
  (1) suffice. The convex-box estimate covers both possible extreme
  point types; it yields a quadratic bounded by `17/40` and `89/90`
  at the endpoints. Thus all intrinsic neighbor heights exceed7/4.
* The bounded Klein image of the **actual entire facet** has a complete
  Euclidean normal fan. It is simple, so the cones give one spherical
  triangle per vertex. Orthogonality gives the claimed Gram determinant;
  its decrease in each squared height and the positive pair products
  justify the half-area branch and bound. Explicitly,
  `tan(Omega/2)<=4sqrt(163)/(53sqrt(65))<32/265`, hence
  `V(F_m)>265*pi/16>49`. No completeness assertion about a finite
  vertex prefix is used.

The inherited capacity and neighboring-excess inequalities apply to this
actual facet with at least fifty vertices. They give `f(P)>=105` and
`W(P)>=23`; their stated volume identity therefore gives
`Vol(P)>=12*pi^2>34*pi^2/3`.

The scope is essential: the regular dodecahedral base, all twenty outgoing
lengths `L=arcosh(phi)`, and the particular actual central wall are
hypotheses. Neither the normal-fan argument nor the prior farthest-facet
theorem forces that configuration in an arbitrary compact RA4 polytope.
This excludes the specified completion branch below the 120-cell volume;
it is not an unrestricted minimum-volume theorem or a nonexistence proof
for larger compact completions.
