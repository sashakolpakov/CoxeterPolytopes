# Independent review: all centre-cutting central walls

20 October 2026. Paper/source review only; no scientific execution.

**PASS for the stated conditional completion theorem.** I read the entire
uniform source
`results/compact_ra4_central_wall_center_cut_uniform_extension_20261020.md`,
SHA256 `388f5996be26240fe6f4eb3d1ef7e5bcbae745d197b3314b7ba15a7078a8c942`,
and the complete unchanged fixed-wall proof
`results/compact_ra4_central_ball_completion_volume_cost_20261020.md`,
SHA256 `df5e75b3cfd8466b9b35f36b44f6c6b37f94b25afcdb276e7a443039b1d6da87`.
No mathematical correction is requested.

The hypotheses remain essential: an actual regular dodecahedral facet,
twenty actual outgoing edges of length `arcosh(phi)`, an actual wall
centred at its centre, and compact right-angled completion. No
recognition of the whole 120-cell or universal minimum is assumed.

I independently checked the following bridges.

* **Complete parameter interval.** The base face planes force `r<h`.
  The orthogonal vertex-wall alternative lies beyond that bound, so
  compatibility with every prescribed vertex wall forces
  `p<a-sqrt(3/2)` strictly. Equality would be forbidden tangency.
  Pairing the central wall with the displayed unit matched centre gives
  exactly `p>=sqrt(2/phi)`. Thus both endpoints and the strict upper
  endpoint of the theorem have the stated meaning.
* **Actual central facet and all its unknown neighbours.** The spanning
  sum of the twenty projected vertex normals makes `O` an intrinsic
  interior point. The positive face-average identity has positive
  remainder throughout the interval, since its limiting difference is
  `sqrt(3/2)-sqrt(2/phi)>0`. Hence every unknown neighbour has strictly
  negative product with every base face wall; the compact normal gap
  gives `y_j>=1`. This does not assume a proposed list of central-facet
  vertices or neighbours.
* **Uniform height contradiction.** Direct expansion in the
  `q_v,b_1,b_2,b_3` frame gives
  `t=kS`, `y_j<=lambda*S`, and
  `1=sum(y_j^2)-phi*lambda^2*S^2` for general `p`. The stated endpoint
  values and monotonicities yield `3<=S<=sqrt(15phi)<5`. I checked the
  constrained convex maximum, strict decrease in `lambda`, both endpoint
  evaluations of the remaining convex quadratic, and the exact squared
  difference `4(60phi-97)>0`. Thus the assumption `t<=sqrt(3)` is
  impossible for every actual intrinsic neighbour, including at the
  centre-reaching lower endpoint.
* **Whole-facet count.** The spherical triangles used are the complete
  Euclidean normal fan of the bounded Klein image of the actual central
  facet. Orthogonality gives the displayed determinant identity, its
  decrease in each squared height, and denominator at least `13/4`.
  Consequently each area is below `sqrt(10)/13`. The sphere-area sum
  gives the claimed integer lower bound of fifty vertices. The
  inherited capacity, neighbouring-excess and volume identities then
  give `f>=105`, `W>=23`, and `Vol(P)>=12*pi^2`.

I also checked the limitation in the last section. The added twelve
normals have the exact regular-dodecahedron Gram and are orthogonal to
the central wall. Their products with every prescribed old wall become
strictly below minus one for sufficiently small radius. A strict common
point near the base facet persists. Thus that finite acute prefix really
blocks the stronger local height inference for all small central balls;
it is not asserted to have finite volume or compact completion.

The inherited inputs are the compact gap, capacity, excess and volume
identities in `manuscripts/compact_ra4_volume_constraints.tex`
(`afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`),
and the actual-wall collar convention in
`results/compact_ra4_exact_dodecahedral_collar_cover_20261020.md`
(`c1a3f94c880d37a6180a7b83840eee05c27883eeb8161d77f63ddaa26b942669`).
This review verifies the written extension relative to those pinned
inputs. It leaves small central balls, off-centre balls, unequal outgoing
lengths and the unrestricted minimum-volume problem open.
