# Outside-author review of the H3/H4 component energies

20 October 2026. **PASS**, full written source and mathematical review.
No scientific program, enumeration, optimizer, Lean or TeX job was run.
This reviewer did not derive the reviewed H3/H4 proof.

Reviewed source:
`results/compact_universal_h3_h4_component_energy_20261020.md`,
SHA256 `3fc546347f0488200e36dbb6407d3a66ebe4e35bc353603936799841747ff0bb`.
I also read its actual-wall/energy inputs in Section1 of
`results/compact_universal_small_a_component_energy_20261020.md`,
SHA256 `801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`,
and the separately reviewed orthogonal-opposite theorem bound there.
This review does not independently certify the small-A source's later
A1/A2 forest budget.

The common-wall argument has the correct direction of inequalities.
Farthest-vertex half-tangents have sum at most
`sqrt(2e(1+k))`; concavity of the logarithm of `2t/(1+t^2)` bounds their
actual tangent product above. The two displayed scalar bounds are
strictly below k at their stated thresholds and increase to a positive
cleared difference thereafter. Chaining adjacent dual directions thus
forces ONE actual opposite. The per-direction tangent bound gives
`R^2=16/25` and `1584/2209` exactly.

I recomputed both inverse Grams. Their determinants are
`1/(4phi^2)` for H3 and `1/(16phi^4)` for H4; all displayed inverse
entries follow. The adjacent dual correlations exceed the thresholds
used in the proof. All endpoint/white-diagram arguments concern actual
finite endpoints of that common wall, rather than formal continuations.

For H3, adjacent contacts create a nonpositive retained triangle.
The endpoint-pair quadratic has exactly the three displayed positive
coefficients and exceeds one already at both contacts1/2. The three
single-contact computations also match. The only survivor is the
native3 contact at node3, with cost `(phi+2)/4<1` and
`h<=4phi^4/75<phi/4`. The spectator argument is complete: no spectator
gives negative h; a sole positive spectator has cost Q>h; at least two
give `h>=phi/4`. The imported `Q_s<h` holds because each other whole
component is retained at the actual finite endpoint, has at most one
contact, and cannot share this opposite wall. No arbitrary finite
spectator list or stronger small-height claim is assumed here.

For H4, any pair of contacts other than the two endpoints gives a
nonpositive cycle after deleting a node outside its connecting path.
The remaining endpoint pair has the retained four-node determinant
`(5-4phi)/16<0`; hence a single contact is forced. I independently
expanded all four final entries: their numerators are respectively
`4938phi+3688`, `7209phi+5334`, `2500phi+2667`, and
`625phi+1219`, with denominator2209. The last is the tightest and is
strictly greater than one from `phi>198/125`. Larger contacts scale
each positive expression upward, while spectators only increase h.

No gap was found in the proof of H3 energy>1/4 or H4 energy>11/36.
The strictness uses actual finite positive-length edges, a strict
interior comparison point, and the stated strict radical inequalities.
The conclusions have no arithmetic hypothesis and allow arbitrary
other finite components. They neither assert these energies at every
vertex nor yield a uniform dimension ceiling by themselves.
