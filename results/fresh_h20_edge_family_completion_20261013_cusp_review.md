# Independent review of the 96-wall H20 cusp repairs

**PASS.** I independently read the [cusp-repair verifier](../scripts/fresh_h20_edge_family_completion_20261013_cusp_repairs.py), its exact arithmetic and projection helpers, and the [lane report](fresh_h20_edge_family_completion_20261013.md). The following replay completed successfully:

```sh
python scripts/fresh_h20_edge_family_completion_20261013_cusp_repairs.py --verify
```

This is a source review and fresh replay of the existing exact implementation, not a second arithmetic implementation or a new search.

The bounded-link argument is sound. For each connected final component the script checks a positive-definite principal minor, a strictly positive kernel vector of the full Gram matrix, and the exact vector identity `Σ a_i r_i = μ delta` with `μ>0`. The minor and kernel establish positive semidefiniteness of corank one. Passing to the Euclidean quotient `delta-perp / span(delta)` gives a positive relation among every component's normals. The full quotient Gram rank is 19. Consequently, a recession vector pairing nonpositively with every normal must pair zero with every normal, and spanning forces that vector to vanish. The link is bounded. The strict timelike center supplies a strictly feasible horospherical point after rescaling and adding a multiple of `delta`, which does not change incident-wall pairings. Thus the proof does not confuse full rank with boundedness; the earlier links had rank 19 but retained spherical factors.

The replay recomputes all 189 new pairs. Wall 94 has 21 right, one parallel and 72 dotted relations; wall 95 has 19 right, one `pi/3` and 75 dotted relations. Both norms are two and both center pairings are strict. Exact radical sign comparisons also check the interactions with earlier caps in different quadratic extensions.

Both new walls are strictly negative on each of the 107 inherited ideal rays, so those rays and their existing incident links are unchanged. Their previously established boundedness is an inherited certificate dependency; this narrow replay does not reproduce the original 107-link construction. The four already certified finite vertices created at walls 88, 89, 90 and 93 are checked against all 96 inequalities and all survive with their original incidence sets. Their original positive-definite Gram and timelike-point proofs remain dependencies of the earlier independently audited 93/94-wall artifacts.

Finally, the replay checks the saved strict spacelike point directly against all 96 walls, with positive norm and negative pairing with the future center. This proves infinite volume. The common strict center and nonpositive wall pairings prove that every wall is an actual facet by orthogonal projection. The report correctly makes no complete cusp or vertex census, finite-volume claim, or unrestricted nonexistence claim.
