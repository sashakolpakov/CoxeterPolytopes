# Independent mixed-source D4 audit

**PASS.** The separate [checker](../scripts/fresh_h9_mixed_d4_coupling_20261014_review.py) and [certificate](fresh_h9_mixed_d4_coupling_20261014_review.json) independently reconstruct the [frozen construction](fresh_h9_mixed_d4_coupling_20261014.json), with no producer or repository imports. Arithmetic uses `Fraction` in `Q(sqrt(5))(sqrt(R))`, where `R=(16+4sqrt(5))/11` and both square roots have their positive real values. Exact signs are reduced by conjugation and rational comparisons; no numerical tolerance enters a predicate.

| Cone | Supports enumerated | Finite extreme rays | Spacelike extreme rays |
|---|---:|---:|---:|
| Compact H5 source | 21 | 10 | 0 |
| Coupled 11-wall H9 seed | 55 | 10 | 20 |
| Final 12-wall H9 chamber | 220 | 12 | 38 |

There are no ideal extreme rays in these cone censuses. This does not assert that their ideal boundaries are empty. Every exact projective vector, squared norm and complete active-wall set matches the primary certificate.

## Independent construction and complete support census

The checker rebuilds the five-node norm-two Gram chain, its inverse, `g=−G⁻¹e4`, `D0=(7+sqrt(5))/4`, the two different prism ends, the D4 Gram block and its outer weight of square one. It verifies the physical coupling `n6=W6/2−omega`, all eleven norm-two identities, the strict future timelike center and signature `(9,1)`.

Unlike the producer's incremental cone clipping, this audit enumerates every support of `dimension−1` walls and computes its nullspace by forward elimination and back-substitution. It tests both orientation and all remaining inequalities exactly, then deduplicates positive projective rays. Explicit independent defining rows establish pointedness. The timelike sum of normals is negative on every feasible nonzero ray and gives a bounded projective slice. Consequently every extreme ray has a rank-`dimension−1` active subset and is encountered. Singular support subsets are not mistaken for rays.

The compact H5 source is reconstructed directly rather than accepted from an archived source census. Its ten extreme rays are all future timelike, proving compactness of this source chamber.

## Every literal polar and every final facet

The first polar bank contains exactly the twenty spacelike seed rays, without omissions or duplicates. Exactly one passes all old Coxeter pair tests. After adding that normal, the complete new cone census supplies exactly thirty-eight spacelike rays; each is tested against every one of the twelve walls and rejected. The audit matches every archived pairing, squared weight, rejection list and accepted index.

The final sixty-six pair relations are **53 right, nine label 3, one label 5 and three dotted**. All twelve walls have independently checked future timelike relative-interior points, obtained by projecting the strict center onto their hyperplanes; every other inequality is strict there. Thus the additions are actual facets rather than redundant inequalities.

The saved rational escape has positive square, negative center pairing and twelve strictly negative wall pairings, all independently reproduced. Its segment to the strict timelike center crosses the light cone through strict inequalities, giving an open ideal-boundary patch and infinite volume. No finite-volume or compact H9 construction is established.

The exhausted procedure is specifically **polar normals of the current spacelike extreme rays**. It does not exhaust arbitrary cap normals, replacement walls or different chambers.

## Completeness of the angle test

The quartic field is not normal. If it were normal, it would contain both square roots of the conjugate values of `R`, and therefore `sqrt(11)`. It would then equal `Q(sqrt(5),sqrt(11))`. A square root of `R` in that field would require either `R` or `R/11` to be a square in `Q(sqrt(5))`. Their rational norms are nonsquares, a contradiction. The code independently checks the relevant norm data; in particular `Norm(R)=16/11`.

Any second quadratic subfield would make the quartic biquadratic and normal. Thus its only nontrivial abelian subfield is `Q(sqrt(5))`. Squared finite Coxeter cosines must lie there, leaving precisely labels `2,3,4,5,6,10`. The two degree-two alternatives involving `sqrt(2)` or `sqrt(3)` are absent. Therefore the literal polar rejections have no untested integer-label tail.

## Conditional quasi-arithmetic obstruction

The separate arithmetic verifies the two Cartan cycles

```text
cycle(0,1) = (3+sqrt(5))/2,
cycle(5,6) = 2+sqrt(5)/2 + (−1+3sqrt(5))sqrt(R)/8.
```

They generate the quartic field. The automorphism fixing `sqrt(5)` and negating `sqrt(R)` changes the second cycle while fixing the entire norm-two five-wall principal Gram block, whose signature is `(4,1)`.

The external theorem dependency is Vinberg's criterion as stated in [Bogachev–Yorov, Theorem 2.2, p. 6](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf): for a finite-volume quasi-arithmetic Coxeter chamber, the Gram field is totally real and embeddings nontrivial on the cyclic field must give positive-semidefinite conjugate Gram matrices. This formulation was checked directly; it is used, not reproved, here.

If a finite-volume completion retaining the seed facets were quasi-arithmetic, this embedding would extend to its totally real Gram field. The extension remains nontrivial on the cyclic field and fixes the same indefinite principal block. That contradicts the criterion. Dividing this norm-two block by two gives the usual unit-normalized block and preserves its signature. Hence such a completion would be non-quasi-arithmetic; the current infinite-volume chamber is not being called a lattice, and the argument does not establish that a completion exists.

Generation and a fresh saved replay pass:

```sh
python scripts/fresh_h9_mixed_d4_coupling_20261014_review.py --verify
```

The manifest binds primary source/JSON, both primary arithmetic dependencies and the independent checker. The producer's prose report is not hashed, so a review link can be added safely.
