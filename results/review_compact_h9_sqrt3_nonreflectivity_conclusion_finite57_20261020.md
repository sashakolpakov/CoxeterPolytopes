# Independent review: executed square-root-of-three host exclusion

20 October 2026. **PASS, within the literal-lattice scope below.** This
is a written mathematical/source review and inspection of saved execution
evidence. I ran no checker, import, enumeration, Lean job, or build.

The reviewed closure is
`results/compact_h9_sqrt3_nonreflectivity_conclusion_20261020.md`, SHA256
`00244030f08babd1e07ae7f9ce4db8f9be2ce9ead6cad70a20b9f12b3b82924d`.
I read this final version in full, including the added compactness bridge.
It correctly concludes that the full integral reflection group of the
literal lattice L=<-sqrt(3)>+<1>^9 has infinite index. It makes no uniform
dimension, K-similarity, or commensurability-class claim.

## Executed interfaces

I checked the three certificate hashes against the actual files and read
their six separate generation/replay ledgers. All six record successful
bounded, wrapper-owned-lock runs, with the reviewed worker hashes:

| Interface | Certificate SHA256 | Actual generation / replay |
|---|---|---|
| Opposite A7 edge | `226a53992f880c06976bbd9b27ae7dd1ad2ef8b5ca7bbf3312201b6b01d40e76` | 0.248 / 0.243 s |
| E8 alpha2 edge | `6cb21f095431b3cf62197e93867a278da3f1ea01f35419fea3c96558ba2bc214` | 0.217 / 0.208 s |
| Integral A8 period | `f2eadb4257bebad6b3e16ef7984eaf6a26015ea91e47dd090e0d9372603f2efd` | 0.314 / 0.276 s |

The A7 certificate selects the E8 simple wall from all three simultaneous
candidates by the positive-root cone. The alpha2 certificate has exactly
one accepted candidate, `bridge3_j1`, and explicitly records that the
proposed upper mirror is the first simple wall. Its returned root is

    r=(23+13s)v-(4+3s)omega2-omega1-omega_beta,

and its endpoint is V1=(73+42s)v-(16+10s)omega2, of square 9(1-s).
I independently checked these returned norm, endpoint and old-wall
products against the E8 inverse identities; they agree with the complete
eighteen-profile source reviewed previously. Cap compatibility is a
post-check and does not truncate that census.

The period certificate checks the integral matrix and inverse, Gram
identity, index-nine glue, both root frames and apices, and trace 5+3s.
Its historical `conditional`/`nonreflectivity_claim: false` wording is
appropriate: that worker proves an algebraic isometry. The separate
executed alpha2 interface now discharges its geometric first-endpoint
premise. The closure does not misrepresent that earlier certificate as
having established the premise itself. Replays reuse their producer
implementations; these are not independent enumeration algorithms.

## Global mathematical composition

The full integral splitting and scalar-facet proof, rather than a finite
wall prefix, identify the entire intrinsic arrangement with that of
M=E8(O)+<1-s>. The original A8 vertex V0 and the returned A8 vertex V1
belong to one full chamber C by the two complete edge interfaces. The
integral future-preserving isometry A sends each outward simple root at
V0 to the reversed outward simple chain at V1. Consequently it maps an
inward chamber germ to an inward germ of the same C. Preservation of the
full integral reflection arrangement therefore gives A(C)=C.

The exact trace 5+3s>9 proves infinite order. The full reflection group
acts freely on chambers, so its intersection with this chamber stabilizer
is trivial. Distinct powers of A therefore give distinct cosets: the
infinite-index conclusion for M does not require treating a finite
prefix as the full chamber. The added compactness argument is also
correct: local finiteness gives finitely many walls meeting a compact
chamber, and its finite spanning vertex set makes its isometry group
finite. Thus the same A precludes the compact M facet required if L
were reflective.

The facet transfer uses the proved complete norm inventory: every ambient
integral mirror crossing the selected scalar mirror is orthogonal, and
every M reflection extends integrally to L. Admissibility and anisotropy
make a hypothetical cofinite full L chamber compact, so its scalar facet
would be a compact full M chamber, a contradiction. This step does not
need the preliminary L first-edge search to be replayed. Finally, a
finite-volume subgroup of the arithmetic lattice O+(L) would have finite
index; hence no subset of these integral reflectors supplies a missing
finite-volume H9 example.

No flaw was found in this composition. The universal norm inventories,
actual chamber/face interfaces and arithmetic compactness remain written
proof inputs; the finite certificates do not formally prove them. This
review supplements, rather than changes, the earlier independent source
reviews of the A7 worker, full E8 prefix/splitting, and period/alpha2
worker. No new hyperbolic polytope or dimension ceiling is obtained.
