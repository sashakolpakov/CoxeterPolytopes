# R5 period execution and literal-family closure

20 October 2026. PASS: reviewed written proof, fixed exact generation,
and separate fresh replay. No new dimension bound or compact example.

The full mathematical proof is
`results/golden_dyadic_r5_chamber_period_candidate_20261020.md`, SHA256
`e2aeff93b42f87673e55399804df28a9eee3713ac2e7cf064c6f885e5bfce671`.
It has complete root and outside-author written reviews. Both reviewers
also read the complete fixed checker before execution.

Root ran exactly two serial bounded jobs with the shared exclusive lock:

* `results/golden_dyadic_r5_chamber_period_generation_20261020.json`:
  PASS,0.525 seconds,58,507,264 bytes peak process-tree RSS.
* `results/golden_dyadic_r5_chamber_period_replay_20261020.json`:
  PASS,0.506 seconds,57,917,440 bytes peak process-tree RSS.

The certificate is
`results/golden_dyadic_r5_chamber_period_20261020.json`,48,288 bytes,
SHA256 `7e15b328cb2af8321d9417c796dfab38b6046a08f005fec588fcc99e473af2b7`.
The checker SHA256 is
`af564d36a57dd29e1ef0b3c17e30e827aa66d5279381c566d76688fbaebe4658`.
The complete108 coefficient tests, six profiles, both integral lattice
basis inclusions, all signed full-germ images, Gram identities, integral
inverse and trace10+11phi all pass. The replay recomputed the arithmetic
from the fixed sources and matched the certificate exactly.

## Geometric and group conclusion

Let W be the full integral reflection group of literal
R5=<-phi>+I5+D4(O_phi), and C its chamber. The pinned paper proof
establishes the actualness and completeness of the two vertex germs.
Their positive matching by the checked future integral isometry D gives
D(C)=C. Since trace(D)=10+11phi>10 in real dimension10, D has infinite
order: any finite-order real matrix has eigenvalues of modulus one and
real trace at most its dimension.

The cosets D^j W in the integral arithmetic orthogonal group are all
distinct. Equality of two would put a nonidentity power of D in W,
stabilizing C, whereas the reflection group acts freely on its chambers.
Thus W has infinite index. The standard arithmetic-lattice theorem is
an explicit external input; the ambient integral orthogonal group has
finite covolume. Consequently W is noncofinite and C has infinite
volume. This also implies noncompactness.

The separately reviewed full coordinate-facet identity in
`results/golden_dyadic_coordinate_facet_uniform_propagation_20261020.md`
propagates noncompactness to every literal R_n, n>=5, of hyperbolic
dimension n+4>=9. The admissible-form arithmetic compactness theorem,
explicitly assumed there, makes these ambient arithmetic lattices
uniform; a cofinite reflection subgroup would have finite index and
compact chamber. Hence every displayed R_n is noncofinite too.

This closes only those literal integral lattices. It does not exclude
other lattices in their similarity or commensurability classes, prove a
dimension ceiling, or identify a new compact polytope. Stop searches of
this family. Frozen predecessor evidence remains unchanged.

## Verification boundary and resources

The arithmetic was checked by scripts; geometric classifications,
actual-wall provenance, chamber matching and the named arithmetic
lattice statements are explicit written inputs. The existing shared
conditional ChamberPeriod Lean kernel records final group logic, but
no R5-specific Lean module or end-to-end geometric Lean is claimed.
No Lean or TeX job was run in this batch.

Both resource gates passed at177.46GiB free disk and at least7.20GiB
available memory. The watchdog retained its59-second/1GiB limits and
one thread. Python -B created no bytecode cache or temporary proof build
files; only the retained certificate and compact ledgers were produced.
