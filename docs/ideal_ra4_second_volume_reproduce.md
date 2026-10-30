# Reproducing the ideal right-angled four-dimensional second volume

For current reading, see the [current second-volume manuscript](../arxiv/ideal_four_dimensional_volume_and_structure/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The theorem is about the **second distinct volume value**, $8\pi^2/3$, in the class of finite-volume entirely ideal right-angled hyperbolic four-polytopes. The unique first minimizer is the ideal24-cell, by Kolpakov. The frozen proof does not classify all equality cases at the second value.

The supporting argument translates the [frozen universal proof](../results/ideal_ra4_second_volume_20261018.md). The [lemma map](ideal_ra4_second_volume_lemma_map.md) separates published inputs, new geometric arguments, exact arithmetic and the partial Lean development. No optimization, large four-polytope census or numerical tolerance enters this proof.

## Exact construction replay

From the repository root, with Python3 and its standard library:

```sh
python scripts/fresh_ideal_ra4_double24_volume_20261018.py --verify
python scripts/check_ideal_ra4_second_volume_saved_20261018.py
```

The first command rebuilds the rational construction and compares it to the saved JSON without writing it. The second is a separate saved-data checker importing no producer. It reconstructs the original active-wall vertex census, reflection, complete face-intersection lattice, cusp incidences and doubled-facet aspect ratios, and checks the antiprism rational identities. Expected results include the double's `f_vector=[42,168,164,38]`, ridge defect12 and volume coefficient `8/3` in units of (pi^2).

The construction JSON's `second_minimum_proved=false` is intentional: that computation alone establishes the example, not the universal lower bound. The latter is the written geometric counting argument. Likewise equality classification is not certified by this replay.

Recorded final runs are producer generation `f42b71`, independent replay `a81264`, and hardened standalone checker `cffdc9`, all exit0/PASS. The integrity/status audit is [ideal_ra4_second_volume_20261018_audit.json](../results/ideal_ra4_second_volume_20261018_audit.json); it is not an additional geometric proof.

## Frozen artifact identities

|Artifact|SHA256|
|---|---|
|[Producer](../scripts/fresh_ideal_ra4_double24_volume_20261018.py)|`16731c71d949e287d4cdce85e0ef98f8e5551beaaef7a6acfdd196c1fb0f5323`|
|[Certificate](../results/fresh_ideal_ra4_double24_volume_20261018.json)|`386d9650b8f669e3edcb3cdee56ed53d1c9ea4a90af2b4ba5c8c1b83f9d7bfa3`|
|[Standalone checker](../scripts/check_ideal_ra4_second_volume_saved_20261018.py)|`05ebb4f12e7413a13c82da8be366fe874933c93c3f6e57e8bc9bb4e932f64c09`|
|[Universal proof](../results/ideal_ra4_second_volume_20261018.md)|`9dff11c375bddcf62f16f3f01ddabd7f5b5c137756cf02ebff275076e6869c25`|

These files are immutable. A later equality argument or stronger theorem belongs in a separate packet.

## Building the manuscript

With a TeX installation containing `latexmk`, `pdflatex`, the AMS packages, `booktabs`, `geometry`, `hyperref` and `xurl`:


This produces `build/ideal_ra4_second_volume/ideal_ra4_second_volume.pdf`. A successful typesetting run is not a mathematical verification.

The manuscript was built successfully with TeX Live2026 and latexmk4.87 on18October2026 (initial build `f74b0a`, five pages). It contains only the frozen second-value theorem, not the later equality-structure candidate.

## Lean scope

The root-owned [IdealRA4Gap module](../formal/CoxeterProofs/IdealRA4Gap.lean) formalizes the integer consequence of the geometric hypotheses. Its `defect_gap_of_partition` takes the partition of facet defects, the number of selected facets, their minimum defect, and the exhaustive facet-case conclusion as hypotheses. These are not filled by declaring the new opposite-cusp lemma an external theorem. The geometric argument and exact construction still require the written proof and the checker.

To check the current formal project using its pinned toolchain, follow `formal/lean-toolchain` and run `lake build` from `formal/`. This guide does not claim that building that project alone establishes the full hyperbolic theorem.

## Review status

Two independent agents read the universal proof, producer and standalone checker, and each inspected the relevant published classification and antiprism premises. Neither review was an independent computation or external human refereeing. The manuscript is a translation of that internally checked proof. New equality-structure research is separate and must not silently strengthen its statement.
