# Ideal right-angled four-polytopes: one paper and its proof dependencies

The current publication source is
[Volumes and octahedral chains of ideal right-angled hyperbolic four-polytopes](../manuscripts/ideal_ra4_volume_and_structure_publication.tex).
It keeps the same theorem scope and five figures, adding the finite-model
algorithms and their termination/completeness arguments. Its
[fresh full source audit](../results/review_ideal_ra4_publication_full_finite57_20261020.md)
passes, with the historical classifications and exact outputs explicit
as inherited inputs. The current [17-page PDF](../manuscripts/build/ideal_ra4_volume_and_structure_publication.pdf)
has a [successful source-bound build](../results/ideal_ra4_volume_and_structure_publication_build_20261018.json)
with zero warnings. A fresh combined package replay for this publication
successor remains pending; the predecessor's checks retain their original scope.

The verified [16-page predecessor PDF](../manuscripts/build/ideal_ra4_volume_and_structure.pdf)
and [source](../manuscripts/ideal_ra4_volume_and_structure.tex) remain
byte-identical at their original paths. Their existing exact/conditional
Lean/build records retain that historical scope. The three earlier
standalone papers are also preserved support. There is still one current
ideal-volume/structure paper.

|Result in this paper|Status|Role in the argument|
|---|---|---|
|The regular ideal 24-cell uniquely minimizes volume, at 4π²/3 (`thm:minimum`).|Published theorem of Kolpakov; explicitly attributed. The paper also gives a direct all-octahedral recognition argument.|Supplies the known minimum and normalization. Recognition of pieces in the cutting proof is proved directly.|
|The next distinct volume is 8π²/3, attained uniquely by the reflected double (`thm:gap`).|Internal written proof, two consolidation reviews, exact supporting checks and conditional Lean verification. No external refereeing or end-to-end geometric formalization claimed.|A defect argument proves the gap. A separate equality argument forces exactly eight doubled-octahedral facets and otherwise octahedral facets. The general chain theorem then gives exactly two cells.|
|Octahedral/doubled-octahedral facets force a chain of 24-cells (`thm:chain`).|Internal restricted structure theorem, with written geometric proof and exact finite support. Conditional Lean verifies the stated composition from explicit geometric premises.|Constructs actual octahedral cuts and proves recognition, opposite intermediate cuts and canonical gluing. It applies without a volume hypothesis, but requires the stated facet spectrum.|

The shared proof proceeds as follows:

1. **Defects and cusp injection** (Sections 2–5): the volume formula and
   disjoint lists of positive-defect neighboring/opposite facets exclude
   every volume strictly between 4π²/3 and 8π²/3.
2. **Equality spectrum** (Sections 6–7): at 8π²/3, the small-facet and
   cusp-shape arguments force precisely eight doubled octahedra. This
   step is independent of the later decomposition theorem.
3. **General structure** (Sections 8–15): rectangular cusp boxes produce
   closed internal surfaces, a six-cusp graph and eight wall normals.
   Lorentz identities place the triangles in an actual octahedral
   section. Cutting, direct all-octahedral recognition and the opposition
   argument give a chain. With m cells its face vector is
   (18m+6, 72m+24, 68m+28, 14m+10).
4. **Application and attainment** (Section 16): eight doubled facets mean
   two cells. The gluing argument identifies their placement as a
   reflection double; an exact halfspace model verifies attainment and
   the face vector (42,168,164,38).

The cube-normal and cutting argument is proved once, in the general
structure theorem. The second-volume theorem invokes it only after
establishing its facet hypothesis. Nothing here classifies arbitrary
ideal four-polytopes as chains or proves universal divisibility of their
volumes by the 24-cell volume.

## Assumptions and reproducibility

The published inputs are stated in Section 1: Kolpakov's volume and
incidence formulas, antiprism bound/equality, the specified complete
nine- and ten-vertex three-dimensional classifications, and
three-dimensional right-angled rigidity. No four-dimensional census is
assumed. See the [primary-source scope audit](../results/ideal_ra4_primary_literature_scope_audit_20261018.md)
and [assumption register](manuscript_assumptions.md).

From the repository root, inspect the plan or run the declared checks:

```sh
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --plan
python3 -B scripts/verify_manuscript.py ideal_ra4_volume_and_structure --mode core --pdf
```

The driver runs checks serially with the repository's resource guards.
`--pdf` builds in a temporary directory and preserves the existing PDF.
Add `--ledger results/ideal4_new_verification.json` to retain a fresh
execution record; that path must not already exist. `--mode full` also
requests the separately marked historical double-model replay. It does
not enlarge the scope of the proof.

The [completed core ledger](../results/ideal_ra4_unified_core_verification_final_20261018.json)
records four exact checks and the two declared Lean modules. The checks
cover the saved double face lattice, small equality-model data and cube
algebra, and the D4 root/opposition model. They do not replace the
universal geometric arguments. The [build ledger](../results/ideal_ra4_volume_and_structure_build_20261018.json)
records 16 pages, five figures and zero TeX warnings. Its source hash
matches the core ledger and the preserved16-page predecessor source.
The publication successor has the separate 17-page build linked above.
It does not have a fresh combined core record; the earlier core ledger
does not acquire the successor's source binding from that build.

The [conditional Lean interface](../formal/CoxeterProofs/ManuscriptIdeal4.lean)
separates the published minimum, internal gap partition, equality-spectrum
stages, cuts, recognition, path structure and gluing. Its compiled axiom
audit has no `sorryAx`; actual hyperbolic geometry and the named script
premises remain explicit hypotheses. The [root review](../results/review_ideal_ra4_volume_and_structure_root_20261018.md)
and [independent consolidation review](../results/review_ideal_ra4_unified_manuscript_resource_20261018.md)
cover the complete written manuscript.

The subsequent all-antiprism and Figure7 exclusion reports are separately
reviewed research extensions. They are not in this checked revision;
their eventual integration should share this paper's geometry rather
than create another standalone manuscript.

The [six-paper audit receipt](../results/six_papers_audit_edition_20261021.json)
pins this guide before the build-status correction above. That exact
6285-byte predecessor is [preserved here](../results/source_archive/ideal_ra4_proof_map_before_build_status_correction_20261021.md.gz),
with gzip round-trip verification. The correction changes no manuscript,
PDF, mathematical result, selected check or historical verification claim.
The [documentation correction receipt](../results/ideal_ra4_proof_map_build_status_correction_20261021.json)
binds the predecessor and the corrected guide.
