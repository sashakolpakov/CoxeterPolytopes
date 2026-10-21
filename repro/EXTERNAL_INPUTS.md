# Third-party inputs for historical research

Third-party papers, source archives and computational supplements are **not
included in Git or the project bulk-data release**. Their local originals
are preserved. The retained provenance files below identify upstream sources;
[external-inputs-manifest.json](external-inputs-manifest.json) records the
original paths, byte lengths and SHA-256 hashes of 37 excluded files, totalling
62,245,573 bytes. This inventory records the local copies, not a redistribution
licence or a claim that the upstream URLs still serve identical bytes.

The six primary papers' selected verification plans do not reference
`data/external/` or the downloaded Chen–Labbé PDF below. Historical scripts
may require these matching
upstream inputs separately. In particular, the catalogue's
`external_enriques` entry points to Martin's downloaded paper and extracted
TeX source; it is an external reference, not a seventh project paper.
The [bulk-data helper](PUBLIC_DATA.md) restores only its own manifest-listed
assets and **does not restore these external inputs**. Neither restoration
route promises a complete historical replay.

## Sources and affected historical routes

| Source | Retrieval and provenance | Historical use |
|---|---|---|
| Shimada, Apéry–Fermi data | [Provenance](../data/external/apery_fermi_shimada_2025/README.md), [Zenodo record](https://zenodo.org/records/14842583) | Apéry–Fermi chamber-fold searches; the extracted `AperyFermiCompData.txt` supplies exact Gram, wall and group data. |
| Brandhorst–Shimada, Enriques surfaces | [Provenance and preprint URL](../data/external/brandhorst_shimada_aut_enriques/README.md) | Published lattice and chamber identifications used by Type-II Enriques certificates. |
| Dolgachev–Shimada, 15-nodal quartics | [Provenance and data URLs](../data/external/dolgachev_shimada_15_nodal_2019/README.md) | Exact wall and symmetry data for 15-nodal quartic chamber searches. |
| Katsura–Kondō–Shimada, characteristic 5 | [Provenance and data URLs](../data/external/katsura_kondo_shimada_char5/README.md) | Supersingular K3 wall/chamber calculations. |
| Kondō, Enriques surfaces | [Provenance](../data/external/kondo_enriques_1986/README.md), [Nagoya record](https://nagoya.repo.nii.ac.jp/records/4843) | Published classification and geometric statements used in historical Enriques audits. |
| Kondō–Shimada, Fermat quartic in characteristic 3 | [Provenance and data URL](../data/external/kondo_shimada_fermat_quartic_char3/README.md) | Exact Fermat-quartic Gram and wall-covector data. |
| Martin, finite Enriques automorphism groups | [Provenance](../data/external/martin_enriques_2019/README.md), [arXiv record](https://arxiv.org/abs/1703.08419) | Published classification, labelled graph figures and the historical `external_enriques` catalogue entry. |
| OSCAR 1.4.0, Type-I Enriques datum | [Provenance and pinned upstream URL](../data/external/oscar_enriques_type_i/README.md) | Serialized lattice/Borcherds context used by historical Type-I and K3-carrier certificates. |
| Shimada, Enriques configuration data | [Provenance and landing page](../data/external/shimada_l10l26note/README.md) | Exact Enriques wall/configuration records. |
| Felikson–Tumarkin, H7 diagram | [Author-maintained diagram](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/7-11.pdf), [local provenance note](../results/fresh_dimension_new_mechanisms_20261014.md) | The historical H7-based construction reconstructs exact geometry from this diagram. The PDF is a third-party reference, not a project figure. |
| Chen–Labbé, Lorentzian Coxeter systems and Boyd–Maxwell ball packings | [arXiv 1310.8608v6](https://arxiv.org/pdf/1310.8608v6), [source-identification note](../results/fresh_all_ideal_seven_simplex_literature_20261017.md) | Figure 14, printed page 24, supplies a human-transcribed historical diagram. No current six-paper verification selection requires this PDF. |

These are dependency examples, not a new completeness audit of every
historical script. Consult the selected script and its source-bound receipts
before reproducing a historical calculation. Published geometric statements
remain separate inputs; retrieving a data file does not prove them.

## Restore the original paths deliberately

1. Follow the relevant provenance link and retrieve the stated version from
   the author, publisher, arXiv or dataset repository. Observe its upstream
   licence and access terms; this project grants no additional rights.
2. Check the download against the manifest before using it. For example,
   `shasum -a 256 FILE` prints its SHA-256 on macOS. Place a matching file at
   the manifest's original relative path. If the upstream bytes differ, retain
   that distinction rather than updating a frozen proof input silently.
3. Retrieve derived members from their matching source archives as described
   below. The manifest also records each extracted member's bytes and hash.

For Martin's paper, use the arXiv record's **Download source** option to obtain
the matching `1703.08419-source.tar`. Inspect its members with `tar -tf` before
extracting it into `data/external/martin_enriques_2019/source/`. The recorded
archive contains `EnriquesFiniteAut.tex`, its `.bbl` and labelled PNG figures.
It also contains upstream `texput.log` and `EnriquesFiniteAut.synctex` build
products. Those two are inventoried for byte provenance but are unnecessary
for scientific replay and need not be restored. Do not compile this external
paper as part of the six primary papers' verification.

For Chen–Labbé, the recorded version is arXiv 1310.8608v6, corresponding to
DOI [10.1007/s10711-014-0004-1](https://doi.org/10.1007/s10711-014-0004-1).
The preserved local copy belongs at
`data/references/fresh_all_ideal_seven_simplex_literature_20261017_chen_labbe.pdf`:
1,619,220 bytes, SHA-256
`1d7a1fd1f20421f9cc7fdeb9b3deacdc05f6ce814e37b7849cd377f72d7f5e99`.
No redistribution permission was recorded, so it was removed from the bulk
payload before any asset upload. The historical checker compares a separately
transcribed graph and hashes the PDF for source identification; downloading
the paper does not verify that transcription.

For the Apéry–Fermi supplement, selectively extract
`AperyFermiCompData.txt` from the matching `compdata.zip`. The unused
369,497,060-byte `AperyFermiFacesData.txt` was intentionally never extracted
into this collection; a full-archive extraction is unnecessary.

The OSCAR datum is the one recorded byte adjustment: the local `.mrdi` copy
has one terminal newline added to the upstream version. Its provenance README
gives the upstream hash, while this manifest hashes the preserved local copy.
An exact historical byte comparison must distinguish these two forms.

No downloads, extraction, scientific checks or manuscript builds are performed
by this document or manifest.
