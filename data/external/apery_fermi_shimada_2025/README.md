# Shimada's Apéry--Fermi K3 chamber data

This directory pins the exact computational supplement used for the
Apéry--Fermi chamber-fold census.

- Author: Ichiro Shimada
- Dataset: *The automorphism group of an Ap'ery-Fermi K3 surface:
  Computational data*
- DOI: <https://doi.org/10.5281/zenodo.14842583>
- Zenodo record: <https://zenodo.org/records/14842583>
- Publication date recorded by Zenodo: 2025-02-16
- License recorded by Zenodo: CC-BY-4.0
- Download date: 2026-08-06

Pinned files:

- `compdata.zip`, 19,520,875 bytes, Zenodo MD5
  `51403629e27a00607e4bed40ceeb1f39`;
- `AperyFermiCompdata.pdf`, 80,037 bytes, Zenodo MD5
  `0ef04c58ec0993c34d660322d5dbf638`;
- `AperyFermiCompData.txt`, selectively extracted from `compdata.zip`.

The archive also contains `AperyFermiFacesData.txt` (369,497,060 bytes).  It
was intentionally not extracted: its codimension-at-most-five face tables are
not needed to parse the Néron--Severi Gram matrix, the 80 walls of `D0`, or the
chamber-stabilizer matrices.  The smaller extracted text contains exactly
those objects as GAP assignments `GramSX`, `a32`, `OSXD0`, `AutXD0`, and
`WallsD0`.  The search uses all of `OSXD0` (order 32), not merely the index-two
geometric subgroup `AutXD0`.  An independent weighted-graph automorphism
enumeration verifies that these 32 permutations are the full abstract
symmetry group of the exact wall Gram.

Run

```text
PYTHONPATH=src:. python3 scripts/search_apery_fermi_chamber_folds.py
```

to produce the canonical stock and the exhaustive subgroup-fold census.
