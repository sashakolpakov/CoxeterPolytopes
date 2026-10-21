# Brandhorst--Shimada Enriques-surface lattice source

`AutEnrVol.pdf` is the authors' open-access preprint of Simon Brandhorst and
Ichiro Shimada, *Automorphism Groups of Certain Enriques Surfaces*,
Foundations of Computational Mathematics **22** (2022), 1463--1512,
doi: [10.1007/s10208-021-09530-y](https://doi.org/10.1007/s10208-021-09530-y),
arXiv:2012.10622.

- authors' source URL:
  <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AutEnrVol/paper/AutEnrVol.pdf>
- downloaded: 2026-08-08
- byte length: `302435`
- SHA-256:
  `e11f06d376996bbb2aa3221c6e67ac82a9fa99e9e56096ebc35279cf206e610b`

The version of record is published open access under CC BY 4.0.  The archived
file is the authors' preprint served from Shimada's institutional page; it is
retained here as a pinned primary source, not relicensed by this repository.

The Type-II compiler uses only the following statements:

- PDF pages 3--4, Definition 1.1 and Propositions 1.2--1.3: for an
  ADE-sublattice `R` of `L10`, the half-sum construction `M_R` has orthogonal
  complement `tilde R(2)`;
- PDF page 6, Table 1.1 row 184: `R=D9`, with primitive closure and orthogonal
  ADE type again `D9`, and embedding/chamber type `12B`; and
- PDF pages 6--7, Definition 1.4(ii): for a `(tau,bar tau)`-generic surface,
  the K3 Neron--Severi lattice with its invariant `L10(2)` is isometric to
  this `M_R` construction; and
- PDF page 39, Section 7.4: a very general Type-II Enriques surface is
  `(D9,D9)`-generic, Number 184, and its nef chamber is the `12B` chamber.

These statements calibrate the abstract surviving anti-lattice as `D9(-2)` in
the negative intersection convention.  They do not themselves label the twelve
chosen lifts curve-by-curve, so the compiler separately proves uniqueness of the
wall anti-Gram up to independent sheet switches.
