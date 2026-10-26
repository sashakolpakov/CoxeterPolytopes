Sources, datasets and availability
==================================

There are three distinct collections: project source and retained verification
records; the scientific bulk payload; and third-party papers/data. Their
availability and licensing must not be conflated. The source repository and
GitHub Pages site are public. **No bulk asset has been uploaded.** The 9 provenance
READMEs below describe historical third-party inputs; they do not redistribute
the original papers or data and are not required by the six primary papers'
selected verification plans.

Project sources and scientific bulk data
-----------------------------------------

The :download:`bulk manifest <../../repro/public-data-manifest.json>` binds 936
scientific files in 33 archives. The :download:`public-data guide <../../repro/PUBLIC_DATA.md>`
explains archive membership, restoration and byte checks. A source checkout
alone is not the full historical working environment.

The helper's plan is read-only and performs no request or extraction:

.. code-block:: sh

   python3 -B scripts/fetch_public_data.py --plan

Restoration commands become useful only when the actual matching assets are
available. A local archive-directory option does not create public availability.
Recreating identical bytes also does not verify the mathematical domain,
coverage or conclusion of a computation; see :doc:`algorithms`.

Third-party originals
----------------------

The :download:`external-input manifest <../../repro/external-inputs-manifest.json>`
records 37 excluded files totalling 62245573 bytes, with original paths, lengths
and SHA-256 hashes. The :download:`retrieval guide <../../repro/EXTERNAL_INPUTS.md>`
gives the restoration procedure. These are records of local versions, not
claims that an upstream URL still serves identical bytes. This documentation
revision has not redownloaded the entire collection.

For every external input, identify its original version, verify its recorded
hash and observe its upstream licence/access terms. Never change a frozen
input hash silently because a newer paper or dataset is available. If upstream
bytes differ, retain the distinction and reassess the affected calculation.
No new redistribution licence is granted by this repository.

.. _apery-fermi:

Shimada's Apéry–Fermi supplement
--------------------------------

The official dataset is `Zenodo record 14842583 <https://zenodo.org/records/14842583>`_,
DOI `10.5281/zenodo.14842583 <https://doi.org/10.5281/zenodo.14842583>`_. The recorded
publication date is 16 February 2025, with CC-BY-4.0 recorded by Zenodo; the local
retrieval was 6 August 2026. The :download:`provenance note <../../data/external/apery_fermi_shimada_2025/README.md>`
retains the archive/PDF byte lengths and recorded MD5 values.

Selectively extract ``AperyFermiCompData.txt`` from the matching ``compdata.zip``.
The assignments ``GramSX``, ``a32``, ``OSXD0``, ``AutXD0`` and ``WallsD0`` supply
the exact Gram, wall and stabilizer data. The recorded fold census uses all 32
members of ``OSXD0``, not only the index-two geometric subgroup ``AutXD0``.
Its independent weighted-graph enumeration concerns the abstract wall-Gram
symmetries. These are separate group statements.

The 369497060-byte ``AperyFermiFacesData.txt`` was intentionally not extracted.
Its codimension-at-most-five face tables are unnecessary for those selected
objects. Restoring a whole archive indiscriminately wastes disk space without
improving the declared replay.

.. _brandhorst-shimada:

Brandhorst–Shimada's Enriques lattice source
-----------------------------------------------

*Automorphism Groups of Certain Enriques Surfaces*, Foundations of Computational
Mathematics 22 (2022), 1463–1512, DOI
`10.1007/s10208-021-09530-y <https://doi.org/10.1007/s10208-021-09530-y>`_,
arXiv 2012.10622. The pinned
`author preprint <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AutEnrVol/paper/AutEnrVol.pdf>`_
is identified by the :download:`local provenance note <../../data/external/brandhorst_shimada_aut_enriques/README.md>`.
The version of record is open access under CC-BY-4.0; the author-preprint file
is retained as a pinned primary source, not relicensed here.

The historical Type-II route uses Definition 1.1 and Propositions 1.2–1.3,
Table 1.1, row 184, Definition 1.4(ii) and Section 7.4. In the negative intersection convention,
these calibrate the anti-lattice as D9(-2) and identify the 12B chamber. They do
not label twelve selected lifts curve by curve; the compiler's sheet-switch
uniqueness argument is a separate project proof obligation.

.. _fifteen-nodal:

Dolgachev–Shimada's fifteen-nodal quartic data
-------------------------------------------------

The paper is *15-nodal quartic surfaces. Part II: The automorphism group*,
`arXiv 1908.05390 <https://arxiv.org/abs/1908.05390>`_, DOI
`10.1007/s12215-019-00464-7 <https://doi.org/10.1007/s12215-019-00464-7>`_. Official downloads:

* `15nodalcompdata.txt <https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/15NodalQuart/15nodal/compdata/15nodalcompdata.txt>`_;
* `five-page data explanation <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdatatex/15nodalCompData.pdf>`_.

The :download:`provenance note <../../data/external/dolgachev_shimada_15_nodal_2019/README.md>`
records the exact local hashes. ``GramS15`` is rank 16, ``D15WallRecs`` contains 314
published induced-chamber walls, and ``OS15D15permutation`` identifies the
published action with S6. A complete published induced-chamber list is not a
classification of every chamber on the surface.

.. _characteristic-five:

Katsura–Kondō–Shimada in characteristic five
-----------------------------------------------

*On the supersingular K3 surface in characteristic 5 with Artin invariant 1*,
Michigan Math. J. 63 (2014), 803–844, DOI 10.1307/mmj/1417799227,
arXiv 1312.0687.
The official material is:

* `the paper <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/KKSssK3char5.pdf>`_;
* `computational data landing page <https://home.hiroshima-u.ac.jp/ichiro-shimada/ComputationData.html>`_;
* `compdataB.txt <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/compdatassK3char5/compdataB.txt>`_;
* `compdataChams.txt <https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/KKSssK3char5/compdatassK3char5/compdataChams.txt>`_.

The :download:`provenance note <../../data/external/katsura_kondo_shimada_char5/README.md>`
pins the exact local files. ``GramSX`` has rank 22; ``walls[i]`` are the supplied
induced-chamber wall lists. The matrices ``SAutXD[i]`` and ``TAutXD[i]`` factor
the period-preserving group :math:`\operatorname{Aut}_X(D_i)`, not necessarily
the entire lattice chamber stabilizer :math:`\operatorname{Aut}(D_i)`.
Published counts of further chambers do not supply their full wall matrices.
The historical certificate is limited to the three supplied groups.

.. _kondo-1986:

Kondō's 1986 Enriques classification source
-------------------------------------------

The pinned ``ko1821.pdf`` is the 147-page Nagoya dissertation copy of *Enriques
Surfaces with Finite Automorphism Groups*. Use the
`Nagoya record <https://nagoya.repo.nii.ac.jp/records/4843>`_ and the
`journal DOI <https://doi.org/10.4099/math1924.12.191>`_. The
:download:`local note <../../data/external/kondo_enriques_1986/README.md>` records
its 9875342-byte hash identity. The historical audit selects the Main Theorem
and curve counts on internal pages 3–7, the specialized Vinberg criterion on
page 12 and the no-hidden-nodal-curves proposition on page 13. Interrupted
partial downloads were not retained as complete sources.

.. _characteristic-three:

Kondō–Shimada's characteristic-three Fermat quartic
---------------------------------------------------

The supplement accompanies `arXiv 1205.6520v2 <https://arxiv.org/abs/1205.6520>`_.
The original
`computational data file <https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/FermatQuartChar3/autFQ/compdataAutFQChar3.txt>`_
and :download:`provenance note <../../data/external/kondo_shimada_fermat_quartic_char3/README.md>`
identify ``compdataAutFQChar3.txt``: 12763 lines and 1085899 bytes in the recorded
retrieval. The selected ``FQWW*`` arrays are wall covectors. Their scalar
products use the inverse Gram matrix ``FQGramS^-1``. Treating them as root
vectors under ``FQGramS`` would change the mathematical data.

.. _martin-2019:

Martin's 2019 Enriques source and TeX archive
---------------------------------------------

*Enriques surfaces with finite automorphism group in positive characteristic*,
Algebraic Geometry 6 (2019), 592–649, DOI
`10.14231/AG-2019-027 <https://doi.org/10.14231/AG-2019-027>`_,
`arXiv 1703.08419 <https://arxiv.org/abs/1703.08419>`_. The
:download:`local note <../../data/external/martin_enriques_2019/README.md>`
pins both the paper and source archive. The Main Theorem is on PDF page 3;
Table 1 on pages 4–5 supplies the seven types and associated data. The retained
labelled source images help compare scanned graphs, while the exact matrix
audit uses Shimada's arrays.

For a historical source restoration, inspect archive member names before
extracting into the original source directory. Keep the required ``.bbl``:
the corresponding bibliography input is unavailable. The inventoried upstream
``texput.log`` and ``EnriquesFiniteAut.synctex`` are disposable build sidecars
and need not be restored. Do not compile the third-party paper as part of the
six primary papers' selected verification.

.. _oscar-type-i:

OSCAR's pinned Type-I Enriques context
---------------------------------------

The source is OSCAR 1.4.0's
`No. 172 datum <https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/data/TauTaubarGenericEnriquesSurfaces/TauTaubarGenericEnriquesSurfaceNo172.mrdi>`_
and `Enriques loader <https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/src/AlgebraicGeometry/Surfaces/EnriquesAut.jl>`_.
The :download:`provenance note <../../data/external/oscar_enriques_type_i/README.md>`
gives both upstream hashes. The loader's assignment ``SY,SX,L26,w,u`` fixes
which serialized object is the Weyl vector. The local datum adds one terminal
newline; distinguish its manifest hash from the upstream hash.

The OSCAR description gives rank 19 and determinant 4. Geometric existence and
the 12A chamber identification require Brandhorst–Shimada's Section 7.4. A
serialized lattice/Borcherds context is not a set of surface equations.
The upstream package is
`GPL-3.0-or-later <https://github.com/oscar-system/Oscar.jl/blob/v1.4.0/LICENSE.md>`_.
This project does not grant new rights to the third-party payload.

.. _shimada-configurations:

Shimada's L10/L26 Enriques configurations
-----------------------------------------

The `official page <https://www.math.sci.hiroshima-u.ac.jp/~shimada/K3andEnriques.html>`_
and :download:`local note <../../data/external/shimada_l10l26note/README.md>` pin
the note, data explanation and ``L10L26compdata2.txt``. The type dictionary is

.. code-block:: text

   I=12A, II=12B, III=20B, IV=20F, V=20A, VI=20E, VII=20C=20D.

The raw records contain ``GramL10``, ten-coordinate ``walls`` and the complete
``configuration`` matrix. The project parser reads balanced integer arrays
and quoted names; it does not execute the supplied GAP source.

Other literature inputs and restoration limits
-----------------------------------------------

The historical H7 construction uses Felikson–Tumarkin's author-maintained
`H7 diagram <https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/7-11.pdf>`_.
The Chen–Labbé reference is `arXiv1310.8608v6 <https://arxiv.org/pdf/1310.8608v6>`_,
with DOI 10.1007/s10711-014-0004-1. Its figure 14 was human-transcribed for a
historical calculation; a hash identifies the source PDF, not the correctness
of the transcription. That PDF was excluded from the bulk payload before
any upload. The retrieval guide records its precise original path and hash.

Local legacy Lean objects and their trace/hash provenance are retained in a
byte-verified ignored archive. Source checkout alone does not restore it.
:doc:`formal_logic` explains compilation, separate audits, timestamps and the
limits of reusing historical objects. No source/data restoration route here
promises complete regeneration of every historical computation.
