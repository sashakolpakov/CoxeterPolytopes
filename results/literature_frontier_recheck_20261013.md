# Literature recheck: bounds, algorithms, and recent classifications

Checked on 2026-10-13 before continuing high-dimensional construction.
The dimension bounds and Allcock algorithm scope were checked again on
2026-10-14 against the linked full texts, following the user's arithmetic
dimension warning. The conclusions below are unchanged.

* **Noncompact arithmetic:** dimension at most 19, or exactly 21.
  Esselmann excludes dimension 20 and every dimension above 21; Borcherds
  attains 21. See [Esselmann's original paper](https://doi.org/10.1006/jnth.1996.0141)
  and [Belolipetsky, §1](https://arxiv.org/pdf/1506.03111).
* **Compact arithmetic:** the preceding 21 bound does not apply. Vinberg's
  upper bound is 29. His compact bound also applies without arithmeticity;
  his arithmetic bound also applies without compactness. See
  [Vinberg's original paper](https://www.mathnet.ru/eng/mmo444).
* **Known dimension records:** the introduction of
  [Guo–Ma–Zang–Zheng, July 2026](https://arxiv.org/pdf/2607.14715)
  still gives compact examples only through dimension 8 and noncompact
  finite-volume examples only through dimension 21. Their new infinitely
  many incommensurable families occur in dimensions 4, 5, 6, 7, and 9;
  this does not raise either dimension record. The
  [Felikson–Tumarkin overview](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html)
  lists finite-volume examples through 19 and in 21. In the sources checked,
  dimension 20 remains a missing dimension, distinct from exceeding 21.
* **Vinberg versus Allcock:** Vinberg enumerates reflecting hyperplanes by
  distance from a base point. Allcock traverses the boundary from a vertex
  using efficient indefinite binary-lattice computations. The published
  [Allcock algorithm](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf),
  introduction pp.1–2, explicitly works for lattices over Z. It is not an
  unrestricted nonarithmetic combinatorial search. A geometric search can
  adapt its edge-traversal organization, but must separately solve metric
  realization, all wall-pair constraints, and global closure. This is the
  distinction used in the repository's
  [algorithm audit](vinberg_allcock_search_strategy_20261013.md).

A recent classification update should also prevent obsolete search targets:
[Mcleod, revised 5 October 2026](https://arxiv.org/html/2608.02894v2),
Theorem 1.1, reports that the compact six-dimensional Coxeter polytope with
ten facets is unique, namely the known Bugaenko example. Together with the
earlier dimension-specific results this completes the compact d+4-facet
classification. This is a reported computer-assisted result; its external
certificates have not been rerun in this repository. It neither supplies a
new dimension nor excludes compact H9 polytopes with more facets.

Consequences for the current search: do not enumerate arithmetic noncompact
hosts in dimensions 20 or above 21. Keep nonarithmetic H20/H22 construction
and compact H9 construction distinct. A coherent infinite-volume chamber,
a completed local edge star, or numerical enclosure of an incompatible
wall bank is not a new finite-volume example.
