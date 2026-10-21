# Root review: C26 geometry and fixed exact checker

20 October2026. Complete independent written/source PASS before execution.

* Proof: `compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md`,
  SHA256 `1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f`.
* Worker: `scripts/check_compact_ra3_c26_graph_metric_20261020.py`,
  SHA256 `47108a5ff3f732bf4b54d7893cc925f8e2a9423106c0ac812f690d355a9976e4`.
* Preparation: `compact_ra3_c26_graph_metric_preparation_20261020.md`,
  SHA256 `09560b0d19d188bff51b5acd26c852256cbd34f6483da92a303ea4ea6aa3b35a`.

I independently checked the primary BGM2012v2 Table2, PDFp15, row26:
https://arxiv.org/pdf/1207.7010v2 . Its fullerene count is exactly one.
The face-capacity argument excludes every face of size>=7 before that
classification is used. The supplied dual has15 vertices,39 edges and26
triangles; its displayed links are cycles of the indicated sizes, so it
is a sphere and exhibits the required unique fullerene type.

The geometry starts from an ACTUAL compact RA3 polyhedron, not merely
the formal model. Andreev uniqueness realizes its graph automorphisms.
The averaged interior point, the faithful order-three rotation and the
horizontal/vertical involutions give the displayed actual normal orbits.
Distinct H faces ensure nonzero horizontal projections, and distinct
primed A faces ensure nonzero axial components. Orthogonality fixes the
remaining radial signs and phase. Thus the symmetry argument does not
assume an order-five extension or infer compactness from a partial Gram.

I rederived all equations. The A triangle gives U^2=2/3. B with its
horizontal reflection gives B^2=1/2. The H adjacencies fix a^2 and b^2;
the A--B orthogonality gives AB=ab/(1+X)>0 and then9X^2+3X-4=0.
The endpoint-third normals of A_i--A_(i+1) give squared cosh
3X/(3X-1)=(7+sqrt17)/4. For B_i--B'_i they give cosh
(1+3X)/2=(1+sqrt17)/4. The elementary strict comparisons separate
these from phi^2 and phi, respectively. Both actual pentagon orbits
therefore contain a nonregular edge. No numerical realization is used.

The RA4 consequence is valid relative to the pinned earlier geometry:
a positive neighbor of D cannot have excess2 by L6 or excess3 by this
theorem. It has at least28 vertices; capacity gives f>=61. The all-D
case has120 facets. Together with W>=17 outside that case and the
exact volume formula, this gives Vol>=47pi^2/6 universally. The
120-cell case itself satisfies this inequality.

I read the entire worker without importing it. Its Q(sqrt17) Fraction
arithmetic, exact sign test and norm division check. The fixed graph,
triangle incidence, all cycle links, twelve automorphisms and witness
endpoint checks have bounded finite loops. The field identities check
squared orthogonality data and the two edge formulas. The external
count-one classification and actual uniqueness/symmetry geometry are
explicitly outside its scope. It neither enumerates fullerenes nor
checks a new hyperbolic realization. All output is deterministic and
bounded; input/control/self hashes are checked before and after, an
owned parent is required, and the certificate is fresh-only and<=64KiB.
The existing parent must supply the exclusive lock, resource gates,
59-second/1GiB process-tree watchdog and one-thread environment.

This review authorizes the two prescribed tiny root-only bounded jobs.
It records no execution outcome or Lean coverage. Unrestricted120-cell
minimality, the second volume and worldwide priority remain open.
