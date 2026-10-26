How the arguments fit together
==============================

The six paper families use local geometric information to force a global
contradiction or a rigidity statement. Their reductions have different
domains. A scalar inequality, finite list or root model becomes relevant
only after the written argument connects it to every actual case required
by the theorem. The detailed chapters below explain these connections;
the :doc:`verification` page separates the evidence for each connection.

.. list-table:: The principal reductions
   :header-rows: 1
   :widths: 23 39 38

   * - Family
     - Object retained throughout the proof
     - Global conclusion
   * - General finite volume
     - One simultaneous payment rule on actual corners; spherical finite links and affine product links.
     - A local charge bound contradicts the dimension-58 incidence inequality; separate higher-dimensional arguments give the ceiling 57.
   * - General compact
     - One original farthest vertex and one spatial vector, with individually funded references.
     - Complete exclusions at actual ranks 18 through 29, together with the published ceiling, give 17.
   * - Compact arithmetic
     - The same geometric vector, plus actual-field root lines, integral cycles and conjugate positivity.
     - Separate exact-17 and exact-16 exclusions, after the general ceiling 17, give 15 at the stated inputs.
   * - Sharp right-angled dimensions
     - Actual facet projections and a strict height descent.
     - A direct spherical-code endpoint gives 8; named code and metric-equality inputs give the entirely ideal ceiling 4.
   * - Ideal right-angled four-polytopes
     - Facet and ridge defects, actual opposite cusp walls, then recognized cutting surfaces.
     - The next volume and its unique reflected double; a chain theorem under the stated facet-type hypothesis.
   * - Compact right-angled four-polytopes
     - Facet capacity and a distinct excess invariant, including actual metric constraints.
     - Recognition at zero excess and a gap to excess at least 17; unrestricted volume minimality remains open.

One corner rule before either vertex calculation
------------------------------------------------

In the finite-volume argument an angle is an actual incidence of a vertex
:math:`v` and a two-face :math:`F`. Its nonnegative payment is
:math:`w(v,F)`. If :math:`s(F)` is the number of sides, the face obligation
is :math:`\sum_{v\in F}w(v,F)\geq5-s(F)`. Triangles therefore need two,
quadrilaterals one, and nonnegativity handles faces with at least five
sides. Double counting gives

.. math::

   D_5(P):=\sum_F(5-s(F))
      \leq\sum_v C(v),\qquad C(v):=\sum_{F\ni v}w(v,F).

The order matters: the final rule is fixed before maximizing any vertex
charge. Context bonuses, priority choices, repairs and transfers are part
of this same function. A useful corner discovered twice cannot pay both
obligations independently. Whole components, their multiplicities and
their spectators remain visible when a local pattern is priced.

At dimension 58 the finite and ideal branches establish, respectively,
:math:`C(v)\leq377` and :math:`\kappa(v)-C(v)\geq36`. With
:math:`v_f,v_i` the two vertex counts, the incidence lower bound is

.. math::

   D_5(P)\geq377v_f+\sum_{v\text{ ideal}}\kappa(v)+4205.

Subtracting the charge upper bound yields the strict contradiction

.. math::

   0\geq4205+\sum_{v\text{ ideal}}(\kappa(v)-C(v))
     \geq4205+36v_i>0.

This also covers :math:`v_i=0`; an ideal vertex is not required to exist.
The constants have written geometric and finite-data interfaces. Checking
the last subtraction alone does not establish either interface.

Finite vertices have spherical forests. Ideal links are products of
Euclidean Coxeter simplices and require their actual corner-fiber
multiplicities. For an affine rank profile :math:`(r_1,\ldots,r_k)`, the
polynomial :math:`\prod_i(1+t+\cdots+t^{r_i})` records the relevant
product-simplex combinatorics. The concentration argument proves a finite
terminal domain and checks positivity before multiplying a margin by its
profile-dependent factor. These steps are developed in
:ref:`fv-finite-affine` and :ref:`fv-positive-concentration`.

The exclusion of 58 is only one part of the dimension theorem. The
bounded range 59--163 and the unbounded range at least 164 have their own
coverage and hypotheses. The :ref:`higher-tail argument <fv-higher-tail>`
does not assume that an arbitrary facet is itself a Coxeter polytope.

One original vector and one norm budget
---------------------------------------

For a compact acute polytope, choose a strict interior point :math:`z`
and an actual vertex :math:`p` farthest from it. With Lorentz signature
:math:`(-,+,\ldots,+)`, write

.. math::

   z=E(p-x),\qquad E>0,\quad x\in p^\perp,
   \qquad \|x\|^2=\sum_C\|x_C\|^2<1.

The sum is over disjoint whole spherical components at the original
vertex. Their energies add. Overlapping chains, endpoints and heads
within one component share a budget. Both the general compact and the
arithmetic papers retain this same original :math:`p,z,x`.

An actual adjacent vertex
:math:`q=\cosh L\,p-\sinh L\,u`, with :math:`\|u\|=1`, gives

.. math::

   -\langle z,q\rangle
      =E(\cosh L-\sinh L\,x\cdot u)\leq E,
   \qquad x\cdot u\geq\tanh(L/2).

Other feasible future unit points give analogous pressures. Each point
must satisfy every wall of the original polytope, in addition to having
the correct norm and future orientation. A construction on selected
walls does not supply the omitted-wall checks. A local argument at an
adjacent vertex also retains its comparison with the original maximum;
it does not choose a second farthest vertex.

For a reference :math:`W`, completing the square gives

.. math::

   \|x\|^2\geq2x\cdot W-\|W\|^2.

This is useful only when the scalar product is funded by proved
pressures. For example, if

.. math::

   w=\sum_i b_iP_i,\qquad Q=Z-\sum_i t_iP_i,
   \qquad W=w+\mu Q,

then the remaining old coefficients are :math:`b_i-\mu t_i`. The proof
checks :math:`\mu\geq0` and :math:`b_i-\mu t_i\geq0` at every required
anchor before using the new pressure. The full reference norm retains its cross terms.
Zero-sum corrections have an orthogonality benefit only where the
specific reference proves that identity. Two separately optimized gains
cannot simply be added inside the same component.

The ordinary component estimates, A-chain windows and B/D-prefix
arguments have separate inventories. Genuine D-fork directions retain
both signs. Weak knots, native angles, foreign contacts, zero-own
contacts and complete endpoint alternatives remain in their stated
cases. A theorem for a whole diagram is different from a theorem for a
proper component with arbitrary spectators. See
:ref:`compact-paid-references` and :ref:`compact-classical-branches`.

The final general ceiling combines the exact-18 exclusion with a
separate cover at the actual ranks 19--29 and the published upper bound
29. It uses neither a rank decrement nor inherited arithmeticity of a
face. Necessary forest covers contain all possible cases needed for a
contradiction; their rows are not claimed to be realized examples.

Arithmetic restrictions after the common geometry
-------------------------------------------------

Arithmeticity adds conditions on the actual ground field and admissible
quadratic form. The argument uses root lines defined over that field,
cyclic Gram products and positivity at nonphysical embeddings. Physical
unit normalization need not preserve field-valued coordinates, and
disconnected components need not acquire equal rescaling factors. The
field-to-coordinate interpretation is therefore an explicit proof step.

These restrictions filter the common geometric possibilities. They do
not provide a second norm budget. Exact B- and D-type configurations,
including their different field inclusions, retain their own proofs.
Four named computed outputs and the manual whole-:math:`D_{16}` dyadic
argument remain distinct inputs to the final dispatch; see
:ref:`arithmetic-four-outputs` and :ref:`arithmetic-dyadic-finish`.

The final implication has the form

.. math::

   n\leq17,\quad n\ne17,\quad n\ne16
       \quad\Longrightarrow\quad n\leq15.

The first premise comes from the general compact paper. The other two
require their complete arithmetic cases. The literal golden-lattice
cutoff is Bugaenko's published theorem, now authenticated against printed
page 6 of the original paper. Its infinite-volume conclusion holds for
the displayed integral form in every dimension at least eight; the
current proof needs only dimension seventeen. The passage from the actual
configuration to that integral group remains an internal reviewed proof,
rather than an appeal to commensurability. See
:ref:`arithmetic-golden-boundary` and the
:download:`current authentication note <../bugaenko_golden_lattice_cutoff_authentication_20261026.md>`.
The :download:`independent primary-source review <../../results/review_bugaenko_golden_cutoff_primary_authentication_20261026.md>`
confirms the exact statement; it does not reconstruct the entire external
proof. The small Lean implication does not prove Bugaenko's theorem, the
four computed outputs or the geometric exclusions.

Actual projections in the right-angled proof
--------------------------------------------

For an interior unit point :math:`q` and an outward unit wall normal
:math:`n` of maximum height :math:`t=-\langle q,n\rangle`, use

.. math::

   q'=\frac{q+tn}{\sqrt{1+t^2}}.

The proof checks its norm, future sheet and every original wall. For a
distinct wall :math:`m`, acuteness and strict interiority give

.. math::

   \langle q',m\rangle
   =\frac{\langle q,m\rangle+t\langle n,m\rangle}{\sqrt{1+t^2}}<0.

At an orthogonal ridge wall the new height is at most
:math:`t/\sqrt{1+t^2}`. Repeated projections therefore satisfy

.. math::

   A_j\leq\frac{A_0}{\sqrt{1+jA_0^2}}<\frac1{\sqrt j}
      \qquad(j\geq1).

Finite starting height supplies the strict endpoint. In this specific
right-angled setting the written face argument proves the needed
inheritance of right angles, finite volume and entire ideality in face
dimensions at least three. It supplies no general face-heredity principle.

In the general branch, three descents from dimension nine give a code in
:math:`S^5` with pairwise products strictly below :math:`1/4`. The direct
cubic tensor argument gives at most 26 facets, contradicting the named
six-dimensional minimum of 27. In the entirely ideal branch, one descent
from dimension five gives a code in :math:`S^3` with products strictly
below :math:`1/2`. Musin's bound and Kolpakov's metric facet-equality
theorem identify the regular ideal 24-cell; its opposite-wall heights
contradict the strict height bound. Cardinality alone would not supply
metric equality. The :doc:`right_angled` chapter states the exact named
inputs and shows the strict tensor endpoint.

Defect, recognition and the scope of the volume results
----------------------------------------------------------

For an entirely ideal right-angled four-polytope, write :math:`v,f` for
vertices and facets. The ridge defect is

.. math::

   W=3(v-f),\qquad
   \operatorname{Vol}(P)=\frac{4\pi^2}{3}+\frac{\pi^2}{9}W.

Facet defects sum to :math:`2W`. Counting neighboring and opposite cusp
facets requires an actual-wall injection so that a positive facet is
not counted twice. The resulting gap is sharpened by
:math:`W\in3\mathbb Z` to :math:`W\geq12` for a noncell. At equality,
the paper first proves the eight-double-facet spectrum. It then invokes
the cutting and recognition argument under the resulting Oct/D facet
hypothesis. Actual sections, omitted walls, opposite cuts and the
no-twist gluing statement all belong to that proof. This order avoids
assuming that every ideal polytope decomposes into 24-cells.

The compact invariant is different:
:math:`W=\sum_R(k(R)-5)=v-5f`, with the sum over ridges. For a particular
facet with :math:`V` vertices, capacity and excess propagation give
:math:`f\geq2V+5` and :math:`W\geq\lceil3(V-20)/4\rceil`. Further
metric and graph arguments recognize the 120-cell when :math:`W=0`
and prove :math:`W\geq17` otherwise. They do not prove the unrestricted
120-cell minimum or a unique second-volume double. The
:doc:`volume_problems` chapter keeps the two invariants and the two
conclusions separate.

Algorithms and formal kernels have stated interfaces
----------------------------------------------------

A proof-relevant algorithm needs a finite domain with its parameters,
symmetries and guards; exact operations; a termination argument;
soundness and completeness for the checked proposition; and an exact
output with a justified geometric interpretation. These obligations
apply equally to a saved certificate and a freshly generated bounded
list, but the two have different replay scopes.

The :doc:`algorithms` chapter works through the forest certificate,
context partition, affine-profile and application-map arguments. The
:download:`computational proof contract <../computational_proof_contract.md>`
records the shared requirements. An optimizer may suggest a certificate;
the acceptance argument uses exact inequalities or independent exact
checks. A byte hash identifies an input and output. It cannot prove
coverage or geometric realization.

Selected Lean kernels prove internal implications at the stated
geometric, analytic, field and finite-output premises. A successful
compilation and a separate actual axiom audit are useful evidence for
that internal logic. Ordinary theorem parameters remain premises even
when they do not appear among the printed foundational axioms. See
:doc:`formal_logic` for the modules, exact receipts and external boundaries.

Illustrations preserve the reviewed mathematical text
-----------------------------------------------------

The new illustrated reading copies add two explanatory TikZ figures to
each of the three long dimension papers:

.. list-table:: Illustrated editions and retained bases
   :header-rows: 1
   :widths: 31 14 14 41

   * - Paper family
     - Base pages
     - Illustrated pages
     - Added figures
   * - General finite volume57
     - 50
     - 50
     - Simultaneous corner accounting; independent higher-tail ranges.
   * - General compact17
     - 285
     - 286
     - Original-vector accounting; exact-rank conjunction.
   * - Compact arithmetic15
     - 216
     - 217
     - Actual-field filters; exact-17/exact-16 dispatch.

Removing the marked 26 October illustration and package blocks restores
each complete base source byte for byte. The
:download:`independent diagram review <../../results/review_illustrated_manuscripts_diagram_semantics_20261026.md>`
checks the labels, captions, dependency directions and restoration against
the :download:`source proposal <../../results/illustrated_manuscript_source_only_proposal_20261026.json>`.
The diagrams depict logical dependencies or abstract accounting, not
metric realizations or new polytope examples. Their separate matching
build receipts are listed on :doc:`verification`.

The three accepted right-angled editions retain their existing one, five
and four figures respectively. Their mathematical text and PDFs are
unchanged. None of these page counts establishes a minimum essential
proof length: the long papers' case appendices discharge actual cases,
and the arithmetic conclusion also uses the general compact theorem.
