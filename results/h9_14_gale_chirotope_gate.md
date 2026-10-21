# The Uniform Rank-Four Gale Gate for the H9/14 Nerve

## Result

Let \(K\) be the boundary complex of a simplicial \(9\)-polytope with
14 labelled vertices. After an arbitrarily small perturbation preserving
the face lattice, \(K\) has a uniform rank-four Gale oriented matroid
\(\mathcal G\) on the same ground set. Its chirotope is specified by

\[
\binom{14}{4}=1001
\]

signs. For every two-set \(S\) and every four-set
\(\{a,b,c,d\}\) disjoint from \(S\), the three signs

\[
\chi(Sab)\chi(Scd),\qquad
-\chi(Sac)\chi(Sbd),\qquad
\chi(Sad)\chi(Sbc)
\]

must contain both \(+\) and \(-\). There are

\[
\binom{14}{2}\binom{12}{4}=45{,}045
\]

such checked relations.

If \(B=\{b_0<\cdots<b_4\}\), define

\[
c_i=(-1)^i\chi(B\setminus\{b_i\}).
\]

Then \(B\) is the complement of a facet of \(K\) if and only if all five
\(c_i\) have the same sign. Thus a proposed complete facet list supplies an
exact finite oriented-matroid feasibility problem: its five-set complements
must be *exactly* the positive circuit supports of \(\mathcal G\). Moreover,
every ground-set element must lie in a positive circuit, because the Gale
oriented matroid is totally cyclic, and every element must also be avoided by
a positive circuit, because every primal label is an extreme vertex.

This is a rigorous necessary polytopality gate. It is not sufficient: an
abstract oriented matroid satisfying these sign axioms need not be realizable
by vectors over \(\mathbb R\).

The gate has now been applied to the explicit abstract completions of H9/14
records 43, 44, and 47. For each record, a frozen 1,001-sign word passes every
relation, and its positive circuits are exactly the completion's 90 facet
complements. Separate exact multiwedge factorizations subsequently prove that
all three nerves are ordinary convex-polytopal. This does not assert that the
particular archived oriented matroids are realizable; it supplies an
independent realization of each face lattice.  A later independent rational
multiwedge realizes record 35, and a later exact certificate excludes record
34. The current theorem-facing partition is (121) impossible, (6)
ordinary-polytopal, and (15) undecided. Coxeter labels and
Lorentzian Gram feasibility remain separate gates; record 35 has since been
excluded specifically from the compact Coxeter track.

## 1. Why the rank is four

Write the fourteen vertices as the columns of their homogenized matrix

\[
A=\begin{pmatrix}1&\cdots&1\\v_1&\cdots&v_{14}\end{pmatrix}.
\]

It has rank \(9+1=10\). A Gale matrix is obtained from a basis of
\(\ker A\), so its vector rank is

\[
14-10=4.
\]

For a 3-dimensional drawing, choose a generic affine chart and reorient each
Gale vector that meets it with the opposite sign.  This produces a
**signed/colored** affine Gale diagram in \(\mathbb R^3\).  The sign or color
cannot be discarded: the Gale oriented matroid is totally cyclic, so it is
not itself the homogenization of an ordinary affine point configuration.
The earlier rank-five formulation was nevertheless off by one; the vector
and oriented-matroid rank is four.

## 2. Why uniformity loses no face lattices

The realization of a simplicial polytope is stable under sufficiently small
vertex perturbations: its finitely many strict supporting and nonsupporting
inequalities remain strict. Inside that open neighbourhood, avoid the
finitely many determinant-zero hypersurfaces indexed by ten vertices. The
resulting homogenized vertex configuration is uniform of rank ten and has
the same face lattice. Oriented-matroid duality then makes its Gale dual
uniform of rank four.

This perturbation changes metric data but not the abstract nerve. It is
therefore legitimate for a combinatorial nonexistence filter, even though a
Coxeter realization itself will usually lie on special algebraic loci.

## 3. Facets are positive Gale circuits

The Gale face criterion says that \(F\) is a face precisely when the origin
lies in the relative interior of the convex hull of the Gale vectors indexed
by the complement of \(F\). A facet of a simplicial \(9\)-polytope contains
nine vertices, so its complement has five elements. In a uniform rank-four
configuration every five-set is a circuit. The origin lies in the relative
interior of its convex hull exactly when the unique circuit dependence can
be oriented with all coefficients positive.

For an increasing five-set \(B=\{b_0,\ldots,b_4\}\), Cramer's rule gives the
circuit coefficient signs, up to simultaneous reversal, as

\[
\operatorname{sgn}(\lambda_i)=(-1)^i
\chi(B\setminus\{b_i\}).
\]

This proves the facet/positive-circuit biconditional used by the checker.
It is important to impose both directions: requiring selected facets to be
positive circuits while permitting additional positive circuits would
silently change the proposed nerve.

The homogenized primal configuration is acyclic: the affine coordinate gives
an all-positive covector. Its dual is therefore totally cyclic. In the
present language this says that every vertex label occurs in at least one
positive circuit. Geometrically, every vertex of the simplicial polytope is
omitted by some facet.

There is a second, distinct condition because all fourteen labels must be
actual extreme vertices, not ghost points. By the Gale face criterion,
\(\{e\}\) is a face exactly when some positive circuit avoids \(e\);
equivalently, \(e\) occurs in at least one selected facet. The standalone
checker imposes both conditions. They are automatic once a genuine closed
pseudomanifold facet list on all fourteen vertices has been certified.

## 4. Exact implementation

The implementation is in
`src/coxeter_polytopes/combinatorics/gale_chirotope.py`.

It has two deliberately separate layers.

1. `add_uniform_rank_four_gale_constraints` and
   `solve_uniform_rank_four_gale` construct the 1001 Boolean signs, impose all
   45,045 uniform three-term Grassmann--Pluecker relations, identify all
   \(\binom{14}{5}=2002\) positive-circuit predicates with the proposed facet
   literals, impose total cyclicity, and require every primal label to be an
   extreme vertex. Z3 is only the search backend.
2. `audit_uniform_rank_four_gale` takes a concrete list of 1001 signs and
   recomputes every relation and every positive circuit directly in Python.
   It neither imports nor trusts the solver.

Global sign reversal does not change a circuit or a Grassmann--Pluecker
product, so the search fixes the first sign to remove this twofold symmetry.

## 5. Polytopal calibrations

The two test fixtures are the dual boundaries of

\[
(\Delta_1)^3\times(\Delta_3)^2
\quad\text{and}\quad
C_5\times\Delta_3\times\Delta_4.
\]

Both are simple \(9\)-polytopes with fourteen facets. Their duals have
fourteen vertices and respectively 128 and 100 facets. The solver finds a
uniform rank-four certificate for each, and the independent checker recovers
exactly the proposed facet complements as its positive circuits. Deleting
one claimed facet makes the *same sign certificate* fail the biconditional;
this tampering check does not assert that every possible certificate for the
altered list is impossible.

For the first product, the tests also include an explicit integral
\(14\times4\) Gale matrix. All 1001 of its maximal minors are nonzero, and
their determinant signs recover exactly the 128 facet complements without
calling Z3. This independently pins the cofactor-sign and positive-circuit
conventions used by both the search and the audit.

Reproduce the calibration with

```bash
PYTHONPATH=.:src pytest -q tests/test_gale_chirotope.py
```

The current exact result is `4 passed`.

## 6. Record 43 and the later completions

Record 43 has missing edges

\[
\{0,11\},\{0,12\},\{1,12\},\{0,13\},\{1,13\},\{2,13\}
\]

and an explicit four-member missing-five family.  Its complete blocker gives
90 five-set facet complements and face vector

\[
(14,85,300,690,1088,1190,880,405,90).
\]

The exact Gale search returned SAT with 92,121 assertions.  Solver-independent
replay verifies all 45,045 Grassmann--Pluecker relations, recovers exactly 90
positive circuits, and checks both total cyclicity and primal extremality.
The packed sign digest is

```text
01f26fdd592d1e6cb6ebd88a3075b2f7e8b2ba442e3cf345a335abdbb1e5dfa7
```

The certificate file SHA-256 and internal critical digest are respectively

```text
2c6463598ab4ae2073b54fc8a9345ab099acba0492c786ba163d0ad5d29d3c66
af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3
```

The independent sign audit is archived in
`data/enumerations/h9_14_record43_selection.json` and replayed by
`scripts/independent_h9_14_record43_audit.py`.  Reproduce the combined
incidence and Gale certificate with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record43.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record43_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record43.py
```

The exact record-43 result is `3 passed`. No real vector realization of this
oriented matroid is supplied. Records 44 and 47 have analogous frozen
certificates and backend-independent audits in
`results/h9_14_record44_selection.md` and
`results/h9_14_record47_selection.md`. All three face lattices are realized
independently by the multiwedge theorem in
`results/h9_14_multiwedge_polytopality.md`.

## 7. Scope and next gate

The global facet-selection model and this Gale model are complementary.
The former enforces ridge degree, Dehn--Sommerville, local strict
\(g\)-inequalities, Macaulay bounds, and the Lanner missing-face order bound.
The latter tests whether a completed facet list can even be the face lattice
of a uniform rank-four oriented matroid dual.

For a final nonexistence theorem, an UNSAT conclusion must be exported in an
independently checkable proof format or replaced by a complete finite
propagation argument. A SAT sign certificate alone passes only the
oriented-matroid gate. For a future completion, one must still establish real
realizability or polytopality by another theorem. For records 43, 44, and 47,
the separate multiwedge theorem now supplies ordinary polytopality without
realizing the archived sign words. Their remaining gates are Coxeter labels
and exact Lorentzian Gram feasibility.
