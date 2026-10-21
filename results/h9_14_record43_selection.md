# An explicit abstract completion of H9/14 record 43

> **Subsequent whole-shadow theorem.**  The complete compact-candidate census
> leaves this selected nerve and one second incidence orbit; both now have
> unrestricted compact-label obstructions.  See
> [`h9_14_record43_whole_shadow_compact_no_go.md`](h9_14_record43_whole_shadow_compact_no_go.md).

## Result

H9/14 second-layer survivor record 43 has an explicit abstract facet
completion.  Its graph6 string is

```text
M?????????O?o?w??
```

and its missing edges are

\[
\{0,11\},\quad \{0,12\},\quad \{1,12\},\quad
\{0,13\},\quad \{1,13\},\quad \{2,13\}.
\]

Adjoin the four missing five-faces

\[
\begin{aligned}
Q_0&=\{3,4,5,6,7\},\\
Q_1&=\{1,2,8,9,10\},\\
Q_2&=\{2,8,9,10,11\},\\
Q_3&=\{8,9,10,11,12\}.
\end{aligned}
\]

The complete blocker of these ten minimal nonfaces consists of exactly
90 sets, all of order five.  Taking their complements gives 90 nine-vertex
facets.  The generated complex has face vector

\[
(14,85,300,690,1088,1190,880,405,90).
\]

Thus record 43 moves from a bounded timeout to an **explicit globally valid
abstract completion**. A subsequent exact multiwedge factorization proves
that this same face lattice is the nerve of an ordinary convex simple
9-polytope; see `results/h9_14_multiwedge_polytopality.md`.

The blocker/Gale certificate in this report alone establishes the incidence
and oriented-matroid statements listed below. The independent multiwedge
theorem supplies ordinary polytopality but no Coxeter labels or Lorentzian
Gram matrix. Hence this is still not a compact hyperbolic Coxeter example.
The current global partition is (121) impossible, (6) explicit
ordinary-polytopal selections, and (15) undecided; record 35 is the later
sixth realization, and record 34 is the later exact exclusion.

## Exact blocker construction

Let \(\mathcal N\) be the clutter formed by the six missing edges and the four
sets \(Q_i\).  A facet complement must be an inclusion-minimal transversal of
\(\mathcal N\).  The certifier tests every one of the \(2^{14}\) vertex masks,
without a SAT backend, and obtains

\[
\#\{T:T\text{ is a minimal transversal},\ |T|=5\}=90
\]

and zero minimal transversals of every other order.  Consequently every
maximal face has nine vertices and the complex is pure of dimension eight.

The frozen list of all 90 masks has SHA-256

```text
e93bd0a0353dc03b18ba499454c08f7f267637a18ccb0e953d22d8faec66effd
```

when each mask is encoded as a two-byte little-endian integer in increasing
order.

## Complete face-table audit

The independent facet auditor regenerates all \(2^{14}\) face flags from the
90 facets.  It verifies:

- the selected complements all lie in the frozen ridge-support core;
- all 405 ridges have degree exactly two;
- the dual facet graph is connected;
- every independent triple is covered;
- the displayed face vector is exactly the record-43
  Dehn--Sommerville schedule with no missing four-face;
- the only inclusion-minimal nonfaces are the declared six edges and four
  five-sets;
- every vertex-, edge-, and triangle-link strict \(g\)-inequality holds;
- every applicable local Macaulay inequality holds;
- every face link has spherical Euler characteristic; and
- the complex and every face link have the mod-two homology of the appropriate
  sphere.

The face table has 4,743 true entries, including the empty face.  Packing the
mask-indexed Boolean table from mask 0 through \(2^{14}-1\), least-significant
bit first in each byte, gives SHA-256

```text
32cb562b346841aa4c4a850a911a04da96d907cf6503ed6947f4960bd0212b51
```

Passing the link-homology audit does not by itself prove that the complex is a
PL sphere. The later multiwedge realization independently proves that this
face lattice is a polytopal sphere.

## Uniform rank-four Gale gate

The 90 facet complements were also submitted to the exact uniform rank-four
Gale chirotope gate.  Z3 found 1,001 basis signs.  The archived signs are then
audited without trusting or importing the solver:

- all 45,045 uniform three-term Grassmann--Pluecker relations hold;
- the positive five-element circuits are exactly the 90 facet complements;
- every element lies in a positive circuit, so total cyclicity holds; and
- every element is avoided by a positive circuit, as required for all 14
  primal labels to be extreme vertices.

The packed sign word has SHA-256

```text
01f26fdd592d1e6cb6ebd88a3075b2f7e8b2ba442e3cf345a335abdbb1e5dfa7
```

This proves feasibility as a uniform rank-four oriented matroid. It does
**not** prove that this oriented matroid is realizable over \(\mathbb R\), and
therefore does not by itself prove convex polytopality. The separate
multiwedge construction proves polytopality of the face lattice without
claiming to realize this archived chirotope.

## Reproduction and trust boundary

The construction and its frozen certificate are:

- `src/coxeter_polytopes/combinatorics/h9_14_record43.py`;
- `scripts/certify_h9_14_record43.py`;
- `data/enumerations/h9_14_record43_selection.json`; and
- `scripts/independent_h9_14_record43_audit.py`.

The critical certificate digest is

```text
af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3
```

Replay the solver-independent certificate and hostile audit with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record43.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record43_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record43.py
```

Only `--write-certificate` invokes Z3 to search for a Gale sign word. Ordinary
verification reads the archived signs and checks every chirotope relation and
positive circuit directly. Ordinary polytopality is now supplied by the
separate multiwedge theorem. The remaining target gates are Coxeter labels
and an exact Lorentzian Gram realization; realizability of this particular
archived chirotope is not claimed.
