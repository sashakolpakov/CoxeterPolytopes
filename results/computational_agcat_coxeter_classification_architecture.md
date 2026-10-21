# Computational classification architecture from the exact AGCat theorem

**Status (2026-08-12): implemented effective middle layer, with exact small-rank
censuses.**  This report turns the proof of
[`agcat_global_maximal_arithmetic_coxeter_finiteness.md`](agcat_global_maximal_arithmetic_coxeter_finiteness.md)
into finite data structures and algorithms.  The new prototype is not a
theorem-text checker.  It enumerates closed Cartan packets, reconstructs
reflection representations up to diagonal root gauge, constructs an integral
Cartan gauge, tests Lorentzian and finite-volume conditions exactly, quotients
facet permutations, and sends every survivor through the repository's
independent chamber verifier.

The global theorem is not yet a global census.  Its field theorem is partly
qualitative, and the order/lattice stage has not yet been enumerated for every
field and rank.  Section 6 states the exact additional finite lists required
for a complete classification.

## 1. How much of the proof is geometry, AG, and arithmetic?

A percentage would be misleading: the proof is a serial composition, so
removing any one of its three main blocks breaks the result.  Its roles split
as follows.

| block | indispensable output | techniques actually used |
|---|---|---|
| AG carrier | exact functor and cycle-coordinate quotient; **no finiteness by itself** | real-framed projective quadric--hyperplane arrangements; diagonal-torus invariant regular functions; real-DM orbifold curves |
| hybrid field gate | a **finite set of minimal Vinberg fields** | polytope geometry forces small packets; Vinberg integrality makes cycle coordinates arithmetic; selected/nonselected boxes, relative Fekete, and cyclotomic norms make them finite |
| Coxeter/hyperbolic geometry | finitely many ranks and a bounded full-rank frame | narrow facets, simple-polytope face averages/Dehn--Sommerville and finite Lanner/transition supports; Vinberg--Khovanskii--Prokhorov \(n\leq995\); Nikulin's connected \(62\)-frame lemma |
| algebra/order theory | finitely many native integral hosts for fixed \((k,m)\) | cycle reconstruction; Burnside; full orders, conductors, class groups, and finite quotient modules |
| reflection geometry | one chamber orbit for a maximal host | \(W(P)=R(L)\), mirror galleries, and exact chamber/Vinberg verification |

Thus the argument is genuinely AG-mediated, but no finiteness implication is
cleanly AG-only.  The field gate is a serial hybrid: the algebraic carrier
exposes integral regular coordinates, polytope geometry finds a small
full-field packet, and arithmetic Fekete/norm estimates make its points
finite.  In particular, the face-average/Hodge forcing is polytope geometry
(with an optional combinatorial-Hodge interpretation), not a claimed AG
boundedness theorem.  The computational core after the field list is produced
is mostly exact number theory, finite graph search, linear algebra, and
reflection algorithms.  It is not a KSBA/BAB surface enumeration: the direct
`AG-DEG` bound on the facet divisor remains open.

The attribution firewall is also exact.

* No conclusion of Agol--Belolipetsky--Storm--Whyte is used.  ABSW is only a
  convenient place where Vinberg's structural reflection lemma is recorded.
* Nikulin's fixed-\((k,m)\) finiteness conclusion and his published global
  ground-field bound are not used.  His geometric connected-frame lemma, with
  constant \(62\), **is** used.  The later conductor proof reconstructs the
  arithmetic core of his fixed-data method independently.
* Vinberg techniques are used substantially: the quadratic-form reduction,
  cyclic arithmeticity criterion, reflection representation, compact
  dimension theorem, and, in a computational realization, the Vinberg root
  algorithm.  Vinberg's compact nonexistence theorem \(n<30\) **is** combined
  with Khovanskii--Prokhorov to obtain the universal \(n\leq995\) ceiling.
  The stronger separate all-arithmetic \(n\leq29\) corollary/cutoff is
  deliberately not used.
* The noncompact dimension ceiling is Khovanskii--Prokhorov; the polygon branch
  imports Long--Maclachlan--Reid; relative Fekete is the terminal field tool.

This is therefore an interesting new composition, not a technique-free proof:
it avoids the three pre-existing *global finiteness conclusions*, while using
their underlying structural geometry where explicitly listed.

## 2. The finite classification pipeline

The proof becomes the following certificate pipeline.

```text
primitive arrangement coordinate / exact field certificate
                         |
                 finite field list K
                         |
      Minkowski enumeration of bounded closed cycles
                         |
 connected graph + complete cycle packet, modulo S_m
                         |
 spanning-tree reconstruction modulo diagonal root gauge
                         |
 exact signatures at every real embedding
                         |
       reflection matrices and full order Lambda
                         |
 class-group representatives + Lambda-stable submodules
                         |
       native projective lattices and full hosts R(L)
                         |
        Vinberg chamber transcript + exact geometry
                         |
 maximal chamber/group, canonicalized up to isometry
```

Here are the finite records needed at each stage.

### 2.1 Field records

A field record consists of

```text
(minimal polynomial, integral basis, multiplication table,
 real-root isolating intervals, distinguished embedding, discriminant)
```

The arrangement proof supplies primitive algebraic integers such as the
three-wall determinant coordinate \(\beta\), plus intervals for every real
conjugate.  Once its degree bound is made numerical, coefficient bounds or a
Minkowski/Fincke--Pohst search enumerate all possible minimal polynomials.
Each polynomial is checked for total reality, the displayed embedding
conditions, and primitiveness \(k=\mathbf Q(\beta)\).  Transition charts,
Lanner charts, and polygon charts each emit the same field-record interface.

### 2.2 Cycle-packet records

For a field \(k\), rank \(m\), and cycle length \(r\leq m\), define the
explicit finite Minkowski box

\[
 \mathcal S_{k,r}=
 \left\{x\in\mathcal O_k:
   |\sigma_0(x)|<62^r,
   \ |\sigma(x)|\leq2^r\ \text{for }\sigma\ne\sigma_0
 \right\}.                                                \tag{2.1}
\]

Enumeration in an integral basis of \(\mathcal O_k\) is a bounded lattice
point problem.  A frame record is

```text
(connected graph G on m vertices,
 b_(i,j,i) in S_(k,2) for every edge,
 b_C in S_(k,|C|) for every simple cycle C).
```

It must satisfy

\[
 b_C=b_{C^{-1}},\qquad
 b_C^2=\prod_{\{i,j\}\in C}b_{iji},                       \tag{2.2}
\]

the acute sign convention at the selected place, and all exact embedding
inequalities.  In practice one assigns a cycle basis, propagates the remaining
cycles, and uses canonical graph augmentation rather than enumerating all
labelled graphs blindly.

### 2.3 Cartan reconstruction

Choose a root and spanning tree.  Put \(a_{ij}=-1\) on every outward tree
edge and \(a_{ji}=-b_{iji}\) on its reverse.  For a non-tree directed edge
\(u\to v\), let \(T_{v,u}\) be the directed tree path from \(v\) to \(u\).
The fundamental cycle gives

\[
 a_{uv}=
 \frac{b_{(u,v,T_{v,u})}}
      {\prod_{e\in T_{v,u}}a_e},qquad
 a_{vu}=\frac{b_{uvu}}{a_{uv}}.                           \tag{2.3}
\]

All remaining cycle equations are then exact consistency checks.  The proof
of the native-host theorem shows that (2.3) produces the unique Cartan matrix
up to

\[
                   a_{ij}\mapsto\lambda_j a_{ij}/\lambda_i.          \tag{2.4}
\]

Solve \(d_i a_{ij}=d_j a_{ji}\) and form
\(Q=\operatorname{diag}(d_i)A\). At the distinguished embedding \(Q\) must
have inertia \((m-1,1,0)\); every
other embedding must be positive definite.  Exact real-algebraic sign oracles
make this a certificate, not a floating-point filter.

Over \(\mathbf Q\), cyclic integrality constructively gives an integral
Cartan gauge.  At every prime \(p\), put

\[
 x_j-x_i\geq-v_p(a_{ij}).                                 \tag{2.5}
\]

A positive constraint cycle is precisely a closed product with negative
\(p\)-valuation.  When all \(b_C\in\mathbf Z\), longest-path relaxation finds
the potentials \(x_i=v_p(\lambda_i)\), and (2.4) becomes integral.  Over a
general field, the replacement is an ideal-valued gauge with local valuation
potentials and an explicit ideal-class record.

### 2.4 Orders and projective lattices

From an integral/projective realization, build

\[
 s_i(e_j)=e_j-a_{ij}e_i,qquad
 \Lambda=\mathcal O_k[s_1,\ldots,s_m]\subseteq
 B=\operatorname {End}_{\mathcal O_k}(L_0).               \tag{2.6}
\]

Compute an \(\mathcal O_k\)-basis of \(\Lambda\) by saturation, then its
conductor

\[
 \mathfrak c=\{a\in\mathcal O_k:aB\subseteq\Lambda\}.     \tag{2.7}
\]

For each ideal-class representative \(\mathfrak a\), set
\(M=\mathfrak aL_0\).  Enumerate the finite submodule lattice of
\(M/\mathfrak cM\), retain inverse images \(L\) that are projective, full,
and \(\Lambda\)-stable, and quotient \(k^\times\)-homothety and normalizer
actions.  This is the effective form of

\[
                         \mathfrak cM\subseteq L\subseteq M.          \tag{2.8}
\]

### 2.5 Full reflection hosts and chambers

For each retained \((q,L)\), enumerate integral reflections by the Vinberg
height algorithm.  A transcript records, height by height, every root and
every rejection; the existing verifier architecture can replay wall
completeness independently.  Exact half-space enumeration then certifies
nonemptiness, Coxeter relations, finite/ideal vertices, absence of hyperideal
rays, and finite volume.  The resulting group is \(R(L)\).  A candidate input
is reflection-maximal exactly when its side group equals this full host.

Finally canonicalize

```text
(field, distinguished place, quadratic similarity class,
 projective lattice class, ordered/unordered chamber Gram data)
```

under field automorphisms and lattice isometries.  The mirror-gallery lemma
means one unmarked chamber orbit per full host; no separate unbounded chamber
enumeration is needed.

## 3. Implemented prototype

The implementation is in
[`cartan_packets.py`](../src/coxeter_polytopes/classification/cartan_packets.py),
with CLI
[`classify_bounded_cartan_packets.py`](../scripts/classify_bounded_cartan_packets.py)
and exact regression tests in
[`test_computational_cartan_classifier.py`](../tests/test_computational_cartan_classifier.py).

Its present scope is deliberately closed and unconditional:

* field \(k=\mathbf Q\);
* a user-supplied finite set of positive rational length-two products;
* ranks at most eight, with a search-size guard;
* connected acute symmetrizable frames;
* exact quotient by diagonal root gauge and every facet permutation; and
* finite-volume Coxeter **simplices**, where the frame itself is the complete
  wall set.

For each pair of vertices it enumerates zero (orthogonality) or one supplied
edge product.  A longer acute cycle is retained over \(\mathbf Q\) exactly
when the product of its edge products is a rational square and the resulting
signed square root is integral.  It then performs (2.3)--(2.5), exact inertia
tests, and the spherical/affine link criterion.  Every survivor is serialized
as a genuine chamber certificate and independently replayed by
`verify_chamber_geometry`; that verifier recomputes the wall Gram matrix,
Coxeter relations, interior point, extreme rays, finite vertices, and cusps.

Run, from the repository root,

```bash
python3 scripts/classify_bounded_cartan_packets.py --rank 3 --summary-only
python3 scripts/classify_bounded_cartan_packets.py --rank 4 --summary-only
python3 -m pytest -q tests/test_computational_cartan_classifier.py
```

## 4. Certified pilot censuses

Take the rational Coxeter edge box

\[
 b_{iji}\in\{0,1,2,3,4\},                              \tag{4.1}
\]

corresponding respectively to labels
\(2,3,4,6,\infty\), with \(4\) interpreted as parallel.  The frozen outputs
are

* [`rational_cartan_cycle_simplex_rank3.json`](../data/enumerations/rational_cartan_cycle_simplex_rank3.json),
* [`rational_cartan_cycle_simplex_rank4.json`](../data/enumerations/rational_cartan_cycle_simplex_rank4.json).

| stage | rank 3 / \(\mathbf H^2\) | rank 4 / \(\mathbf H^3\) |
|---|---:|---:|
| labelled pair assignments | 125 | 15,625 |
| connected assignments | 112 | 15,104 |
| labelled integral \(\mathbf Q\)-cycle packets | 68 | 2,976 |
| unlabelled packets modulo facet permutation | 18 | 207 |
| Lorentzian packets | 13 | 200 |
| finite-volume simplex chambers | **13** | **18** |
| compact | 4 | 1 |
| noncompact | 9 | 17 |
| independent exact geometry replays | 13/13 | 18/18 |

The frozen executable evidence is:

| artifact | SHA-256 |
|---|---|
| classifier module | `f712f87aa334d06554a56804b8ad2e0d475727c9bb0eb1a199cb9a3eb8691ee0` |
| CLI | `c801b40e4699d9b8667b16e55908f8ee603211e16873b1abe5636d2cad32c914` |
| exact tests | `4f0b83fe3c6ad8fccdb61564a433aa9a2b983c6fce06908668bdf57111ef822b` |
| rank-three JSON | `a1cb9923af6636f0e6396e55d07b065e326cac67d19a3473f82683112429cca9` |
| rank-four JSON | `80514dc961565c004d0f3fd34bc0e9defe921322becbe85a5cbf0e5898913528` |

The classifier suite has **10/10 passing tests**; combined with the owner and
independent native-host suites, **21/21 tests pass**.  The checks include an
exact obstruction from a nonintegral closed cycle, constructive integral
gauge recovery, gauge/permutation invariance, Euclidean-versus-hyperbolic
separation, exact rank-three and rank-four census counts, and independent
finite/ideal vertex reconstruction.

The triangle result has a second transparent check: its thirteen unordered
label triples are

\[
\begin{gathered}
(2,3,\infty),(2,4,6),(2,4,\infty),(2,6,6),(2,6,\infty),
(2,\infty,\infty),\\
(3,3,\infty),(3,4,4),(3,6,6),(3,\infty,\infty),
(4,4,\infty),(6,6,\infty),(\infty,\infty,\infty).
\end{gathered}                                             \tag{4.2}
\]

They are exactly the triples in the label alphabet with hyperbolic angle sum,
connected nonorthogonality graph, and rational integral cycle field.  This is
an exhaustive census in the displayed box, not a random search.

The eighteen tetrahedra are candidate chambers, not claims of new examples:
the pilot intentionally does not yet compare a literature database or test
whether each side group is the full reflection group of its native lattice.
Those are the host/maximality stages in Sections 2.4--2.5.

## 5. Completeness theorem for the prototype

> **Proposition 5.1 (implemented finite-box completeness).**  Fix a rank
> \(3\leq m\leq8\) and a finite set \(S\subset\mathbf Q_{>0}\) of integral
> length-two Cartan products recognized as Coxeter relations.  The prototype
> returns exactly the isomorphism classes, modulo facet permutation and
> diagonal root gauge, of finite-volume rational Coxeter simplices whose
> nonzero length-two products lie in \(S\), whose Cartan graph is connected,
> and whose closed Cartan products are integral.

Indeed, the program exhausts \((\{0\}\cup S)^{\binom m2}\).  Equation (2.2)
forces the unique acute longer-cycle value whenever it lies in \(\mathbf Q\).
Equation (2.3) reconstructs the unique gauge class.  Exact symmetrization and
inertia keep precisely Lorentzian forms.  A Lorentzian Coxeter simplex has
finite volume precisely when every deleted principal Gram matrix is positive
definite or positive semidefinite of corank one.  Finally all \(m!\) facet
permutations are canonically quotiented, and the separate chamber engine
replays every accepted geometry.

The proposition does **not** say that the finite label box contains every
arithmetic simplex, nor that an accepted simplex group is reflection-maximal.

## 6. What is still needed for a complete global classification?

The global finiteness proof becomes a terminating global enumeration once the
following inputs are emitted as explicit files rather than existential finite
sets.

1. **Field manifest.**  Enumerate and certify the full finite set
   \(\mathscr K_{\mathrm{Cox}}\).  The transition/fork label bounds are
   explicit, but the Long--Maclachlan--Reid polygon lists and all fixed-label
   Fekete integral points have not yet been expanded into a common manifest.
2. **Minkowski packet manifests.**  For every
   \((k,m)\in\mathscr K_{\mathrm{Cox}}\times\{3,\ldots,996\}\), enumerate the
   boxes (2.1), preferably lazily through cycle-basis constraints rather than
   materializing their Cartesian products.
3. **General-field gauge backend.**  Replace rational prime potentials by
   prime-ideal potentials and ideal-class data; add exact signs at every real
   embedding.
4. **Order/lattice manifests.**  Compute one \(L_0\), the full order
   \(\Lambda\), conductor \(\mathfrak c\), class-group representatives, and
   every retained submodule of \(M/\mathfrak cM\).
5. **Full-host transcripts.**  Run and replay Vinberg enumeration for every
   retained lattice, then canonicalize the full chambers.  Termination is
   guaranteed for hosts actually arising from the theorem, but a practical
   queue needs height/memory partitioning and independently checkable partial
   transcripts.

Non-simplex chambers do not require guessing an a priori facet bound.  Their
connected (m)-wall frame passes through the same packet compiler; after the
order/lattice stage, Vinberg enumeration reconstructs the *entire* wall set of
the full host.  The simplex-only restriction belongs to the present chamber
validator loop, not to the theoretical classification architecture.

Conditional on those five explicit manifests, completeness is not heuristic:
every maximal chamber chooses a bounded connected frame, its cycle packet
appears in the Minkowski enumeration, reconstruction returns its
representation, the conductor sandwich returns its native lattice, and
maximality identifies its side group with the enumerated full host.

The largest present obstacle to a usable worldwide census is therefore not
the Cartan reconstruction; that part is now executable.  It is expanding the
qualitative field union and the conductor intervals into manageable explicit
manifests.

## 7. Finding examples and improving dimension information

The architecture supports two different computational programs.

**New examples.**  Start with low-degree fields from the field manifest and
sparse connected frame graphs.  Use local spherical/affine/Lanner constraints,
cycle-basis propagation, and nonselected definiteness before constructing an
order.  Surviving projective lattices enter a Vinberg search.  A new example
requires a complete chamber transcript plus failure to match the canonical
host/chamber database; a raw Lorentzian packet is only a lead.

**Dimension bounds.**  To raise the known lower bound for the largest
dimension, search high-rank admissible packets and complete one to a certified
finite-volume host.  To lower the upper ceiling \(995\), prove that the packet
CSP is unsatisfiable above a rank: combine local Lanner obstructions,
nonselected definiteness, arithmetic cycle boxes, and cusp-link constraints,
then emit a DRAT/LRAT or independently replayable case-partition certificate.
The current prototype improves no global dimension bound; it supplies the
exact packet compiler needed for those searches.

The most productive next implementation step is the quadratic-field backend
for \(\mathbf Q(\sqrt5)\), because the repository already has exact arithmetic,
embedding signs, chamber verification, and Vinberg searches for that field.
It exercises every new issue--nontrivial embeddings, algebraic Coxeter labels,
and ideal-valued gauges--without immediately confronting the full field list.
