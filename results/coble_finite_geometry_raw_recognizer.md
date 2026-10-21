# Arbitrary-input recognition of five Coble finite-geometry cores

**Status (2026-08-09): the missing generic recognition gate is closed.**
There is now an executable, surface-free class `C_FG` whose membership is
decided from a normalized pair

```text
{lattice_gram, primitive_wall_vectors, orientation.positive_half_space}
```

alone.  The recognizer accepts arbitrary integral unimodular coordinate
changes and arbitrary wall permutations.  It does not read an archive name,
a source digest, a surface or Kondō label, an expected type, or a stored
output Gram.

This is a real improvement over a five-row replay, but it is **not L4 under
the strict five-gate contract**.  Its honest rating is `L3 + closed generic-
recognition gate`: a theorem-shaped scoped inverse class whose payoff gate is
still open.  Even on `C_FG ∩ D_CobVin`, where the independent finite-volume,
Vinberg, cusp, and exact-return hypotheses are certified, the present grammar
contains five already known abstract finite geometries.  Kondō pushback therefore
produces no new abstract polytope and no sharper global dimension bound.

The executable is
[`certify_coble_finite_geometry_raw_recognizer.py`](../scripts/certify_coble_finite_geometry_raw_recognizer.py),
with neutral generated
[`fixtures`](../data/fixtures/coble_finite_geometry_raw/manifest.json),
the frozen
[`certificate`](../data/enumerations/coble_finite_geometry_raw_recognizer.json),
and
[`tests`](../tests/test_coble_finite_geometry_raw_recognizer.py).
The separate
[`fixture generator`](../scripts/generate_coble_finite_geometry_raw_fixtures.py)
knows the five calibration archives only to extract their physical Gram and
walls; no archive information enters a fixture or the recognition path.

## 1. Input contract and normalization

The input is an abstract integral lattice presentation, not a named surface.
The current normal form is the odd-unimodular physical Picard presentation:

1. `lattice_gram` is integral, symmetric, unimodular, and has exact inertia
   `(1,r-1,0)`.
2. `primitive_wall_vectors` is a list of distinct primitive integral row
   vectors spanning the lattice rationally.
3. `orientation.positive_half_space` is literally
   `pairing(x, wall) >= 0`.  With this orientation, pairwise intersections of
   distinct effective wall classes must be nonnegative.

H14 and H17 also have arithmetic reflection-lattice presentations in their
upstream archives.  Those must first be de-scaled to the odd-unimodular
presentation.  The present theorem begins after that arithmetic
normalization; it does not silently pair raw arithmetic walls with a physical
Gram.

The recognizer then reconstructs, rather than receives, the boundary frame.
Square-`-4` walls give the multi-boundary frame.  In the one-boundary case it
finds the unique characteristic primitive square-`-1` wall
`b=B/2=-K`.  It verifies

```text
B_i.B_j = 0,
sum B_i = -2K,
K^2 = -b,
K characteristic,
```

and partitions every remaining wall into

```text
root C:       C^2=-2, K.C=0, C.B_i=0;
crossing E:   E^2=-1, K.E=-1, boundary profile (1,1);
cusp lift E:  E^2=-1, K.E=-1, boundary profile (2).
```

All subsequent decisions use exact pairings and exact integral contraction
bases.  There is no Gram-template comparison.

## 2. The theorem-ready class `C_FG`

Define `C_FG` to be the normalized oriented wall data passing the preceding
physical Coble alphabet and exactly one of the following five structural
recognizers.

| dimension / walls | reconstructed core | intrinsic certificate |
|---|---|---|
| `H14/60` | polar/chord | five different unimodular `P2` contractions; each has four boundary lines and one conic, four contracted tangency roots, and all six pair/complement residual lines |
| `H17/48` | grid/projectivity | eight boundaries form `K4,4`; the sixteen disjoint crossings are the complete `4×4` grid; the 24 roots are exactly all `S4` permutation matchings of bidegree `(1,1)` |
| `H11/206` | Möbius `S(3,4,10)` | ten disjoint crossing centers contract to a determinant-one `P1×P1`; the 30 roots are four-center `(1,1)` blocks and every one of the 120 triples occurs exactly once |
| `H10/78` | Desargues triangle/cut | the unique ten-exceptional plane frame gives ten degree-one triangle triples and ten degree-two six-edge cuts; the colored incidence graph has exactly 120 identifications with the edge design of `K5` |
| `H10/72` | tetrahedral `4+6` | exactly five unimodular `P2` frames share six ordinary centers and carry the ten edges of `K5`; every frame has 12 lines, four conics, four quartics, and exactly 24 colored design markings |

For `H10/72`, the recognizer first computes the maximal-clique histogram
`7:80, 9:25, 10:5` in the disjointness graph of the 51 nonboundary
square-minus-one walls.  Across the five ten-frames, six walls occur five
times, ten occur twice, and 35 occur zero times; every pair of frames meets
in seven walls, so the ten twice-occurring walls are intrinsically the ten
edges of `K5`.  In each frame it reconstructs
`H=(b+sum(E_i))/3`, checks the full unimodular `I_(1,10)` marking and
integrality of all 72 transported walls, then reads the root coefficients.
An exhaustive 720-permutation check—not a graph hash—finds exactly 24
identifications with the abstract `4+6` colored design.  The returned record
contains only invariant counts and profiles, so it is literally unchanged by
the tested coordinate shear and wall reversal.

The dimensions and wall counts are routing invariants only.  Acceptance also
requires the full endpoint basis, divisor-class identities, and finite design.
The tests apply a nontrivial determinant-one integral shear to every lattice
basis and reverse every wall list; all five outputs remain unchanged at the
structural level.

## 3. Four logically separate layers

### (i) Recognized grammar

Membership in `C_FG` is now a genuine P-first computation.  A caller can
supply a new normalized Gram and wall list through `--input`; the program
does not need the fixture manifest or any known-row metadata.

### (ii) Uniform explicit realization

For every accepted datum, the recovered unimodular `P2` or `P1×P1` marking
transports the corresponding standard finite-field construction back to the
input lattice.  The polar/chord, grid/projectivity, and Möbius cores are
realized in characteristic three; the Desargues triangle/cut and tetrahedral
cores are realized in characteristic five.  Thus every member of `C_FG` has a marked strong-
Coble **root-seed** realization.

The relevant construction audits are the
[`four-point/Möbius hostile audit`](coble_four_point_constructible_incidence_hostile_audit.md)
and the
[`Desargues inverse theorem`](coble_type_vi_char5_desargues_inverse.md), together
with the
[`tetrahedral inverse theorem`](coble_type_vii_char5_tetrahedral_inverse.md).
This layer realizes the recovered boundary, crossing, and root core.  It does
not, by itself, say that all supplied cusp walls are the complete negative-
curve stock.

### (iii) Exact same-`P` return

For an input which independently lies in `D_CobVin`, the
[`uniform inverse-transfer theorem`](coble_vinberg_incidence_inverse_transfer.md)
adds exactly the missing statements:

```text
finite volume;
Vinberg finite index and root exhaustion;
complete primitive-cusp/single-boundary lift regeneration;
hidden-wall exclusion;
literal equality of the reconstructed full nef chamber with the input P.
```

Consequently every member of `C_FG ∩ D_CobVin` has a surface carrier and
returns literally the same oriented polytope.  The recognizer deliberately
does not infer these hypotheses merely from finite-sidedness or from the
finite-geometry core.

### (iv) Kondō pushback

After the carrier and exact return exist, Kondō's classification applies as
an independent algebro-geometric theorem.  On this thin class it recovers the
four distinct dimensions already visible in the grammar, namely
`{10,11,14,17}`.  Therefore the formal implication is valid, but its present
pushback content is thin: it does not yet lower a worldwide dimension bound,
exclude a broader P-only family, or create a new abstract Coxeter polytope.

## 4. Hostile mutations

Each of the five positive fixtures has a paired mutation which reflects one
root through an incident root while leaving every other wall fixed.  If
`a^2=b^2=-2` and `a.b=1`, replace only `b` by `b+a`.  Then

```text
(b+a)^2 = -2,
rank unchanged,
facet count unchanged,
complete wall-norm histogram unchanged,
```

but `a.(b+a)=-1`.  Thus these are accurately classified as **acute/orientation
mutations**: they are rejected before the deeper contraction-design gate,
even though all crude rank/facet/norm invariants agree with the positive
fixture.  They are not advertised as deep incidence counterexamples.  A
separate nonprimitive-wall mutation is also rejected.

No orientation-preserving deeper mutation is frozen here.  In the natural
single-root block model this is already rigid for the Möbius core: every
triple of ten centers lies in exactly one `S(3,4,10)` block, so any new
four-set shares three centers with an unchanged old block and would create a
negative root pairing.  A multiwall orientation-preserving mutation would
need a separate exhaustive construction and audit; inventing one is not
needed for the present positive recognition theorem.

## 5. Five-gate rating

| gate | verdict | reason |
|---|---|---|
| independent P-only class `C` | `PASS_SCOPED` | arbitrary-input membership in `C_FG` is executable and invariant under basis/wall reorderings |
| explicit carrier/root seed `Car` | `PASS_SCOPED_ROOT_SEED` | all five recovered designs have uniform explicit finite-field realizations |
| exact same-`P` return | `CONDITIONAL_PASS` | unconditional on `C_FG ∩ D_CobVin`, but not certified from the two arrays by this executable alone |
| independent AG classification | `CONDITIONAL_PASS` | applies after carrier and return, not inside recognition |
| new/nonthin pushback | `FAIL_NONTHIN_NOVELTY` | current grammar already fixes five known designs and four dimensions |

Hence the honest label is:

```text
L3 plus closed generic-recognition gate;
theorem-shaped scoped inverse class, payoff gate still open;
not L4, flagship, or global.
```

## 6. Replay

Recognize one genuinely supplied input:

```bash
python3 scripts/certify_coble_finite_geometry_raw_recognizer.py \
  --input data/fixtures/coble_finite_geometry_raw/positive_05.json
```

Replay the full positive/negative suite against the frozen certificate:

```bash
python3 scripts/certify_coble_finite_geometry_raw_recognizer.py \
  --input data/fixtures/coble_finite_geometry_raw/manifest.json \
  --check-output data/enumerations/coble_finite_geometry_raw_recognizer.json
```

The focused test suite has ten tests covering fixture regeneration, all five
positive recognitions, the five-frame tetrahedral reconstruction,
coordinate/wall-order invariance, decoy metadata, norm-preserving hostile
mutations, schema gates, genuine CLI input, frozen archive replay, and the
five-gate scope.
