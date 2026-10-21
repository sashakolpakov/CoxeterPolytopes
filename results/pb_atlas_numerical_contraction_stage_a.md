# PB-Atlas numerical Stage A: the arrays-only all-contraction compiler

## Result

The first genuinely arbitrary-input piece of `PB-Atlas` is now implemented.
The recognizer accepts only

\[
  (G_{\operatorname{Pic}},\; (w_i)_{i=1}^m,\;
    (b_j,q_j)_j,\; K),
\]

namely the integral Picard Gram, primitive effectively oriented wall vectors,
marked boundary indices and anti-bicanonical weights, and the canonical
class.  The relation checked at input is

\[
  \sum_j q_jw_{b_j}=-2K.
\]

No wall label, surface name, Kondō type, stored contraction frame, stored
endpoint Gram, incidence coordinate, or output witness is admitted by the
parser.  The executable rejects extra fields.

For every such finite wall array, the compiler enumerates **every sequential
boundary-preserving contraction frame supported by the supplied walls**.  If
the supplied array is the complete negative-curve/facet array of a full-nef
carrier, this is every geometric `(-1)`-contraction frame.  At ranks one and
two it emits, respectively, all `P2` and all Hirzebruch `F_e` endpoints.
`F_1` is retained and marked nonminimal: it is a valid ruled contraction
target even though its negative section may be contracted once more to a
plane endpoint.

This closes **numerical Stage A only**.  It does not prove that any universal
incidence scheme is nonempty and it is not an `L4` result.

## Exact contraction rule

Let `S=(w_{i_1},...,w_{i_s})` be the strict classes already selected and let
`A=S^tGS`.  The current image of a candidate wall `c`, pulled back to the
original Picard lattice, is

\[
 e_c=c-SA^{-1}S^tGc.
\]

The next contraction is allowed exactly when

\[
 e_c^2=c^2-(c^tGS)A^{-1}(S^tGc)=-1,
 \qquad K\mathbin{\cdot}e_c=-1.                 \tag{1}
\]

The inverse update is the integral rank-one Schur update already used by the
older `H10/30` and `H11/26` plane certifiers.  Thus no floating point and no
geometric label enters a transition.  Every accepted pivot has determinant
`-1`, so all reachable exceptional Grams are unimodular and their inverses
remain integral.

The implementation extracts and unifies three existing project-native
pieces rather than introducing a parallel convention:

1. the Schur-inverse transition from
   `certify_coble_type_ii_b3_plane_inverse.py`;
2. the exact pullback update `D -> D+(D.E)E` from the rational-surface MMP
   certifiers;
3. the rational congruence inertia routine used in the exact K3 lattice
   verifiers.

## Why the symmetry quotient is exhaustive

Column Hermite normal form records the index of the wall span in `Pic`.  When
that index is one, every permutation preserving the marked weighted wall Gram
lifts automatically to an integral Picard isometry: it defines a map on the
integral wall generators, Gram equality preserves every relation, the inverse
permutation gives the inverse lattice map, and the canonical pairings force
`K` to be fixed.  A nonsaturated full-rank wall array is **not rejected**.
Instead, each of its marked-Gram permutations is solved for its unique
rational Picard lift, and only integral determinant-one isometries preserving
`K` and every wall are retained.  The identity always survives, so arbitrary
arrays remain enumerable even when no symmetry compression is available.

The full **integrally liftable** marked wall-Gram automorphism group is
enumerated.  At a reachable mask `M`, the remaining walls are partitioned
into orbits under the full stabilizer of `M`.  Formula (1) is invariant under
that stabilizer, so testing one representative of each candidate orbit tests
every labelled extension.
The new mask is then replaced by its full-group canonical representative.
Induction on the contraction depth proves that the orbit traversal omits no
labelled frame.  Every terminal orbit is expanded back to its complete list
of labelled masks in the archive.

If an automorphism list is supplied explicitly rather than discovered by the
core, the validator now also proves that it is a permutation **group**: the
identity must occur, duplicate rows are forbidden, and the generated subgroup
must equal the supplied set.  Checking individual Gram-preserving
permutations is not enough, because canonical masks and stabilizers would not
define group orbits for a list missing a product or inverse.

This compression is substantial.  For example, the `H17/48` calibration has
`311238` labelled frames at depth eight but only `329` marked-Gram orbits.

## Endpoint reconstruction

For a terminal mask with strict-class matrix `S`, the integral orthogonal
projector is

\[
 P_S=1-S(S^tGS)^{-1}S^tG.
\]

Column HNF of `P_S` gives the entire orthogonal-complement lattice, not merely
its rational span.

* In complement rank one, the primitive nef generator `H` is required to
  satisfy `H^2=1` and `K.H=-3`; this certifies `P2`.
* In complement rank two, the two nef rays are recomputed from all supplied
  wall inequalities.  An ordering `(f,h)` is accepted precisely when

  \[
     f^2=0,\quad f.h=1,\quad K.f=-2,
     \quad h^2=e\ge0,\quad K.h=-(e+2).
  \]

  Then `s=h-ef` has `s^2=-e`, `s.f=1`, and
  `K=-2s-(e+2)f`, so the endpoint is `F_e`.  The two rulings of `F_0` are
  canonically treated modulo their swap.

For every terminal orbit representative, the archive contains:

* every labelled member mask in the orbit;
* a valid sequential contraction order;
* the strict contracted Gram and every total exceptional transform;
* the integral `P2` or `(s,f)` model basis;
* the full unimodular endpoint-plus-exceptional basis and its Gram;
* canonical coordinates in that basis;
* the coordinates of **every input wall**, both before and after forgetting
  the exceptional coordinates;
* boundary indices, anti-bicanonical weights, and endpoint images;
* an exact per-step contraction trace.

This is the numerical packet needed by the next universal
cluster/incidence-equation generator.  It contains no chosen coordinate
point or Hilbert witness.

## Four full calibrations

All four requested frozen physical branches were exhausted, including the
two high-rank cases; the qualifier “if feasible” is no longer needed.

| calibration metadata (not recognizer input) | marked-Gram group | labelled ruled endpoints | labelled `P2` endpoints |
|---|---:|---:|---:|
| `H10/30` | `24` | `F0: 837`, `F1: 2064`, `F2: 390` | `573` |
| `H11/26` | `8` | `F0: 86`, `F1: 212`, `F2: 54`, `F3: 24`, `F4: 12` | `52` |
| `H14/60` | `120` | `F0: 105`, `F1: 310`, `F2: 120`, `F3: 20` | `5` |
| `H17/48` | `1152` | `F0: 1` | `0` |

The complete labelled-frame trails, starting at depth zero, are:

```text
H10/30: 1,17,87,329,948,2094,3472,4417,4455,3291,573
H11/26: 1,14,39,72,125,218,335,424,486,512,388,52
H14/60: 1,45,355,1620,5200,13032,26880,45100,58640,
         56270,38136,17210,4615,555,5
H17/48: 1,16,216,2288,13628,50256,126664,230896,311238,
         309232,221416,109584,34940,6320,504,16,1,0
```

The first two terminal plane counts independently reproduce the older
`573` and `52` sequential-endpoint enumerations.  The latter two are new full
all-contraction calibrations.  Aggregated over the four controls there are
`4235` ruled frames and `630` plane frames.

## Mutation and hostile-input audits

Each of the four branches is replayed after:

1. a nontrivial determinant-one Picard basis shear; and
2. reversal of the entire wall array, with boundary indices and the full
   automorphism group transported.

The first mutation preserves every terminal endpoint-coordinate packet
entrywise, modulo the unavoidable swap of the two `F0` rulings.  The second
transports every labelled terminal mask and preserves the complete abstract
endpoint-orbit multiset.  A hostile full-rank index-two wall span is accepted,
while a synthetic marked-Gram swap with a half-integral Picard lift is
discarded and the identity is retained.  Thus nonsaturation reduces symmetry
compression but never corrupts or narrows the contraction census.

## Exact boundary of the result

Stage A does **not** do any of the following:

* generate universal point/proximity/incidence equations;
* generate collision, Jacobian, irreducibility, or Hilbert avoidance factors;
* prove a primitive-binomial component survives those factors;
* prove the resulting incidence stack is nonempty;
* construct an arbitrary strong Coble carrier;
* exclude hidden negative curves or return literally the same `P`;
* establish a non-singleton `L4` class or an ABSW-type dimension theorem.

The next exact task is therefore **PB-Atlas Stage B**: consume the frozen
endpoint packets, construct the universal cluster and curve-evaluation
matrices without row-specific routing, and generate every bad factor from
those matrices.  Only after a regular open survives is there a geometric
carrier to which the K3-cover rank bound can be applied.

## Reproducibility

* arrays-only core: `scripts/pb_atlas_contraction_core.py`;
* fixture extractor: `scripts/freeze_pb_atlas_stage_a_fixtures.py`;
* calibration/mutation certifier:
  `scripts/certify_pb_atlas_numerical_stage_a.py`;
* input schema: `schemas/pb-atlas-physical-branch.schema.json`;
* output schema: `schemas/pb-atlas-numerical-stage-a.schema.json`;
* frozen arrays: `data/fixtures/pb_atlas_numerical_stage_a_inputs.json`;
* frozen certificate:
  `data/enumerations/pb_atlas_numerical_contraction_stage_a.json`;
* tests: `tests/test_pb_atlas_numerical_contraction_stage_a.py`.
