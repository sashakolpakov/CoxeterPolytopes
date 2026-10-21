# All mixed indicators preserve the forced-base obstruction

The forced-base dual also satisfies **every mixed three-port indicator**.
No vertex mass needs to increase. The exact bound remains

\[
 \boxed{\frac{31154121494284588800}{2547155390157199561}
 \approx12.2309465746>12.}
\]

The [source](../scripts/fresh_sqrt5_mixed_forced_base_transfers_20261017.py)
and [feature export with exact certificate](fresh_sqrt5_mixed_forced_base_transfers_20261017.json)
close the previous union-of-constraints gap for this partial model. They
combine the [forced-base analytic repair](fresh_sqrt5_forced_base_dual_repair_20261017.md)
with the [mixed-port incidence theorem](fresh_mixed_face_port_incidence_20261017.md).
No optimizer was needed or run.

## What is reconstructed and checked

The earlier [mixed/shared certificate](fresh_sqrt5_mixed_shared_transfers_20261017.md)
already reconstructed every provenance of the 2,807 old positive rows.
The forced-base repair retains those masses exactly. This source reconstructs
the old deterministic orbit numbering from the complete vertex forests and
their marked triple codes, then remaps every frozen coefficient. It verifies
the number of provenances against the unchanged merged quad rows. All 1,858
positive quad provenances are retained; their complete mixed vectors agree
in the frozen certificate.

For **all 4,195 new base stars**, including zero-mass ones, the source builds
the literal graph from the saved connected core and its actual orthogonal
spectators. At every finite corner it reconstructs the whole rank-13 forest,
the omitted root and the complete simplex port block. It independently
recounts both the pair and one-port coefficients and checks them against the
frozen [base-domain patch](fresh_sqrt5_forced_base_simplex_domain_20261017.md).
There are 9,496 reconstructed corner occurrences.

At a corner of a d-simplex the occupied pair footprint is the complete graph
on its d ports. A three-port subset contains at least two occupied pairs
exactly when all three ports lie in that block. Thus each corner has
`C(d,3)` mixed indicators, and a full simplex row has
`(d+1) C(d,3)` occurrences.

In particular, **maximal edges and triangles have zero mixed indicators**.
An edge still retains its two endpoint one-port charges and its signed
base-transfer constraint. The new tetrahedron and four-simplex rows receive
their complete mixed footprints. No mixed payment is fabricated for an
edge or for a single triangle corner.

| Exact verification | Count |
|---|---:|
| Positive-cost pair/clique/port columns checked | 211,339 |
| Signed shared-transfer columns checked | 5,693 |
| Positive face masses | 5,344 |
| Positive vertex forests | 752 |
| Triple-orbit columns on positive forests | 56,045 |
| Extra zero-load forests in full new-star export | 27 |
| Total exported triple-orbit columns | 57,956 |
| Old positive-row plus new full-star feature records | 7,002 |
| Mixed columns exceeding the existing vertex mass | **0** |

The export explicitly reconstructs all 286 triples at each forest, groups
them by the complete marked forest code, and assigns each orbit its literal
orbit size as the vertex cost. The 27 extra forests occur only in zero-mass
new stars; their loads are zero. Forests outside this export have no positive
face load, so their mixed-column inequalities also have zero left side.

## Combined exact dual

Keep every face mass from the forced-base repair. Recompute the load of
every ordinary positive-cost column and every mixed indicator column. For
each vertex type v the required mass is the maximum of its old mass and
every mixed load divided by that orbit's size. These maxima equal the old
masses at all 752 positive vertex types.

The source then checks all positive-cost loads against these vertex masses
and recomputes all 5,693 signed transfer loads, which remain nonpositive.
The numerator and normalization therefore remain

\[
 \sum_f b_f\alpha_f=\frac{1620989504989}{123750000000},\qquad
 \sum_v\beta_v=\frac{2547155390157199561}{2378376000000000000}.
\]

The JSON saves the face masses, augmented vertex masses, every signed load,
all nonzero mixed loads, the old-to-new orbit map and all 7,002 feature
records. `corner_features` on the new stars retains the literal omitted root
and vertex type for each contribution.

## Provenance and scope

The inherited quad coefficients come from the complete frozen
all-provenance calculation. Signed base-transfer coefficients still use
conservative coordinatewise minima over all merged provenances. These
minimum coefficients are valid for every represented actual star, but this
certificate does **not** assert a simultaneous split of a merged row mass
into actual full-star multiplicities satisfying stronger provenance-specific
incidence constraints.

The new result is an exact obstruction to the union of the current pair,
22 clique, one-port, mixed-three-port and represented/new shared-transfer
weighting families on the saved partial domain. It is not a compact
realization, a proof of an actual global incidence distribution, or a
dimension theorem. Other missing base types, fixed arithmetic-form classes,
uncovered simplexes, all-finite quad aggregates and other face domains remain
outside.

All new-star features are exported, but only **positive old-row** mixed
features are inherited here. That suffices to verify this dual: omitted old
rows have zero multiplier. A future primal optimization using these mixed
columns would first need their features on the remaining old rows; the
present export must not be mistaken for that complete optimization matrix.

**FINAL/FROZEN:** generation and fresh saved reconstruction replay passed.
Direct source manifests, Python syntax and local report links were checked.
The old feature reconstruction is explicitly inherited from its frozen
certificate; no independent second arithmetic implementation is claimed.
