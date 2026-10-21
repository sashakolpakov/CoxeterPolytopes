# Mixed right-angled hyperbolic dimension is at most eleven

19 October2026. **Internally proved, with independent complete written
reviews, exact generation and replay, and conditional Lean verification
of the final count logic.** The published ceiling is12; examples reach8.
This theorem improves the upper bound to11, leaving a gap of three.
It makes no new-example, sharpness, external-refereeing or worldwide
priority assertion.

**Theorem.** A finite-volume right-angled hyperbolic polytope with both
finite and ideal vertices has dimension at most11.

## Proof and the new step

The published Dufour ceiling12 reduces the problem to dimension12.
Let f be the number of finite vertices and B the number of canonical
three-dimensional triangular bipyramids, each with two finite and
three ideal vertices. At each finite vertex the triangles in its
bad-triangle graph correspond exactly to these three-faces. Thus the
actual finite opposition graph has B edges and degree sum2B.

The complete earlier polygon/three-face coupling gives

    2B > (184460/1619) f > 113.93452 f.

The new [positive-part opposition proof](right_angled_mixed_positive_part_opposition_cap_research_20261019.md)
uses the actual invariant h=3q-sum n_i on each opposition component.
Its square is3. A fixed future null vector perpendicular to h gives
positive energies E(q), and every formal triple has ratio
r=E(q_S)/E(q)=5/4+d with |d|<=3/4. The exact uniform220-triple
moments and a full-interval quartic majorant

    P(d)>=max(0,r^-2-1/2)

give average P<=1420561/5500000 at every finite vertex.
Nonnegativity permits inserting absent formal triples. Summing actual
ratios over ordinary finite vertices pairs each actual undirected edge
as r^-2+r^2-2>=0, and yields

    2B <= (1420561/12500) f = 113.64488 f.

This is an ordinary average bound, with no maximum-height or weighted
vertex substitution. The positive cross-product difference is

    184460*12500-1420561*1619=5861741.

The two inequalities contradict each other. The same argument with
the accepted compact4 and entirely ideal6 ceilings also gives the
overall finite-volume right-angled ceiling11.

## Complete verification record

* New proof source SHA256:
  `19ec3d814759cd2d9b9b1df8ef128f1df7d7afee68590aa3310479c954a65c9f`.
* [Root full mathematical review](review_mixed_positive_part_opposition_cap_root_20261019.md)
  and [independent complete review](review_mixed_positive_part_opposition_cap_resource_20261019.md),
  SHA256 `358969924c4ce3681084ee38755468e0917cf3fc2cdda5b9747dcd33758abf56`.
  Both reconstruct the actual-object geometry, moments, majorant and
  entire opposing polygon/Q8/exceptional-incidence lower count.
* [Exact coefficient certificate](right_angled_mixed_positive_part_opposition_cap_research_20261019.json),
  SHA256 `2b1e08dd3b2ce25efec25e5d3947ba581ed8ef3f06582454d8776296fb13c7e1`.
  A standard-library Fraction checker independently expands the degree-six
  polynomial and all220 subset moments through degree four, and records
  every Bernstein coefficient, central bound and final rational comparison.
  Separate [generation](right_angled_mixed_positive_part_cap_generation_20261019.json)
  and [replay](right_angled_mixed_positive_part_cap_replay_20261019.json)
  pass in0.086/0.084s, with peak process-tree RSS below56MB.
* `CoxeterProofs.RightAngledMixed11` compiles and passes its separate
  six-theorem printed-axiom audit in0.56/0.34s. Source SHA256
  `f0b1cffede6210dbddf7cc7d95f29d14baa5984bcd2e19641e029afb1e7f864d`.
  It derives the coupling and strict lower bound from explicit incidence
  interfaces, proves the numerical contradiction and composes the
  dimension12 tail. Only propext and Quot.sound occur; the constant
  identity is axiom-free. No sorryAx or unrecorded dependency build.

The scalar checker does not enumerate polytopes. Standard acute/link
geometry, Dehn--Sommerville, the published82-facet eight-dimensional
minimum and Dufour's ceiling12 are external inputs. The actual opposition,
Q8 attachment and exceptional-face arguments are new/earlier internal
written geometry. The Lean module assumes the cleared aggregate analytic
transport comparison, not its real-valued proof or the polynomial's
geometric realization. These boundaries remain explicit.

The [proof and reproduction map](../docs/right_angled_dimension11_proof_map_20261019.md)
records the single focused manuscript in the mixed family. Its
[11-page PDF](../manuscripts/build/right_angled_dimension11.pdf) now
contains the complete geometric,scalar and counting proof,with three
visually reviewed vector figures and zero build warnings. The
[current source-bound package ledger](mixed_dimension11_package_verification_20261019.json)
collects all exact,conditional-kernel and final-PDF records. All temporary
TeX files and raster previews were removed.
The previous transport manuscript and its conditional8 results remain
supporting sources; no additional publication is inferred from that history.

## Literature comparison after the proof

Dufour's published12 ceiling is recorded in
[Alexandrov2023](https://arxiv.org/abs/2303.09533) and still used in
[Egorov's31March2026 preprint,Introduction](https://arxiv.org/html/2603.29925v1).
The latter also records examples only through8 and the published compact4
and entirely ideal6 ceilings. This check was made after deriving our proof;
it establishes the comparison used here, not a complete priority search.
