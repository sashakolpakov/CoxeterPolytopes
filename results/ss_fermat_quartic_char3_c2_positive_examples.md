# Two cyclic finite-volume sections of the characteristic-3 Fermat chamber

## Outcome

The full lattice stabilizer of the complete 5,944-wall Kondo--Shimada
chamber has two cyclic fixed sections which survive exact global facet
reduction.  Both are finite-volume, noncompact hyperbolic Coxeter polytopes:

| involution class size | dimension | facets | finite vertices | ideal vertices | hyperideal rays | decision |
|---:|---:|---:|---:|---:|---:|---|
| 9,072 | 10 | 196 | 212,824 | 10,560 | 0 | finite-volume, noncompact |
| 4,536 | 11 | 206 | 327,566 | 18,354 | not needed for the literal proof | finite-volume, noncompact |

These claims do not use surface-automorphism reflectivity.  The algebraic
surface supplies an exact Neron--Severi lattice, a complete chamber wall
stock and its finite lattice stabilizer.  The restricted cones are then
certified by rational polyhedral and Coxeter computations.

## Complete Coxeter reductions

For the `H10` class, 5,944 walls form 2,982 involution orbits and 2,031
distinct primitive restricted rays.  Exact positive conic certificates
remove 1,835 inequalities, leaving 196 strict facets.  All 19,110 facet pairs
are Coxeter-compatible.  In the positive-wall convention their square
histogram is `2^80,(4/3)^80,4^36`.

For the `H11` class, the walls form 2,978 orbits and 1,892 distinct restricted
rays.  Exact positive conic certificates remove 1,686 inequalities, leaving
206 strict facets.  All 21,115 facet pairs are Coxeter-compatible.  Their
positive squares are `2^176,4^30`.

Every retained inequality has an exact timelike strict relative-interior
witness.  Every removed inequality is given as a positive rational
combination of retained normals.  Thus both displayed facet lists are
complete and irredundant, rather than locally selected wall samples.

## Literal Vinberg finite-volume proofs

For `H10`, the complete connected elliptic/affine census produces 10,560
distinct full parabolic zero-neighborhoods, all valid.  Vinberg's literal
edge/endpoint audit gives 1,427,720 elliptic rank-nine edges:

- 779,360 finite--finite;
- 569,520 finite--ideal;
- 78,840 ideal--ideal.

Every edge has exactly two endpoints.  There are 212,824 finite vertices and
10,560 ideal vertices, with no bad or unrecognized edge.

An independent GMP rational double-description run gives 223,384 projective
extreme rays: exactly the same 212,824 timelike and 10,560 lightlike rays,
and no spacelike/hyperideal ray.

For `H11`, the connected component census produces 18,354 distinct full
parabolic zero-neighborhoods, again all valid.  The literal audit gives
2,688,957 elliptic rank-ten edges:

- 1,177,760 finite--finite;
- 1,247,706 finite--ideal;
- 263,491 ideal--ideal.

Every edge has exactly two endpoints.  There are 327,566 finite vertices and
18,354 ideal vertices, with no bad or unrecognized edge.  This complete
Vinberg test is already a decisive finite-volume proof; a second `H11` cdd
run was deliberately stopped under system memory pressure and is not needed
for the theorem.

## Fixed lattices and lineage

The saturated `H10` fixed lattice has rank 11, determinant 6,144 and Smith
diagonal `2^10,6`.  Dividing its Gram by two and applying the archived
unimodular change of basis gives `<3> + <-1>^10`, the determinant-three
lattice of McLeod's known 14-wall `f_3^10` chamber.  All 196 mirrors remain
integral roots, with half-scaled types
`80*(-1,div1) + 36*(-2,div1) + 80*(-6,div3)`.  Solomon's exact Euler sums identify the new chamber as a
reflection-subgroup chamber of index 75,202,560, with covolume
`pi^5 * 46376/893025`.  Thus this is not a new lattice, maximal reflection
group or commensurability class; priority of the literal 196-wall subgroup
and its divisor extraction remains open.

The saturated `H11` fixed lattice is even of signature
`(1,11)`, determinant 1,024 and Smith diagonal `1^2,2^10`.  Its inverse Gram
has half-integral diagonal entries, hence its 2-elementary parity is
`delta=1`.  Nikulin's indefinite 2-elementary classification therefore
identifies it as

`U + E8(-2) + A1(-2)^2`.

This identifies the lattice genus; it does not identify the 206-wall chamber
with a previously published chamber.  The pinned CoxIter corpus has at most
50 facets per record, so neither the 196- nor 206-facet diagram occurs
literally in that snapshot.  This is a bounded exact-match statement, not a
claim of worldwide publication priority or a new arithmetic
commensurability class.

All 206 facet reflections act integrally on this fixed divisor lattice.  In
primitive fixed-lattice coordinates their square/divisibility types are
`164*(-2,div2) + 30*(-4,div2) + 12*(-2,div1)`.  The exact reflection audit is
archived together with the companion `V4-H11` stable-overlattice computation
in `data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json`.

## Reproducible artifacts

- Full cyclic class census and representatives:
  `data/enumerations/ss_fermat_quartic_char3_c2_element_class_screen.json`.
- Complete facet reductions and conic certificates:
  `data/enumerations/ss_fermat_quartic_char3_c2_survivor_screen.json`.
- Positive Coxeter Grams:
  `data/enumerations/ss_fermat_quartic_char3_c2_h10_coxeter_screen.json` and
  `..._c2_h11_coxeter_screen.json`.
- Connected component audits:
  `data/enumerations/ss_fermat_quartic_char3_c2_h10_component_audit.json` and
  `..._c2_h11_component_audit.json`.
- Literal Vinberg certificates:
  `data/enumerations/ss_fermat_quartic_char3_c2_h10_vinberg_certificate.json`
  and `..._c2_h11_vinberg_certificate.json`.
- `H10` independent cone certificate and prepared GMP input:
  `data/enumerations/ss_fermat_quartic_char3_c2_h10_exact_cone_certificate.json`
  and `..._c2_h10_exact_cone.ine`.
- Fixed-lattice bases and Grams:
  `data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json`.
- H11 mirror integrality and stable-overlattice audit:
  `data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json`.
- `H10` McLeod lineage, Euler index and full face vector:
  `data/enumerations/ss_fermat_quartic_char3_c2_h10_mcleod_index.json`.

The principal reproducers are
`scripts/certify_ss_fermat_quartic_char3_c2_element_classes.py`,
`scripts/screen_ss_fermat_quartic_char3_c2_survivors.py`,
`scripts/certify_ss_fermat_quartic_char3_c2_literal_vinberg.py` and
`scripts/certify_ss_fermat_quartic_char3_c2_exact_cones.py`.  The exact
lineage/index reproducer is
`scripts/certify_ss_fermat_char3_c2_h10_mcleod_index.py`; the H11 mirror
reproducer is `scripts/certify_ss_fermat_char3_h11_mirror_lattices.py`.
