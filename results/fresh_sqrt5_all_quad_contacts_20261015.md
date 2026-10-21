# One-white contact domains for all Q(sqrt5) quadrilateral cores

The complete necessary one-white contact census has finished for all **8,298 marked cores** in the [all-union census](fresh_sqrt5_all_quad_unions_20261015.md). The [terminal certificate](fresh_sqrt5_all_quad_contacts_20261015.json) records 81,772 full contact trials from 448,983 domain-search states. Unlike the earlier unpaid-core test, there is no distance-two restriction and no H4-payment shortcut. Every processed domain accounts for all finite labels 2, 3, 4, 5, 6 and 10; dotted distances remain parameters.

| Core disposition | Count |
|---|---:|
| Certified affine principal subdiagram; excluded | 319 |
| Sealed nonsingular finite core; full-star rank excluded | 1,224 |
| Sealed singular finite core; retained | 636 |
| Finite core with accepted white contacts; retained | 4,044 |
| Symbolic dotted core; retained | 2,075 |

All 319 affine-core exclusions occur among the dotted classes. Of the 2,075 retained dotted cores, 175 have empty nonzero contact domains and are structurally sealed, with distance and rank feasibility still unresolved. Across all classes, 18,964 finite and 15,532 dotted profiles are retained. The other trial dispositions are 13,269 fixed finite-cycle failures, 21,513 physical negative-index failures, 11,682 conjugate-PSD failures, and 812 certified affine extensions.

Among the 4,044 finite cores with extensions, 24 are singular and 4,020 are nonsingular. The frozen source/data and complete saved-domain replay pass; syntax, direct source bindings and report links pass as well. Independent reconstruction of the domains and exact filters is being reviewed separately.

The [producer](../scripts/fresh_sqrt5_all_quad_contacts_20261015.py) tries every contact vector for one added white node, requiring all four enlarged corner diagrams to be spherical. Unassigned contacts are right during the search. Failure of positive definiteness cannot be repaired by adding further acute contacts, so this prunes complete branches safely. The all-right vector is omitted from the recorded nonzero domain and remains possible. No lower bound on a geometric edge length or finite label cutoff is imposed beyond the complete arithmetic-field alphabet.

Every full profile is tested for fixed finite-cycle field obstructions. Cycles containing a dotted edge are ignored by that test, since its squareclass has not been specified. For a finite-contact core, the whole extension must have Lorentz index at most one and positive-semidefinite conjugate Gram matrix. Certified affine principal subdiagrams exclude compact realizations. Failure to find one of the implemented affine witnesses is not asserted to establish absence of all affine subdiagrams.

The exact arithmetic uses the already reviewed positive radical rescaling into Q(phi). Let C be the old spherical corner, B its contacts to the two remaining black nodes, and S their two-by-two Schur complement. For a new white node with old-corner contact column w, remaining-black column v, and diagonal d, the cached calculation gives

```
rho = d - w^T C^(-1) w,
q   = v - B^T C^(-1) w,
S_new = S - q q^T / rho.
```

The enlarged corner is finite, so rho is positive in both embeddings. The physical S_new has two negative diagonal entries and passes precisely when its determinant is nonpositive. Its conjugate passes precisely when its two diagonal entries and determinant are nonnegative. During generation, `--check-full` compared all 52,159 finite, field-compatible trials with the frozen complete elimination, with exact agreement. The optimization does not change the necessary matrix tests.

A finite core with no accepted nonzero profile is sealed against all external white contacts, not just a bounded number of additions. Any external white node touching the old core would itself give one of the tested profiles. A path through other external whites cannot avoid its first direct contact with the core. Thus all external whites are orthogonal to every core node, and all marked corner components of the core are entire components in the full corner.

If such a sealed core is nonsingular, it cannot occur in a full quadrilateral star. Its Gram matrix and the positive-definite Gram matrix of the additional white support form a nonsingular block diagonal matrix, whereas the full two-face star has one more normal than the ambient Lorentz dimension. This is a full-star rank exclusion. A sealed singular core remains as a possible local star with fixed entire marked components. Cores with accepted profiles and symbolic dotted cores need further compatibility, distance, rank or geometric analysis; they are not silently excluded.

If a symbolic dotted core has no accepted profile, the same direct-contact argument seals its marked corner components uniformly in the unknown distances: the implemented dotted tests are necessary relaxations. Its rank and distance feasibility remain undetermined, so it is retained under the dotted status rather than classified as a nonsingular rank exclusion.

Each complete per-core record and the next class index were committed atomically to [SQLite](../data/enumerations/fresh_sqrt5_all_quad_contacts_20261015_running.sqlite). The terminal next-class index is 8,298; the database and source/data files are now immutable, despite the historical `running` filename. Resumption used that index, so no completed class was rerun or skipped. Saved replay checks the database, the contiguous completed prefix and every recorded domain.

The contact lists can also impose necessary restrictions on a finite-component embedding in a weight model: each proposed outside vertex's contacts to the embedded old corner component must be the restriction of an accepted full-core profile, or of the omitted all-right profile. Separate compatibility of all outside vertices is still required for a realization; using only these individual restrictions is a necessary relaxation.

No new dimension bound follows from this package alone. All compactness exclusions here concern compact arithmetic Coxeter polytopes over Q(sqrt5), with the inherited ambient finite-corner and field hypotheses.

The [independent review](fresh_sqrt5_all_quad_contacts_20261015_review.md) passed generation and saved replay. It separately reconstructs all 81,772 contact vectors, all 52,159 finite arithmetic filters and 1,131 positive affine kernels, and audits the sealing argument and terminal database. It does not claim a new dimension bound.

```sh
python3 scripts/fresh_sqrt5_all_quad_contacts_20261015.py --verify
```
