# Independent review of the local edges and port balances

The [independent checker](../scripts/fresh_sqrt5_universal_self_edges_20261016_review.py)
and [certificate](fresh_sqrt5_universal_self_edges_20261016_review.json) verify
the [universal local edge construction](fresh_sqrt5_universal_self_edges_20261016.md)
using the project's separate Golden arithmetic. No producer edge code or
producer number-field implementation is imported. Generation and fresh saved
replay pass.

The finite component catalogue is an explicit frozen dependency. The reviewer
matches all 68 types, verifies each literal Gram's labels, integrality and
positive definiteness at both embeddings, and independently inverts it.
For all 724 root positions it constructs the duplicated-root Gram, checks
the physical and conjugate Schur signs, and verifies the saved endpoint
vectors directly against that full Gram. The endpoints are future-compatible
timelike rays with the required zero sets. Their distance and the integral
reflection coefficients are checked. The largest conjugate threshold is
independently recovered as 30.

A generating-function calculation grouped by rank multiplicities, using the
coefficients of `(1-z^r)^(-m_r)`, independently recovers all 311,574 rank-21
finite forests. Orthogonal spectators extend the component witnesses to
each forest. For D5+2E8 the reviewer checks all 21 full edge witnesses and
their strict real-radical null escapes. It also reproduces the negative
definite symmetric blocks of all 210 canonical copied squares.

The literal rational dual inequalities, the exact lower bound above 40, all
10,626 root-port balances and the simple abstract port graph pass separate
checks. Reconstructing the original LP's face coefficients remains the scope
of its frozen producer and independent earlier reviews.

These checks establish local edges and the stated incidence relaxation.
They do not realize simultaneous edge stars, arbitrary quadrilateral
closures, a common arithmetic lattice or a finite-volume polytope. No
dimension bound is improved by this diagnostic.

```sh
python3 scripts/fresh_sqrt5_universal_self_edges_20261016_review.py
python3 scripts/fresh_sqrt5_universal_self_edges_20261016_review.py --verify
```
