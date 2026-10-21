# Independent audit of H4-block searches

The five saved H9, H12, and H13 finite prefixes all have valid exact
spacelike escape certificates. They do not give compact polytopes.
The three H9 presentations are integrally isometric to previously
screened E8 hosts; the audit independently verifies all three transports.
The H12/H13 searches use a third H4 block and are separate bounded
searches, without any assertion of worldwide lattice or group novelty.

The verifier is
[fresh_h4_blocks_independent_verify.py](../scripts/fresh_h4_blocks_independent_verify.py).
It imports only the Python standard library. Its arithmetic uses pairs
of `Fraction` coefficients, polynomial reduction by `phi^2=phi+1`,
exact sign decisions using rational squares, and its own matrix reduction.
It does not import the search engine, the project's Golden arithmetic,
or its cone verifier. The
[audit JSON](fresh_h4_blocks_independent_audit_20261013.json)
records all replayed artifacts and calibrations.

## Certificate and exhaustion checks

For every artifact, the standalone audit reconstructs the claimed block
Gram matrix, checks all root coordinates and reflection coefficients are
integral, and verifies the norms, pairwise nonpositive products, and
strict timelike interior point. It checks every level's actual height
`alpha*x^2/q`, the norm and time coordinate of each root assigned to the
level, chronological ordering, and the completed height bound.

The control-coordinate list is also independently exhaustive. Integer
bounds for `|x|` and `|x'|` give an enclosing box for the coefficients
`x=a+b*phi`; all elements in that box are filtered exactly using the
two norm bounds and time reflection integrality. The resulting set
must equal the saved control list, including levels with no accepted root.

Height-zero completeness is independently re-enumerated. For each
candidate norm `q`, every primary-dominant H4 weight coordinate satisfies
`w_i^2 <= 2q` at both real places. The audit takes the full Cartesian
product of that finite alphabet, reconstructs vectors using its own H4
inverse, checks their exact energies, joins the blocks and scalar, and
compares the entire dominant-root census with the saved data. This uses
no Schur-complement or wall pruning from the search implementation.

It finds one norm-two orbit for each H4 block, and the scalar norm-one
orbit when the scalar is present. Independent closure under the four
simple reflections gives exactly 120 integral H4 roots, all positive or
negative in the simple basis, and exactly one dominant norm-two root.
Thus the simple H4 walls account for the full height-zero system.

For each escape, the audit recomputes its positive norm, checks every
wall inequality, confirms the time coordinate is positive, and checks
the claimed active normals have rank equal to the hyperbolic dimension.
A strict timelike interior and a feasible spacelike direction certify
infinite volume of the finite cone. They do not establish nonreflectivity
of the full lattice: later mirrors could remove the direction.

| Dimension | alpha | H4 blocks | Scalar | Walls | Controls | Exact escape norm |
|---:|---|---:|---|---:|---:|---|
| 9 | phi | 2 | yes | 22 | 14 | `(246-152phi)/5` |
| 9 | 2phi | 2 | yes | 23 | 7 | `(466-287phi)/5` |
| 9 | sqrt5 | 2 | yes | 21 | 19 | `56-69phi/2` |
| 12 | 2phi | 3 | no | 36 | 7 | `(7576-4682phi)/5` |
| 13 | 2phi | 3 | yes | 65 | 7 | `290-179phi` |

All five searches complete height 30. Every displayed norm is exactly
positive at the chosen embedding. This audit verifies escapes for the
H9 artifacts; it does not independently reproduce their full vertex
and extreme-ray counts.

## Source-level enumeration review

No correctness error was found in
[fresh_h4_blocks_h9_search_20261013.py](../scripts/fresh_h4_blocks_h9_search_20261013.py)
for these screened parameters. The bounds have the following exact
justification.

The H4 block has diagonal two, chain products `(-phi,-1,-1)`,
determinant `5-3phi=phi^-4`, and an integral inverse. Its leading
principal minors are positive at both embeddings. For dominant weights
`w`, the energy is `w^T H4^-1 w`. After fixing the first `k` weights,
the minimum energy over the remaining unconstrained coordinates is
`w_prefix^T (H4_prefix)^-1 w_prefix`. This is precisely the source's
`SCHUR` array. It is a valid lower bound separately at both embeddings.

The extra truncation bound using the leading part of `H4^-1` is applied
only at the chosen embedding. Every entry of that inverse and every
dominant weight is nonnegative there, so omitted terms can only increase
the full energy. The source correctly does not apply this argument at
the conjugate embedding, where the signs need not agree.

The coordinate bound `w_i^2 <= 2*target` is Cauchy--Schwarz against a
simple root of norm two. Integral weights reconstruct integral vectors
because H4 is unimodular over `Z[phi]`. The product join retains every
block combination with the required total energy; the scalar-square
lookup is injective because only nonnegative primary scalar coordinates
are retained. Earlier-wall inequalities are tested exactly, and each
equal-height packet is checked for pairwise acuteness.

The complete norm alphabets are `{1,2}` for `alpha=phi` and `2phi`, and
`{1,2,2+phi,4+2phi}` for `alpha=sqrt5`. Unimodularity of the spatial
block gives primitive `q | 2alpha`; the ring is a PID and totally
positive units are squares. In the `2phi` case the only additional
candidate is norm four. Such a root would have every spatial coordinate
divisible by two, by reflection integrality and spatial unimodularity.
The norm equation modulo four then forces its time coordinate to be
divisible by two, since two is inert. The vector is therefore not
primitive. This argument works with or without the scalar summand.

The norm classes used here cannot represent parallel integral roots
related by a nonunit scale: no quotient of distinct listed norms is a
square in the field. Thus the absence of general parallel-root
deduplication in this engine does not affect the screened parameters.

The independent audit re-enumerates height-zero vectors and control
coordinates. It does **not** independently re-enumerate every positive
height vector; that part is supported by this source review, saved
certificates, and the following calibrations.

## Exact positive calibrations and duplicate detection

The one-block host `<-phi>+H4` produces five walls. Independent computation
of their Gram matrix gives the compact simplex chain `[5,3,3,4]`.
Enumerating every possible extreme-ray support gives exactly five rays,
all strictly timelike on the future sheet. The calibration artifact is
[fresh_h4_blocks_independent_h4_phi_h30.json](fresh_h4_blocks_independent_h4_phi_h30.json).
It is a known-diagram calibration, not a proposed new example.

The two-block H8 host `<-2phi>+H4+H4` produces eleven walls.
The independently replayed matrix `T` has determinant `5-3phi`, integral
entries and integral inverse, and satisfies

```
T^t B_E8 T = B_H4+H4.
```

It maps every one of the eleven saved normals exactly to the documented
Bugaenko normals under the saved permutation, so their whole Gram
matrices agree. Independent exhaustive support enumeration finds 41
future timelike extreme rays and no other feasible rays. Thus the H8
calibration reproduces the known compact Bugaenko chamber, with a
different Vinberg basepoint.

The audit also reconstructs the three E8 target forms and replays each
saved H9 isometry for `alpha=phi,2phi,sqrt5`. In every case it verifies
Gram congruence, determinant `5-3phi`, integral forward matrix and
integral inverse. Consequently all three H9 H4 presentations have
exactly the same integral orthogonal and reflection groups as the
corresponding earlier `<-alpha>+E8+<1>` hosts. Their different finite
wall lists do not establish new groups or new polytopes.

## Reproduction

```sh
python3 scripts/fresh_h4_blocks_independent_verify.py \
  results/fresh_h4_blocks_h9_phi_h30_20261013.json \
  results/fresh_h4_blocks_h9_2phi_h30_20261013.json \
  results/fresh_h4_blocks_h9_sqrt5_h30_20261013.json \
  results/fresh_h4_blocks_h12_2phi_h30_20261013.json \
  results/fresh_h4_blocks_h13_2phi_h30_20261013.json \
  results/fresh_h4_blocks_h8_2phi_h30_20261013.json \
  results/fresh_h4_blocks_independent_h4_phi_h30.json \
  --calibrations --output results/fresh_h4_blocks_independent_audit_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_h4_blocks_independent.py
```

All 13 regression tests pass. The suite also rejects a forged escape norm, a timelike
interior substituted for the escape, a missing zero-root control level,
and a false height record.
