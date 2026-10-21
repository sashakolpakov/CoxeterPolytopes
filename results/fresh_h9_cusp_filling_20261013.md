# Cusp filling of two known H9 chambers (13 October 2026)

No compact H9 example was found. Two exact obstructions explain why fairly
large classes of angle changes cannot work; a broader numerical scout found
no candidate. These conclusions do not exclude finite-volume chambers or
general compact H9 constructions.

## Construction and scope

The inputs are the two 18-facet noncompact chambers
`T1_D16_W4_7_12_14` and `T2_D16_W1_5_9_14` in the
[existing exact archive](../data/enumerations/alexeev_h9_doubling_cubes.json).
All indices below are zero based. Reverse and divide the archived Gram matrix
by two, giving diagonal 2, off-diagonal 0 for right angles, −1 for π/3,
and −2 for parallel walls.

The proposed operation deletes some old walls and changes old parallel pairs
to intersecting or disjoint pairs. Additional angle changes and additional
acute-compatible walls are allowed only as specified below. The aim is to
remove cusps while preserving Lorentz rank ten and Coxeter angles.
Throughout, all retained and additional outward normals must have nonpositive
pairings with one another and admit a strict timelike interior point. Arbitrary
changes below concern magnitudes of allowed pairings, not positive pairings.

## Exact obstruction for T1

**After deleting at most five old walls, changing only the old parallel
pairings cannot give a compact chamber, even after adding arbitrary
acute-compatible walls.** All retained nonparallel old pairings must remain
fixed for this statement.

There are 54 displayed affine principal subdiagrams unaffected by these
changes: 36 copies of affine A3 and 18 copies of affine D4. The certificate
lists one surviving subdiagram for each of the 8,568 five-wall deletion sets.
Deleting fewer walls cannot remove more subdiagrams. Every affine A3 has
positive kernel vector (1,1,1,1); every affine D4 has kernel vector
(2,1,1,1,1), with the center first. Exact rational elimination checks that
each block is positive semidefinite of corank one.

For completeness, the π/3 graph splits into K4,4 and a double star with
adjacent centers of degree five. The K4,4 contributes 36 affine A3 circuits
and eight affine D4 stars; the double star contributes ten affine D4 stars.
The smallest set hitting all 54 displayed circuits has six vertices. An
explicit transversal is `{0,1,2,9,10,17}`. This is a minimum for this circuit
family, not a claim that deleting those six walls permits a compact filling.

The geometric implication also covers additional walls. Let the positive
kernel coefficients be a_i and put δ=Σa_i r_i. For a strict interior
timelike vector c, ⟨c,δ⟩<0, so δ is nonzero and lies on the future null cone.
Every retained or additional acute-compatible normal w has ⟨δ,w⟩≤0.
The vectors c+sδ remain timelike and inside all half-spaces for every s≥0,
and approach the ideal point [δ]. Thus the chamber cannot be compact.

## Exact obstruction for T2

**After deleting at most five old walls, changing any nonzero old pairings
while preserving all old orthogonalities cannot give a compact H9 chamber,
even with additional acute-compatible walls.** This allows the π/3 edges,
as well as parallel edges, to change without a bound on their labels.

The π/3 graph consists of two disjoint stars with centers 1 and 3 and eight
leaves each. Among the 16 leaves, the only nonzero pairings form eight
disjoint parallel pairs. After five deletions, the retained leaf count is
11, 12, or 13, occurring in 4,368, 3,640, and 560 deletion sets respectively.

Normalize every retained normal to squared length two. The leaf Gram matrix
remains block diagonal, with singleton blocks [2] and pair blocks

```
[ 2  −w ]
[−w   2 ],   w ≥ 0.
```

Such a pair has determinant 4−w². Its only singular value on this range is
w=2, which is an affine A1 pair and prevents compactness by the same positive
null-vector argument. Otherwise every block has full rank. At least eleven
retained leaves would therefore have Gram rank at least eleven, impossible
in the ten-dimensional Lorentz space for H9. Additional walls do not change
this rank contradiction.

## Numerical scout beyond the fixed-angle obstruction

The separate [scout](../scripts/fresh_h9_cusp_filling_20261013.py) uses T1.
It chooses a ten-wall Lorentz basis containing no old parallel pair. The
remaining eight walls each have exactly one parallel partner in that basis.
After mutating basis angles, a Schur-complement norm equation gives the
possible new weight at each old parallel pairing. Cross equations then
test whether at least four extra walls can coexist, retaining at least
14 of the 18 original walls.

| Pass | Bases | Basis-angle mutations | Rank-ten proposals |
|---|---:|---:|---:|
| Initial | All 212 eligible bases | 171,272 | 0 |
| Broader | 13 symmetry representatives covering all 212 bases | 1,383,498 | 0 |

The initial pass changes one or two old π/3 edges, or changes one old right
angle. The broader pass allows up to four old π/3 edges to change. Labels
2, 4, 5, and 6 are used for multiple changes; single changes also include
8 and 10. A changed old right angle uses label 3, 4, 5, or 6. Every changed
edge is inside the declared basis. All other nonparallel pairings stay fixed.
New intersecting pair weights are recognized numerically among labels 2–60;
disjoint weights are allowed. Identically free diagonal equations receive
only a finite sample of weights. Parallel retained pairs are rejected
because the target is compactness.

This is floating-point scouting, not an exhaustive nonexistence proof.
Basis signatures, norm equations, and cross equations use tolerances. The
absence of candidates here does not extend either exact obstruction to
arbitrary angle changes. The two pass counts overlap and are not counts of
distinct chambers.
An independent source review checked this scope and the saved counts,
including the unique exterior-wall parallel partner for all 212 bases and
the 13 orbit sizes summing to 212. It did not rerun the broad numerical scout.

## Reproduction and independent replay

The obstruction archive and verifier are separate from the numerical scout.
The verifier imports no producer, numerical library, or optimization solver.
It uses standard-library rational arithmetic to check the source digest,
every positive affine kernel and corank, the complete five-deletion cover,
the six-vertex transversal, the T2 matching, and all leaf counts. It rejects
three corrupted certificates: a wrong kernel, a missing deletion case, and
a wrong matching. Result: **PASS**.
An independent agent reran the verifier and reviewed both geometric
implications, including the acuteness and strict-interior assumptions.

```sh
python scripts/fresh_h9_cusp_filling_obstructions_20261013.py
python scripts/fresh_h9_cusp_filling_obstructions_verify_20261013.py
python scripts/fresh_h9_cusp_filling_20261013.py
python scripts/fresh_h9_cusp_filling_20261013.py --broad
```

Artifacts:

- [Exact producer](../scripts/fresh_h9_cusp_filling_obstructions_20261013.py)
- [Exact certificate](../data/enumerations/fresh_h9_cusp_filling_obstructions_20261013.json)
- [Independent verifier](../scripts/fresh_h9_cusp_filling_obstructions_verify_20261013.py)
- [Replay result](fresh_h9_cusp_filling_obstructions_verify_20261013.json)
- [Initial numerical scout](../data/enumerations/fresh_h9_cusp_filling_20261013.json)
- [Broader numerical scout](../data/enumerations/fresh_h9_cusp_filling_20261013_broad.json)

Further work would need to alter the affine circuits of T1, change some old
leaf orthogonalities of T2, or use substantially different initial chambers.
Neither obstruction addresses all such constructions.
