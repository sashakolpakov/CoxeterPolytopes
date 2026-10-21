# Independent review: the scalar facet lattice

20 October 2026. **PASS, within the scope below; no mathematical flaw found.**
This is a full written read and independent algebraic review, not a scientific
execution or a global compactness certificate. No worker, import, replay,
Lean, or typesetting job was run for this review.

## Bound sources and scope

Reviewed in full:

* `results/compact_h9_sqrt3_scalar_facet_next_edge_20261020.md`, SHA256
  `f6be0faf51bdc8188a3f8363fe2c5e3ee35e16c356ba49dd1546574b5195567d`.
* The earlier complete-host/first-edge derivation is covered by the separate
  review `results/review_compact_h9_sqrt3_diagonal_route_resource_20261020.md`,
  SHA256 `49274ae03f72772c9ee95980fc181e6b013f1e00855469655a19332e3618fa4c`.

The saved first-edge certificate is SHA256
`10c8afa329dd9328c924c17562ba7fd33500dde2b4686f5e6eccb6f876d05eb3`.
Its existence and identity do not turn this written review into a replay.
The successor accurately distinguishes the recorded generation from its
then-deferred replay.

## Independent checks

Write s=√3, ε=2+s, and O=Z[s]. Dividing the selected vector by ε gives
the stated primitive scalar root ρ=(a;b⁹), with square 3−s. The selected
square is ε²(3−s)=9+5s. The displayed relative-unit identity identifies
this mirror with the earlier analytic upper mirror; it is not a new
additional cap.

The excluded norm class 3+s is genuinely empty in the literal host lattice.
Reflection integrality forces all spatial coordinates to be divisible by s.
After dividing the norm equation by s, reduction modulo (s) gives A²=−1
in F₃. For the remaining norm-6 class the same calculation instead gives
A²=1. Thus two such roots have product divided by s congruent to −AA′,
which cannot vanish. In particular distinct such mirrors cannot be
orthogonal. The finite-angle norm and integrality tests exclude every
oblique intersection with the other possible root classes. For arbitrary
mirrors, rather than chamber-neighbor labels, the angle 5π/12 must also be
allowed in the intermediate list: its squared-cosine norm is the same as
for π/12, so the exclusion is unchanged. Discreteness makes any intersecting
pair generate a finite dihedral subgroup. Parallelism is excluded by the
positive-definite conjugate embedding. Consequently every actual neighbor
of ρ is orthogonal to it.

I independently derived the integral change of coordinates and its inverse.
The coprimality of a and c gives A=ct and Σy=at; the integral shift
yᵢ=zᵢ+ht gives Σz=(4+s)t and Q=Σzᵢ²−εt². In the free basis t,z₁,…,z₈,
the Gram matrix is exactly the displayed one. Its Schur complement is
(1−s)/9 and its determinant is 1−s. This is an integral identification,
not a diagonalization merely over K.

The primitive-reflector divisibility q∣2 det(M), total positivity and unit
normalization leave the necessary classes 1,ε,2,2ε. Modulo 2, the norm
is (1+s)t² and takes only residues 0 and 1+s. The classes 1 and ε have
different residues and are impossible. Thus every primitive M reflector
has norm 2 times a unit. The surviving classes are necessary possibilities;
the note does not assert that each is represented at every contact profile.

The extension argument is sufficient without a separate gluing theorem.
For every such reflector r∈M⊂L, the coefficient 2/(r,r) is an O-unit.
Since L is integral, its reflection preserves all of L and fixes ρ.
Saturation M=L∩ρ⊥ also preserves primitivity. Conversely an ambient
orthogonal reflection restricts to M. Nonorthogonal ambient mirrors are
disjoint from ρ and impose constant-sign inequalities on it. This justifies
the full M-chamber identification, and hence the conditional implication
that a compact full host chamber would give a compact M chamber.

Finally the A₇ extension inventory and all six right sides in (3) check.
The two isolated norm classes give 144 and 144ε; the A₈, D₈, E₈ and B₈
rows give 81,36,9,18ε. The middle branching position is affine. Directly
substituting z₁=(4+s)t−S₈ and the two-level eight-coordinate vector yields
Z²−8(s−1)t²=72q−9k(8−k)g² with Z=9S₈−8(4+s)t.

## Limits

No exact identification with a known nonreflective integral lattice was
obtained. The visible A₈ orthogonal splitting has nontrivial index-nine
gluing; similarity over K does not settle lattice reflectivity. No complete
orbit interval or finite upper mirror for the next relative equation is
proved here. The note therefore supplies a valid smaller conditional
stopping test and a concrete next edge, not a new compact H⁹ polytope,
an all-height root census, or a nonreflectivity theorem.
