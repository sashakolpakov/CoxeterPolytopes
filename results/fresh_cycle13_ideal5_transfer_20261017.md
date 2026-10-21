# The missing cycle13 distance-six transfer

An actual quadrilateral with a varying affine cycle of order13 and marks
at distance6 has a retained payment at another corner. This permits paying
only the65 cycle pairs of distance at most5. Its13 distance-six pairs may
be removed. No additional finite coupon is needed, and no restriction on
arithmeticity, cusp simplicity or spectator rank is imposed.

The [source](../scripts/fresh_cycle13_ideal5_transfer_20261017.py) and
[certificate](fresh_cycle13_ideal5_transfer_20261017.json) adapt the complete
[cycle14 transfer](fresh_cycle14_ideal5_transfer_20261017.md) to arc lengths6
and7. The decisive new input is the
[strict spectator spectrum below one](fresh_spectator_energy_below_one_20261017.md).
The earlier order13 support diagnostic explicitly left this universal
short-path obligation open; absence from its preferred maps alone did
not remove these cycle pairs.

Use the actual factor-fiber face star W+a+b+c+d, with opposite pairs a,c
and b,d. At the chosen ideal corner, the cycle F is orthogonal to a
positive-definite spectator diagram Z. White-plus-one-black edge diagrams
are strictly positive definite. Opposite heights h,k at its positive cusp
kernel are positive. Adjacent black contacts or a shared nonzero spectator
contact already supply a retained short path, so assume neither occurs.

Deleting a from F leaves A12. Its complete one-root positive-semidefinite
extensions are isolated, a bond3 or4 at either endpoint, a simple bond at
either second node, or simple bonds at both endpoints. Exceptional finite
and affine types cannot contain A12. This is the same all-label path
extension classification used in the preceding cycle theorem.

Each of the six single-contact profiles gives a finite marked corner:
the opposite fixed simple leaf and its fixed degree-two neighbor cannot
be part of a classical affine completion. Additional spectators are
orthogonal to the old path, so those degrees remain fixed. Affine cycles
have no leaf, affine C has doubled terminal bonds, affine B has a doubled
terminal at one end and a fork at the other, and affine D has forks at
both ends. Exceptional affine diagrams are too small. Hence the complete
connected corner component is finite, even if spectators attach at the
new root. Its marked distance is at most7 and receives the finite8 unit.
The checker reconstructs all six positive-definite cores, distances, and
fixed leaf/neighbor degrees. Every other corner component is a principal
part of the positive-definite spectator diagram, so the actual corner
is finite as claimed.

Only isolated and both-end3 profiles remain. The inherited normalized
Schur determinant is separately concave in the two reciprocal heights.
On the full allowed rectangles its four ordered profile-pair minima are

|Profiles|Minimum|
|---|---:|
|isolated, isolated|-5/21|
|isolated, both-end3|1/14|
|both-end3, isolated|1/14|
|both-end3, both-end3|8/21|

All16 rectangle values are checked exactly. Spectators add a nonnegative
quadratic form. Thus any case containing both-end3 is impossible under
the unpaid-face hypothesis.

For two isolated profiles, put alpha=1/6+1/7=13/42. The determinant identity,
reverified formally in the source, is

    det(S)/(h²k²)=1-2*alpha/h²-2*alpha/k²
                 +alpha*(x/h-y/k)^T G_Z^(-1)(x/h-y/k).

The isolated positive heights are native opposite-contact magnitudes:
each equals1 or is at least sqrt2. If either is at least sqrt2, the first
three terms are at least1-3*alpha=1/14>0, a contradiction. Consequently
h=k=1 before using any unscaled contact energy. The determinant equation
then forces

    D=(x-y)^T G_Z^(-1)(x-y)=4-1/alpha=10/13.

The supports of x,y are disjoint, and Z+c,Z+d are strictly positive
definite principal edge diagrams. The complete strict spectrum therefore
applies: below one the possible energies are0, r/(r+1) for integral r>=1,
2-sqrt2, and1-sqrt5/5. The rational value10/13 is neither irrational
exception and would require r=D/(1-D)=10/3, which is not integral.
This contradiction exhausts the final profile. Rotation and reflection
cover all13 distance-six pairs. The proof uses neither a finite cutoff
on spectator rank nor a finite list of native angles.

Generation and a fresh separate-process saved replay terminated exit0/PASS
on the six single-profile diagrams,16 exact rectangle values, formal
isolated determinant, native-height gate, spectrum exclusion and65/13 pair
partition. A separate full internal source and mathematical read passed
the complete inherited profile scope, fixed-leaf exclusion of affine
completion, strict edge hypotheses and height-before-energy argument.
That read did not repeat the6/16 arithmetic checks. Direct bindings,
syntax and report links pass. Source, data and report are final/frozen.
This is a local transfer, with no new ambient dimension or chamber claim;
external refereeing and proof-assistant formalization are not asserted.

    python3 scripts/fresh_cycle13_ideal5_transfer_20261017.py
    python3 scripts/fresh_cycle13_ideal5_transfer_20261017.py --verify
