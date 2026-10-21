# Root review: repeated A2 corners and the global drift

20 October 2026. Complete independent written review of
`compact_a2_corner_short_edge_capacity_20261020.md`, SHA256
`db7016f34053492c4b35b15b181823c7cb952d741654d709fe7a85daa8982048`.
PASS. No arithmetic hypothesis, execution, build or Lean verification.

The A2 inverse Gram, height formula and exact threshold31/50 check.
I rederived the entire short-wall inventory, including a release at an
uncontacted mate. A double contact cannot give a short edge: the minimal
3,3 profile needs another block, and the next3,4 profile already exceeds
the threshold. Three single-contact blocks are also excluded. Two blocks
leave precisely3,5 after the affine3--6 exclusion and the exact7-angle
bound. All other short walls are isolated single contacts, with both
finite and dotted classes retained.

I checked the actual height products and spatial contributions in every
pair case. Two3,5 walls cannot be orthogonal; neither can one coexist
with an isolated wall. For isolated walls in the same block the positive
spatial term is too large for acuteness. In different blocks, the only
allowed negative product is-1/2, so three would give an affine triangle.
A finite isolated wall cannot participate even in a pair. Thus the
two-wall case contains only dotted walls, and the only possible finite
two-direction case is the isolated7 wall. The count concerns directions,
not merely the number of closing walls.

The7 constant is the sensitive part of the endpoint estimate. I checked
the cubic for cos(2pi/7), both rational brackets, and the two half-angle
comparisons. The second has exactly the stated positive integer margin
884. The first rational comparison also has the required strict sign.
The one-short-direction cases retain an additional long direction, so
the common lower bound on two selected directions is justified.

All released white diagrams coincide, giving one shared f. The norm of
the dual sum is sqrt(2n); the edge factor is sqrt(4/3). At n=14 the
squared rational sum574/125 exceeds21, and the ratio increases thereafter.
The argument therefore excludes A2^r from the drift witness for n>=14
for ANY positive shared white weighting. Mixed diagrams and higher-rank
components are still missing; this is not a dimension-13 theorem.
