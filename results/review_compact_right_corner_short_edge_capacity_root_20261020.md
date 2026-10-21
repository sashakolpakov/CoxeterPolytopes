# Root review: two short directions at an all-right corner

20 October 2026. Independent full written review of
`compact_right_corner_short_edge_capacity_20261020.md`, SHA256
`aad3f1d9509d716cab98f228061b88783f7e0c848aa849e3f09351fc49b05771`.
PASS. No execution, external classification or dimension theorem.

I independently derived the actual-wall coordinates and edge formula.
An equivalent shortness condition is beta_i^2>3*a^2. The endpoint
positive-definiteness condition gives a^2<1/2. The missing parallel case
is correctly excluded by the positive null sum of the two unit normals,
which would satisfy every other wall halfspace and give an ideal point.

The support census is exhaustive. Four contacts would force three
retained squares equal to1/4, and shortness would force the fourth equal
to1/4 too, making a=0. With three contacts the released square can only
be1/4,1/2 or(3+sqrt5)/8. The first two leave no possible sum of two native
squares; the third forces the other two to be1/4. With two contacts the
smaller square must be1/4 or1/2. The inequalities give exactly the two
finite families and the dotted-plus3 family in the source. A single
contact must be dotted. No formal profile is asserted to be realizable.

Two actual short closing walls have product greater than-1/2, hence
must be orthogonal by all-wall acuteness and the Coxeter gap. Their
supports share exactly one coordinate. I checked each exceptional
profile against the resulting equation beta_j*gamma_j=a*b. The
(3,3,5) case requires both shared contacts to be3 but its squared height
h satisfies3*h/8<1/16, excluding equality. The remaining two-contact
families can share only their3 coordinate. Their height product must
then be1/4; the finite heights are below1/2 and the dotted heights above
1/2, so at most one of each occurs.

The direction count is stronger than a distinct-wall count and is
correctly proved. A finite wall closing two short directions has
a^2<1/12, while compatibility with a dotted-plus3 wall would require
a^2>1/6. Thus the two-wall case closes at most one short direction per
wall; every one-wall case closes at most two.

At an A1^n corner all rank-(n-1) white diagrams coincide, so ANY globally
shared positive white weight cancels from the drift comparison. I
checked the exact n=14 inequality (the squared difference is193) and
its monotonic extension. This excludes A1^n as a drift witness for
n>=14, without any arithmetic restriction. Other vertex diagrams
remain uncontrolled. It is not a dimension-13 theorem.
