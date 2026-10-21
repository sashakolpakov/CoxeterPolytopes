# Root review: quartic order-exchange obstruction

20 October 2026. Full paper review; no worker, Lean or TeX execution.

Source: `results/compact_biquadratic_unit_order_exchange_obstruction_20261020.md`,
SHA256 `043f1f1f6605f6e0b80c1051fd6d288f119e6fbb65487af45d092e1669320299`.
Its two actual-edge inputs are explicitly pinned and reviewed separately.

**PASS for the literal N9 obstruction and its stated n>=8 extension.**
This excludes a construction family; it gives no uniform dimension bound,
new example, field-wide or commensurability-class exclusion.

I checked the alternative actual vertex using positive-definite
completion. For retained old normals n_J, put beta=-c.n_J and
w=B_J^-1 beta>=0. If delta=1-beta.w>0, the unnormalized point
`p+(a_c/delta)*(c+sum w_j*n_j)` is future timelike and pairs
nonpositively with every actual wall. Here the released old wall is
strict, so the finite completion is an actual new vertex. The explicit
solution is Z=(1;1,1), with unit-scaled negative time vector T_Q.

The Q lattice is the full orthogonal direct summand, not the index-s
span of its two simple roots. Its displayed E,F Gram determinant and
the time-vector determinant are units. For a nonorthogonal crossing
mirror, the rank-three reflection subgroup is finite B3. In its
ordinary root model a conjugate of one old norm-two root has the
stated contacts. Integral conjugation preserves both its exact norm
and primitivity. The two possible Q projections are both nonintegral,
contradicting the integral direct-summand projection. This argument
does not assume that the original B2 was parabolic. The converse
extension of intrinsic reflectors is also justified by the direct sum.

Thus both orders P then Q and Q then P are actual constructions in
one chamber, and sealing prevents mixed additional mirrors at their
vertices. I checked the A/B frame Gram equality, their unit determinants,
the time-sheet orientation and every positive simple-root image.
These give an integral full-arrangement isometry stabilizing the same
chamber, which is stronger than an abstract equality of Gram matrices.

Independently, using the alternative Q basis
`(beta+phi*c)/s, beta`, I derived the trace as
`K^2*H^2-K*(K*H+1)*(I+J)/a-K*(2+C)/k+C*(J-I)/s`.
The identities I+J=a(H-1), I-J=sK, H=alpha K and C=K-alpha reduce
it to the same polynomial and exact coefficient expansion as the
source. The trace exceeds the active dimension, so the isometry has
infinite order. Free chamber action then gives infinitely many
distinct cosets modulo the full reflection group. The named arithmetic
lattice theorem supplies the infinite-covolume conclusion.

For n>=8 the recursive residual dimension is at least6, exactly where
both complete edge inventories remain valid. No small-B3/F4 exception
is extrapolated to dimension four. I contributed the independent trace
derivation and adversarial checks; this is root review, not external
refereeing or completed Lean verification of this geometric instance.
