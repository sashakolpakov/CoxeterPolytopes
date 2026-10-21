# Full-source review: off-center wall and conditional volume cost

20 October 2026. **PASS at the written-proof level.** I read the complete
11,953-byte source and rederived its geometry and inequalities. No mathematical
gap or required source correction was found. No scientific execution, build,
formal verification or external refereeing is claimed.

Reviewed source, unchanged:

* `results/compact_ra4_uniform_outgoing_center_cut_volume_cost_20261020.md`
* SHA256 `78135c4cca5d64c583727d612d5f22cf8f99f50a6a9b0aedbae17780e01c2b2c`.

## Contribution and independence boundary

This is a full check of the final source, not a claim of independent discovery.
I developed the preceding central-wall argument and contributed to the
off-center height closure. Before the combined source was assembled, I
independently rederived the local-frame identities, checked the adjacent-roof
gap case split, rational quadratic envelope, strict product contradiction and
fan arithmetic in discussion with the finite57 lane. The latter supplied the
combined off-center cone and height proof; root supplied the initial
no-orthogonal-roof lead. This review explicitly acknowledges that collaboration.

For the final source I also checked every earlier off-center step, including
the normalization, excluded orthogonality, support estimates, cone coefficients,
boost constants and uniqueness of the cut wall. No verdict is inferred from
other reviews.

## Checks of the actual geometric bridge

An additional wall is disjoint from the actual whole facet D, so its timelike
projection gives the stated future unit z. Acuteness against D's twelve walls
also puts z in D, justifying x_j>=0 in Section2. The q0-cut condition is exactly
p>=p0 Z. The orthogonal-roof branch leads to x_j>s and
p^2>=3/sqrt(5), after which the adjacent-roof gap leaves only x_j=p/phi.
The resulting norm equation contradicts x_j>s. Thus all twenty roof products
are strictly below minus one, not merely nonpositive.

I rederived the support bound and monotonicity of g, p^2<4/3 and r<1/20.
The rational comparison is valid:
2323/2400<897/926. The two-cut-wall distance estimate indeed makes their
normal product positive, establishing uniqueness without assuming their
centers coincide.

Both face-average identities in Section4 are exact. With
mu/s=phi*h/(p*sqrt(phi)*C-phi*Z), the residual is in the roof span;
its total coefficient and strict convex-hull condition reduce precisely to
the displayed inequality (4). Its two rational sides exceed351/2000 and
are below7/64 respectively. The strict interior convex-hull representation
can use all twenty coefficients positively. The denominator also excludes
m from the four-dimensional roof span, so the projected roofs span the
whole intrinsic four-space. This establishes the strict sign for EVERY
unknown intrinsic neighbor, not merely a selected list of future corners.
The resulting O_z is the actual intrinsic interior point used in the proof.

The boost bounds in Section5 are correctly oriented. Their exact constants
give the stated margin
-1058717/1080000>-1, forcing a genuine orthogonal roof by the compact gap.

## Height closure and complete facet count

In the vertex frame t=s*sqrt(phi)*a has no off-center error. The three
adjacent-roof alternatives are zero product or product strictly below
minus one; this gap is essential. Two top values violate y_j>=1, and
zero top values violate the displayed convex quadratic. Exactly one top
value gives a>=phi^(3/2).

I checked the S/a>=9/5 polynomial coefficient by coefficient. Its rational
envelope has exact endpoint values

    -0.46130238 at a=2,
    -0.020390924 at a=267/100,

where these terminating decimals denote rationals. They imply the claimed
weaker endpoint bounds. The replacement of T by its proposed upper bound
is legitimate because T>=2>U, so U^2+(T-U)^2 increases in T.

A is positive before replacing Z by hp. The spatial estimate follows from
||W_ell||^2/a^2<=J^2+phi-3<9/16. The explicit final margin is valid:
75649/2000000>3/80. Thus orthogonality to m contradicts the assumed low
height; no saved corner inventory is substituted for actual geometry.

Finally I checked the entire Klein normal-fan argument. Its triangles are
the cones of ALL actual vertices of the compact facet. The determinant
formula and exact comparison

    2062687500<2105030569

give area<4/15 per triangle. Consequently V>15pi>46, and trivalence makes
V even. The resulting V>=48, f>=101, W>=21 and Vol>=23pi^2/2 follow from
the stated inherited inequalities.

## Dependencies and conclusion

The two inherited source pins were freshly rehashed and match:

* Fixed-wall proof: `df5e75b3cfd8466b9b35f36b44f6c6b37f94b25afcdb276e7a443039b1d6da87`.
* CompactRA4 primary source: `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`.

Those sources supply the compact wall gap, normal-fan identity, capacity,
neighboring excess and volume formula. This check does not rerun their
historical formal or build evidence.

The proved consequence is conditional on an actual marked regular D and
all twenty outgoing lengths L. At or below the 120-cell volume, its matched
q0 is interior to P. This does not establish whole-120-cell containment,
universal minimality, a second-volume theorem, or a restriction on every
nonuniform collar. The reviewed source states these limitations accurately.
