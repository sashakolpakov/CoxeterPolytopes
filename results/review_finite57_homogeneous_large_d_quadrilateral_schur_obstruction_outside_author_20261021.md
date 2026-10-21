# Full outside-author review: homogeneous large-D quadrilateral obstruction

21 October 2026. **PASS at the exact restricted scope below; no correction requested.**
This is a manual mathematical source review. No scientific worker, symbolic
package, enumeration, test, Lean or TeX process was run or imported. Only
source reads and streaming file hashes were used. No bound file was changed.

Reviewed source: `finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md`,
8,485 bytes, SHA-256
`3ef4130fc30f557425fdb06144342054e019c6eedda3e54a197c532b158b10ac`.

I did not author this source or its Schur argument. I authored the cited
A17 endpoint-producer scope audit and earlier related finite57 incidence/
financing work, and reviewed the Context662 composition. Those inherited
contributions are disclosed; this review independently derives the new
restricted geometric obstruction and does not re-referee my own contextual
foundations.

## Precise theorem reviewed

Source lines11–28 exclude an **actual quadrilateral two-face** whose complete
positive supporting-white diagram is `A_(r-4)+2A1+Z`, with `r>=7`,
`rank Z=n-r`, with Z orthogonal to every black side, and with all four
active corners equal to whole D_r components marked at their long endpoint
and fork center. The complete white support, common Lorentz realization,
and unchanged orthogonal Z are hypotheses. An unsealed D_r subdiagram at
one corner does not establish them. The conclusion is independent of
arithmeticity and of the types, labels and multiplicities inside Z.

## Actual marks and native gap

The marked-corner identification in lines49–68 is valid. Deleting the long
endpoint and fork center of D_r leaves a unique A_m path and two isolated
leaves, where `m=r-4>=3`. The long mark has exactly one white contact; the
center has exactly three. A particular actual black wall therefore cannot
switch these two roles between its two corners, and the roles alternate
around the quadrilateral. Fixing a center completion fixes which end of the
white path is adjacent to the center; both long completions attach to the
other end. The two center completions join that same path end and both
isolated leaves. This proves the claimed simultaneous white identifications,
rather than separately relabeling four diagrams. The long mark and center
are nonadjacent in every D_r under consideration, so all adjacent black
products are zero. Neither D4 triality nor an isolated-white/path ambiguity
occurs here.

All roots may be normalized to square2, because every nonorthogonal bond
inside these D_r corners is simple. For distinct actual chamber facets,
the source's gap `-(a,c) in {0} union [1,infinity)` is correct in this
normalization. If their two-root Gram is positive definite, the positive
projection argument in the pinned TeX at lines175–184 gives a common finite
point strict against all other walls. Thus their angle is an actual native
facet angle; a nonorthogonal finite product has magnitude `2 cos(pi/l)>=1`.
The singular and indefinite two-root cases have magnitude at least2.
This does not assume that an arbitrary face inherits Coxeter angles.

## Independent Schur and rank derivation

Let `s=m+1`. For the norm-two A_m Cartan block,
`(B^-1)11=(B^-1)mm=m/s` and `(B^-1)1m=1/s`; each isolated
leaf contributes inverse diagonal1/2. Subtracting the complete white
projection therefore gives

* long diagonal `(m+2)/s`, center diagonal `1/s`;
* every long/center entry `-1/s`;
* opposite long entry `-g-m/s`;
* opposite center entry `-h-(2m+1)/s`.

The last entry includes a total1 from the **two** isolated leaves. These
are exactly all entries asserted at lines81–97.

In normalized opposite-difference coordinates the eigenvalues are `2+g`
and `2+h`. The complementary opposite-sum block is

    T = [[2/s-g, -2/s], [-2/s, -2m/s-h]].

Expansion yields

    det(T) = gh+(2mg-2h-4)/s
           = [(m+1)gh+2mg-2h-4]/s.

For `g=0`, this is strictly negative for every `h>=0`. For `g>=1`, its
numerator is `h[(m+1)g-2]+2mg-4>=2m-4>0`. These exhaust the native gap;
no opposite-wall label or dotted distance is omitted. In the second branch
both diagonal entries of T are negative, so it would already give two
negative directions. In either branch the four-dimensional Schur block is
nonsingular.

The core has `(m+2)+4=r+2` roots, and its Gram consequently has rank
`r+2`. Z is a positive block of rank `n-r`, orthogonal to the entire core,
so the full Gram has rank `n+2`. Actual vectors in the original
`(n+1)`-dimensional Lorentz space cannot have this rank. This proves the
theorem, including the orthogonal opposite-long case `g=0`. There is no
assumed cusp radical, added energy at different vertices, or missing
realizability implication in this contradiction.

## Budget and dependency boundaries

The contextual calculations in lines138–168 agree with the pinned current
rule. D17 has long end0 and center14; the pair `(0,14)` retains weight1/2.
The ordinary component costs are D17=`227/2`, A17=`453/4`, and
D23=`8*23-32-1-3-2-1=145`. There are no applicable cross/context
additions in the forests displayed here. Hence `2D17+D23` costs372,
while `A17+D17+D23` costs371.75. The finite-simplex shape coefficient is
`[5*784-56*57]/2=364`; fixing inert charge258.5 leaves allowance105.5.
Deleting only the A17 endpoint half from the latter forest still leaves
deficit7.25. These are formal budget profiles, not claimed polytopes.

The distinct earlier node3 guard debits pairs `(3,13),(3,14)` with its
own retained/new coupons; it does not provide another saving for the
new `(0,14)` obstruction. Context662 finances selected actual face
recipients with once-only incidence accounting; it does not bound every
other finite vertex. Thus the source correctly claims **no improvement
of372, no universal finite-vertex allowance, no whole-coupon deletion and
no dimension56 consequence**. Split, ideal, role-changing, and
spectator-attached corners remain outside the homogeneous theorem.

All four direct source bindings were streamed and matched:

| Path | Bytes | SHA-256 |
|---|---:|---|
| `manuscripts/finite_volume57.tex` |156120|`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`|
| `results/finite57_complete_context662_source_composition_20261021.md` |7476|`7e53dc180baf41bce6d5cab815426a3391ca2de40fc2aa0849c4a6b58226664f`|
| `results/finite57_current_a17_endpoint_producer_scope_audit_20261021.md` |12283|`c9e2107004641fd3bcab1fc53e1b797aa264b38cf110e6e83f7bc0e65f8f1d71`|
| `results/fresh_d17_node3_large_d_guard_transfer_20261018.md` |12994|`eac18fe44e8e1d76175adfbf7446a73ff16c3f104bf92d7f4170edec4138ffa4`|

The last three are contextual scope inputs, not substitutes for the Schur
proof. This PASS is a complete written-source review of the restricted
theorem, not runtime verification, global financing, a new example, or a
change to any of the six frozen papers.
