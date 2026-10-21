# Coordinate facets propagate the literal golden dyadic obstruction

20 October 2026. Paper-only conditional implication; no computation.
The dimension-nine base period remains subject to its independent reviews
and fixed exact check. Nothing here promotes that base result.

Let `O=Z[phi]`, `phi^2=phi+1`, and

    R_n=<-phi>+I_n+D4(O), n>=5.

Import the UNIFORM primitive norm pool and full zero-height germ from
`results/golden_dyadic_d4_extension_first_edge_20261020.md`, SHA256
`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`.
That unchanged source proves, in every stated rank, that a primitive
reflecting normal can be unit-rescaled to norm `q=1,2,4`; norm4 normals
pair the whole lattice into `2O`; and the full germ at `(1;0)` is
`B_n+F4`. In particular the final coordinate wall `e_n^perp` of its
chosen B_n chamber is an ACTUAL facet of the full chamber `C_n`.

## Exact full-facet arrangement

Set `H=e_n^perp`, so the exact orthogonal lattice is `R_(n-1)`.
Consider an ambient mirror, distinct from H, crossing H in hyperbolic
space. Write its normalized integral primitive normal as
`r=x e_n+r_0`, with `r_0` integral in `R_(n-1)` and `r^2=q`.
The crossing gives `x^2<q` at the physical embedding. Positive
definiteness at the other real embedding gives `x'^2<=q`.

For `x=A+B phi`,

    Tr(x^2)=2A^2+2AB+3B^2
           =2(A+B/2)^2+(5/2)B^2.                 (1)

This supplies the COMPLETE coordinate table, without a native-angle
assumption:

|normal norm|necessary coordinate strip|all surviving x|
|---|---|---|
|1|`x^2<1, x'^2<=1`|0|
|2|`x^2<2, x'^2<=2`|0, +1, -1|
|4|`x in2O, (x/2)^2<1, (x'/2)^2<=1`|0|

For the first and last rows, (1) is an integer strictly below2, while
its smallest nonzero value is2. In the middle row it is strictly
below4. Thus B is0 or+/-1: B0 gives x0,+/-1, and B+/-1 gives only
`+/-phi,+/-phi^-1`. The former violate the physical bound and the
latter violate the conjugate bound, since `phi^2>2`. This proves
the displayed list exhaustively. The norm4 divisibility used in its
row is exactly `(r,e_n)=x in2O` from the imported pool.

Every induced mirror is therefore a full `R_(n-1)` mirror:

* If x0, r is an integral primitive reflecting normal in the exact
  perpendicular lattice; its ambient reflection restricts to it.
* If x=+/-1, necessarily q2. The projected normal `r_0=r-xe_n`
  is integral and has unit norm1, hence is primitive and reflective
  on `R_(n-1)`. Its hyperplane is exactly the intrinsic intersection
  of the original mirror with H.

Conversely EVERY reflecting normal of `R_(n-1)` extends to R_n by
zero e_n-coordinate. The reflection preserves the orthogonal summand
`O e_n`, and primitivity is preserved by this direct sum. Thus all
intrinsic mirrors occur, not just a selected finite prefix. We have
proved the equality of FULL arrangements on H.

Since H is an actual coordinate facet, its intrinsic relative
interior is a chamber of this full induced arrangement. Mirrors
disjoint from H impose no extra boundary within hyperbolic H.
Consequently the actual facet of `C_n` is exactly a full chamber of
`R_(n-1)`, with the inward B_(n-1)+F4 germ fixing its chamber choice.

## Conditional uniform conclusion and its boundary

Assume the full literal `R_5` chamber is not compact. If some C_n,
n>=5, were compact, successively taking these actual compact
coordinate facets would give a compact full R_5 chamber, a
contradiction. Hence the base noncompactness implies noncompactness
of the full literal chambers R_n for EVERY n>=5 (hyperbolic
dimensions n+4>=9).

For the stronger NONCOFINITE conclusion the external arithmetic
input is explicit: the integral orthogonal group of each R_n is a
uniform arithmetic lattice. The form is admissible, since phi'<0
and its other real form is positive definite; anisotropy and the
standard Borel--Harish-Chandra compactness theorem apply. A cofinite
reflection subgroup would have finite index in this uniform lattice,
and hence a compact reflection chamber. Thus the same base
noncompactness implies noncofiniteness of every literal R_n, n>=5.
This last arithmetic implication does not follow from facet geometry
alone.

The prospective base source is
`results/golden_dyadic_r5_chamber_period_candidate_20261020.md`, SHA256
`e2aeff93b42f87673e55399804df28a9eee3713ac2e7cf064c6f885e5bfce671`.
Its hypothesis is deliberately left conditional here until its
review/check status is settled. The propagation concerns ONLY these
displayed integral lattices and their FULL arrangements. It neither
excludes other lattices in their K-similarity/commensurability classes
nor bounds arbitrary compact Coxeter polytopes.
