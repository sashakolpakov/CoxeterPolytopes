# Complete all-retained single-factor cusp sections

An interior Coxeter H20 section through any of the 18 archived cusp representatives cannot be active in one affine factor while retaining **all of that factor's original inequalities as distinct actual facets of the Euclidean cusp link**. This closes the all-retained case, with arbitrary real section normals, arbitrary null shifts, and every Coxeter angle label. The case where one original factor facet is omitted remains open in higher rank.

The [exact checker](../scripts/fresh_h20_single_factor_classification_20261014.py) and [certificate](fresh_h20_single_factor_classification_20261014.json) give a complete finite reduction over all 53 affine factors. There are 287 possible nontrivial rank-one cuts. Of these, 119 fail a local angle test; the remaining 168 give bounded Euclidean Coxeter product links. Every one of those 168 normal planes has an actual obtuse global ridge for every interior shift. Seventy-seven coincide with previously checked all-retained leaf-sum planes; 91 are additional directions. None coincides with a previous twin-leaf difference plane.

The classification and all continuous rejections use standard-library rational arithmetic. Twenty-eight rational midpoint witnesses were found by numerical linear programming, then reconstructed and checked exactly; the checker does not import an optimizer or trust numerical coordinates. No new Coxeter polytope is established.

## Why every all-retained direction appears in the census

Work in the positive-definite quotient at the chosen cusp. Let the active irreducible affine factor have rank `n`, original roots `r_0,...,r_n`, Gram matrix `G`, and positive affine relation

`D = sum m_i r_i`, with `m_i>0`, `D²=0`.

Its Euclidean simplex is bounded. An interior hyperplane section has dimension `n-1`; under the hypothesis it has `n+1` distinct facets. A bounded Euclidean Coxeter chamber is a product of affine simplexes. The difference between its number of facets and its dimension is the number of factors, so this section has exactly two factors. Their facet sets give a partition `S,T`, with both sizes at least two. The rank-one and rank-two original factors consequently have no all-retained case.

Normalize the spatial section normal to square one and write `a_i=(e,r_i)`. The projected Gram is

`H = G-aa^T`.

The two product factors are orthogonal, so `G_ST=a_S a_T^T`. This cross matrix is nonzero because the original affine graph is connected. It has rank one. Since every nonzero old off-diagonal entry is negative, its nonzero entries form a complete bipartite graph, and all nonzero entries of `a_S` have one sign while those of `a_T` have the opposite sign.

The checker proves directly that none of the 53 literal affine graphs has a four-cycle: no two nodes have two common neighbors. A complete bipartite graph without a four-cycle has one side of size one. Thus every possible crossing graph is a star. The exhaustive enumeration chooses each possible center and each nonempty subset of its neighbors, marks exactly those edges as crossing, and propagates a two-coloring through the connected graph. Inconsistent choices are rejected; consistent partitions are identified up to swapping the colors. It retains only partitions with at least two nodes on each side. This enumerates every rank-one cut, not just selected Dynkin shapes or coefficient patterns.

The cut also determines the spatial normal uniquely. Set

`v = sum_{i in S} m_i r_i`, and `c = sum_{i in S} m_i a_i`.

Since `m^T a=0`, the cross-block identity gives `(v,r_i)=c a_i` on **both** sides of the cut. Here `c` is nonzero because the nonzero entries of `a_S` have one sign. Moreover `v²=c²`. The original roots span the Euclidean factor, so its unit normal is proportional to `v`, modulo the null direction in the original Lorentz space. Consequently every real normal in this direction has the form, up to nonzero overall scale,

`e(t)=v+tD`, with `q=e(t)²=v²>0`.

No additional continuous spatial parameter remains.

## Exact local angle and facet-survival gates

For each of the 287 cuts the checker reconstructs `v` and the entire projected Gram

`H_ij = G_ij - (v,r_i)(v,r_j)/q`.

All entries are rational. It verifies positive diagonal entries, zero cross-block entries, nonpositive off-diagonal entries, a positive kernel on each block, and a positive-definite principal minor of the required size. Thus no projected inequality is discarded or mistaken for a zero normal. The two blocks are genuine affine simplexes of total rank `n-1`.

For rational squared cosine, the complete Euclidean Coxeter list is

`0, 1/4, 1/2, 3/4, 1`.

Indeed, `4 cos²(pi/m)=2+2 cos(2pi/m)` is a rational algebraic integer if the squared cosine is rational. Its range leaves exactly these values, including the parallel limit. This excludes every other integer label without an arbitrary angle cutoff. The checker records a forbidden rational squared cosine for each of the 119 rejected cuts and checks every pair in the 168 survivors.

The two positive affine relations after projection are

`sum_{S} m_i r_i' = -tD`,

`sum_{T} m_i r_i' = (1+t)D`.

Therefore the bounded product link exists exactly for `-1<t<0`. In this interval both block constants have the required strict sign on the future horosphere. Outside it, the original normal is a one-sign combination of old roots, so the hyperplane misses the strict interior. All facets survive throughout the interval; it is not necessary to assume this from a raw Gram calculation alone.

Distinctness is checked using the unique null line of `span(v,D)` and the proportionality class of the normal's pairings with all incident roots. The 168 planes are distinct, 77 match the previous leaf-sum census, and none matches an archived difference plane. The two earlier C2 sum planes with a redundant local inequality correctly do not appear here.

## Actual global ridges for all 168 planes

The projection calculations alone establish local links, not global Coxeter sections. The following certificates supply actual global ridges and rule out every surviving plane.

An ambient face point is constructed by projecting the strict source center `h` onto a positive-definite old wall set. The checker verifies its timelike norm, future orientation, all 210 inequalities, and its exact zero set. In particular, a positive combination of two such points is timelike, and its zero walls are precisely the intersection of the two zero sets.

Choose either side of the cut as the active set `A`, call the other side `B`, and use `v_A=sum_A m_i r_i`. The new parameter is still denoted `t`; swapping the sides replaces the original parameter by `-1-t`.

For all 168 planes the certificate finds an old root `i` in `A` and an auxiliary old root `j` outside the factor, such that

`u=(v_A,r_i)>0`, `(i,j)=0`,

`(v_A,r_j)=(D,r_j)=-b<0`,

and `j` is orthogonal to all roots in `B`. The projected pairing is

`(i',j') = ub(1+t)/q > 0`.

This becomes an obstruction only after the actual-ridge argument below.

For 154 planes, both ambient wall sets `B union {i,j}` and `A union {j}` are positive definite. Normalize their exact face points `L,R` to `(D,L)=(D,R)=-1`. They have

`(v_A,L)=-1`, `(v_A,R)=0`,

and common zero set exactly `{i,j}`. The point

`x(t)=(-t)L+(1+t)R`

is future timelike, lies in the section, and has exactly `i,j` zero for every `-1<t<0`. The checker separately verifies a strict lower bound for the projected two-root determinant throughout the whole interval. These are therefore actual obtuse ridges, with both inequalities actual facets.

The fourteen remaining planes occur at cusps 12, 13, and 14. A single pair of full-interval face points is unavailable, but two overlapping certificates cover each plane. For each orientation of the cut, the checker verifies:

- the finite face point `L` above at `t=-1`;
- an exact rational future timelike point `M` at `t=-1/2`, strict against every wall except `i,j`;
- the same strictly positive projected pairing and determinant bound for `-1<t<=-1/2`.

After normalizing `(D,M)=-1`, the interpolation is

`x(t)=(-2t-1)L+2(1+t)M`.

Its coefficient of `M` is positive throughout that half interval, so every nonridge inequality is strict, including at `t=-1/2`. Applying the other orientation covers `-1/2<=t<0` in the original parameter. Thus there is no uncovered midpoint or endpoint limit. The [28 rational midpoint recipes](fresh_h20_single_factor_classification_20261014_recipes.json) are all checked directly against the original 210 walls; the proof does not depend on the linear program that suggested them. Exact positive-definiteness tests reject singular faces that a preliminary floating-point factorization had accepted; none of those faces is used as a finite face in this certificate.

For either interpolation, the determinant bound is elementary. With `N_i=G_ii-u²/q`, the determinant at `s=1+t` is

`N_i G_jj - b² G_ii s²/q`.

The checker proves it remains strictly positive by substituting the largest `s` in the declared interval: `1` for a full interval, `1/2` for a half interval. This also proves independence of the two projected normals. The strict zero-set checks then establish a genuine codimension-two ridge, not merely two intersecting supporting hyperplanes.

## Verification and remaining constructive scope

```sh
python scripts/fresh_h20_single_factor_classification_20261014.py --verify
```

The certificate exhausts the graph cuts and rational local angle tests, then certifies the actual global obstructions. The reduction from arbitrary real normals to graph cuts is the mathematical argument above, not a claim that a finite coefficient grid was exhaustive.

A [separate SymPy audit](fresh_h20_single_factor_classification_20261014_root_review.md)
rechecks all168 local acceptances,119 local rejections,182 continuous
ridge proofs and353 distinct endpoint points without producer helpers.
A second [independent midpoint check](fresh_h20_single_factor_classification_20261014_midpoint_review.md)
focuses on the28 reconstructed rational points in the fourteen delicate
two-half proofs.

The result assumes that every original facet of the active affine factor survives as a distinct actual **local cusp-link facet**. It does not exclude a factor slice omitting one inequality, including cases where two projected inequalities coincide and only one is retained. It also does not classify sections avoiding these cusp representatives, changed ambient chambers, or subdivisions of the resulting non-Coxeter sections. Those remain separate constructive routes.
