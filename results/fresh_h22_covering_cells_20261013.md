# Finite covering tests and the missing bridges

There is a finite exact test for the proposed H21-to-H22 ball/horoball covering. It requires the vertices of a subdivision defined by the cap heights, not merely the original chamber vertices. A cheaper necessary test already detects **19 disconnected finite-contact components** in the 325 upper caps of the frozen 860-facet branch. A completion retaining those caps must connect the components using new finite-angle contacts.

The accompanying small example also disproves the tempting shortcut: **covering every original finite vertex and bounding every inherited cusp link does not imply finite volume**, even when every pair of walls satisfies the Coxeter angle condition. This is an exact H3 counterexample, with all its ideal vertices checked, rather than a numerical objection.

[Standalone checker](../scripts/fresh_h22_covering_cells_20261013.py) · [exact output](fresh_h22_covering_cells_20261013.json) · [previous all-point covering criterion](fresh_h22_incidence_mixed_height_20261013_covering.md) · [one-side reduction](fresh_h22_extension_feasibility_20261013.md) · [860-facet construction](fresh_h22_multiquadratic_20261013.md).

The [source/proof and fresh replay review](fresh_h22_covering_cells_20261013_review.md) passes and independently reproduces the 860-branch contact graph.

## 1. Which covering regions can meet?

Use the physical form `Q ⊕ <1>`. A positive-height cap is `r=(u,a)`, with `a>0`, positive norm, and nonzero future-causal old projection u. It defines

`B_r = { x in H^n : H_r(x)=−(u,x)/a ≤ 1 }`.

For future-timelike u, this is a closed ball. After normalizing r to norm one, write

`r=(c/sinh R, coth R)`, where `c²=−1`.

For two such caps, with center distance d,

`(r_i,r_j) = [−cosh d + cosh R_i cosh R_j] / [sinh R_i sinh R_j]`.

Consequently a finite angle π/m gives

`cosh d = cosh R_i cosh R_j + cos(π/m) sinh R_i sinh R_j`.

This is strictly smaller than `cosh(R_i+R_j)`: the balls overlap. Acute Coxeter angles therefore do **not** force a packing by disjoint balls. Dotted pairs have disjoint balls, and parallel pairs have a single tangency point. Horoballs satisfy the same last two statements.

Here is a proof that includes horoballs without a limiting argument. Normalize both cap norms to one, and suppose a unit old point x belongs to both covering regions. Write `u_i=a_i H_i x+v_i`, with `v_i ⟂ x`. Then

`v_i²=1−a_i²(1−H_i²)≤1`, and

`(r_i,r_j)=(v_i,v_j)+a_i a_j(1−H_i H_j)≥−1`.

All H values are nonnegative. Thus a dotted pair cannot meet. Equality −1 forces both H values to equal one and `v_i=−v_j`, with unit length. The intersection is the unique old point

`x=(u_i+u_j)/(a_i+a_j)`.

These statements concern the covering regions directly. They do not infer actual hyperbolic ridges merely from a positive-definite pair Gram matrix.

## 2. A contact-graph obstruction

Form a graph whose vertices are the upper caps, joining i and j when their normalized Gram entry is strictly greater than −1. In an acute Coxeter family these are exactly the finite-angle pairs, including right angles. Parallel pairs are not edges.

Assume old dimension n≥2, a finite cap family, and that every cap covers a nonempty open part of P. **Full coverage of P requires this graph to be connected.**

Indeed, regions in distinct graph components are either disjoint or touch at a single point. Remove the finitely many possible cross-component tangencies from `P°`. An open convex n-dimensional region minus finitely many points remains connected for n≥2. If the caps covered it, the finitely many component unions would give a partition into disjoint nonempty relatively closed sets, hence a disconnection. Each component remains nonempty because it contains a relatively open covering region.

The nonempty-region premise is automatic for the present retained-old-wall construction. A timelike old projection lies in the old chamber cone, and its center satisfies `H<1`; perturbing into `P°` preserves that inequality. For a null old projection δ, points `h+tδ`, normalized to square −1, have H tending to zero as t grows, where h is any strict old timelike center. The pure vertical exception u=0 is excluded here: it has a positive pairing with a retained cap of its own height sign. Both signs are already present in the 860 branch.

This filter is necessary, not sufficient. A connected graph can still leave large gaps. It is nevertheless useful before any cell enumeration: graph components specify which retained families need finite-angle bridges. One new cap may join several components, so the component count is not a lower bound on the number of additional caps.

### The frozen 860-facet branch

The checker independently reconstructs the 301 rational parent null projections and their same-height Gram matrix. Their 45,150 unordered pairs consist of 23,520 entries −1, 20,370 entries 0, and 1,260 entries −2, at norm two. It then uses the hash-bound 860 certificate's finite-contact identities and checks the rigorous saved intervals for every omitted new contact.

There are **325 upper caps and 43,900 finite-contact edges**. Their components are:

| Component | Upper caps |
|---|---|
| Large component, size 305 | All 301 parent null caps, plus recipe types 0, 2, 6, 22 |
| Triangle, size 3 | Recipe types 17, 23, 24 |
| Seventeen singleton components | Types 1, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21 |

Thus no finite acute-compatible completion retaining these facets can succeed by adding only caps dotted or parallel to every retained upper cap. Even a construction that repairs every inherited old cusp still needs finite-angle bridges between these components. This does not exclude changing the core or adding caps with such contacts.

The complementary [inherited-cusp footprint audit](fresh_h22_inherited_cusp_footprint_20261013_review.md) verifies that only 303 of the 1,095,990 old cusp rays currently have bounded links. Retaining this core forces at least 2,191,374 additional facets to address the others. That necessary boundary obligation and the contact-graph obstruction are separate: satisfying the former alone does not satisfy the latter or prove a cover.

The graph calculation inherits the primary certificate's exact new-cap normals, norms and finite-angle identities. It is not a second full reconstruction of those normals. Parent-parent arithmetic and all omitted graph edges are checked independently in this script, and its output binds all source hashes.

## 3. A finite exact covering criterion

Let `C(P)` be the old polyhedral chamber cone, contained in the future causal cone. This containment is available because the old finite-sided chamber has finite volume. For each upper cap set

`b_i=u_i/a_i`, `H_i(v)=−(b_i,v)`.

Define the closed minimum-height cell

`C_i={v in C(P): H_i(v)≤H_j(v) for every j}`.

It is a polyhedral cone: the additional inequalities are `(b_j−b_i,v)≤0`. On this cell the covering height is the one linear function H_i. These cells are closely related to weighted hyperbolic Voronoi diagrams; the affine nature of hyperbolic diagrams in the Klein model is treated by [Nielsen and Nock](https://arxiv.org/abs/0903.3287). The criterion below follows directly from Lorentzian convexity and does not require identifying these particular weights with a named diagram convention.

**Criterion.** The upper caps cover every old unit timelike point of P if and only if, for every extreme ray v of every nonempty cell C_i,

`H_i(v) ≤ sqrt(−v²)`.

For a timelike ray, check `H_i(v)²+v²≤0`; H_i is already nonnegative. For a null ray, the condition is exactly `H_i(v)=0`. A null ray must not be discarded as an unbounded direction of an affine chart.

To prove sufficiency, express any cell point as `v=Σ λ_k v_k` using its extreme rays and nonnegative coefficients. The reverse Lorentz triangle inequality gives

`H_i(v)=Σ λ_k H_i(v_k) ≤ Σ λ_k sqrt(−v_k²) ≤ sqrt(−v²)`.

The reverse triangle inequality itself follows by expanding the square and applying the reversed Cauchy inequality to pairs of future-causal vectors. Necessity at timelike rays is immediate. At a null ray, if the minimum height were positive, continuity would leave it positive on nearby timelike points while `sqrt(−v²)` tends to zero. That would violate full coverage. This also covers a cell consisting only of a null ray, without assuming that the cell contains timelike interior points.

Every cell is pointed, since it is contained in the future causal cone. Thus lineality cannot occur under these hypotheses. Empty cells can be omitted. Lower-dimensional cells are retained and treated in their actual linear spans; their extreme-ray decompositions obey the same proof. If a proposed implementation sees nonzero lineality, it has either omitted the future component/old inequalities or violated the hypotheses.

Combined with the previous all-point covering equivalence and one-side symmetrization, passing this finite test proves finite volume of the symmetric extension. Failure supplies a spacelike escape: at a failed timelike ray choose the added height strictly between its proper time and all cap heights. At a failed null ray, first perturb toward a strict old timelike point. A further small perturbation makes old boundary inequalities strict. This is an exact obstruction certificate, not merely a missed numerical sample.

### A practical certificate protocol

1. Freeze the old cone, cap normals, positive height signs, field embeddings and strict common center. Check or bind exact positive norms, acuteness, and the old causal-cone premise. Pure vertical caps require their separate treatment from the one-side reduction.
2. Run the finite-contact graph filter. For a disconnected graph, report the components and necessary bridge contacts; no cell search is needed to reject the current family.
3. For each nonempty C_i, enumerate all extreme rays exactly. A useful projective slice is `(h,v)=−1` for a fixed strict old future timelike h. Its intersection with the future causal cone is compact, so null rays appear as ordinary finite vertices in this slice. Degenerate cells must be included.
4. Certify completeness of that enumeration, not just feasibility of the saved rays. An exact double-description replay, or an independently verified equality between the H and V descriptions, suffices. A numerical list of candidate vertices does not.
5. Check the displayed quadratic inequality at every ray; require exact zero for null-ray height. Record a failed ray together with its active old facets and tying caps. Those constraints identify a concrete gap face for constructing a new cap.

At an ordinary cell vertex inside an old face of codimension d, the independent constraints comprise d old boundary equations and `n−d` cap-height equalities. Thus a generic interior cell vertex involves n+1 tying caps. Original P corners are only the case d=n. This explains precisely what a corner-only covering test omits.

The [complete 395-orbit old corner walk](fresh_h21_source_corner_walk_20261013.md) provides a reliable initial combinatorial skeleton. It does not certify the new cell vertices. Old symmetries reduce the subdivision only when they also preserve the chosen cap family. This finite criterion can still be expensive in H21; it is a proof protocol and a way to target gaps, not a claim that the current complete subdivision has been computed.

## 4. Exact counterexample to the corner-and-cusp shortcut

Start in H2 with form `diag(−1,1,1)` and norm-one old walls

`(0,1,0), (0,0,1), (1,−1,−1)`.

They bound a finite-volume triangle of angles `(π/2,0,0)`. Its finite vertex is `(1,0,0)`, and its ideal vertices are `δ1=(1,−1,0)` and `δ2=(1,0,−1)`.

Extend to `diag(−1,1,1,1)` and add six norm-one caps:

`(3,0,0,±sqrt(10)), (3,−3,0,±1), (3,0,−3,±1)`.

All nine walls share the strict timelike center `(3,−1,−1,0)`. The complete pair census is 13 right pairs, 4 parallel pairs and 19 dotted pairs. Exact projections of the center certify that every inequality is a true facet.

The old finite vertex is covered, and its two vertical endpoints `(1,0,0,±3/sqrt(10))` are actual finite vertices, each of square −1/10. At each old cusp, the two incident old walls and the two new caps form two orthogonal affine A1 blocks. Their positive kernels sum respectively to δ and 6δ; physical rank three and quotient rank two prove the H3 link is bounded. A complete enumeration of projective cone vertices finds exactly two finite, two ideal and eight spacelike vertices. Both ideal vertices are precisely these bounded cusps.

Nevertheless

`x=(3,−1,−1,11/4)`

is a strict future spacelike escape, with `x²=9/16`. Its three old pairings are −1. Its largest null-cap pairing is −13/4, and its largest finite-cap pairing is `−9+11sqrt(10)/4<0`, certified by `1210<1296`. Thus the chamber has infinite volume.

The three upper height functions have three polyhedral minimum cells. An exhaustive exact enumeration gives cell ray counts 5, 3 and 3, totaling seven distinct rays. The three original vertex rays pass coverage; all four newly introduced rays fail it. For example,

`v=(1,−1+1/sqrt(10),−1/sqrt(10))`

lies on the old third wall and ties the finite cap with the first cusp cap. Its square is `(1−sqrt(10))/5`, its minimum height squared is 9/10, and its covering gap is `(11−2sqrt(10))/10>0`.

This is a reproducible low-dimensional failure of the proposed shortcut. It is not a new finite-volume Coxeter example or a claim of literature priority.

## 5. The failure persists in every old dimension

There is also a general construction explaining why closing all old corners and cusps is insufficient. For every distinct future vertex ray v_i of a finite-volume Coxeter P, put `c_i=−v_i²≥0` and add signed norm-one caps

`r_i^±=(T v_i, ±sqrt(1+T² c_i))`.

All incident old-wall pairings are zero. Nonincident old-wall pairs become dotted for large T. For distinct rays, set `A_ij=−(v_i,v_j)>sqrt(c_i c_j)`. A sufficient strict dotted bound for same-height caps is

`T² > (c_i+c_j+2A_ij)/(A_ij²−c_i c_j)`.

It follows by squaring `sqrt((1+T²c_i)(1+T²c_j)) < T²A_ij−1`; the displayed bound also makes the right side positive. Opposite-height pairs are automatically dotted or parallel. Thus sufficiently large T yields a coherent Coxeter chamber.

Every old finite vertex is strictly covered by its own cap, since `T sqrt(c_i)/sqrt(1+T²c_i)<1`. Every old cusp gains an orthogonal vertical affine A1, so its previously bounded old link becomes bounded in the next dimension.

Choose a strict old timelike h not proportional to a vertex ray and write `c_h=−h²`, `A_i=−(v_i,h)`. Increasing T further so that

`T² > c_h/(A_i²−c_h c_i)`

for every i makes all cap heights at `h/sqrt(c_h)` exceed one. The central gap therefore remains. This argument applies to the old H21 chamber as well, without constructing a huge cap list.

The useful next search is consequently **coherent finite-angle bridge construction followed by a complete covering-cell test**, with the complete old cusp census supplying necessary boundary obligations. Neither more distant caps nor the cusp count alone establishes progress toward finite volume.

## Reproduction and scope

```sh
python scripts/fresh_h22_covering_cells_20261013.py --verify
```

The script uses only the Python standard library. It independently checks all nine-wall example identities, both bounded cusp links, the complete projective vertex census, every minimum-cell ray, and the strict escape. Its H22 graph is an exact computation relative to the explicitly bound 860 certificate, with no new volume or unrestricted nonexistence claim. The all-dimensional constructions and cell criterion are proved above; no enumeration of the H21 covering cells is claimed.
