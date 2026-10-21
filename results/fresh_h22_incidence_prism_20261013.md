# An exact H22 incidence chamber with π/7 caps

The construction has **835 actual facets**, includes **23 compatible new caps with π/7 contacts**, and has **23 certified finite vertices** involving those caps. It still has **infinite volume**: an exact rational vector is strictly feasible and spacelike. There is no new finite-volume example, dimension record, or literature-novelty claim.

This is a changed construction, not an extension retaining the preceding multi-million-facet cusp bank. Its useful feature is the departure from integral-lattice completion: any finite-volume completion retaining its rational old wall block and a π/7 contact would have to be **non-quasi-arithmetic**.

## Construction and exact checks

Start with the old 210 Borcherds walls in a rational space with form Q of signature (21,1). The 42 norm-two walls split into the points and lines of the archived incidence structure, with 21 walls in each part. Write P_i for the point normals, v′ for the archived null-bank center projection, and h for the strict timelike center, h²=−14.

For a binary vector u of weight w, use the old projection

`x=(A v′−Σ u_i P_i)/2`, with `x²=(w−A²)/2`.

The selected null bank has 280 profiles with (A,w)=(3,9) and 21 profiles with (A,w)=(4,16). Add both signed height-√2 walls for each profile. In coordinates where the last entry multiplies √2, the form is **M=diag(Q,2)** and these normals are `(x,±1)`. This gives 812 walls. The literal 301-projection bank is used; no full-symmetry or exhaustive-profile claim is needed for its validity.

For each new cap, 21 retained walls prescribe right angles and an anchor prescribes angle π/7. Let T=2cos(π/7), the root in (1.8,1.9) of

`T³−T²−2T+1=0`.

The [rational recipes](../data/enumerations/fresh_h22_incidence_prism_20261013_recipes.json) specify u,v, a sign ε, and q∈{2,4}, satisfying u⊥v, v²>0, u·anchor=1, v·anchor=0, and anchor²=q. Then

`r=−(qT/2)u + ε sqrt((q−q²T²u²/4)/v²) v`

has norm q and the required anchor angle. The 23 selected caps are pairwise strictly dotted. Every cap cuts its archived feasible positive-square parent ray. All finite contacts with parent walls have labels 2,4,7; every remaining cap–parent relation is strictly dotted.

The [primary checker](../scripts/fresh_h22_incidence_prism_20261013_verify.py) reconstructs the 812 parent normals, checks their complete pair table and Lorentzian signature, proves all finite-angle identities by polynomial reduction, and proves strict signs using outward-rounded integer intervals. A common strict timelike center and pairwise nonpositive products imply that every inequality is an actual facet: projecting the center orthogonally onto that wall gives a timelike point strict on all other walls. See the [primary certificate](fresh_h22_incidence_prism_20261013_verify.json).

The [independent review](fresh_h22_incidence_prism_20261013_independent.md) checks the formulas and interval gates, and separately reconstructs two representative caps with Fraction arithmetic and rational bisection. It verifies all 1,624 of their parent relations. This is independent arithmetic for two caps, not a second implementation of all 23.

## Finite vertices and remaining volume obstruction

The [geometry checker](../scripts/fresh_h22_incidence_prism_20261013_geometry.py) certifies 23 positive-definite 22-wall Gram matrices by interval elimination. Each consists of one new cap and 21 parent walls. Acuteness implies that the inverse of this positive-definite Gram matrix is entrywise nonnegative; orthogonal projection of the strict center therefore gives an actual finite vertex, strict on every other wall. Numerical prefilters only select candidates for proof: **23 is not an exhaustive vertex count**.

The [geometry source review](fresh_h22_incidence_prism_20261013_geometry_review.md) checks this projection argument and the scope of the candidate prefilters. The same [geometry certificate](fresh_h22_incidence_prism_20261013_geometry.json) records a rational strict escape with

`x²=3956990142533/12500000000000 >0`, `(x,h)<0`, and `(x,r)<0` for every one of the 835 walls.

Joining x to the strict timelike center reaches a null point still strict on every wall, yielding an open set of ideal directions and infinite hyperbolic volume. All comparisons are exact rational or outward-rounded algebraic sign checks. Adding walls that cut this one vector would require a fresh volume check.

## Conditional obstruction to quasi-arithmetic completion

The retained old 42 norm-two walls have a rational Gram matrix G₀ with row sums −3. Thus `1ᵀG₀1=−126`; its unit-normalized version remains rational and indefinite. A π/7 pair contributes the length-two cyclic product `4cos²(π/7)=T²` of twice the normalized Gram matrix. This generates a cubic field, since its polynomial

`X³−5X²+6X−1`

has no rational root. These assertions have a [reproducible checker](../scripts/fresh_h22_incidence_prism_20261013_arithmetic_obstruction.py) and [exact artifact](fresh_h22_incidence_prism_20261013_arithmetic_obstruction.json).

Here is the conditional implication. If a finite-volume completion were quasi-arithmetic, its full Gram field would be totally real. An embedding moving T² extends to that field and is nonidentity on the cyclic field. Vinberg's criterion would require the conjugate Gram matrix to be positive semidefinite. Its rational principal block G₀/2 is unchanged and indefinite, a contradiction. If total reality fails, that already violates the criterion. This uses the criterion stated in [Bogachev–Yorov, Theorem 2.2](https://nvbogachev.github.io/Quasi_arithmetic_Coxeter_prisms-2024.pdf).

This excludes **quasi-arithmetic finite-volume completions retaining these facets**, not arbitrary completions. It neither proves that a nonarithmetic completion exists nor rules one out.

## Replay

All four commands have passed on the frozen inputs:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h22_incidence_prism_20261013_verify.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h22_incidence_prism_20261013_geometry.py --verify
python3 scripts/fresh_h22_incidence_prism_20261013_independent.py --verify
python3 scripts/fresh_h22_incidence_prism_20261013_arithmetic_obstruction.py --verify
```

The [subsequent exact edge frontier](fresh_h22_pi7_edgewalk_20261013.md) classifies all 506 directed releases from these 23 vertices: 497 unclosed rays, seven finite endpoint directions, and two actual ideal endpoints with unbounded links. It has a [source review](fresh_h22_pi7_edgewalk_20261013_review.md). The next construction task is to solve those missing endpoint and link obligations, allowing new metric parameters and field extensions. The current pair checks and local finite vertices do not supply global closure.

The [837-facet extension](fresh_h22_pi7_cusp_pair_20261013.md) subsequently bounds both of those ideal links by reflected parallel mates. An independent [missed-cusp certificate](fresh_h22_incidence_prism_20261013_uncapped.json) identifies another inherited cusp untouched by all 23 π/7 caps and their reflected mates. It gives a simpler strict escape of square 1,586. The [covering reformulation](fresh_h22_incidence_mixed_height_20261013_covering.md) and [review](fresh_h22_incidence_covering_20261013_review.md) explain why finitely many further timelike-projection caps cannot repair that particular cusp.
