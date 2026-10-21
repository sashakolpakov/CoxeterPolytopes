# Four mixed-face normals and two infinite-volume fourteen-wall intersections

The exact search finds four new norm-2 normals beyond the frozen spacelike-generated cone classification: two for each same-kind raw23 union, and none for the mixed-kind union in the finite-anchor domain tested here. Each same-kind pair can coexist. Both pairs cut all 32 original spacelike extreme rays and preserve all 18 original finite vertices, but their complete fourteen-wall cones have **72 new spacelike rays each**. Explicit strict escape certificates prove that both resulting hyperbolic 9-dimensional intersections have infinite volume. They are not compact or finite-volume examples.

This is a new geometric domain and a constructive result for four normals. It does not exclude arbitrary completion, classify all mixed directions, install a capacity or face-mass inequality, change the finite arithmetic matrix, or invoke an optimizer.

The producer is [fresh_raw23_lorentz_two_face_native_roots_20261018.py](../scripts/fresh_raw23_lorentz_two_face_native_roots_20261018.py); the certificate is [the saved JSON](fresh_raw23_lorentz_two_face_native_roots_20261018.json). The [standalone saved checker](../scripts/check_raw23_lorentz_two_face_roots_saved_20261018.py) independently checks the new finite-anchor domain, roots, changed cones and escape arithmetic. Frozen geometry is inherited from the [original union/polar packet](fresh_raw23_union_polar_completion_preflight_20261018.md), and the previous [spacelike pole classification](fresh_raw23_spacelike_pole_classification_20261018.md) remains unchanged.

## Complete original two-face inventory

The known twelve-wall Lorentz space has dimension ten, signature `(9,1)`, and 50 extreme rays in each of the three frozen cases `(1,1)`, `(1,6)`, `(6,6)`. The producer inspects all 1,225 unordered ray pairs per case. It gathers their full common active-wall set, checks its exact rank eight, and then gathers **every** frozen ray on those walls. Every accepted face has precisely its two displayed generators. Thus no nonsimplicial face is omitted: a pointed two-dimensional polyhedral cone has two extreme rays.

Each case has the same complete inventory:

| Generators | Faces | Restricted Gram inertia |
| --- | ---: | --- |
| Two spacelike rays | 112 | `(2,0,0)` |
| One timelike and one spacelike ray | 64 | `(1,1,0)` |
| Two timelike rays | 49 | `(1,1,0)` |

There are 675 two-faces in total, including 192 mixed faces and 147 timelike-pair faces. This computation does not reuse a positive-definite bound on a Lorentzian span.

All timelike generators have negative product with the same strict future-timelike center. A positive combination of two future-timelike generators is future timelike. Such a projected normal has negative product with every future causal point of the known hyperbolic section; a future-null projection is nonpositive there. A positive passive part can still make the **full** normal have norm two, so this is not a nonexistence claim for full normals. It says these causal projections cannot cut the known section or remove its open ideal patches. The finite search below consequently treats a nonpositive projected norm as a noncutting case, not as an impossible norm-2 normal.

## Finite-anchor reduction on the mixed faces

Write a projected normal in a mixed face's relative interior as

\[
x=aT+bS,\qquad a,b>0,
\]

where the frozen generators are timelike `T` and spacelike `S`. Each mixed face has exactly four known-wall contact forms that are not identically zero. Their coefficients `p_j=-<T,n_j>`, `q_j=-<S,n_j>` are nonnegative. Consequently every such contact is strictly negative on the relative interior. Exactly five of their six unordered pairs are linearly independent.

This packet covers the necessary candidate domain in which **two linearly independent nonzero contacts are finite Coxeter angles** in the fixed arithmetic model over `Q(sqrt(5))`. For norm-2 roots the allowed nonzero finite squared products in that model are

\[
1,\quad2,\quad(3+\sqrt5)/2,\quad3,\quad(5+\sqrt5)/2,
\]

corresponding to labels `3,4,5,6,10`. This finite native catalogue is an inherited hypothesis of the compact-arithmetic model; the new search does not extend it to arbitrary ground fields.

If the known scaled root has squared scale `s_j`, an anchor fixes the positive raw contact to `sqrt(s_j*f_j)`. Two independent anchors therefore determine `a,b` uniquely by an invertible two-by-two system. Every one of the 25 ordered square choices is tested for every independent pair: `192*5*25=24,000` exact tuples. The multiquadratic square roots are ordered positive; no numerical grid or approximation occurs.

The tests are:

1. Both coefficients must be positive.
2. The Lorentz projection must be spacelike and have norm at most two. A full passive component lies in the orthogonal positive-definite complement, so it can only increase the norm. No passive contact is assumed zero.
3. Every known-wall physical two-cycle must belong to the ground field. Every finite contact must be in the native catalogue; every two-cycle must be integral and its conjugate lie in `[0,4]`.

The last gate is a necessary condition of the stated arithmetic model. No conjugate-projection-norm condition is used to exclude a candidate. All four survivors happen to have coefficients in the base field and projected norm exactly two, which forces the passive part to vanish. The causal mixed candidates are future causal: along `T+tS`, the quadratic norm starts negative and has positive leading coefficient, hence exactly one positive zero. Its causal positive-parameter portion is connected to the future-timelike `T`; it cannot change time orientation there.

| Exact outcome | Tuples |
| --- | ---: |
| Outside the positive relative interior | 7,160 |
| Causal projection, hence noncutting on the known section | 2,210 |
| Projected norm greater than two | 11,156 |
| A physical two-cycle outside the ground field | 2,834 |
| A native, integrality or conjugate pair gate fails | 636 |
| Surviving norm-2 roots | 4 |

These counts sum to 24,000. There are no duplicate surviving roots. The mixed-kind case has no survivor in this domain. Directions with no finite nonright contact, only one such contact, or only proportional finite contact forms remain outside the finite-anchor classification. In particular, no statement that all mixed-face normals have been exhausted is made.

A useful necessary consequence is that every spacelike cutting candidate on these mixed faces must have at least two parallel or dotted known-wall contacts, under the same native arithmetic assumptions. If three of the four nonzero forms were finite, some pair would be independent; the complete search would force one of the displayed roots, each of which has exactly two finite contacts. Cases with at most two finite contacts already have at least two nonfinite ones.

## Four exact roots

Ray indices refer to the frozen original cones. Coefficients in this table multiply the timelike and spacelike generators in that order. The listed simple contacts have physical product `-1`; the listed dotted contacts have physical product `-2-sqrt(5)` and square `9+4sqrt(5)`. All remaining old-wall products are zero.

| Case | `(T,S)` rays | `a` | `b` | Simple walls | Dotted walls |
| --- | --- | --- | --- | --- | --- |
| `(1,1)` | `(28,38)` | `(3+sqrt(5))/4` | `(5+3sqrt(5))/4` | `3,6` | `0,10` |
| `(1,1)` | `(40,42)` | `(3+sqrt(5))/4` | `(5+3sqrt(5))/4` | `2,7` | `0,11` |
| `(6,6)` | `(3,13)` | `(-1+sqrt(5))/4` | `1` | `3,6` | `9,10` |
| `(6,6)` | `(15,17)` | `(-1+sqrt(5))/4` | `(5-sqrt(5))/4` | `2,7` | `9,11` |

The certificate stores their full ten coordinates and all twelve raw products. Direct evaluation proves norm two, nonpositive products with every old wall, strictly negative products with all eighteen original finite vertices, and positive products with each of the thirty-two original spacelike rays. Thus they preserve those finite vertices and cut every original escape ray.

Within either same-kind case, the new normals have product `-1-sqrt(5)` with each other, so their hyperplanes are ultraparallel and their outward halfspaces are compatible. Each new root instead has **positive product `+1`** with that case's old normalized pure pole. The two-root completion is an alternative to the old pole; appending all three would fail the acute-wall condition.

## Complete changed cones and arithmetic checks

For each same-kind case, the twelve old walls plus both new roots define a new fourteen-wall cone. The producer enumerates every one of its `choose(14,9)=2,002` active nine-wall subsets, computes exact kernels and orientations, and deduplicates rays by the full fourteen products. The standalone checker repeats this with reverse-column elimination. Both cones have:

- 106 extreme rays: 34 future-timelike and 72 spacelike, with no null extreme ray;
- a strict future-timelike interior point, the sum of the finite rays;
- fourteen actual facets, each with a saved timelike point having zero product with precisely that wall and strictly negative products with the other thirteen;
- physical scaled Gram inertia `(9,1,4)` and conjugate inertia `(10,0,4)`;
- all 91 pair gates valid and all 93 unoriented simple cycles integral.

For the facet certificate, the sum of finite rays on a new wall can lie on other walls as well. The code therefore adds a sufficiently small positive multiple of the sum of **all** rays on that wall, retaining negative norm and making the other thirteen inequalities strict. Merely checking rank nine of its zero-ray span would not establish an actual hyperbolic facet.

The scaled Gram is symmetric with diagonal `2s_i`; the two new scales are one. For a simple cycle the physical cyclic product is exactly the product of its raw entries divided by the product of the vertex scales. Integrality in `Q(sqrt(5))` is checked by integral trace and norm. The checker enumerates the full simple-cycle domain independently. These are exact native and arithmetic checks on the fourteen walls, not a finite-volume arithmetic-lattice claim.

## Exact open ideal escapes

The changed cones have new spacelike rays even though the two normals cut all the old ones. Let `c` be the saved strict timelike center and `r` the new ray with index six. The certificate gives:

| Case | `<r,r>` | `epsilon` in `y=r+epsilon*c` | `<y,y>` |
| --- | --- | --- | --- |
| `(1,1)` | `(117-sqrt(5))/55` | `1/128` | `253259/180224-1150623sqrt(5)/4505600` |
| `(6,6)` | `(7-sqrt(5))/5` | `1/256` | `3636393/3604480-5724387sqrt(5)/14417920` |

Both final norms are strictly positive by exact ordered-field signs. All fourteen products of `y` are strictly negative, and the full coordinates/products of `r`, `y`, and `c` are saved. The quadratic norm polynomial of the segment `c+t(y-c)` is negative at zero and positive at one. At its first zero, the segment therefore reaches the future light cone at a nonzero point while **every wall inequality is still strict**. An open neighborhood of that ideal point also satisfies all inequalities. The intersection consequently has an open ideal patch and infinite hyperbolic volume.

This escape proof is field-independent and allows any additional positive passive coordinates to remain zero. It excludes finite volume for these particular fourteen-wall intersections. It does not rule out adding further compatible walls, or completing the original union in another way. The newly created 72 spacelike rays are the concrete next completion target.

## Verification and limits

Generation `75854/7618cd` and separate producer replay `75869/756f71` exited zero with the stated complete domains. The standalone saved check `36223/345bed` also exited zero and passed all 24,000 tuples and 4,004 changed active sets. The final helper-only guard explicitly checks that each facet interior point is on the center's future sheet; its separate final run `15891/02ef40` exited zero and passed.

Construction's independent full producer/report/helper mathematical read passed, including the immediate inherited cone/kernel/native-square interfaces. It checked complete face and anchor domains, passive-norm scope, changed cones, arithmetic gates and strict infinite-volume escapes by reading. It did not execute tuples, native arithmetic, cones or historical source incidence. Root's full review remains pending; no such review is claimed here.

The checker imports no producer. It uses the frozen standalone polar helper's exact `Q(sqrt(5))` and reverse-column cone implementation, plus a separately coded, reversed-radicand multiquadratic construction for all 24,000 tuples. That construction shares the mathematical ordered-extension design, not a claim of unrelated arithmetic. It rechecks complete pair/anchor domains, all candidate decisions, full roots, all changed active sets, pair/cycle arithmetic, facet points and strict escapes. The old fifty-ray catalogues and actual raw23 source incidence, and the ground-field/native assumptions, remain explicit frozen inputs.

No old unknown compatibility, owner, face mass, charge, multi-owner cost or matrix entry is changed. There is no new compact example, target-12 certificate, capacity theorem or global arithmetic conclusion.
