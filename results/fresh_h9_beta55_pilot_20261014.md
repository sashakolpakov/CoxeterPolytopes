# A different discriminant55 host: a complete small Vinberg prefix

The new lattice

`L55 = E8 ⊥ [[−2φ,−1],[−1,4]]`, with `φ=(1+√5)/2`,

has a complete height20 Vinberg prefix with twelve actual walls,23 finite vertices and77 finite edges. A strict spacelike escape proves this prefix has infinite volume. The experiment does not decide reflectivity of the lattice. It differs from the now-excluded discriminant11 host even as a quadratic space up to similarity over `Q(√5)`.

The [enumeration source](../scripts/fresh_h9_beta55_pilot_20261014.py) and [certificate](fresh_h9_beta55_pilot_20261014.json) cover58 positive control levels and2042 search nodes. The separate [geometry checker](../scripts/fresh_h9_beta55_pilot_20261014_geometry.py) and [geometry certificate](fresh_h9_beta55_pilot_20261014_geometry.json) test all220 nine-wall supports. The arithmetic uses the frozen standalone Fraction helper from the earlier independent Vinberg replay; this is a new producer, not an independent audit of itself.

## Host and complete reflection norms

Put `α=2+φ`, `β=4−φ` and `δ=1+8φ`. Then

`δ=φαβ`, `Norm(α)=5`, `Norm(β)=11`, `αβ=7+φ`.

The determinant is `−δ`, whose field norm is `−55`. The physical form has signature `(9,1)`. Its conjugate is positive definite: the conjugate first binary diagonal is positive and the conjugate binary determinant `−δ′` is positive. Both α and β are totally positive, and their ideals are distinct primes. The determinant ideal is squarefree, so there is no proper integral overlattice of this lattice.

The [previous ring and norm-alphabet review](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) supplies the ring-only facts: `Z[φ]` is a norm-Euclidean PID, its units are `±φ^k`, and every totally positive unit is a square. For a primitive integral reflection vector of norm `q=2b`, evenness gives `b∈Z[φ]`. Reflection integrality and a primitive-coordinate Bezout identity imply `Gr∈bO^10`. Applying the adjugate and Bezout again gives `b|det(G)`. Total positivity of q and squarefreeness therefore give exactly four possible norm representatives up to unit rescaling:

`2`, `2α`, `2β`, `2αβ`.

A vector of one of these norms is automatically primitive: a common nonunit divisor would have its square dividing the squarefree b. The search additionally verifies every reflection-integrality entry directly.

The ratio of determinant norms for this host and the discriminant11 host is5, which is not a rational square. Similarity in rank10 changes determinants by a square in the field, so these forms are not similar over `Q(√5)`. The same norm-square argument separates this form from unit-determinant forms. These calculations do not assert novelty or a literature-wide classification. A focused repository check found no literal prior discriminant55 search.

## The zero sector has only eight simple roots

Use coordinates `(x,y1,…,y8,t)` and `H=(4,0,…,0,1)`. Then

`H²=−4δ`, `−(H,r)=δx`.

The last coordinate vector has norm4 but is **not** a lattice reflector: its reflection coefficient against the first basis vector is `−1/2`. Hence imposing its halfspace would discard valid roots.

For a zero-height root `x=0` of norm `2b`, the reflection filter and the unimodular E8 block force b to divide every y-coordinate and t. A nonunit b would then give an even lattice vector of norm `2/b`, impossible. For b=1, write the coordinates as `u+φv` in the rational positive lattice `E8⊕<4>`. The norm equation gives `u²+v²=2` and `2(u,v)+v²=0`. Both rational norms are even nonnegative integers. These equations force v=0 and u an ordinary E8 root, with t=0.

Thus there are exactly240 zero roots, and the simple sector consists of the eight negative E8 simple roots. The certificate checks their literal ordering against the frozen complete ordinary E8 shell. H is not an initial finite vertex.

## Complete bounded enumeration

Let `C` be the E8 Cartan matrix, `J=C⁻¹`, `λ=Cy`, and `l=x−4t`. Zero-sector acuteness means `λ≥0`; **l has either sign**. The norm equation is

`4λᵀJλ+l² = 4q+δx²`,

and the geometric height is `δx²/(4q)`. Reflection integrality is equivalent to `b|λ_i,l` together with `l≡x (mod 4O)`. The remaining first binary coefficient follows from `4(Gr)_0=l−δx` and coprimality of4 and b.

All four norms satisfy `q<18`, `q′<13`; also `δ>13` and `−δ′>3`. At height20 this bounds `x<11`, `|x′|<5`. Every positive-coordinate trial satisfies `|l|<39`, `|l′|<8`, while the λ bounds are stronger by a factor of2 in the squared norm. Rational bounds on φ convert these to explicit finite coefficient boxes. Every candidate is then filtered exactly at both embeddings. The output records the complete control and coordinate alphabets.

The recursion fixes signed l first and then visits `λ8,…,λ1`. Schur bounds are recomputed from the literal positive form. For an earlier wall `(X,Y,T)`, with `L=X−4T`, acuteness is

`4λ·Y+lL ≤ δxX`.

After l is fixed, every unassigned contribution is nonnegative because `Y=JΛ≥0`. This justifies partial-wall pruning without an incorrect sign restriction on l. All equal-height survivors are checked mutually acute before acceptance.

Four positive walls appear. The first two have norms `2αβ` and `2α`, with vectors `(φ−1,0^8,−2)` and `(φ−1,0^8,φ)`. The next two have norm2. In particular a short/prime5 pair has Cartan two-cycle `2+φ`, giving an actual π/10 relation; the old beta11 mixed-norm obstruction does not apply here. Every finite pair angle in the twelve-wall output is checked exactly.

## Geometry and stopping point

The twelve normals span the ten-dimensional ambient space and admit a strict future timelike center. The usual acute projection formula supplies an explicit relative-interior point on every facet.

The complete220-support cone census has23 finite vertices, no ideal vertices, and27 spacelike extreme rays;170 supports give infeasible orientations. All finite vertices have exactly nine incident walls and positive-definite incident Grams. There are77 finite edges and53 remaining finite-corner germs. A displayed ray plus a small positive multiple of the strict center yields a strict spacelike escape, checked against all twelve walls.

Generation and saved replays pass for both the complete prefix and its geometry. This result is a new arithmetic-host experiment with a reproducible finite stopping point. It supplies no finite-volume example and no host-wide obstruction. A fresh height100 extension and finite-corner isometry check are the next bounded tests; they must retain all four norm classes and the signed residual coordinate.
