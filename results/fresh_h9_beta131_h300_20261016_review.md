# Independent norm131 continuation and chamber symmetry review

The [independent bounded-height checker](../scripts/fresh_h9_beta131_h300_20261016_review.py)
reconstructs the complete height300 Vinberg prefix on the norm131 lattice.
The [two-corner reviewer](../scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py)
then reconstructs an infinite-order integral symmetry of its full reflection
chamber. The resulting host-wide conclusion is nonreflectivity, rather than
an inference from the infinite volume of a finite wall prefix. No new
finite-volume example is claimed.

The bounded census uses the frozen independent Golden-arithmetic engine,
with natural covector order, a literal E8 Cartan matrix, independently larger
coordinate boxes, and a weaker coordinate-energy bound than the producer.
It imports neither the producer nor its arithmetic helper. All200 height
controls, the282 roots, their pairings, a strict interior point and all actual
facet witnesses are checked. The independently enumerated height100 control
prefix is preserved exactly. The source manifests bind both the frozen
height100 review and the independent census engine.

## Complete primitive norm alphabet and zero sector

Use O=Z[φ], φ²=φ+1, and coordinates (x,y₁,…,y₈,t). The literal form is

`G = E8 ⊕ [[−2φ,−1],[−1,6]]`,

with the binary block on x,t and the E8 numbering specified in the checker.
Write δ=1+12φ and β=φδ=12+13φ. Then detG=−δ, N(δ)=−131 and N(β)=131.
The physical signature is (9,1), while the conjugate form is positive
definite. The basepoint is H=(6,0⁸,1), with GH=(−δ,0⁹).

The form is even, so a primitive reflector has norm2b. Reflection integrality
gives Gr/b integral. Multiplying by the adjugate and using primitivity over
the principal ideal domain O gives b dividing detG. Its prime norm131 leaves
only the unit and δ ideals. The totally positive units are the squares φ²ᵏ;
normalizing by a unit therefore leaves exactly norms2 and2β. This proves
completeness at arbitrary height, not just at the current search bound.

At height zero x=0. A norm2β root would have y and t divisible by β, using
the unimodularity of E8 and the x-covector −t. Dividing the entire root by β
would give norm2/β, contradicting even integrality. For norm2, the bounds
6t²≤2 at both embeddings force the integral norm N(t) to vanish. Thus t=0.
Writing y=u+φv in the E8 root lattice gives u²+v²=2 and 2(u,v)+v²=0.
Evenness and positivity of E8 force v=0. The zero roots are exactly the240
ordinary E8 roots, with the eight checked simple sector walls.

These arguments are inherited from the frozen pilot and were checked again
for their use in the full-chamber conclusion. The height300 search enumerates
both norm classes, signed t-related coordinates, and every admissible control
within its proven boxes. There is no omitted long-root or zero-height sector.

## Independent two-corner witness

The symmetry review imports only the project's Golden arithmetic. It rebuilds
the literal form and reads the two specified root frames, whose source indices
are `[0,1,2,3,4,5,6,7,10]` and whose ordered target indices are
`[5,2,4,3,1,244,42,11,7]`. Their exact positive-definite Gram matrices agree.

It computes each corner apex directly from H by orthogonal projection off
its nine normals. If S is their Gram matrix and
`b = S⁻¹ ((H,Nᵢ))`, then `V = H − Σ bᵢNᵢ`.
The checker verifies bᵢ≤0, V²<0, (H,V)=V²<0 and every one of the282 wall
inequalities. Exactly the nine indicated walls contain each apex. This does
not depend on the producer's assertion that its70-corner census is complete.

The map is reconstructed from the two root-plus-apex frames. Its positive
apex scale is determined by their determinant ratio and checked against the
two apex norms. All nine oriented root images agree, the matrix and inverse
are integral, and `TᵀGT=G`. The matrix has determinant−1 and trace

`32 + 53φ > 10`.

This proves that T has infinite order: a finite-order real10-by10 matrix has
eigenvalues on the unit circle and real trace of absolute value at most10.
Positive apex scaling between two future points proves preservation of the
future component.

## Why the witness concerns the full reflection chamber

Let r be any later simple root in the complete Vinberg process. It satisfies
(H,r)<0 and (Nᵢ,r)≤0 for every existing wall. At either verified corner the
independent projection gives `H = V + Σ bᵢNᵢ`, with bᵢ≤0. Hence

`(V,r) = (H,r) − Σ bᵢ(Nᵢ,r) < 0`.

Thus no later simple wall reaches or cuts either corner. The full lattice
mirror arrangement is locally finite in hyperbolic space, so their exact
finite corner sectors persist in the full chamber. The integral isometry T
preserves all lattice mirrors and maps one oriented corner sector to the
other. The image chamber therefore agrees with the original chamber near a
corner and is the same global chamber.

An infinite-order chamber stabilizer makes the full reflection subgroup
infinite index in the arithmetic lattice O⁺(L). The admissible integral form
gives O⁺(L) finite covolume, so its reflection chamber has infinite volume.
Any reflection group all of whose wall reflections preserve this literal
lattice is contained in the full reflection subgroup and cannot have finite
covolume. This excludes the host, without claiming anything about other
lattices or arbitrary geometric walls.

The distinction of determinant norm131 from the earlier norm11 and norm55
hosts is only a distinction from those recorded similarity classes; this is
not a classification or priority claim. The full corner census and map search
helped find the witness, but their exhaustiveness is not a premise of this
independent nonreflectivity proof.

The final machine-readable outputs are the
[prefix review](fresh_h9_beta131_h300_20261016_review.json) and
[symmetry review](fresh_h9_beta131_h300_symmetry_20261016_review.json).
The full symmetry status is gated on successful independent reconstruction
of the complete bounded prefix. A `--witness-only` run writes no certificate
and deliberately reports that this dependency is still required.

```sh
python3 scripts/fresh_h9_beta131_h300_20261016_review.py
python3 scripts/fresh_h9_beta131_h300_20261016_review.py --verify
python3 scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py
python3 scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py --verify
```
