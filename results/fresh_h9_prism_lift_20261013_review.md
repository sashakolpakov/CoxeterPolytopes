# Independent review of the single-top prism obstruction

Source and mathematical review: **PASS**. A fresh invocation of the [standalone Fraction checker](../scripts/fresh_h9_prism_lift_20261013.py) reproduced [the saved certificate](fresh_h9_prism_lift_20261013.json):

```sh
python scripts/fresh_h9_prism_lift_20261013.py --verify
```

The reviewed [theorem report](fresh_h9_prism_lift_20261013.md) excludes a compact or finite-volume H5 prism obtained from the certified H4 prism by retaining all six side hyperplanes and adding a single top with finite Coxeter angles to all six, closing the edge above every old vertex. It does not exclude several caps, dotted side relations, changed side walls, other H5 faces, or general H9 completion.

The field argument is valid. For `a=(9+3phi)/11`, the norm `9/11` proves that `a` is not a square in `Q(sqrt(5))`. The degree-four field `Q(sqrt(a))` contains `Q(sqrt(5))`. If it were normal, the conjugate square root would force `sqrt(11)` into it, making it `Q(sqrt(5),sqrt(11))`. Squaring an element of this quadratic extension shows that either `a` or `a/11` would then be a square in `Q(sqrt(5))`; their nonsquare rational norms exclude both cases. Thus `sqrt(a)` is not in any finite abelian extension. Every finite collection of finite-angle cosine pairings lies in an abelian cyclotomic extension, so the exact six-side relation forces the top's pairing with `wall17` to vanish. This uses all finite labels, with no label cutoff.

The geometric step is also valid. Normalize the proposed spacelike top to norm two. The saved actual corner has the positive-definite H4 Gram matrix of type `5–3–3`. Its inverse has strictly positive entries, and every diagonal entry exceeds two. A nonzero finite Coxeter pairing has magnitude at least one, so any nonzero coupling to this corner makes the five-wall Schur complement negative. This is incompatible with both a finite endpoint, where the incident Gram matrix is positive definite, and an ideal endpoint, where it is positive semidefinite. The exact side relation then forces all remaining side pairings to vanish.

The resulting top is orthogonal to the old H4 ambient space. Inside the proposed intrinsic H5 lift only the freed common A1 direction remains. The top therefore repeats the old base hyperplane or vanishes on the whole face; it cannot close a full-dimensional prism. Components along the four common A1 walls retained in defining the H5 face disappear on that face and do not change the conclusion.

This review checked the field and endpoint arguments independently and inspected the checker's source bindings, norm conventions, inverse Gram calculations and extraction of the actual H4 corner. Its fresh replay uses the author's separate standard-library arithmetic. Actual prism incidence is inherited from the hashed [96-wall geometry certificate](../data/enumerations/fresh_h9_face_bridge_20261013_verify.json), which already has a separate [symbolic and higher-precision review](fresh_h9_face_bridge_20261013_review.md). This is not a second enumeration of that chamber.
