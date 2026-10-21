# The compact H4 prism has no single-top H5 prism lift with finite side angles

The [certified compact H4 prism](fresh_h9_face_bridge_20261013.md) cannot be extended to a full-dimensional compact H5 prism by retaining its six side hyperplanes and adding one top with finite Coxeter angles to all six. The obstruction covers **every finite label**, not a bounded label list. It also rules out finite-volume prisms in this same construction, allowing ideal top vertices.

The scope includes the prism requirement that the new top closes the edge above every old H4 vertex. It does **not** exclude several new caps, dotted top/side relations, facet replacements, or other completions of the H9 chamber. No new polytope is claimed here.

## The radical forces one right angle

Use norm-two normals throughout. Rescaling a proposed spacelike top normal to norm two is harmless. Let `phi=(1+sqrt(5))/2` and

```text
T = [seed6, crosspair0, nextH4, compact0],
L = Gram(T), the chain 3–5–3,
b = wall17,
r = facePrism = -d_3 - sqrt(a)b,
a = (9+3phi)/11,
(d_i,T_j) = delta_ij.
```

Here `b` is perpendicular to all four `T_j`. The six retained prism sides are `T_0,T_1,T_2,T_3,b,r`. For a proposed top `w`, write `p_j=(w,T_j)`, `p_b=(w,b)` and `p_r=(w,r)`. Their linear relation is

```text
p_r = -sum_j (L^-1)[3,j] p_j - sqrt(a) p_b.
```

Each finite Coxeter pairing is `-2cos(pi/m)` for an integer `m>=2`, so it belongs to the maximal abelian extension `Q^ab`. The coefficients of `L^-1` belong to `Q(sqrt(5))`, also a subfield of `Q^ab`.

In contrast, `sqrt(a)` does not belong to `Q^ab`. Its polynomial is

```text
11 X^4 - 21 X^2 + 9.
```

Indeed, `Norm(a)=9/11` is not a rational square, so `a` is not a square in `K=Q(sqrt(5))`, and `F=Q(sqrt(a))` has degree four. This field contains `K`, since `phi=(11a-9)/3`. If `F` were normal over `Q`, it would contain the conjugate `sqrt(a')`, hence also `sqrt(11)` because `sqrt(a)sqrt(a')=3/sqrt(11)`. Since `sqrt(11)` is not in `K`, this would force `F=K(sqrt(11))`. An element of this biquadratic field whose square is in `K` lies either in `K` or in `sqrt(11)K`; thus either `a` or `a/11` would be a square in `K`. Their norms, `9/11` and `9/1331`, exclude both alternatives. Therefore `F` is not normal, while every finite subfield of `Q^ab` is Galois and abelian over `Q`.

The pairing relation now forces **`p_b=0`**: otherwise it would express `sqrt(a)` as a quotient of elements of `Q^ab`.

## One existing H4 corner forces all the other right angles

The existing prism vertex with sides `[T_1,T_2,T_3,r]` has the positive-definite H4 Gram matrix

```text
C = [[ 2, -phi,  0,  0],
     [-phi,  2, -1,  0],
     [ 0,   -1,  2, -1],
     [ 0,    0, -1,  2]].
```

Every entry of `C^-1` is strictly positive. Its diagonal is

```text
[8+12phi, 12+18phi, 6+8phi, 2+2phi],
```

whose entries are all greater than two. Let `c` be the nonnegative magnitudes of the top's four finite pairings to this corner. If any is nonzero, it is at least one, so `c^t C^-1 c>2`. Consequently the five-normal Gram matrix has negative Schur complement `2-c^t C^-1 c` and is indefinite. It cannot describe a finite or ideal endpoint of the edge above this old vertex: normals incident to such an endpoint have positive-definite or positive-semidefinite Gram matrix, respectively.

Thus closing this edge forces `p_1=p_2=p_3=p_r=0`. Together with `p_b=0`, the linear relation and

```text
(L^-1)[3,0] = -(4+5phi)/11 != 0
```

force `p_0=0` as well. The top is therefore perpendicular to the entire old H4 ambient space `span(T,b)`, of signature `(4,1)`.

In an H5 lift along one of the remaining common A1 directions, the orthogonal complement of this space is the single positive line defining the old base. The proposed top hyperplane repeats that base; its two orientations either duplicate the old inequality or collapse the lift onto the old H4 hyperplane. Neither yields a full-dimensional finite-volume prism. Components along the other four common A1 directions do not help: they vanish on the intrinsic H5 face, leaving the same repeated hyperplane.

## Exact replay and independence

The [standalone certificate](fresh_h9_prism_lift_20261013.json) binds the existing prism recipe and verified geometry by SHA-256. Its [checker](../scripts/fresh_h9_prism_lift_20261013.py) uses separate standard-library `Fraction` arithmetic in `Q(phi)`. It reconstructs the inverse tetrahedral Gram, the cap relation and norm, the rational nonsquare checks, the actual H4 corner from the saved vertex Gram, its positive pivots, and every inverse-entry sign. It imports no cap producer or existing field arithmetic. The general abelian-field and endpoint arguments above supply the all-label theorem; this is not a finite-label enumeration.

```sh
python scripts/fresh_h9_prism_lift_20261013.py --verify
```

Result: `PASS_ALL_LABEL_SINGLE_TOP_H5_PRISM_OBSTRUCTION`. The actual prism/face incidence is inherited from the linked, hashed geometry certificate; this replay does not re-enumerate the entire 96-wall chamber.

The [independent mathematical/source review and fresh replay](fresh_h9_prism_lift_20261013_review.md) pass. The root agent also reviewed the proof and checker and reproduced the certificate.
