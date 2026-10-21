# A new search field: conductor-nine dihedral hosts

Four arithmetic H9 hosts over the real cubic field of discriminant 81 were searched through Vinberg height 300. Their prefixes have **472, 204, 410 and 204 Coxeter walls**, and all retain independently certified infinite-volume escapes. Two H8 variants fail the same volume gate. This changes the ground field from the earlier quadratic and cubic-seven searches; it does not establish a new finite-volume example or nonreflectivity of any full host.

**Subsequent full-host result:** the [infinite-order chamber certificate](fresh_cubic9_global_nonreflectivity_20261013.md) now proves nonreflectivity of the α=Tβ H3 host and rules out finite-volume integral reflection completions of its H9 counterpart by the necessary-face argument. The other three coefficients are not closed by that result. The bounded artifacts below remain valid with their original scope.

## Field and complete primitive norm coverage

Put `T=2cos(π/9)`, so `T³−3T−1=0`, with distinguished embedding `T∈(1,2)` and other roots in `(-1,0),(-2,-1)`. Define

```
β = 4−T²,       D = 2+T−T²,
B = [[2,−T], [−T,2]],
Lα = <-α> + B⁴ + <1>,       α ∈ {T, 2T, Tβ, 2Tβ}.
```

Each α has signs `(+,-,-)`, while B is positive definite at all three embeddings. Thus Lα has signature `(9,1)` at the distinguished place and is positive definite at the other two.

The polynomial discriminant is 81. A proper integral-basis index would give field discriminant 9 or 1, contradicting the totally real cubic Minkowski lower bound. Hence `O=Z[T]`. Its ideal-class Minkowski bound is 2, and 2 is inert because `X³+X+1` is irreducible modulo 2. There is no ideal of norm two, so the class number is one. The units T, D and −1 have signatures `(+--),(++-),(---)`, generating all eight patterns. Since the unit group has eight classes modulo squares, every totally positive unit is a square. Here `D²=2−T`, `N(D)=−1`, and the totally positive β has norm 3.

The block Smith factors are `(1,β)`. The dual exponent e is β for α=T,Tβ and 2β for α=2T,2Tβ. For a primitive integral reflection vector, its norm q divides 2e: reflection integrality gives `2e*r/q∈L`, and primitivity removes the vector coordinates. Together with class number one and the unit-sign calculation, this gives the safe norm pool `{1,2,β,2β}`, with `{4,4β}` added for even α. These last two classes are imprimitive. Since B is invertible modulo 2, reflection integrality makes all spatial coordinates even; the norm equation then makes the time coordinate even too. Keeping these conservative classes in the producer does not omit a mirror; exact ray deduplication removes parallel copies.

## The full zero-height chamber

The complete zero-height chamber is **`I₂(18)^4×A1`**, with nine simple walls and 146 oriented primitive roots. In particular the coordinate swap in each binary block must be included: retaining only an `I₂(9)` chamber would miss a reflection.

Write `F(x,y)=x²−Txy+y²`. A nonzero integral F is totally positive, with trace at least three. At norm one this eliminates every binary block, leaving the two scalar roots. At norm two, the scalar coordinate is even by reduction modulo the inert prime 2; trace positivity eliminates it and leaves exactly one block with F=1.

For β-divisible norms, reflection integrality writes each binary vector as `adj(B)t` and the scalar as βs. After dividing by β the norm equation is `2ΣF₊(t)+βs²=1` or 2, where `F₊(x,y)=x²+Txy+y²`. The first has no solution, by trace positivity and β-valuation. In the second, reduction modulo 2 makes s even; the trace bound eliminates it and leaves exactly one F₊=1 block. Norms four and 4β supply no additional primitive mirrors.

An independent unpruned coefficient-box census finds all eighteen oriented F=1 roots and eighteen dual roots of norm 2β in each block. Its coordinate bounds follow from `x²,y²≤4/β` at every embedding and the exact trace-dual basis. The simple block roots are `−e₁` and `D(e₁−e₂)`, of norms 2 and 2β. Their angle is π/18, with `4cos²(π/18)=2+T`. The latter reflection swaps the coordinates, giving the dominant covector inequalities `0≤u≤v` used in enumeration.

## Complete bounded searches

For height `h=αx²/q≤H`, the positive-definite conjugate forms imply `x'²≤−q'/α'` at both other embeddings. Exact trace-dual bounds enclose a finite integral alphabet. The producer uses rational isolating intervals to prune it, then exact field signs for membership. Each block covector satisfies `u²,v²≤2(q+αx²)` at all embeddings. Inverse-block reconstruction enforces coordinate and reflection integrality, and energy joins enumerate every four-block multiset and every placement. Equal-height packets are selected together.

Every accepted pair is separately checked for a Coxeter angle or disjointness. The finite-angle alphabet is `4cos²(π/m)∈{0,1,2,3,T²,2+T}`, namely labels `m∈{2,3,4,6,9,18}`. Arithmetic integrality alone is not used to assert a Coxeter angle.

| Time coefficient α | H30 controls | H30 walls | H300 controls | H300 walls | Volume gate |
|---|---:|---:|---:|---:|---|
| T | 12 | 18 | 42 | 472 | Exact infinite-volume prefix |
| 2T | 12 | 20 | 45 | 204 | Exact infinite-volume prefix |
| Tβ | 12 | 18 | 42 | 410 | Exact infinite-volume prefix |
| 2Tβ | 12 | 19 | 45 | 204 | Exact infinite-volume prefix |

For example the α=2T, H300 cone contains the spacelike vector

```
x=1;
B1=B2=B4=(0,0), B3=(2T,2T);
scalar=T²−T−1;
q(v)=−9−29T+18T² > 0.
```

Every wall pairing is nonpositive. The independent checker adds a sufficiently small positive multiple of the certified strict timelike interior, keeping the square positive and making every wall inequality strict. This proves an open patch at infinity and infinite hyperbolic volume, rather than merely reporting a failed compactness test.

Removing the scalar coordinate gives H8 hosts. At height 300 the α=2T and Tβ cases have **182 and 148 walls**, respectively, and both have independently verified strict spacelike escapes. These are bounded negative outcomes, not classifications of all H8 chambers over this field.

A subsequent [necessary-face screen](fresh_cubic9_high_label_faces_20261013.md) tests the lower hosts with one and two binary blocks plus the scalar. The four coefficients in the table give H3 prefixes with respectively **7,7,10,7 walls**, and H5 prefixes with **57,21,44,21 walls**, all through height 300. All eight retain independently replayed strict spacelike escapes. These lower-dimensional bounded failures do not exclude a later finite-volume completion or the full H9 host.

## Independent checks and calibration

The standalone [Fraction verifier](../scripts/fresh_cubic9_dihedral_independent_audit_20261013.py) imports no repository arithmetic or geometry. It uses generic polynomial reduction and rational root isolation, reconstructs the entire zero-height census, and checks saved roots, integral reflections, chronology, every pair angle, strict interiors and escapes. It rejects forged control heights and timelike fake escapes.

The separate [enumeration verifier](../scripts/fresh_cubic9_enumeration_audit_20261013.py) uses unpruned coefficient boxes and a Cartesian spatial join. It independently repeats every time control, spatial candidate and acceptance decision through height 30. At height 300 it independently repeats **all time/norm controls**, but does **not** repeat the full spatial census. The saved H300 walls and geometry are checked separately. The H8 cases receive saved-root/geometry replay; their positive spatial enumeration is not independently repeated.

The same implementation exactly recovers the known compact **(2,4,18) triangle** using one binary block and α=T. All three vertices are timelike, and a positive dependence of spanning spatial normals proves affine boundedness. This is a positive calibration, not a new example.

Audit artifacts:

- [Original three H9 hosts and triangle](fresh_cubic9_dihedral_independent_audit_20261013.json).
- [Complete H30 enumeration and α=2T H300 controls](fresh_cubic9_enumeration_audit_20261013.json).
- [Other original H300 controls](fresh_cubic9_additional_controls_audit_20261013.json).
- [Fourth-host root and geometry audit](fresh_cubic9_fourth_host_audit_20261013.json) and [enumeration audit](fresh_cubic9_fourth_host_enumeration_audit_20261013.json).
- [H8 root and geometry audit](fresh_cubic9_h8_independent_audit_20261013.json).
- [Eight H3/H5 root and geometry audits](fresh_cubic9_lower_faces_audit_20261013.json).
- [Independent source review](fresh_cubic9_source_review_20261013.md).
- [Arithmetic and certificate regression tests](../tests/test_fresh_cubic9_20261013.py).

Example replay:

```sh
PYTHONPATH=src:. python scripts/fresh_cubic9_dihedral_h9_20261013.py --alpha 2theta_beta --height 300 --output results/fresh_cubic9_dihedral_h9_2theta_beta_h300_20261013.json
python scripts/fresh_cubic9_dihedral_independent_audit_20261013.py results/fresh_cubic9_dihedral_h9_2theta_beta_h300_20261013.json --output /tmp/cubic9_geometry_audit.json
PYTHONPATH=src:. python -m pytest -q tests/test_fresh_cubic9_20261013.py
```

The field contains the π/18 cyclic product, so these wall sets are not simply a presentation over the earlier quadratic or cubic-seven fields. No worldwide novelty follows from that observation. Higher mirrors for the remaining coefficients, other basepoints or wall selections, and other lattices over this field remain open; the subsequent full-host theorem closes the α=Tβ direct-sum family.
