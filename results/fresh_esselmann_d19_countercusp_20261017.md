# An explicit bad cusp in the genus of D19

Status: **exact generation and fresh saved replay PASS**, 2026-10-17. The [standalone certificate](../scripts/fresh_esselmann_d19_countercusp_20261017.py) and [saved data](fresh_esselmann_d19_countercusp_20261017.json) construct an even positive lattice M of rank 19 in the genus of D19 whose integral reflection roots span only 18 dimensions. Consequently U⊕D19 is nonreflective. Here U is the integral hyperbolic plane with Gram matrix [[0,1],[1,0]].

This is one explicit obstruction illustrating the exceptional dimension gap. It is **not** a reconstruction of Esselmann's classification of every rational Lorentz lattice of rank 21. In particular, it makes no exclusion claim for nonarithmetic hyperbolic 20-dimensional Coxeter polytopes. No novelty claim is made.

## The lattice and its primitive D5

Use positive E6 with simple roots indexed 0 through 5 and bonds

```
(0,1), (0,2), (2,3), (0,4), (4,5).
```

Every root has square 2; a bond has inner product −1. If G is this Gram matrix and ω=G⁻¹e3, then

```
ω = (2,1,5/3,4/3,4/3,2/3),     ω² = 4/3,
w = 3ω = (6,3,5,4,4,2),        w² = 12.
```

Deleting root 3 leaves D5 on indices 0,1,2,4,5, and w is orthogonal to it. The coordinates of w have gcd 1, so E6∩Rω=Zw.

The inverse matrix G⁻¹ is entrywise positive and its smallest diagonal entry is 4/3. A finite Weyl group moves any dual vector into the dominant weight cone, where its fundamental-weight coefficients are nonnegative integers. It follows that every nonzero vector of E6* has square at least 4/3; ω attains this bound. The exact root orbit has 72 vectors. For completeness, a dominant norm-2 vector cannot have two nonzero weight coefficients, since their diagonal contribution alone is at least 8/3. The only integral single fundamental weight of square 2 is G⁻¹e1, which occurs in that orbit. Thus the orbit is the entire norm-2 set, not a sampled subset.

In E6*⁴, take the ternary glue code generated modulo 3 by

```
(1,0,1,1),   (0,1,1,−1).
```

Its nine words comprise zero and eight words of weight 3. Let

```
N = E6⁴ + Z(ω,0,ω,ω) + Z(0,ω,ω,−ω).
```

Both glue generators have square 4 and their mutual product is zero. N is even, has index 9 over E6⁴ and determinant 3⁴/9²=1. This reconstructs the E6⁴ Niemeier lattice directly; the familiar glue description is also recorded in [Bertin–Lecacheux, Table 1, p.603](https://www.jstage.jst.go.jp/article/jmath/72/2/72_599/_pdf).

No nonzero codeword is supported on only the first component. Therefore **N∩span(first E6)=E6**. The five chosen simple roots span a primitive sublattice of E6, hence their D5 is primitive in N. The certificate also checks this directly by Smith normal form: the five embedding invariants are all 1. A nonzero glue coset has at least three nonzero E6* components, so its vectors have square at least 3·4/3=4. Consequently the norm-2 roots of N are exactly E6⁴.

## The rank-19 complement and the genus comparison

Let M=D5⊥∩N. A concrete full-rank sublattice is R=Zw⊕E6³, with determinant 12·3³=324. In this decomposition, adjoining

```
(w/3,0,ω,ω),   (0,ω,ω,−ω)
```

gives an even overlattice of index 9 and determinant 4. It lies in D5⊥∩N. Since D5 is primitive in the unimodular N, its full complement also has determinant 4, so this displayed overlattice is the full M. The certificate saves a literal 19×19 integral Gram matrix, its positive LDL pivots, and its integral embedding into N.

The discriminant group of M is cyclic of order 4, with generator square 3/4 modulo 2Z. This is both checked from the displayed Gram inverse and consistent with the opposite discriminant form of D5: a D5 generator has square 5/4, and −5/4=3/4 modulo 2Z. A D19 generator has square 19/4=3/4 modulo 2Z. Thus M and D19 have the same signature and discriminant form and belong to the same genus.

The last implication uses Nikulin's Corollary 1.9.4 on even lattices ([original article](https://www.mathnet.ru/eng/im1677)). The original English PDF was not successfully retrieved in this session; the precise statement and attribution were checked in [Arima–Shimada, Proposition 1.1, p.1](https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/A10A9/A10A9.pdf). Equality of determinants alone would not justify the genus assertion.

## All integral reflection directions, including long roots

Let r be a primitive vector of M whose reflection preserves M, and write q=r² and d=gcd{(r,x):x∈M}. Because r is primitive, the class r/d in M*/M has exact order d. The discriminant exponent is 4, so d divides 4. Reflection integrality says q divides 2d. Since M is even, the only possible positive squares are

```
q ∈ {2,4,8}.
```

This bound concerns every primitive reflecting vector, without a height cutoff. The same elementary divisibility condition appears in [Scharlau–Walhorn, p.279](https://www.numdam.org/item/AST_1992__209__279_0.pdf).

For q=4 or 8, project r to any of the three surviving E6 real factors. Pairing with the E6 sublattice and reflection integrality force that projection to lie in (q/2)E6*. If nonzero, its square is at least

```
(q/2)²·4/3 = 16/3 > 4,     or     64/3 > 8,
```

already larger than the entire vector's square. All three projections must therefore vanish. Then r lies on the first-component axis. The no-single-factor glue fact proved above gives M∩Rω=Zw, whose least positive square is 12. This rules out both q=4 and q=8.

For q=2, the norm-2 description of N shows that the roots in M are precisely E6³. No first E6 root lies on the perpendicular axis of minimum square 12. Hence M has 216 norm-2 roots, spanning exactly 18 dimensions, and no other reflecting directions.

Here **nonreflective positive lattice** means that its integral reflection roots fail to span the full Euclidean space. Its reflection group is of course finite; group infinitude is not the criterion in positive signature. This is the definition used in [Scharlau–Walhorn, Definition 1.1, p.280](https://www.numdam.org/item/AST_1992__209__279_0.pdf).

## Why this obstructs the hyperbolic lattice

Vinberg's cusp lemma requires the Euclidean quotient at every primitive isotropic vector of a reflective Lorentz lattice to be reflective. Its genus version, [Scharlau–Walhorn, Lemma 1.3, p.281](https://www.numdam.org/item/AST_1992__209__279_0.pdf), states that reflectivity of U⊕L for positive L of rank at least 2 requires every lattice in genus(L) to be reflective. Applying it to L=D19 and the explicit M proves that **U⊕D19, of signature (20,1), is nonreflective**. Geometrically, its bad cusp has a Euclidean direction not bounded by any integral reflection wall.

This is why inspecting only the obvious D19 cusp would miss the obstruction: D19 itself has a full root system, while another lattice in its genus does not. The argument uses the full cusp/genus condition, not just the visible roots at one chosen cusp.

## A concrete contrast one dimension higher

The same Niemeier source also makes clear why root-span failure need not descend monotonically with dimension. Delete only the D4 on indices 0,1,2,4 of the first E6. Its perpendicular plane contains

```
x = 2ω3−ω5 = (2,1,2,2,1,0),
y = −ω3+2ω5 = (2,1,1,0,2,2),
Gram(x,y) = [[4,−2],[−2,4]].
```

The complement is obtained from Zx⊕Zy⊕E6³ using the same two glue generators, with ω3=(2x+y)/3. Both x and y pair evenly with this whole complement, so their norm-4 reflections are integral. The certificate constructs its 20×20 Gram matrix and verifies all 20 spanning reflections directly. Its discriminant form is (Z/2)² with each nonzero square equal to 1 modulo 2Z, as for D20.

Thus this particular D4 complement has reflection roots spanning all 20 directions; the D5 complement loses its last direction. This comparison **does not prove that the entire genus of D20 is reflective**, nor does it supply the complete classification behind Borcherds' dimension-21 example.

## Published statement and limits of this reconstruction

[Esselmann's 1996 paper](https://doi.org/10.1006/jnth.1996.0141) is the classification source. Its full text was not accessible in this session, so no claim is made to have checked all of its rank-21 elimination cases. [Allcock, p.3](https://web.ma.utexas.edu/users/allcock/research/rk3monolithic-17.pdf), explicitly attributes to Esselmann that the only reflective lattices in ranks at least 21 are two rank-22 cases: the even sublattice of I21,1 and its dual. This is the precise published exception statement used to distinguish hyperbolic dimensions 20 and 21 in the rational/integral arithmetic setting.

The present certificate supplies a complete, independently reproducible instance of the positive-lattice/cusp mechanism, including the long-root check that a mere norm-2 root count would omit. It does not replace the classification needed to exclude every other rational Lorentz lattice of signature (20,1).

Reproduction:

```
python scripts/fresh_esselmann_d19_countercusp_20261017.py
python scripts/fresh_esselmann_d19_countercusp_20261017.py --verify
```

Both completed with `PASS_EXACT_ILLUSTRATIVE_D19_GENUS_NONREFLECTIVE_CUSP 19 18`. The saved payload binds the exact source hash. Arithmetic uses SymPy rational matrices and integer normal forms; no repository lattice tables or producer arithmetic are imported.
