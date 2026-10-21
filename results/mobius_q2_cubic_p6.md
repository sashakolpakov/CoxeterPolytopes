# The `q=2` Möbius endpoint is the Hermitian cubic and the classical `P^6`

## Result

Let `Q=P1 x P1` over `F4`, and let

```text
Gamma_21 : y=x^2,              Gamma_12 : x=y^2.
```

Their five intersections are the points

```text
((s:t),(s^2:t^2)),       (s:t) in P1(F4).
```

Blow them up and call the resulting surface `X`.  The exact certificate proves:

1. the five points are in del-Pezzo general position;
2. `X` is a degree-three del Pezzo surface, and its anticanonical model is the
   smooth Hermitian cubic

   ```text
   z0^2 z3 + z0 z3^2 + z1^2 z2 + z1 z2^2 = 0
   ```

   over `F4`;
3. the five exceptionals, ten ruling fibres, ten `(1,1)` projectivity graphs,
   and the two Frobenius graphs are exactly all 27 `(-1)`-curves;
4. their dual nef chamber in `NS(X)=I_(1,6)` is the classical
   Potyagailo--Vinberg right-angled polytope `P^6`, Euclidean-dual to Gosset's
   `2_21`;
5. its exact projective cone has

   ```text
   72 finite vertices + 27 ideal vertices + 0 hyperideal rays.
   ```

Thus this is a positive, exact algebraic-surface presentation of a genuine
finite-volume hyperbolic Coxeter polytope in `H^6`.  It is **not** a new
abstract polytope, and the Hermitian cubic is not a new surface.

The new project-level point is the frozen uniform passage

```text
F4 Frobenius graphs -> Hermitian cubic -> 27 divisor walls -> exact P^6,
```

together with its placement as the `q=2` companion of the `q=3` Möbius/Coble
`H11/206` construction.  Worldwide priority for this particular presentation
is left open.

## 1. The five Frobenius intersections

Write `F4=F2[w]/(w^2+w+1)`.  In homogeneous coordinates
`x=(x0:x1)` and `y=(y0:y1)`, the two graph equations are

```text
y0 x1^2 + y1 x0^2 = 0,         x0 y1^2 + x1 y0^2 = 0.
```

Every point of `P1(F4)` satisfies `x^4=x`, so the five displayed points lie on
both graphs.  The intersection product of bidegrees `(2,1)` and `(1,2)` is

```text
(2,1).(1,2)=2*2+1*1=5.
```

They are transverse: in the affine charts the equations are
`y+x^2` and `x+y^2`, whose characteristic-two Jacobian has determinant one;
the same calculation in reciprocal coordinates handles infinity.  Bézout
therefore proves that these are all geometric intersections, each with
multiplicity one.

The two coordinate projections of the five centers are injective.  Exact
`F4` row reduction gives:

```text
every triple:       rank 3 in H^0(O(1,1))^*,
every four-subset:  rank 4 in H^0(O(1,1))^*.
```

Consequently each triple determines one `(1,1)` curve and no such curve
contains four centers.  For all ten triples, the corresponding `2 x 2`
bilinear coefficient matrix has nonzero determinant.  Hence every curve is a
smooth irreducible graph of a projectivity.  Their incidence is the complete
Steiner design `S(3,3,5)`:

```text
10 blocks; point replication 6; pair replication 3; triple replication 1.
```

In particular, the five centers satisfy the standard degree-three del Pezzo
blowup criterion on `P1 x P1`: no two lie on a ruling fibre and no four lie on
one `(1,1)` curve.

## 2. The anticanonical cubic is explicitly Hermitian

Use the divisor basis

```text
F, G, E1, E2, E3, E4, E5,
```

with `F.G=1`, `E_i^2=-1`, and all other basis pairings zero.  Then

```text
K_X = -2F-2G+sum E_i,            K_X^2=3.
```

The five-center evaluation matrix on `H^0(Q,O(2,2))` has exact rank five.
Thus `H^0(X,-K_X)` has dimension four.  With the deterministic kernel basis
frozen in the JSON certificate, substitution into homogeneous polynomials in
`z0,...,z3` gives relation-kernel dimensions

```text
degree 1: 0,       degree 2: 0,       degree 3: 1.
```

The unique cubic relation is

```text
z0^2 z3 + z0 z3^2 + z1^2 z2 + z1 z2^2 = 0.
```

This is a nondegenerate Hermitian form of rank four.  Its gradient is

```text
(z3^2, z2^2, z1^2, z0^2),
```

which has no projective zero even over the algebraic closure.  The cubic is
therefore geometrically smooth.  It is the classical Hermitian cubic surface
`H(3,4)`.  As a numerical cross-check,

```text
#X(F4) = #(P1 x P1)(F4) + 5*(#P1(F4)-1) = 25+20=45.
```

## 3. The 27 curves are exact, not merely an obvious substock

The geometrically realized curves and their divisor classes are

| species | count | class |
|---|---:|---|
| exceptionals | 5 | `E_i` |
| ruling fibres | 10 | `F-E_i`, `G-E_i` |
| projectivity graphs | 10 | `F+G-E_i-E_j-E_k` |
| Frobenius graphs | 2 | `F+2G-sum E_i`, `2F+G-sum E_i` |

Every class has square `-1`, canonical pairing `-1`, and is represented by a
smooth rational curve.

There is also a short exhaustive proof that no other `(-1)`-curve class can
exist.  Contract the six pairwise-disjoint curves

```text
L_i=F-E_i  (1<=i<=5),       L_6=2F+G-sum E_i,
```

and put

```text
H=3F+G-sum E_i.
```

The seven columns `H,L_1,...,L_6` have determinant `+-1` and Gram matrix
`diag(1,-1^6)`, while

```text
K_X=-3H+sum L_i.
```

Write any exceptional lattice class as

```text
D=dH-sum_(i=1)^6 m_i L_i.
```

The equations `D^2=K.D=-1` become

```text
sum m_i^2=d^2+1,              sum m_i=3d-1.
```

Cauchy gives

```text
(3d-1)^2 <= 6(d^2+1),
```

so `d` is only `0`, `1`, or `2`.  Exact integer enumeration gives

```text
d=0: 6,             d=1: 15,             d=2: 6.
```

These are respectively `L_i`, `H-L_i-L_j`, and
`2H-sum_(j != i)L_j`, and the 27 classes agree literally with the 27 curves
above.  This proves exhaustiveness without importing the answer “a cubic has
27 lines.”

Since `-K_X` is ample, the standard del Pezzo cone theorem says that the Mori
cone is generated by these exceptional curves.  Therefore their 27
nonnegative half-spaces give the complete nef cone, not an outer
approximation.

## 4. The complete nef chamber

The Néron--Severi Gram matrix is

```text
U + <-1>^5,
```

of determinant one and signature `(1,6)`.  It is the odd unimodular lattice
`I_(1,6)`.  Between distinct wall curves, the complete pairing histogram is

```text
pairing 0: 216 pairs,          pairing 1: 135 pairs.
```

After negating the lattice form to the `I_(6,1)` wall-normal convention,
pairing zero means orthogonal intersecting mirrors and pairing one means
parallel mirrors.  Hence every actual dihedral angle is `pi/2`: the chamber is
right-angled.

For each wall class `C`,

```text
s_C(x)=x+2(x.C)C
```

is an integral reflection.  The certificate checks all 27 matrices
entrywise: each preserves the Gram matrix, has determinant `-1`, and is the
identity modulo two.

### Exact double description

The 27 integral inequalities were sent directly to `cddexec_gmp`.  Every
returned primitive ray was replayed against all inequalities and classified
by its exact Lorentz square:

| ray kind | count | active facets |
|---|---:|---:|
| finite | 72 | 6 each |
| ideal | 27 | 10 each |
| hyperideal | 0 | -- |

Every one of the 27 inequalities occurs as an active facet.  The ample
witness `-K_X=2F+2G-sum E_i` has square three and pairing one with every wall.
The complete cone therefore proves finite volume and noncompactness; the
absence of hyperideal rays is essential here.

### Literal Vinberg criterion

Independently of cddlib, the complete `27 x 27` Coxeter Gram matrix passes
Vinberg's elliptic-edge endpoint and full affine zero-neighborhood test:

```text
elliptic rank-5 edge diagrams:              648
edges with (finite,ideal) endpoints (1,1):  432
edges with (finite,ideal) endpoints (0,2):  216
bad finite-volume edges:                      0
finite vertices:                             72
ideal vertices:                              27
```

Each ideal vertex has ten incident facets split into five affine two-node
components.  All 135 connected affine-component zero-neighborhood audits are
valid parabolic neighborhoods of rank five.  Thus the literal verdict is

```text
finite-volume noncompact.
```

All 648 edges fail the compact-only endpoint test because every edge has an
ideal endpoint.  This explicitly separates noncompact finite volume from
both compactness and infinite volume.

## 5. Literal identification with `P^6` and `2_21`

In the unimodular marking above, send

```text
D=dH-sum m_i L_i       to       u(D)=(m_1,...,m_6,d)
```

in the sign-reversed Lorentz lattice `I_(6,1)`.  The 27 normals become exactly

```text
(-e_i,0)                                  6 vectors,
(e_i+e_j,1), i<j                         15 vectors,
(sum_(j != i)e_j,2)                       6 vectors.
```

This is an equality of integral normal-vector multisets with Table 2 of
Everitt--Ratcliffe--Tschantz, not just a match of face counts or graph hashes.
Their theorem identifies these mirrors with the sides of the
Potyagailo--Vinberg `P^6`, the Coxeter chamber of the level-two congruence
reflection subgroup of `PO_(6,1)(Z)`.  Their `P^6` census is exactly the
independently recovered `72+27`, and its Euclidean dual is Gosset's `2_21`.

The line-disjointness graph has 216 edges and is the `2_21` one-skeleton; its
complement is the 135-edge intersection graph of the 27 cubic lines.

## 6. The integral-reflection dichotomy in the whole Möbius factory

For a prime power `q`, the Frobenius factory blows up `q^2+1` centers.  A
Möbius circle/block `A` has `q+1` centers and strict-transform class

```text
C_A=F+G-sum_(i in A) E_i.
```

Therefore

```text
C_A^2=2-(q+1)=1-q=-(q-1).
```

Moreover `C_A.E_i=1` for every incident exceptional.  Thus `C_A` is primitive
of divisibility one in the full blowup Néron--Severi lattice.  Its reflection
is

```text
s_C(x)=x + 2(x.C_A) C_A/(q-1),
```

so it is integral on the full Néron--Severi lattice if and only if

```text
q-1 divides 2.
```

Among prime powers, this gives precisely

```text
q=2: norm -1 circles, the Hermitian cubic and P^6;
q=3: norm -2 circles, the Type-MI Coble seed and H11/206.
```

This is a genuine necessary-and-sufficient criterion for the **full-NS
integral-reflection Möbius factory**.  It does not prove that abstract
`q>3` Coxeter polytopes do not exist, nor exclude other lattices, surfaces, or
wall constructions.

## 7. Lineage and novelty boundary

The following are primary lineage sources:

- Potyagailo and Vinberg, *On right-angled reflection groups in hyperbolic
  spaces*, Comment. Math. Helv. 80 (2005), 63--73,
  [doi:10.4171/CMH/4](https://doi.org/10.4171/CMH/4): the finite-volume
  right-angled `P^n` series.
- Everitt, Ratcliffe, and Tschantz, *Right-angled Coxeter polytopes,
  hyperbolic 6-manifolds, and a problem of Siegel*, Math. Ann. 354 (2012),
  871--905, [arXiv:1009.3231](https://arxiv.org/abs/1009.3231): the exact 27
  normal vectors, `72+27` census, congruence-two group, and `2_21` duality.
- Coxeter, *The polytope `2_21` whose twenty-seven vertices correspond to the
  lines on the general cubic surface*, Amer. J. Math. 62 (1940), 457--486,
  [doi:10.2307/2371466](https://doi.org/10.2307/2371466).
- Hirschfeld, *The double-six of lines over `PG(3,4)`*, J. Austral. Math.
  Soc. 4 (1964), 83--89,
  [doi:10.1017/S144678870002276X](https://doi.org/10.1017/S144678870002276X):
  the classical `F4` double-six/27-line surface.

Accordingly:

- new abstract hyperbolic polytope: **no**;
- new algebraic surface: **no**;
- new classical cubic-line/Gosset correspondence: **no**;
- project-new frozen Frobenius/Möbius presentation and `q=2,3` factory
  dichotomy: **yes**;
- worldwide priority for that presentation: **not claimed**.

The result is an exact structural `L3` calibration.  It does not by itself
supply the desired class-wide `POLYTOPE -> SURFACE -> AG pushback` theorem.

## 8. Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_mobius_q2_cubic_p6.py --check-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_mobius_q2_cubic_p6.py
```

Frozen artifacts:

- `scripts/certify_mobius_q2_cubic_p6.py`
- `data/enumerations/mobius_q2_cubic_p6.json`
- `tests/test_mobius_q2_cubic_p6.py`

An independent implementation is frozen in the
[hostile audit](mobius_q2_cubic_p6_hostile_audit.md), with script
`scripts/audit_mobius_q2_cubic_p6_hostile.py`, certificate
`data/enumerations/mobius_q2_cubic_p6_hostile_audit.json`, and focused tests
`tests/test_mobius_q2_cubic_p6_hostile_audit.py`.
