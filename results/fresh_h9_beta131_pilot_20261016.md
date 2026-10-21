# A complete norm131 arithmetic H9 pilot

The even lattice

\[
 L=E_8\perp\begin{pmatrix}-2\phi&-1\\-1&6\end{pmatrix},
 \qquad \phi=(1+\sqrt5)/2,
\]

has a complete all-root Vinberg prefix through height20 with12 actual
walls. Its complete cone has22 finite vertices, no ideal vertices, and28
hyperideal extreme rays. An exact strict spacelike escape proves this
prefix has infinite volume. Reflectivity of the full lattice remains open
at this checkpoint; no finite-volume example is claimed.

Generation and saved replay pass for the
[prefix producer](../scripts/fresh_h9_beta131_pilot_20261016.py) and
[data](fresh_h9_beta131_pilot_20261016.json),
[complete cone checker](../scripts/fresh_h9_beta131_pilot_20261016_geometry.py)
and [data](fresh_h9_beta131_pilot_20261016_geometry.json), and
[corner-map screen](../scripts/fresh_h9_beta131_pilot_20261016_maps.py)
and [data](fresh_h9_beta131_pilot_20261016_maps.json).

## Host and duplicate screen

Coordinates are \((x,y_1,\ldots,y_8,t)\). The E8 block has diagonal2 and
product−1 on the zero-based edges
\((0,2),(1,3),(2,3),(3,4),(4,5),(5,6),(6,7)\).
Set \(\delta=1+12\phi\), \(\beta=\phi\delta=12+13\phi\), and
\(H=(6,0^8,1)\). Then \(\det G=-\delta\), \(N(\delta)=-131\),
\(N(\beta)=131\), \(H^2=-6\delta\), and \(GH=(-\delta,0^9)\).
The physical signature is(9,1); the conjugate form is positive definite.
The latter also prevents a nonzero K-rational isotropic vector.

The prime determinant norm131 distinguishes this lattice from the earlier
norm11 and norm55 Q(√5) lattices even up to rational similarity: the ratio
of determinants has nonsquare rational norm, whereas a square in K has
square rational norm. In dimension10 a similarity scalar changes the
determinant by a tenth power. The norm31 and silver hosts are over Q(√2),
and are not identified with this K-lattice. A literal repository search
found no prior norm131 host. This is a scoped duplicate screen, not an
exhaustive classification or originality claim.

## Primitive norms and the zero sector

O=Z[φ] is norm-Euclidean: rounding coefficients leaves |a|,|b|≤1/2 and
|a²+ab−b²|≤3/4<1. Its totally positive units are the squares φ^(2k).
These facts are also bound through the earlier independently reviewed
norm-pool package. For a primitive integral reflector r, evenness gives
r²=2b and reflection integrality gives Gr∈bO^10. Multiplication by the
adjugate and a primitive-coordinate Bezout relation gives b|detG.
The ideal(δ) is prime because its norm is the rational prime131. Positivity
at both embeddings and rescaling by units therefore leave exactly two
norm classes, 2 and2β. Both are searched. Squarefreeness also prevents a
nonprimitive vector in either norm class.

At x=0 a long root would have t and all E8 coordinates divisible by β:
the time covector is−t and E8 is unimodular. Then r/β would be in the even
lattice with norm2/β, impossible. For norm2, positivity at both embeddings
gives |t|,|t′|≤1/√3; its integer field norm then forces t=0. Writing the
E8 part u+φv gives u²+v²=2 and2(u,v)+v²=0. Positive evenness forces v=0
and u an ordinary E8 root. Thus the zero sector is exactly E8 with240
roots and eight simple walls. The scalar coordinate vector has norm6 and
is not a reflector: its reflection coefficient against the time basis
vector is−1/3. No ninth zero wall is inserted.

## Complete bounded enumeration

Put J=G_E8^-1, λ=G_E8 y, and l=x−6t. Then

\[
 l^2+6\lambda^TJ\lambda=6q+\delta x^2,
 \qquad h=\delta x^2/(6q).
\]

l has both signs; dominant λ coordinates are nonnegative. Complete
congruences are l≡x mod6O and all λ,l divisible by b=q/2. Coprimality of β
and6 makes the remaining time reflection congruence automatic; each
resulting vector is also checked directly.

The exact bounds q<67, q′<8, δ>20 and−δ′>6 imply
x²<402T/20, |x′|<3, l²<402(1+T), |l′|<7. Dual Cauchy gives the stronger
coordinate bound3λ_i²≤6q+δx² at both embeddings. Explicit rational
coefficient boxes are followed by exact two-place tests. Schur pruning
forms are recomputed from the literal inverse-energy matrix. Previous-wall
pruning occurs only after signed l is fixed, when all remaining λ
coefficients are nonnegative. Norms, reflection integrality and chronological
acuteness are checked exactly.

T=20 gives52 complete controls and1,573 recursion visits. Four positive
walls join eight zero walls. All66 pairs are Coxeter:54 right, sevenπ/3,
twoπ/5 and three dotted. A strict timelike center and an actual facet
witness for every wall are saved.

## Geometry and the finite symmetry screen

Exact cone enumeration gives all22 finite and28 hyperideal rays. Mixing
a hyperideal ray slightly with the center gives a strict spacelike escape
of norm1136595/1024−(43841/64)φ>0, negative against every wall.

All78 norm-colored outward-sector maps between known finite corners were
checked, using unordered source/target pairs and including self-maps.
Of these,23 are integral isometries; inverses cover reverse orientations.
A common future timelike fixed vector is recorded, of norm(13−9φ)/29<0.
Thus the group generated by these saved maps is finite: it lies in a
compact point stabilizer and a discrete integral isometry group. This does
not exclude maps at later corners or other chamber symmetries.

```sh
python3 scripts/fresh_h9_beta131_pilot_20261016.py --verify
python3 scripts/fresh_h9_beta131_pilot_20261016_geometry.py --verify
python3 scripts/fresh_h9_beta131_pilot_20261016_maps.py --verify
```

All three saved replays pass. Independent natural-order Golden-arithmetic
prefix verification is separate. A later complete height can introduce
new corner sectors and cut this escape; this bounded pilot does not decide
that. All prior frozen obstructions remain unchanged.
