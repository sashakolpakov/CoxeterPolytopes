# Complete finite-angle alphabets for the two octic cap recipes

For the two [nine-contact recipes](fresh_h9_after_fourth_rescore_20261015.md), the complete field-theoretically possible finite labels are:

| Pair of normals | Complete possible finite labels |
|---|---|
| Old quartic field only | `2,3,4,5,6,10` |
| First octic cap versus any old wall | `2,3,4,5,6,10,20` |
| Second octic cap versus any old wall | `2,3,4,5,6,10` |
| The two octic caps together | `2,3,4,5,6,10,20` |

In particular, **label 20 must be allowed for the first cap and their mutual pair**. Labels 8, 12, 15, 16, 24 and 30 are excluded here; a degree bound alone would not exclude all of them. This theorem is a necessary field gate. It neither asserts that a particular pair realizes a permitted label nor certifies a cap or finite-volume chamber.

The [stdlib Fraction checker](../scripts/fresh_h9_octic_angle_20261015.py) and [certificate](fresh_h9_octic_angle_20261015.json) verify the exact radicand identities, signs, squareclass tests and bounded totient list used below. Generation and saved replay pass. They import no project arithmetic. The subfield proof is given here rather than represented as a claim that the finite arithmetic checks alone prove it.

Put `F=Q(sqrt(5))`, `phi=(1+sqrt(5))/2`, and

```
D      = 6 + (16/5)sqrt(5),
delta3 = 50 + 22sqrt(5),
delta5 = 245/4 + (111/4)sqrt(5).
```

The inherited field is `K0=F(sqrt(D))`; the cap fields are `L3=K0(sqrt(delta3))` and `L5=K0(sqrt(delta5))`. All physical square roots are positive. For normals `r,s` with positive squares `q_r,q_s`, use the invariant

`A = 4(r.s)^2/(q_r q_s)`.

At a Coxeter intersection with label `m`, this equals `4cos²(pi/m)`. This quotient belongs to the coordinate/form field even when the old wall norms vary and their square roots do not lie in it. One must also check `r.s<=0`; squared equality alone loses the angle's sign. With the sign and a common strict interior, `A=4` is the parallel case and `A>4` is the dotted case. For the finite labels common to both fields, the values of `A` are respectively

```
m = 2,3,4,5,6,10:
A = 0,1,2,(3+sqrt(5))/2,3,(5+sqrt(5))/2.
```

For label 20, `A=2+sqrt(delta3)/(2phi²)`.

To identify every possible finite label, recall that for `m>=3`,
`Q(4cos²(pi/m))=Q(zeta_m+zeta_m^-1)` is totally real and abelian of degree `phi(m)/2`. These facts follow from cyclotomic irreducibility and the fixed field of complex conjugation; see [Conrad, Cyclotomic Extensions](https://kconrad.math.uconn.edu/blurbs/galoistheory/cyclotomic.pdf). The arguments below identify all abelian subfields needed, rather than merely using the degree of the coordinate field.

At the nonidentity embedding of `F`, the signs of `(D,delta3,delta5)` are `(-,+,-)`. Complex conjugation at this embedding flips precisely the imaginary radical generators. Every totally real element is fixed by that conjugation. Conversely, the resulting fixed fields below are totally real because their defining radicands are positive at both embeddings of `F`. Thus the maximal totally real subfields are exactly

```
K0: F,
L3: E3=F(sqrt(delta3)),
L5: E5=F(sqrt(D*delta5)).
```

All degree assertions are nondegenerate. In fact every nonempty product of `D,delta3,delta5` has nonsquare rational `F`-norm, as recorded in the certificate; hence the three squareclasses are independent over `F`. In particular both cap fields have degree eight over `Q`, and their compositum has degree sixteen.

The identity `delta3=phi^4(10+2sqrt(5))` shows
`E3=Q(cos(pi/10))=Q(zeta20+zeta20^-1)`. Its Galois group is cyclic of order four, since `(Z/20Z)^*/{±1}` is cyclic of order four. Its unique quadratic subfield is `F`. For `E5`, set

`R5=D*delta5=(1623+725sqrt(5))/2`, with `Norm_F(R5)=1501=19*79`.

If `F(sqrt(R5))` were Galois over `Q`, the conjugate radicands would define the same quadratic extension of `F`. Thus `R5_conjugate/R5` would be a square in `F`, equivalently `Norm_F(R5)` would be a square in `F`. A rational number is a square in `F` only if it is a rational square or five times a rational square. Neither holds for 1501. Thus `E5` is not Galois and cannot be an abelian quartic field. It also has no quadratic subfield other than `F`: another quadratic subfield would make it a biquadratic, hence Galois, quartic field. Therefore its possible real cyclotomic subfields are only `Q` and `F`.

For completeness, a real cyclotomic field of degree dividing four has label in

```
degree 1: 2,3,4,6;
degree 2: 5,8,10,12;
degree 4: 15,16,20,24,30.
```

The degree-two fields for labels 8 and 12 are `Q(sqrt(2))` and `Q(sqrt(3))`, so neither occurs. The conductor-16 quartic contains `sqrt(2)`; the conductor-24 quartic is biquadratic. Neither can equal `E3`. Finally the conductor-15/30 field is `F(sqrt(30-6sqrt(5)))`. Its radicand divided by `10+2sqrt(5)` is `3/phi²`, which is not a square in `F`; hence it differs from `E3`. Label 20 is present, giving the first two octic alphabets in the table.

The same conclusion covers the mutual pair, whose pairing may lie in the degree-sixteen compositum. Its maximal totally real subfield is

`E=F(sqrt(delta3),sqrt(R5))`, of degree eight over `Q`.

The three quadratic extensions of `F` inside `E` are `E3`, `F(sqrt(R5))`, and `F(sqrt(delta3*R5))`. The latter radicand has norm 120080; neither this number nor 1501 is a square in `F`. Those two quartic fields are therefore not Galois over `Q`, while `E3` is cyclic, not biquadratic. It follows that `F` is the only quadratic subfield of `E`: a different quadratic subfield would generate a biquadratic quartic together with `F`. Moreover, `E` itself is not Galois: an extension of conjugation on `F` would force `sqrt(1501)` into `E`, contradicting that uniqueness. Every abelian quartic field has a quadratic subfield; if contained in `E`, that quadratic subfield must be `F`. Thus the only abelian quartic subfield of `E` is `E3`. An abelian degree-eight subfield would equal `E`, which is impossible. This proves the mutual-pair alphabet without a finite-label cutoff. The potential degree-eight real cyclotomic labels `17,32,34,40,48,60` are excluded by this Galois argument.

The finite totient enumeration is rigorously bounded: `phi(m)^2>=m/2`, since in its prime-factor product every odd-prime factor contributes at least one and the possible factor from 2 contributes at least one half. Thus `phi(m)<=16` implies `m<=512`. The checker enumerates that finite range and verifies the displayed degree lists. No larger label can evade the field-degree and subfield arguments.

```sh
python scripts/fresh_h9_octic_angle_20261015.py
python scripts/fresh_h9_octic_angle_20261015.py --verify
```
