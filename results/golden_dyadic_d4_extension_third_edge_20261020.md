# The third edge: an extra B coordinate, rather than an isolated short wall

20 October 2026. Complete paper interval; no worker, numerical or root
search. This proves an actual finite edge, not a compact example.

Import the lattice and first edge from
`results/golden_dyadic_d4_extension_first_edge_20261020.md`, SHA256
`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`,
and the actual second edge from
`results/golden_dyadic_d4_extension_second_edge_20261020.md`, SHA256
`a894a30bb3b17e63a434049195931e8d6bd51a160b44d85141a820a08dd2e72e`.
Thus `R_n=<-phi>+I_n+D4(O)`, n>=5, with the complete primitive
normalized norm pool1,2,4. Keep the roots a1,a2,a3,r of those notes,
and write s=a4. Put `beta=3phi-4>0`.

## 1. Ray and a finite upper

At the actual vertex S of type `B_n+A2+I2(5)`, release a1. The
whole white diagram is `B_n+A1(a2)+I2(5)(r,s)`. With time coordinate1
the outgoing ray is

    P(t)=(1;0^n;(-1-2t)/3,(-1+t)/3,(-1+t)/3,-1), t>=0,
    P(t)^2=(-beta+2t^2)/3.

The vector

    c=phi*(1;0^n;-1,0,0,-1)

is integral, primitive and of norm1: its D4 sum is -2phi, and
its squared norm is `phi^2(2-phi)=1`. Therefore it is an integral
reflector. It is orthogonal to all white roots, and meets the ray at

    t0=beta/2,
    P(t0)^2=-beta/phi^2<0.                         (1)

This upper is strictly before light because beta<2. It does not
yet assert that c is a SIMPLE wall of the full arrangement.

## 2. Complete finite first-wall profiles

Let a potential first new simple wall have normalized norm q and
coordinates `(a;I;z1,z2,z3,z4)`. The old vertex p belongs to the
chamber closure, so `a>=0`. Let

    epsilon=I^2, b=-w.a2>=0,
    f=-w.r>=0, g=-w.s>=0.

Both old norm4 roots r,s have divisibility2. A finite extension of
their I2(5) component can therefore only be q4, with contact2 to
exactly one endpoint: this gives H3. A q2 contact would have label4
and hence a nonfinite `(4,5)` path. A q4 contact2phi gives a
nonfinite `(5,5)` path; q1 has no integral finite contact. Multiple
contacts create a cycle. Joining such an H3 extension also to a2
would give `(4,3,5)`, which is not finite.

The whole B_n, n>=5, has only its norm2 long-end extension. It may
also connect by a3-bond to a2, but cannot join the I2(5) component
or introduce a5-bond. The complete inventory is therefore:

|role|q|epsilon|b|f,g|
|---|---:|---:|---:|---|
|isolated|1,2,4|0|0|0,0|
|a2 only, label4|1|0|1|0,0|
|a2 only, label3|2|0|1|0,0|
|a2 only, label5|2|0|phi|0,0|
|a2 only, label4|4|0|2|0,0|
|B_n only|2|1|0|0,0|
|B_n through new root to a2|2|1|1|0,0|
|H3 at r|4|0|0|2,0|
|H3 at s|4|0|0|0,2|

There are11 literal profiles. All three white components and all
finite joins are included. The epsilon1 profiles have the specific
long-end coordinate `I=e1` after choosing the usual B_n chamber.

## 3. The I2(5)-orthogonal profiles

First take f=g0. Then `z4=-a`, `z1+z2+z3=-a`. Define

    V=-2z1+z2+z3.

Since `z3-z2=b`, direct reconstruction and the norm equation give

    z1=(-a-V)/3,
    z2=(-2a+V-3b)/6, z3=(-2a+V+3b)/6,
    V^2-2beta*a^2=R,    R=6(q-epsilon)-3b^2.       (2)

The ray product is `(-beta*a+tV)/3`. A new wall has a>0, V>0,
and collision at or before(1) gives V>=2a. Applying both embeddings
to(2) yields the complete necessary strip

    a^2<=R*phi^2/6,        a'^2<=beta*R'/10.       (3)

The R-values and full positive O-lists are

|R|roles|complete a list|
|---|---|---|
|3|q1/a2; B_n/a2 bridge|empty|
|6|isolated q1; B_n extension|phi|
|9|q2/a2, label3|phi|
|9-3phi|q2/a2, label5|empty|
|12|isolated q2; q4/a2|1,phi|
|24|isolated q4|1,phi,phi^2,2phi|

For completeness write `a=A+Bphi`; then
`B=(a-a')/sqrt(5)` is integral. Using `8/5<phi<81/50`:

* R3 gives `a<6/5`, `|a'|<3/5`; hence B0, but integer a1
  fails the conjugate bound.
* R6 gives `a<=phi`, `|a'|<3/4`; hence B0 or1. Integer1 and
  phi-1 fail their conjugate bounds, leaving phi.
* R9 gives `a<2`, `|a'|<9/10`; the same coefficient argument
  leaves only phi.
* R=9-3phi gives `a^2<2`, `a'^2=<(3/(2phi))<1` in(3).
  Thus B0 or1; 1 fails the conjugate bound, phi-1 has conjugate
  absolute value phi, and phi is outside the physical bound.
* R12 gives `a<23/10`, `|a'|<21/20`, hence B0 or1. Only
  1 and phi remain.
* R24 gives `a<=2phi`, `|a'|<3/2`, hence B0,1,2. B0 gives1;
  B1 gives phi or phi^2 (phi-1 fails); B2 gives only2phi, since
  the next smaller possibility2phi-1 has conjugate absolute value
  sqrt(5)>3/2.

All displayed surviving values satisfy the exact bounds(3). This is
a complete coefficient strip, not an imposed height cutoff.

Since `beta*phi^2=2phi-1`, R9 at a=phi gives
`V^2=7+4phi`, of rational norm61, a nonsquare. For R12 the two
values give `V^2=4+6phi=2phi^4` and `10+4phi`. The former is
nonsquare because2 is inert; the latter has nonsquare rational norm124.

For R24, reflection integrality of q4 requires a in2O. Of its four
values only2phi remains. It reconstructs `V=4phi` and the vector
2c, which is imprimitive and therefore supplies no new primitive
normal.

For R6, a=phi reconstructs `V=2phi`, b0, and

    (z1,z2,z3,z4)=(-phi,0,0,-phi).

This gives precisely c in the isolated q1 row, and `v=c+e1` in
the q2 long-end B_n row. Both pass the full D4 sum condition.
No other I2-orthogonal profile survives.

## 4. The two H3 profiles cannot reach this interval

Now put `j=f/(2phi)`, so `z4=j-a`, and use `g=-w.s` to obtain
`z1+z2+z3=g+j-a`. Keep the same V and define

    U=beta*a+4j+g.

The ray product is `(-U+tV)/3`. Completing the square in the norm
equation, for either of the two H3 rows, gives exactly

    beta*V^2-2U^2=24(2phi-3).                     (4)

Indeed the rows are `(j,g)=(phi-1,0)` and `(0,2)`, both q4 with
b=epsilon0. Collision before t0 implies `V>=2U/beta`; hence

    U^2<=28-16phi<4.                             (5)

On the other hand a>=0 and beta>0 give
`U>=4(phi-1)>2` in the first row, or `U>=2` in the second.
This contradicts(5). Thus no H3 attachment was lost in the census.

## 5. Actual simple wall and full endpoint

The norm1 root c and the norm2 root v=c+e1 have the same collision
parameter on the I=0 ray. They are not two new simple walls.
Choose the B_n simples

    beta_i=-e_i+e_(i+1) (i<n), beta_n=-e_n.

Then v has product -1 only with beta1 and extends the whole B_n
to B_(n+1). Moreover

    c=v+beta1+...+beta_n.

Thus c is a nonsimple positive root at the new corner. The complete
finite profile census leaves v as the unique new simple wall, with
all full-lattice coordinates integral and primitive. The actual
endpoint T=P(beta/2) has the whole finite diagram

    B_(n+1) + A1(a2) + I2(5)(r,s).                 (6)

There is no finite enlargement joining the B_(n+1) four-bond to
the5-component. Nor can the orthogonal norm2 root a2 complete the
norm4 I2(5) to H3: a connected3/5 component has a single root
length squareclass, whereas2 and4 differ by a nonsquare in K.
This also directly checks saturation of(6), rather than identifying
the nonsimple c as an extra isolated component.

The result is one complete new actual edge. Although c is a unit
normal, its existence does not make it an actual facet of this full
chamber. Its mirror cuts the full zero arrangement at T and is
subordinate there to v. No full-facet shortcut, period, compactness,
or dimension improvement is inferred.
