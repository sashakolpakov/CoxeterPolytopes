# Silver p7 T26: exact full-shell polar leaves a hyperideal patch

For the arithmetic host

`<-(2sqrt(2)-1)> + T26`

in hyperbolic dimension 26, the complete first Vinberg packet has infinite
hyperbolic volume.  This replaces the earlier, weak coordinate-ray estimate
by an exact vertex of the complete 117,936-wall shell polar.

## Complete first packet

The primary catalogue Gram matrix is even of rank 26, determinant 3, Smith
form `1^25,3`, minimum 4, and kissing number 117,936.  Its pairing radical
modulo three has minimum 24 and kissing number 1,638.

The complete primitive mirror-norm alphabet has eight classes: independent
optional dyadic, coefficient-prime-seven, and discriminant-three factors.
There are ten arithmetic controls at seven heights through
`H=5+4sqrt(2)`; all eight earlier controls are empty, and there are no
height-zero walls.  The unique realized component is

`r_b=(2+sqrt(2),(1+sqrt(2))*b),  b^2=4 in T26`,

of norm 2.  Hence the first packet is the complete 117,936-vector minimal
shell, with normalized walls `(sqrt(2),b)`.

## Exact 66-contact polar vertex

Let `y` be `1/30` times

```
(43,119,50,-11,1,9,60,5,-27,-3,-20,30,0,-18,-43,-35,
 6,30,4,-11,2,-1,-16,30,-19,14).
```

Exact arithmetic gives `y^2=97/150`.  A full PARI enumeration of all 58,968
minimal representatives modulo sign proves

`max_{b^2=4}|(y,b)|=1`,

with exactly 66 contacting sign pairs.  The canonical representative stream
has SHA-256
`69cdb97ab9d5b89e24db1b1fe94b83f5488e68f14e71e65cf04a8d991b8c4870`.

The boundary polar vector has spatial part `(4-sqrt(2))*y` and Lorentz norm

`316/25-(538/75)sqrt(2)>0`.

After shrinking by `999/1000`, the strict witness

`z=(1,(999/1000)*(4-sqrt(2))*y)`

still has norm

`315418291/25000000-(44768699/6250000)sqrt(2)>0`,

while its maximum pairing with every first-packet wall is

`-1/500-(3/1000)sqrt(2)<0`.

Thus a strict open hyperideal patch survives the complete first packet.  This
is a bounded-prefix statement only: later mirrors may cut it.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_t26_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_t26_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb4.html
