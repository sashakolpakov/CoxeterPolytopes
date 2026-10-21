# Silver p7 Lambda22: exact first packet is infinite-volume

For the arithmetic host

`<-(2sqrt(2)-1)> + Lambda22`

in hyperbolic dimension 22, the complete first nonempty Vinberg packet does
**not** close a finite-volume chamber.  An exact strict hyperideal witness
survives every wall in the packet.

## Complete local and chronological screen

The Smith form of the authoritative Nebe--Sloane `Lambda22` Gram matrix is
`1^20,2,6`.  Combining its dyadic valuations `0,2,3,4`, the optional inert
prime-three factor, and the optional prime-seven factor `3+sqrt(2)` gives
exactly 16 primitive mirror-norm classes.  For a prime-seven class the
spatial vector is first divided by `3+sqrt(2)`; the dyadic and ternary
pairing gates are then applied to the quotient.

There are 16 reflective controls at 10 distinct heights through
`H=5+4sqrt(2)`.  All 13 earlier controls are empty after exact target
division, minimum-four Gram enumeration, and the pairing-radical minima
`min(S2)=12`, `min(S3)=24`, `min(S6)=96`.  There are no height-zero walls.

At `H=5+4sqrt(2)`, the simultaneous norm-1 and norm-4 controls are empty.
The unique primitive component is

`r_d=(2+sqrt(2),(1+sqrt(2))*d),  d^2=4 in Lambda22`,

of norm 2.  PARI `qfminim` enumerates the full signed minimal shell, so this
packet contains exactly **49,896 walls**.  Distinct minimal vectors satisfy
`d.e<=2`, which makes all distinct same-height wall pairings strictly
negative.

## Strict hyperideal certificate

Let

```
k=(-4490,27287,17909,13568,8412,18464,12340,6324,3634,1942,
   2575,1337,544,272,408,214,32,64,38,22,12,8).
```

For the vendored Gram matrix, exact integer arithmetic gives
`G*k=3*e_22` and `k^2=24`.  The complete minimal-shell enumeration gives
`|d_22|<=2`.  Therefore

`z=(1,(4-sqrt(2))*k/8)`

has

`z^2=31/4-5sqrt(2)>0`,

while its largest pairing with any of the 49,896 walls is

`-(2+3sqrt(2))/4<0`.

Thus the complete first-packet prefix has infinite hyperbolic volume.  This
is not an eventual nonreflectivity result: a later mirror may cut the
displayed witness, and no global compactness or reflectivity claim is made.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_p7_lambda22_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_p7_lambda22_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA22.html
