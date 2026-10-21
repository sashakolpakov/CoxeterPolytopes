# Silver p7 Lambda25: the exact first packet is infinite-volume

For the arithmetic host

`<-(2sqrt(2)-1)> + Lambda25`

in hyperbolic dimension 25, the complete prefix through the first positive
packet does **not** close a finite-volume chamber.

## Exact local and chronological screen

The catalogue Gram matrix has determinant `2^26`, Smith form
`1,1,2^22,16`, minimum 8, and 196,656 minimal vectors.  Its complete
primitive mirror alphabet has the ten dyadic valuations
`0,2,3,4,5,6,7,8,9,10`, with an independent optional factor `3+sqrt(2)`
above seven: twenty classes total.  Exact pairing-radical minima for moduli
`1,2,4,8,16` are respectively `8,8,16,16,64`.

There is exactly one height-zero wall.  With the chosen orientation it is
represented by

```
w0=(-2,1,0,1,2,1,1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2),
```

and exact arithmetic gives `w0^2=16` and `G*w0=-8*e_24`.

There are 17 arithmetic controls at nine heights through
`H=5+4sqrt(2)`.  All 14 earlier controls are empty.  At the final height the
norm-one and norm-two controls are empty; the unique realized component is

`r_b=(2+2sqrt(2),(1+sqrt(2))*b)`

with `b^2=8` and `G*b=0 mod 2`.

## Why the full-shell compactness signal was false

The norm-four reflection gate is essential.  The full Lambda25 minimal shell
has 196,656 vectors, but PARI `qfminim` finds only **98,256** minimal vectors
in `S2={b:G*b=0 mod 2}`.  Therefore the first packet has 98,256 walls, not
196,656.  A contact polar computed from the full shell is not a polar of the
Vinberg cone.

Moreover, every reflecting minimal vector has `b_24=0`.  Since
`G*w0=-8*e_24`, all 98,256 first-packet walls are exactly orthogonal to the
height-zero wall.

## Strict hyperideal certificate

Choose the height-zero side `(z,w0)<=0` and put `z=(1,-w0)`.  Then

`z^2=17-2sqrt(2)>0`,

its pairing with the height-zero wall is `-16`, and its pairing with every
first-packet root is `-6-2sqrt(2)<0`.  Thus a strict open hyperideal patch
survives the complete prefix, proving that prefix has infinite hyperbolic
volume.  Later mirrors may still cut it; this is not a global
nonreflectivity result.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_p7_lambda25_first_packet.py
PYTHONPATH=src pytest -q tests/test_silver_p7_lambda25_first_packet.py
```

Source: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/LAMBDA25.html
