# Silver Leech norm-223: the H114998 limiter is an accepted root

The complete necessary audit through `H=250000` had a unique limiting
component at

`H=114998+81316*sqrt(2), q=2, x=130+92*sqrt(2)`.

Its 216 binary Gram types have exact norm-four-direction maximum
`679+479*sqrt(2)`. Exactly two types attain it:

- `(115270, 57364, 81316)`, with augmented determinant
  `14`;
- `(115266, 57366, 81316)`, with augmented determinant
  `0`.

Both maxima uniquely have `(d.u,d.v)=(679,479)`. In the first type, putting
`e=u-170d` and `f=v-120d` gives `(e^2,f^2,e.f)=(10,4,6)`, so
`(e-f)^2=2`; this type is impossible in the Leech lattice. In the second,
`(e^2,f^2,e.f)=(6,6,6)`, hence `e=f`.

## Explicit realization

In standard MOG coordinates scaled by `sqrt(8)`, take

```text
d = (4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
e = (-4, -2, 2, 0, 0, 0, 2, 2, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2)
u = 170*d+e
v = 120*d+e
```

Then `d^2=4`, `e^2=6`, `d.e=-1`, and
`(u^2,v^2,u.v)=(115266,57366,81316)`. Thus

`r=(130+92*sqrt(2), u+sqrt(2)*v)`

has square two and realizes the limiting direction pairing. The vector is
primitive: divisibility by the only possible nonunit `sqrt(2)` would force
`u in 2*Lambda`, whereas `d.u=679` is odd. Norm two makes reflection
integrality automatic.

## Chronological acceptance

There are `34` earlier nonempty necessary
components. Writing the spatial part of `r` as

`((679+479*sqrt(2))/4)d + (1+sqrt(2))(e+d/4)`

gives `(e+d/4)^2=23/4`. Exact Schur/Cauchy bounds are strictly negative for
every earlier component. The tight rejection-square certificate is control
ordinal `565` and equals
`57+41*sqrt(2)>0`. The simultaneous norm-one component is empty. The only
possible positive pairings with a distinct same-height norm-two root are
`1` and `sqrt(2)`; their spatial difference targets `2` and
`4-2*sqrt(2)` have empty Leech Gram censuses. Hence the displayed root is
Vinberg-accepted.

## Scope

This upgrades the old necessary-superset threat to a genuine accepted wall,
but it does not produce a finite-volume example. The archived strict shrink
still survives every necessary wall through `H=250000`, with exact Lorentz
excess `34745715955068223/1324163109523323 + (-73706794102313518/3972489328569969)*sqrt(2)>0`.
This is a bounded result, not Vinberg termination or global nonreflectivity.
