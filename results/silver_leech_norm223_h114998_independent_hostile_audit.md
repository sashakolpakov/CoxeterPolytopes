# Independent hostile audit: norm-223 H114998 root

The explicit norm-two root at
`H=114998+81316*sqrt(2)` is **chronologically Vinberg-accepted**.

This replay independently rebuilt all `2855`
surviving binary Leech Gram types in the 35 frozen nonempty components and
conservative norm-four pairing maxima (tightening two frozen bounds). The
endpoint has 216 types. Exactly two
attain `679+479*sqrt(2)`: one is impossible because it creates a Leech
norm-two difference; the other is the explicit MOG realization
`u=170d+e`, `v=120d+e`, with `d^2=4`, `e^2=6`, and `d.e=-1`.

For the `2639` types in all
`34` earlier components, a single exact
Schur--Cauchy envelope is strict. Its tight row is control ordinal
`565`:

```text
gap = (63+45*sqrt(2))/4
gap^2 - radicand = 57+41*sqrt(2) > 0.
```

The simultaneous norm-one component is empty. For distinct norm-two roots at
the endpoint, a positive pairing would give `Delta^2<4` at the selected place
and `Delta'^2<7` at the conjugate place, hence `Tr(Delta^2)<11`. But
`Delta=p+sqrt(2)q` has trace `2(p^2+2q^2)`: it is at least 16 if `q!=0`,
while the only nonzero sub-11 case for `q=0` is `p^2=4`, contradicting
`Delta^2<4`. Thus every distinct same-height wall is nonacute.

The root is primitive (`d.u=679` is odd) and its norm two makes reflection
integrality automatic. This upgrades the necessary-prefix limiter to a
genuine accepted wall, but it is not a new finite-volume or compact example.
