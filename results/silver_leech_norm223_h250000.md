# Silver Leech norm-223 necessary prefix through H=250000

This freezes an exact bounded necessary-wall certificate for

`<-a> + Lambda_24`, where `a=-13+14*sqrt(2)` and `N(a)=-223`.

## Exact census

- coefficient factorization: `a=(15+sqrt(2))*(sqrt(2)-1)`;
- normalized local norm classes: `1`, `2`, `15+sqrt(2)`, and
  `30+2*sqrt(2)`;
- the two odd-prime classes use spatial divisor `15+sqrt(2)` and their
  binary target is divided by its square;
- complete arithmetic controls through rational cutoff `H=250000`: `1129`;
- nonempty necessary components: `35` at `24` distinct heights;
- raw binary Gram triples: `2875`;
- short-combination exclusions: `20`;
- surviving necessary triples: `2855`.

The exact serialized control-skeleton digest is
`07d22789e741b59f956e22829fea5fdd51f743ee1e2c5245c000b94a8bd36a0d`.

## Limiting norm-four direction

The tightest necessary component is

- control ordinal `1084`;
- `q=2`;
- `H=114998+81316*sqrt(2)`;
- `x=130+92*sqrt(2)`;
- reduced target `229998+162632*sqrt(2)`;
- `216` surviving binary triples;
- exact minimum-four/3-by-3-Gram maximum `679+479*sqrt(2)`.

It gives

`c_max = 3802/2159 - (698/2159)*sqrt(2)`

and exact Lorentz excess

`122315101/4661281 - (86488302/4661281)*sqrt(2)`,

which is approximately `4.658215156e-4` and remains positive. The smallest
strict rational shrink denominator is `29193`; the witness

`w = 110987984/63027687 - (20376016/63027687)*sqrt(2)`

has exact excess

`34745715955068223/1324163109523323`
` - (73706794102313518/3972489328569969)*sqrt(2)`,

approximately `1.67945985e-9`, and has strictly positive margin against every
component in the prefix.

## Scope

This is a necessary **super-set** of walls. The retained binary Gram triples
are not asserted to be realized Leech pairs, primitive roots, or
Vinberg-accepted walls. In particular, the limiting component is a
necessary-superset limiter/threat, not a realized cutter. The certificate
proves only a strict hyperideal obstruction for the complete bounded prefix;
it is not a compactness, reflectivity, or global nonreflectivity theorem.

Regenerate and verify with:

```text
PYTHONPATH=src python3 scripts/audit_silver_leech_norm223_h250000.py
PYTHONPATH=src pytest -q tests/test_silver_leech_norm223.py
```
