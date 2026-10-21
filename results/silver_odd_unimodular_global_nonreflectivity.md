# A Pell-boost global nonreflectivity criterion and silver family

## General criterion

Let `M=<-a>+L` over a real quadratic integer ring and let `d in L` have
square `m`.  Suppose integral `A,U` satisfy

`A^2-a*m*U^2=1`, `a*U` is integral, and `(A-1)/m` is integral.

Then

`g(x,y)=(A*x+U(d,y), y+((a*U)x+((A-1)/m)(d,y))d)`

is an integral Lorentz isometry with integral inverse.  Any oriented mirror
separating `(1,0)` from its image has height at most

`H_sep=A^2-1`.

For a spacelike test ray `z=(1,t*d)`, the same endpoint inequality supplies
an explicit positive lower bound for `(z,r)`.  Consequently a complete
strict-ray certificate through `H_sep`, together with a complete mirror-norm
alphabet excluding time-zero mirrors, proves that `g` stabilizes the full
reflection chamber.  Since `A>1`, it is an infinite-order hyperbolic
translation; the chamber has infinite volume and the reflection subgroup
has infinite index.  More explicitly, if `Gamma=O^+(M)` and `W` is generated
by all integral reflections, then `W` is normal because `Gamma` permutes all
reflective mirrors.  The `W`-stabilizer of a chamber is trivial, so the powers
of the infinite-order `g` inject into `Gamma/W`.  Thus `[Gamma:W]=infinity`.
For an admissible host `Gamma` is an arithmetic lattice, and infinite index
forces the `W`-chamber to have infinite volume.

## Universal silver theorem

Take

`a=2sqrt(2)-1`, `m=3`, `A=16+12sqrt(2)`, `U=7+5sqrt(2)`.

Then

- `aU=13+9sqrt(2)`;
- `(A-1)/3=5+4sqrt(2)`;
- separating slope `(d,y)/x >= 8-4sqrt(2)`;
- `x^2/q <= 297+210sqrt(2)`;
- `H_sep=543+384sqrt(2)=A^2-1`.

For `z=(1,(43/55)d)`, every separator would satisfy

`(z,r) >= ((399-282sqrt(2))/55)x > 0`.

The complete six-class reflection arithmetic has no time-zero mirror for
any positive-definite integral unimodular `L` of minimum three.  The exact
rank-independent frontier through the separator cutoff has
139 controls, 297
necessary binary Gram types and 26 nonempty
heights.  Every augmented Gram pairing is strict.  The worst row remains

`H=455+322sqrt(2), q=2, x=16+11sqrt(2), Gram=(460,226,322)`,

with maximum `37+26sqrt(2)` and margin
`(-51+37sqrt(2))/55>0`.

Therefore:

> For every positive-definite integral unimodular lattice `L` with
> `min(L)=3`, the arithmetic host
> `<-(2sqrt(2)-1)>+L` has an infinite-volume complete reflection chamber;
> its full integral reflection subgroup has infinite index.

This statement is independent of rank and makes no shell-transitivity,
kissing-number, or automorphism-orbit assumption.

At the nonidentity embedding, `-(2sqrt(2)-1)` and the spatial form are both
positive definite.  Hence the host is anisotropic, `Gamma` is cocompact, and
the group-theoretic infinite-index conclusion above applies without a cusp
or finite-volume caveat.

## Published-classification corollary

Combining the universal theorem with the published complete classifications
gives a finite-slice consequence.  The numbers of rootless odd unimodular
lattices in ranks `23,24,25,26,27,28,29` are respectively

`1, 1, 0, 1, 3, 38, 10092`,

for a total of **10,136**.  Every associated orthogonal silver host
`<-(2sqrt(2)-1)>+L` is nonreflective by the theorem above.  This eliminates
the complete rootless odd-unimodular orthogonal slice in hyperbolic
dimensions `23,24,26,27,28,29` (rank 25 is vacuous).  It is emphatically
not whole-dimension or field-wide nonexistence.

The count corollary uses [Borcherds, rank 26](https://arxiv.org/abs/math/9912236),
[King/Bacher--Venkov, ranks 27--28](https://arxiv.org/abs/math/0012231), and
[Allombert--Chenevier, ranks 28--29](https://arxiv.org/abs/2410.19569).
The core Pell-boost theorem does not depend on those classification counts.

The explicitly audited catalogue instances are `S26`, `T27`, `U28`, and
`dim29odd`, pinned by `data/enumerations/silver_s26_h26_first_packet.json` and
`data/enumerations/silver_odd_unimodular_h27_h29_h1000.json`.  These names are examples,
not an asserted classification of every minimum-three unimodular lattice.
The theorem is not field-wide nonexistence over `Q(sqrt(2))`.

```sh
PYTHONPATH=src python3 scripts/audit_silver_odd_unimodular_global_nonreflectivity.py --freeze
PYTHONPATH=src python3 scripts/audit_silver_odd_unimodular_global_nonreflectivity.py --verify
PYTHONPATH=src pytest -q tests/test_pell_boost_nonreflectivity.py tests/test_silver_odd_unimodular_global_nonreflectivity.py
```
