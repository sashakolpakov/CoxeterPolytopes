# A complete H12 source with an actual A12 corner

**The standard odd unimodular lattice `I_(12,1)` has the required actual `A12` corner.** A small exact reconstruction gives its complete 14-facet reflection chamber, with 35 finite vertices and two bounded cusp links. This supplies two explicit intrinsic face-completion normals for the new H20 `E8` face. Their compatibility with other H20 facets is a separate gate.

The [standalone rational checker](../scripts/fresh_h12_a12_source_20261014.py) and [certificate](fresh_h12_a12_source_20261014.json) pass generation and saved replay. No earlier chamber catalog, floating-point acceptance, or arithmetic enumeration in dimension 20 is used.

## The explicit corner and missing facets

Use coordinates `(t,x1,...,x12)` with square `-t^2 + sum(xi^2)`. The initial `B12` sector has eleven roots `-ei+e(i+1)` and the last root `-e12`, all with time coordinate zero. Add

```text
r = (1,1,1,1,0,0,0,0,0,0,0,0,0),  r^2=2,
s = (3,1,1,1,1,1,1,1,1,1,1,1,0),  s^2=2.
```

The eleven norm-two chain roots followed by `s` form `A12`. Their common point is

```text
v = (11,3,3,3,3,3,3,3,3,3,3,3,3),
v^2 = -13,   (v,r)=-2,   (v,-e12)=-3.
```

Thus all other facet inequalities are strict at the corner. Its squared radial-sine distance height from `(1,0,...,0)` is `108/13 < 9`.

Put `T=v/13`, so `T^2=-1/13`, and let `omega_i` be the fundamental weights of the ordered `A12` chain. The two missing facets have the particularly simple expressions

```text
r     = -omega3  + 2T,   norm2,
-e12  = -omega11 + 3T,   norm1,
(r,-e12)=0.
```

They pair minus one with chain positions 3 and 11 respectively, and zero with the other chain roots. Reversing the chain replaces those weights by `omega10` and `omega2`. These formulas are directly transferable into any isometric intrinsic form `Cartan(A12) + <-1/13>`. They do not guarantee admissible pairings with ambient walls having nonzero projection to the complementary `E8` span.

## Completeness and volume checks

A primitive integral reflection normal in the unimodular lattice has norm 1 or 2: reflection integrality makes its norm divide twice every coordinate, and the coordinate gcd is one. The height-zero roots give the full `B12` chamber. At positive time `t`, its dominant vectors satisfy `x1>=...>=x12>=0` and `sum(xi^2)=t^2+q`, with `q` equal to 1 or 2.

The checker enumerates all such integer square partitions through `t^2/q <= 9`. There are only **22 candidates in seven levels**. It records a positive-pairing obstruction for every rejection, verifies that each accepted packet is mutually acute, and accepts exactly `r` at height `1/2` and `s` at height `9/2`. This is a complete bounded Vinberg enumeration, already sufficient for the stated corner because `108/13<9`.

The chamber is also complete globally. On `t=1`, its first twelve inequalities give `x1>=...>=x12>=0`; the `r` inequality gives `x1+x2+x3<=1`. Hence the projective slice is bounded. Exhausting all 91 twelve-row active sets yields **37 distinct extreme points: 35 timelike and two null**, with no spacelike point. All fourteen facets have explicitly checked strict timelike facet points. At each null point, every connected incident block has a positive kernel, a positive-definite codimension-one principal minor, and a positive weighted root sum proportional to that null point. The total Euclidean quotient rank is 11. These identities prove bounded cusp links and finite volume.

No later Vinberg wall can occur: it would have positive time and nonpositive pairings with all fourteen accepted walls. After normalizing its time coordinate to one, it would lie in the projective slice just certified to be contained in the closed future light cone, contradicting its positive norm. This upgrades the small bounded enumeration to the full reflection chamber.

## The lattice glue check

There is also a direct explanation for the corner. Start with `A12 + <-13>`. A first fundamental weight satisfies `omega1^2=12/13`. Adjoining `omega1 - 5v/13` gives a norm-minus-one integral glue vector. The checker realizes it as an actual integral vector in the displayed `I_(12,1)` coordinates and verifies that eleven chain roots, `v`, and this glue vector form an integral basis of determinant plus or minus one. Thus no classification of indefinite unimodular lattices is needed to identify this particular overclosure.

The perpendicular integral lattice is exactly `A12`. In a glue representative with residue `k` modulo 13, orthogonality to `v` gives `5k-13m=0`; since `gcd(5,13)=1`, the residue is zero. This perpendicular lattice has no norm-one vectors and has precisely the `A12` norm-two root system. It therefore agrees with the complete finite stabilizer found above. Changing the sign of `v` gives the equivalent glue formula with `+5v/13`.

This is a ready-made finite-volume **H12** source. It is not an H20 construction, a new arithmetic chamber, or a proof that the two missing facets lift compatibly. The [H20 face-transfer certificate](fresh_h20_contact_closure_20261014_h12_transfer.json) has its own ambient compatibility checks.

```sh
python scripts/fresh_h12_a12_source_20261014.py --verify
```
