# Silver S26 host in H26: exact first packet

## Outcome

The genuinely new rootless host

\[
 L=\langle-(2\sqrt2-1)\rangle\oplus S_{26}
 \quad\text{in }\mathbf H^{26}
\]

has no height-zero walls. Its complete first nonempty Vinberg packet is at
`H=5 + 4*sqrt(2)` and contains **102,180** primitive
walls. An exact strict hyperideal point survives the whole packet and, after
replacement, every possible wall through `H<=1000`. This bounded prefix has
infinite volume and is not a compact example. Walls above the frozen cutoff
may cut the point; this is not a nonreflectivity theorem.

## Exact host and local arithmetic

The spatial lattice is Conway--Borcherds' odd unimodular rank-26 lattice
`S26`, of minimum three and kissing number 3,120. The catalogue lower-triangle
Gram matrix is stored literally; exact Bareiss elimination gives determinant
one. PARI independently recovers cumulative signed counts 3,120 and 105,300
through norms three and four, hence the exact norm-four shell has 102,180
vectors.

For `rho=2sqrt(2)-1`, the primitive pairing-ideal argument leaves exactly six
norm classes

`(1, 2+sqrt(2), or 2) * (1 or 3+sqrt(2))`.

The second factor is accompanied by the indispensable spatial divisibility
`y in rho*S26`. All six height-zero classes are empty.

## Complete chronology through the first packet

There are exactly **12** arithmetic controls at
**9** distinct heights through the first packet.
All **10** earlier controls have empty exact
minimum-three Gram censuses.

At `H=5 + 4*sqrt(2)` the sole surviving component is

\[
 q=2,\quad x=2+\sqrt2,\quad
 y=(1+\sqrt2)d,\qquad d^2=4.
\]

The Gram triple `(u^2,v^2,u.v)=(4,4,4)` forces this scalar form, so every
norm-four vector occurs and the packet has exactly 102,180 walls.

## Exact infinite-volume certificate

Fix a minimal vector `e` of `S26`, so `e^2=3`, and put

\[
 z=(1,e).
\]

For every `d^2=4`, integrality and minimum three applied to `e+-d` give
`|e.d|<=2`. After dividing a
packet root by `1+sqrt(2)`, its section inequality is
`z.d <= 4-sqrt(2)`, and

\[
 |e.d|\le2<4-\sqrt2.
\]

The wall margin is `2 - 1*sqrt(2)`. Meanwhile

\[
 z^2=3,\qquad z^2-(2\sqrt2-1)
 =4 - 2*sqrt(2)>0.
\]

Thus `z` is strictly spacelike and strictly inside every first-packet slab.
The segment from the Vinberg basepoint to `z` crosses the light cone away
from every facet, producing an open ideal patch. The complete first-packet
chamber has infinite volume.

## Persistence through every possible wall at H<=1000

The fixed minimal direction can be continued much farther. Put

\[
 z=(1,(43/55)e),\qquad e^2=3.
\]

The exhaustive arithmetic census contains **134** controls
at **104** distinct heights. Exactly
**30** controls at **24**
heights retain a necessary binary Gram type, for
**256** types in total.

For each type `(u^2,v^2,u.v)`, the audit maximises `e.(u+sqrt(2)v)` over
every integral augmented Gram that is positive semidefinite and has no norm
one or two among combinations in `[-3,3]^3`. This local test is explicitly
only a necessary gate, not a proof of true ternary minimum. Retaining the
larger superset makes the resulting upper bound valid for every actual S26
root.

The tight row occurs at

\[
 H=455 + 322*sqrt(2),\quad q=2,\quad
 x=16 + 11*sqrt(2),\quad
 (u^2,v^2,u.v)=(460,226,322).
\]

Its maximum is `37+26sqrt(2)`, attained at `(e.u,e.v)=(37,26)`, while the
normalised wall bound is `28+21sqrt(2)`. The exact remaining margin is

\[
 -51/55 + 37/55*sqrt(2)>0.
\]

Every other row has larger margin. The witness remains spacelike by

\[
 3(43/55)^2-(2\sqrt2-1)
 =8572/3025 - 2*sqrt(2)>0.
\]

Consequently an open hyperideal patch survives every possible wall, hence
every Vinberg-accepted wall, through the complete `H<=1000` prefix.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_s26_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_s26_first_packet.py
```
