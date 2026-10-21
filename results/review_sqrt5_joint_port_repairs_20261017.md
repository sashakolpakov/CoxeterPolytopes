# Independent arithmetic review of the joint-port repair union

**PASS.** The three-, four- and five-port diagnostics retain exactly the
same 2,219 positive face masses, including 355 actual selected dotted-prism
provenances. Their repairs are simultaneous: each stage raises vertex
masses by a maximum over its new constraints, preserving all previous
positive-cost constraints, signed flags and simplex-packing covers.

The [standalone Fraction checker](../scripts/review_sqrt5_joint_port_repairs_20261017.py)
imports no producer or project code and calls no optimizer. Starting from
the [complete restricted local-cone certificate](fresh_sqrt5_restricted_flag_local_cones_20261017.md),
it checks every saved orbit load, reconstructs every per-vertex repair,
and independently checks the positive face coefficients of every exported
violated column. It verifies unchanged face masses and objective numerator,
all successive normalizers, and preservation of all preceding critical
columns.

| Added constraints | Inequalities checked | Violated columns | Repaired ratio |
|---|---:|---:|---:|
| Strengthened three-port | 44,260 | 24 | 12.35526807578 |
| Four-port facet | 90,619 | 47 | 12.35033157101 |
| Fifteen five-port facets | 2,153,595 | 115 | 12.34494473291 |

The final exact ratio is

\[
\frac{66908102699388132000000000000000}
     {5419878674792788039126804769429}>12.
\]

Thus the same restricted obstruction survives the union of these local
constraint families. This is not a numerical optimization result, a complete
joint local cone, primal feasibility, or a dimension theorem.

The geometric proofs and producer sources also received written internal
review. For three ports, triangular adjacency either fills the triple or
forces any two selected squares into one prism. For four and five ports,
the selected dotted-prism footprints are stars; between two triangle blocks,
connected square components are precisely such stars. The five-port domain
has a [separate review](fresh_sqrt5_joint_five_port_domain_20261017_review.md).
The five-port implementation correctly uses exact intersection-size
histograms even for nonmonotone simplex coefficients. All stages average
the whole inequality over complete finite-forest automorphism orbits,
retain actual provenance, and check their binomial coefficient-sum identities.

This independent arithmetic checker does **not** reconstruct forest orbits,
raw prism Grams or literal feature histograms. Those remain explicit inputs
from the separately replayed and source-reviewed
[three-port](fresh_sqrt5_joint_triangle_square_ports_20261017.md),
[four-port](fresh_sqrt5_joint_four_port_diagnostic_20261017.md) and
[five-port](fresh_sqrt5_joint_five_port_diagnostic_20261017.md) producers.
The [saved review](review_sqrt5_joint_port_repairs_20261017.json) records the
scope and source hashes. Generation and fresh saved replay pass.

```sh
python3 scripts/review_sqrt5_joint_port_repairs_20261017.py --verify
```
