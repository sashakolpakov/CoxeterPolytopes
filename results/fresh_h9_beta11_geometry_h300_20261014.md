# Global geometry of the height-300 full reflection prefix

The 35-wall chamber has **880 finite vertices, no ideal vertices and 1,489
hyperideal extreme rays**. A strict future spacelike escape proves infinite
volume. New wall 25 removes the height-100 escape, but the stored replacement
survives all 35 inequalities. Neither escape decides the full lattice's
reflectivity.

The [height-300 root census](fresh_h9_beta11_vinberg_h300_20261014.md) and its
[independent review](fresh_h9_beta11_vinberg_h300_20261014_review.md) certify
both primitive norm classes through the stated bound. A separate exact
double-description calculation reconstructs the complete extreme-ray census
and all 880 finite vertices. It checks their timelike norms, future orientation
and all wall inequalities. The escape is obtained from a hyperideal ray by a
small positive displacement toward a strict interior point; every decisive
inequality is stored explicitly.

See the [geometry script](../scripts/fresh_h9_beta11_geometry_h300_20261014.py)
and [certificate](fresh_h9_beta11_geometry_h300_20261014.json). This computation
uses the project's quadratic arithmetic. It does not independently rerun the
root enumeration, and a complete bounded prefix is not a termination proof.

```sh
python scripts/fresh_h9_beta11_geometry_h300_20261014.py --verify
```
