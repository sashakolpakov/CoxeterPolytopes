# Global geometry of the height-100 full reflection prefix

Exact cone enumeration gives **260 finite vertices, no ideal vertices and 563
hyperideal extreme rays** for the 21-wall chamber. A stored future spacelike
vector satisfies every wall inequality strictly, so this chamber has infinite
volume. This does not decide reflectivity: later walls can remove the escape.

The [bounded Vinberg enumeration](fresh_h9_beta11_vinberg_20261014.md)
includes both primitive norm classes. Geometry is computed by a separate exact
incremental double-description module. The certificate reconstructs every
finite vertex, checks its timelike norm and wall inequalities, and records the
complete extreme-ray incidence census. Its selected hyperideal ray is not
simple; nine independent active walls determine it, and all its additional
active walls are checked explicitly.

The [script](../scripts/fresh_h9_beta11_geometry_20261014.py) and
[certificate](fresh_h9_beta11_geometry_20261014.json) pass generation and saved
replay. Geometry uses the project's exact quadratic arithmetic; this is not a
second independent arithmetic implementation. Completeness of the root prefix
has its own [independent review](fresh_h9_beta11_vinberg_20261014_review.md).

```sh
python scripts/fresh_h9_beta11_geometry_20261014.py --verify
```
