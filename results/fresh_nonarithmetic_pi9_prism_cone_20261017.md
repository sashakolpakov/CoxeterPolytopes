# Complete cone of the62-wall chamber

The [exact prism construction](fresh_nonarithmetic_pi9_prism_completion_20261017.md)
has a complete ambient cone with **193 finite vertices, no ideal vertices,
and2765 exterior projective vertices**. An explicit strict spacelike escape
still proves infinite volume. The compact triangular prism is a whole
three-dimensional face, not a compactness claim for the chamber.

The [source](../scripts/fresh_nonarithmetic_pi9_prism_cone_20261017.py) and
[certificate](fresh_nonarithmetic_pi9_prism_cone_20261017.json) start from the
frozen complete60-wall cone. All2868 preceding rays are converted to the
exact computational root basis and reconstructed from nine independent
active equations. This reduces interval cancellation without changing their
exact oriented rays. The reconstruction stream is hashed in the certificate.

Clipping by60 removes35 exterior rays and creates85 rays, including two
finite vertices. Clipping by61 removes45 exterior rays and creates85 rays,
including three finite vertices. Each crossing is certified by an independent
eight-normal minor. Every resulting ray is reconstructed from its active
equations; all183396 final wall inequalities and every causal sign are checked.
All188 preceding finite vertices persist.

Complete active-ray sets identify the retained compact tetrahedron and
quadrilateral, with four vertices each, and the new compact triangular prism,
with six vertices. The original S6=[0,1,2,3,5,6] face has19 finite and nine
exterior projective vertices. Its complete edge inventory is in the
[boundary ledger](fresh_nonarithmetic_pi9_prism_edge_ledger_20261017.md).

The renewed escape is a positive perturbation of the exterior ray supported
on S6+[7,33,49], with coefficient2^-17 toward the common strict center.
Its norm is strictly positive, its future pairing is negative, and all62
wall products are strictly negative. These inequalities establish infinite
volume directly.

Generation and a fresh saved replay passed. A separate agent's read-only
mathematical/source review passed; the arithmetic uses the inherited exact
interval engine and is not an independent reconstruction of all root recipes.

```sh
python3 scripts/fresh_nonarithmetic_pi9_prism_cone_20261017.py --verify
```
