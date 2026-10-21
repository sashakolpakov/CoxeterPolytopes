# Independent geometry replay of the 79-wall split-seed extension

**PASS.** A separate exact quadratic-tower implementation checks all 79 normal–Gram identities and norms, the common strict timelike interior, the first strict future spacelike escape, and the actual finite vertex defined using `second33`. That vertex satisfies all 79 inequalities, and its nine physical incident normals have the archived positive definite Gram matrix.

The source is [the frozen 79-wall geometry](../data/enumerations/fresh_h9_split_second_frontier_20261013_geometry.json). Pure loading was used once to create a [frozen input](fresh_h9_split_second_frontier_independent_geometry_20261013_input.json); the replay checks its source hashes and imports no construction producer. It uses a separately implemented recursive field multiplication and exact sign comparison, including dependent radicals.

The strict escape has norm

\[
\frac{11764468488123570049-6234856568308221881\phi}
{16000000000000000000}>0
\]

and center pairing $-870214989/250000000<0$. Every wall inequality is strict. This independently proves that the displayed 79-wall cone has infinite volume, despite its certified finite vertices.

The [audit result](fresh_h9_split_second_frontier_independent_geometry_20261013.json) records all checks. Full cross-pair angle gates and completeness of the affine census remain the scope of the construction agent's separate producer-free verifier; this narrow second audit does not duplicate them. Given acuteness, projection of the strict interior to any wall preserves the other strict inequalities and stays timelike, so each retained wall is an actual facet.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_split_second_frontier_independent_geometry_20261013.py --verify
```

Checker: [independent geometry](../scripts/fresh_h9_split_second_frontier_independent_geometry_20261013.py). No numerical geometry is used in replay. This result does not exclude further additions or replacements.
