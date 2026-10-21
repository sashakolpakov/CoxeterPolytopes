# Independent geometry replay of the changed 73-wall H9 chamber

The [separate arithmetic replay](fresh_h9_vertex_walk_independent_geometry_20261013.json) passes case 4 of the [third vertex-walk archive](../data/enumerations/fresh_h9_vertex_walk_20261013_third_geometry.json), also named in the [final checkpoint](../data/enumerations/fresh_h9_vertex_walk_20261013_final.json). This is the chamber with four replacements, `walk0`, `bridge0`, `bridge2`, and `step0`.

The check covers all 73 normal norms and Gram identities, their common strict timelike interior, the first strict future spacelike escape, and the actual finite vertex incident to `step0`. Every one of the 73 inequalities is checked at that vertex. Exactly nine walls vanish; their physical Gram entries equal the saved positive definite rank-nine Gram. In particular, `step0` meets `seed6` at angle pi/3. The other three saved new vertices are outside this narrow second audit and are covered by the construction agent's producer-free replay.

The exact strict escape has norm

$$
rac{156183485405898214+207283610259155447\phi}{10^{19}}>0
$$

and negative pairing with the common timelike center. Consequently this displayed 73-wall chamber still has infinite volume. Actual finite vertices do not establish global closure.

The checker uses a separate recursive exact quadratic-tower implementation, including dependent or zero radicals. The construction's pure loader is used once to create the [frozen input](fresh_h9_vertex_walk_independent_geometry_20261013_input.json). It contains encoded towers, all normal coordinates, the seed Gram, the vertex order, and hashes of the source archives. Ordinary replay imports no construction producer. Full pair-angle and affine-circuit completeness are covered by the construction agent's separate producer-free verifier; this second audit does not duplicate them.

```sh
PYTHONPATH=src:scripts:. python scripts/fresh_h9_vertex_walk_independent_geometry_20261013.py --verify
```

The replay takes about two seconds. The optional `--extract` flag deliberately regenerates the frozen input.
