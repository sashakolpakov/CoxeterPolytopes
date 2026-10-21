# Independent geometry replay of the 73-wall H9 construction

The [separate arithmetic replay](fresh_h9_gap_completion_independent_geometry_20261013.json) passes the final case 0 in the [coupled chamber archive](../data/enumerations/fresh_h9_gap_completion_20261013_coupled_geometry.json). It checks all 73 normal norms and Gram identities, their common strict timelike interior, the first strict future spacelike escape, and the actual joint finite vertex involving `coupled13` and `vertex12` at angle pi/3.

Every one of the 73 inequalities is verified at that vertex. Exactly nine walls vanish there, and their physical Gram entries match the saved positive definite rank-nine Gram. Thus the saved point is an actual finite vertex of the chamber. The strict spacelike escape proves that the entire displayed chamber still has infinite volume.

The checker uses a separate recursive exact quadratic-tower implementation. It handles dependent and zero radicals without assuming field independence. Input extraction reuses the construction's pure loader once; the resulting [frozen input](fresh_h9_gap_completion_independent_geometry_20261013_input.json) contains encoded field towers, all normal coordinates, the seed Gram, and source hashes. Replay imports no construction producer. The all-pair Coxeter-angle and complete affine-circuit checks are covered by the construction agent's separate producer-free verifier, not repeated by this narrow geometry audit.

```sh
PYTHONPATH=src:scripts:. python scripts/fresh_h9_gap_completion_independent_geometry_20261013.py --verify
```

The replay takes about two seconds. Regenerate the frozen input only deliberately with `--extract`; ordinary verification checks the existing source hashes and exact archived claims.
