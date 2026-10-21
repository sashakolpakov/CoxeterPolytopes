# Record 14 ridge-engine capacity hazard

The earlier Record-14 ridge-aware C++ attempt is **not mathematical evidence**.
Its reused candidate bitset had thirteen 64-bit words, hence capacity 832, while
the five compact `a=1` source branches have candidate counts

```text
994, 994, 784, 994, 994.
```

Thus four branches wrote outside the bitset during table initialization.  A
crash, kill, empty transcript, or apparent search termination from that build
must not be read as nonexistence.

The Record-18 engine remains separate.  Record 14 now has an isolated
sixteen-word fork (capacity 1,024) with a runtime assertion
`candidate_count <= capacity` before any candidate bit is set.  A regression
test verifies that the thirteen-word guard rejects exactly source orbits
`0,1,3,4`, admits the 784-candidate orbit `2`, and that the checked
sixteen-word engine admits all five.  It also verifies that the isolated fork
differs from the shared engine only by its provenance header, word count, and
capacity guard.

Frozen files:

- `scripts/explore_h9_14_record14_ridge_aware_engine.cpp`
- `scripts/explore_h9_14_record14_ridge_aware_additions.cpp`
- `scripts/audit_h9_14_record14_ridge_capacity.py`
- `data/enumerations/h9_14_record14_ridge_capacity_hazard.json`
- `tests/test_h9_14_record14_a1_symmetry.py`

This is an implementation-safety certificate only.  It proves neither the
existence nor the nonexistence of a compact or finite-volume hyperbolic
Coxeter polytope.
