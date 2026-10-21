# Structural triage of the proposed E8+A1 host

Continuation: [two complete first edges are now certified](fresh_h9_e8_first_20261013.md). The full indefinite conjugacy question remains unresolved. The historical triage below predates that positive-height search.

The proposed lattice

`M = <-3(T−1)> + E8 + <1>`, where `T^4−4T²+2=0`,

has different lattice-preserving mirrors from the [closed tripled D9 host](fresh_h9_tripled_d9_20261013.md) in their shared rational coordinates. The saved infinite chamber symmetry does **not** preserve M. This prevents transferring that particular proof. An unrelated integral isometry or conjugacy of their full indefinite reflection groups has not been ruled out; no new full-group coverage or positive-height search is claimed here.

The [exact structural certificate](fresh_h9_tripled_d9_20261013_e8_triage.json) and [checker](../scripts/fresh_h9_tripled_d9_20261013_e8_triage.py) give an explicit E8 basis in the first eight Euclidean coordinates. Its Gram matrix is the E8 Cartan matrix, with determinant one. In the ambient form `diag(-3(T−1),1,...,1)`, M has determinant `-3(T−1)` and dual exponent 3. The previously certified PID and totally positive unit-square results therefore give the safe primitive norm pool `q | 6`: ten dyadic/3 classes. This is a safe pool, not a claim that all ten occur.

The full zero-height reflection chamber is E8+A1, with 121 mirrors. Its stabilizer is finite because the spatial form is positive definite at every real embedding and the integral group is discrete in the product of those compact orthogonal groups. It contains W(E8) and the last coordinate reflection. An irreducible finite rank-nine reflection group would have type A9, B9 or D9, each of order smaller than W(E8); hence none can contain it. The reducible alternative is exactly E8+A1. This differs from the B9 zero-height chamber of D9, but only distinguishes these pointed descriptions.

Both M and `M0 = <-3(T−1)> + O^9` are maximal integral lattices: their determinant ideal is the inert prime (3), of valuation one, whereas a proper integral overclosure requires the square of a nontrivial index ideal to divide the determinant. They are distinct neighbors in the shared coordinates. Their intersection is ` <-3(T−1)> + D8 + <1>`, of index ideal (2), cardinality 16, in each. Thus saturation does not identify these two displayed lattices. This argument does not classify all maximal lattices or identify their full automorphism groups.

Two elementary witnesses show that the reflection sets differ in both directions:

- The E8 root `h=(-1/2,...,-1/2,0)` has norm 2 and preserves M. Its reflection sends the D9 vector `e1+e2` to `e1+e2+h`, which has half-integral coordinates and is not in D9.
- The D9 mirror with primitive normal `2e1`, of norm 4, flips one coordinate. It sends h outside E8 and therefore does not preserve M.

In addition, conjugating the saved tripled-D9 infinite transport into the displayed E8 basis produces two nonintegral matrix entries; its inverse also has two. These exact entries are archived. The failed transport is useful structural evidence, not a nonreflectivity or reflectivity verdict for M.

The structural replay passes:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_e8_triage.py --verify
```

No expensive edge search was started. Any continuation should first address possible indefinite lattice or full-group conjugacy, rather than treating a changed zero-height diagram alone as proof of a new arithmetic group.
