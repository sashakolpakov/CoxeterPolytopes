# Cubic-seven unit `D24+`: exact H=100 prefix

For

\[
L=\langle-(\theta-1)\rangle\oplus D_{24}^{+},
\]

the glue-complete Vinberg search through height 100 has 21 control levels
and exactly two nonempty equal-height packets.  The complete canonical
primitive mirror norms are `1,2`.

The spatial enumeration uses the genuine orthonormal `D24+` model.  It
exhausts the D24 root coset and all seven nonzero coefficient residues of the
spinor coset, with dominance `y1>=...>=y23>=abs(y24)`.  No spinor root is
accepted through this bound.

The first packet is the norm-one affine highest root at

```text
x = 2+5theta+2theta^2
h = 5+14theta+6theta^2
y = (2theta+theta^2)(e1+e2).
```

Together with the 24 genuine D24 simple roots it gives the affine alcove:
25 extreme rays, split as 2 finite, 0 ideal, and 23 hyperideal.

The second packet is at

```text
x = 2+8theta+4theta^2
h = 6+18theta+8theta^2,
```

and contains three norm-two, root-coset walls.  Their dominant spatial
multisets are

```text
(1+2theta+theta^2)^2, (2theta+theta^2)^2, 0^20
(2+3theta+theta^2), (1+2theta+theta^2), (1+theta), 0^21
(1+4theta+2theta^2), 1, 0^22.
```

All three survive every earlier wall and are mutually acute.  Exact double
description, with adjacency certified by rank 23 of common active
functionals, gives the packet snapshots:

| walls | rays | finite | ideal | hyperideal |
|---:|---:|---:|---:|---:|
| 25 | 25 | 2 | 0 | 23 |
| 28 | 134 | 4 | 0 | 130 |

Thus the second packet does not approach finite-volume closure; it increases
the hyperideal-ray count.  Every later control through H=100 has zero
dominant mirror candidates.  This remains a bounded-prefix conclusion, not
a global nonreflectivity theorem.

- Backend: `src/coxeter_polytopes/search/cubic7_d24_plus_unit.py`
- Audit: `scripts/audit_cubic7_unit_d24_plus_h100.py`
- Test: `tests/test_cubic7_d24_plus_unit.py`
- Archive: `data/enumerations/cubic7_unit_d24_plus_h100.json`
