# The tripled-coefficient D9 host is also nonreflective

**The full reflection chamber of `<−3(T−1)>+D9` has infinite volume.** Four complete first-edge censuses reach two full B3+B6 corners. A newly constructed integral isometry maps their oriented simple roots, preserves the actual D9 lattice in both directions, and has trace greater than ten. The connected complete corner sectors therefore give an infinite-order stabilizer of the full chamber.

This is a new calculation in a changed quadratic space, not a transfer of the preceding D9 obstruction. No finite-volume or compact H9 example is claimed.

## Host, primitive norms and zero chamber

```text
K=Q(T), T^4−4T²+2=0, T=2cos(pi/8), O_K=Z[T],
alpha=T−1, beta=3alpha,
L=<−beta>+D9,
D9={y in O_K^9 : sum(y_i) in 2O_K}.
```

Repository search found only the preceding checkpoint's untested proposal for this host. Its determinant is `−12alpha`; the ratio to the old D9 form is three. The [field certificate](fresh_h9_glued_host_20261013_structure.json) proves that three is inert, so this ratio has odd valuation and is not a square. In even rank ten, the determinant squareclass is also unchanged by scalar similarity. Thus the new form is not similar over K to the previous D9 form.

The [structure certificate](fresh_h9_tripled_d9_20261013_structure.json) checks dual exponent twelve. For a primitive reflective normal `r` of norm `q`, reflection integrality gives `2r/q` in the dual and hence `24r/q` in L. Since O_K is a PID, primitivity implies `q|24`. All totally positive units are squares, so the safe complete pool has **26 classes**:

```text
{1,pi,T²,2/pi,2,2pi,2T²,4/pi,4,4pi,4T²,8/pi,8} × {1,3},
pi=2+T.
```

The unique dyadic valuation ranges from zero through twelve; the inert prime-three valuation is zero or one. Every root is tested in the actual D9 module, including spatial parity, division by T and by three, and reflection divisibility against a full lattice basis. No norm class is omitted because it was absent in the old host.

The spatial lattice is unchanged, so the complete zero-height reflection chamber remains B9, with 81 mirrors. Signed permutations preserve D9; total definiteness makes the spatial integral group finite, and the rank-nine finite reflection classification excludes a larger irreducible reflection group containing B9.

## Four complete edges

| Step | Finite profiles | Interval points | Simultaneous first mirrors | Complete endpoint |
|---|---:|---:|---:|---|
| Initial coordinate-axis edge | 27 | 14 | 1 | B8+A1 |
| Axis-chain step 2 | 39 | 24 | 2 | B7+B2 |
| Axis-chain step 3 | 29 | 16 | 2 | B6+B3 |
| Initial A2+B6 edge | 41 | 1,107 | 2 | B3+B6 |
| Total | **136** | **1,161** | | |

The coordinate-axis and chain steps use `E3=K(sqrt(3alpha))`. Its positive relative norm-one generator is `Z3+X3 sqrt(3alpha)`, where

```text
Z3=−29689−16068T+50682T²+27429T³,
X3=−18616−10074T+31780T²+17199T³,
Z3²−3alpha X3²=1.
```

The first axis normal is `(2X3,2Z3,0,...,0)`, of norm four. After each coordinate boost, the next two first facets have transformed form `(X3,−e_k+Z3 e_(k+1))`, of norm two. The boosts are coordinate changes in the surrounding integral module; actual D9 lattice and reflection conditions are rechecked after pulling every candidate back.

The initial A2+B6 edge has its first three spatial coordinates equal. Writing `z` for their sum and `r0` for the fixed profile component gives

```text
z²−9alpha x²=3(q−r0²).
```

Thus `z+3x sqrt(alpha)` uses `E1=K(sqrt(alpha))`. If

```text
(−1+2T²+T³+(2T²+T³)sqrt(alpha))^41 = Z41+X41 sqrt(alpha),
```

the selected first normal is

```text
(2X41/3, (2Z41+2)/3, (2Z41+2)/3, (2Z41−4)/3, 0,...,0),
norm four.
```

The fields were already certified, but their **absolute norm targets and root tests were recomputed for this lattice**. The two [initial-edge](fresh_h9_tripled_d9_20261013_edge1.json) [certificates](fresh_h9_tripled_d9_20261013_edge3.json) and two [chain](fresh_h9_tripled_d9_20261013_step2.json) [certificates](fresh_h9_tripled_d9_20261013_step3.json) preserve the literal PARI output and full fundamental-unit maps. Their relative norm-one kernels have rank one and are proved complete from the full unit bases. Correcting representatives by checked square units and checking exact predecessor/successor inequalities exhausts each finite interval.

All finite-angle profiles use the full field-compatible labels `2,3,4,6,8,16`. No candidate is discarded for incompatibility with an old prefix before first contact. Exact positive-cone identities resolve simultaneous mirrors and identify the full local simple facet. A bounded unit scout supplies an upper endpoint only; it is not the completeness argument. Certified PARI norm/unit completeness remains an explicit dependency.

## Full-chamber transport

The [global certificate](fresh_h9_tripled_d9_20261013_symmetry.json) constructs a new ten-dimensional isometry from the two B3+B6 frames. Old incident index `i` maps to new incident index

```text
[8,1,2,3,4,5,6,7,0][i].
```

There are no sign changes of the primitive normals. The corner-vector scale is positive. The isometry and its inverse are integral in the actual time-plus-D9 basis; the spatial sum condition is checked explicitly. Both finite reflection closures have 45 mirrors. The matrix has coefficients of at most 72 decimal digits, and its exact trace exceeds ten.

The source corner is joined to the initial B9 corner by three complete first steps; the target is joined directly by the new A2+B6 step. The checker reconnects the eight common oriented walls at every step and verifies strict arriving/released-wall signs. These paths place both full local sectors in the same chamber of the full reflection group.

A lattice automorphism normalizes that full root arrangement. Mapping the two complete oriented sectors therefore preserves the same chamber. Trace greater than ten proves infinite order. The integral isometry group is discrete, so infinitely many disjoint translates of a small interior ball lie in the chamber, proving infinite volume. In particular, no reflection subgroup using reflections preserving this exact lattice can have finite covolume. Other lattices and unrestricted H9 constructions are outside this conclusion.

The intermediate construction has thirteen compatible actual facets and five reached finite corners. Its exact strict escape, checked in the [separate geometry replay](fresh_h9_tripled_d9_20261013_chain_independent.json), is

```text
(1,145/100,75/100,50/100,40/100,30/100,20/100,15/100,10/100,5/100),
square=156/25−3T>0.
```

That finite-prefix escape alone would not close the full host; the chamber symmetry does.

## Independent replay

The [initial-edge Fraction checker](../scripts/fresh_h9_tripled_d9_20261013_independent.py) and [chain Fraction checker](../scripts/fresh_h9_tripled_d9_20261013_chain_independent.py) independently regenerate all profiles in reversed order and replay all unit intervals, primal and reflection integrality, simultaneous-root choices, endpoints and combined angle gates. The initial [saved audit](fresh_h9_tripled_d9_20261013_independent.json) and chain audit both pass. Their arithmetic is separate from the producers; the certified PARI norm/unit interface is shared explicitly.

The [independent transport checker](../scripts/fresh_h9_tripled_d9_20261013_transport_independent.py) rebuilds the matrix using Fraction arithmetic rather than importing the producer's K arithmetic. Its [saved result](fresh_h9_tripled_d9_20261013_transport_independent.json) checks both D9 module directions, primitive roots, both 45-mirror closures, all four oriented path connections and the trace: **PASS**.

The [independent mathematical review](fresh_h9_tripled_d9_20261013_transport_review.md) found no gap in the full-corner and chamber-preservation premises. All producer and separate saved-result replays passed after the artifacts were frozen.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013.py --edge 1 --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013.py --edge 3 --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_steps.py --stage 2 --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_steps.py --stage 3 --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_chain_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_tripled_d9_20261013_symmetry.py --verify
python scripts/fresh_h9_tripled_d9_20261013_transport_independent.py --verify
```

Earlier D9 host artifacts remain unchanged. Further blind increases of a cutoff or small changes of this coefficient are not justified by these results; the repeated transport mechanism should be screened structurally first.

The [parameterized transport analysis](fresh_h9_tripled_d9_20261013_mechanism.md) explains both completed obstructions with one explicit Pell-matrix factorization and four D9 congruences, and proposes an E8+A1 spatial pivot with the duplication checks still explicit.
