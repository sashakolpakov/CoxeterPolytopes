# RC7 fraction-state regression audit

This archive compares the same pinned AlVin commit
`d92006a842cb4668e7c89e8c8e2d30337e4b0f26` before and after one RC7
state-restoration repair.  The repaired image is
`coxeter-alvin:rc7-fractionfix-d92006a`, immutable image ID
`sha256:45b0180c691d5d07958b9e7c77163ace4aa3db4ba963faf56ad3af3b75b7043c`.

## Source diagnosis

`RCyclotomic7Integer_VFs::computeNextAlVinFractions()` advances `rciE` twice
by `conjugate(2)` to obtain the two nonidentity real embeddings.  The pinned
source restored `rciMinusAlpha0` but did not restore `rciE`.  It consequently
formed and stored the Vinberg fraction with `sigma^2(norm)` and also left the
shared entry of `aiPossibleNorms2` rotated for the next batch.  For the
mandatory example, the normalized missing wall has

```text
x0   = [0,1,-1]
root = ([0,1,-1], 0, [-1,0,-1])
norm = [-6,-4,-4]
```

but the unpatched program never calls `findVector(x0,norm)` for this pair.
The repair in
`../../../../tools/alvin/alvin-rc7-fraction-conjugation.patch` adds the
missing third `rciE->conjugate(2)` before fraction construction.

## Exact source/basis reconciliation

AlVin and the two theses use
`lambda_j=2*cos(2*j*pi/7)`.  With `theta=lambda_1`, conversion to the exact
power basis is

```text
[u,v,w] = (w-2v) + (u-w)*theta + (v-w)*theta^2.
```

Guglielmetti's thesis, page 132, gives the admissible form exactly as
`[-1,0,0],1,[-2,-5,-7]`.  Bottinelli's thesis, page 84, prints `-phi_0` rather
than `-2phi_0` in the last diagonal coefficient; that displayed coefficient
has a negative conjugate and is therefore a typo.  Bottinelli's listed roots
are exact reflecting vectors for Guglielmetti's form and not for the printed
nonadmissible coefficient.

Bottinelli's fourth root

```text
([-3,-1,-3], 0, [-4,-1,-3])
```

is multiplied by the totally positive norm-one unit
`5-theta-2*theta^2=[-2,-3,-1]` to obtain the normalized missing wall above.
Together with the two cone roots and the third Bottinelli root, it gives the
four-wall compact chamber.  The exact scale-invariant product on the last
edge is

```text
4*cos(pi/14)^2 = 1 + theta + theta^2 = [-2,-2,-3],
```

so AlVin's reduced numerator/denominator pair is
`[-2,-2,-3]/4`; this independently confirms the order-14 recognition patch.

## Before/after evidence

The unpatched capture is
`../pinned_rc7_d92006a_20260812/hostile-m14`.  It incorrectly terminates with
five walls and no order-14 edge.  Its stdout SHA-256 is
`5034c4822dac5c7016bc316f697c517fc1a213a86410d893a343c302125c7d7c`.

The repaired `hostile-m14` capture in this directory terminates with four
walls and Coxeter matrix

```text
1  2  1   2
2  1  2   1
1  2  1  14
2  1 14   1
```

The independent exact verifier proves that the supplied polygon is compact,
has four finite vertices and no cusps, and that all reflections preserve the
lattice.  `hostile-m14/result.json` has `acceptance_passes: true`; its stdout
SHA-256 is
`a21e1e4dcce55c23b3a9c9d628ed863cfb98d53d35be10ffd3a4e72a22d377ab`.

The repaired `control` capture also has `acceptance_passes: true`: it retains
the known eight-wall compact chamber, with 16 finite vertices and no cusps.

The exact in-tree regression is
`tests/test_cubic7_alvin_verifier.py::test_vinny_rc7_four_wall_order_14_regression_is_exact_and_compact`.
The complete verifier file passes with `9 passed`.
