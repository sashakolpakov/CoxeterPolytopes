# Cubic-seven Coxeter--Todd-squared host through its first cutter

## Outcome

For

```text
<-2*theta> + K12 + K12  in H24,
```

the point `z=(1,(3/2)h,0)`, `h=G_K12^-1 e1`, is first cut at

```text
H = 33 + 94*theta + 42*theta^2.
```

An explicit primitive accepted norm-two wall pairs with `z` by exactly
`theta>0`.  The replacement `z'=(1,(7/5)h,0)` has positive norm
`196/75 - 2*theta` and is strict against every possible wall
through the same height.  Therefore this complete prefix still has infinite
volume.  This is a bounded obstruction, not a compact example and not a
Vinberg-termination claim.

## Complete chronology

There are 11 locally admissible controls through the
cutting height.  Exactly 4 levels,
containing 5 necessary coefficient
Gram types, survive.  The first three are the previously accepted packets at

```text
5+14theta+6theta^2,
16+45theta+20theta^2,
25+70theta+31theta^2.
```

The last height has two necessary types.  The cutting type is realized by

```text
U = [{'a': '4', 'b': '10', 'c': '4'}, {'a': '0', 'b': '0', 'c': '0'}, {'a': '0', 'b': '0', 'c': '0'}, {'a': '2', 'b': '2', 'c': '0'}, {'a': '0', 'b': '0', 'c': '0'}, {'a': '0', 'b': '0', 'c': '0'}, {'a': '-2', 'b': '-5', 'c': '-2'}, {'a': '1', 'b': '1', 'c': '0'}, {'a': '1', 'b': '4', 'c': '2'}, {'a': '-1', 'b': '-1', 'c': '0'}, {'a': '-1', 'b': '-4', 'c': '-2'}, {'a': '2', 'b': '5', 'c': '2'}]
```

in archived coordinate form, with coefficient Gram

```text
[[12,30,12],[30,84,36],[12,36,16]],
2U-2V+3W=0,
(h.U,h.V,h.W)=(4,10,4).
```

Consequently

```text
r=(1+6theta+3theta^2; U+theta*V+theta^2*W, 0),  r^2=2,
```

and `(z,r)=theta`.

## Chronological acceptance

Complete exact cross-Gram maximization over the minimum-four necessary
supersets gives maximum root pairings

```text
0,
-3theta-2theta^2,
-4-8theta-3theta^2
```

against the three earlier packets.  Thus every earlier wall is nonpositive.
Against the other necessary type at the cutting height, the full relaxed
Schur/divisibility superset has maximum
`-2 - 8*theta - 4*theta^2<0`.
For its own type, a rank-two Schur/divisibility enumeration retains
8173 cross blocks; the
exact integral minimum-four gate retains
6960.  Excluding the identical
root, their maximum is `-4*theta - 2*theta^2<0`.
The cutter is therefore accepted, not merely locally realizable.

## Replacement

The replacement has

```text
(z',z') = 196/75 - 2theta > 0,
(z',r) = -2/5 - (2/5)theta^2 < 0.
```

The exact necessary directional maxima at the three earlier types and the
two cutting-height types are respectively `(0+4theta+2theta^2)`,
`(2+6theta+3theta^2)`, `(2+8theta+4theta^2)`,
`(3+10theta+5theta^2)`, and `(4+10theta+4theta^2)`.  The fourth entry is
the other necessary H33 type, not the explicit cutter type.  Inserting
scale `7/5` gives a strict negative pairing at all five types.  Empty
necessary levels contain no wall.  Hence `z'` survives the complete prefix.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_cubic7_coxeter_todd_squared_h33_cutter.py --verify
PYTHONPATH=src pytest -q tests/test_cubic7_coxeter_todd_squared_h33_cutter.py
```
