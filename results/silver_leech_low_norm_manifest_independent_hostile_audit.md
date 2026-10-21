# Independent hostile audit of the silver Leech norm-128 manifest

## Verdict

**REPAIRED BOUNDED-PREFIX CERTIFICATE.**  The arithmetic census and every
archived exact number pass an implementation-independent replay: 32
prime-ideal atoms give 57 squarefree coefficient-ideal directions, all 57
first necessary packets agree field by field, and all 57 strict
fixed-direction witnesses agree field by field.

The originally frozen theorem did **not** pass.  A hyperideal patch surviving
the first necessary packet proves only that this finite prefix has infinite
volume.  Later-height mirrors can cut that patch.  The earlier claims that
there was “no finite-volume survivor” or that all 57 final chambers were
non-finite-volume were therefore invalid scope inflation.  They have been
removed.  The repaired theorem is:

> For every one of the 57 squarefree coefficient-ideal directions of ideal
> norm at most 128, the complete first necessary packet leaves a strict open
> hyperideal patch.  Thus none closes at its first packet.  No final chamber
> volume or reflectivity conclusion follows without later-packet exhaustion.

Two additional scope/provenance defects were repaired:

1. Multiplication of a coefficient by a nonunit field square identifies the
   diagonal **quadratic space**, but need not give an isometry of integral
   `Z[sqrt(2)]`-lattices.  The 57 rows are search directions modulo field
   squares, not all integral lattice isometry classes with coefficient norm
   at most 128.
2. The old JSON artifact did not contain the advertised
   `all_first_packet_obstructions` array.  The regenerated v3 owner artifact
   contains exactly 57 such entries and states the bounded-prefix scope.

## Clean-room boundary

The independent backend
`src/coxeter_polytopes/search/silver_leech_manifest_cleanroom.py` imports only
the exact `Silver` field arithmetic.  It does not import the owner manifest,
its prime enumeration, its control enumerator, its Leech Gram filter, or its
pairing oracle.  Its independently generated JSON is compared structurally
with the owner JSON only after both computations finish.

## Prime ideals and coefficient directions

For an odd prime `p`, the replay tests whether `x^2=2 mod p` has a root.
This independently gives:

- the ramified atom `p2=(sqrt(2))`, of ideal norm 2;
- inert atoms `(3)`, `(5)`, `(11)`, of ideal norms 9, 25, 121;
- two distinct atoms over each of
  `7,17,23,31,41,47,71,73,79,89,97,103,113,127`.

Thus there are `1+3+2*14=32` atoms within the cutoff.  Exhausting squarefree
products whose ideal norms are at most 128 gives 57 lineages.  The count by
ideal norm is

```text
1:1  2:1  7:2  9:1  14:2  17:2  18:1  23:2  25:1  31:2
34:2  41:2  46:2  47:2  49:1  50:1  62:2  63:2  71:2
73:2  79:2  82:2  89:2  94:2  97:2  98:1  103:2  113:2
119:4  121:1  126:2  127:2
```

The split ideals are retained separately.  In particular, the two norm-seven
generators are `3-sqrt(2)` and `3+sqrt(2)`, while their product gives the
norm-49 lineage with canonical coefficient `-7+7sqrt(2)`.  The inert-content
lineages are also present: `(3)` gives `-3+3sqrt(2)`, and its product with the
dyadic atom gives `3sqrt(2)` of absolute field norm 18.

Every atom generator is chosen totally positive.  Multiplication by the
sign-changing unit `epsilon=1+sqrt(2)` produces sign `(+,-)`, and multiplication
by `epsilon^(2k)` gives the unique half-open Pell-domain representative.  All
57 reconstructed coefficients have exact field norm equal to minus their
ideal norm and have the required two embedding signs.

## Local mirror ideals and spatial gates

For a primitive root `r=(x,y)`, put

```text
I(r) = (a*x, (y,Lambda)),       q(r) | 2 I(r).
```

At an odd prime ideal `P` occurring once in the squarefree coefficient ideal,
the local alternatives are exhaustive:

- if `y` is primitive modulo `P`, then `I(r)` is a unit and `v_P(q)=0`;
- if `y` lies in `P Lambda`, primitivity forces `x` to be a unit and
  `v_P(q)=1`.

The second case therefore carries the spatial gate `y in P Lambda`.  Choices
at distinct odd prime ideals are independent.  At the ramified dyadic prime,
evenness of the Leech lattice and `q|2I` give valuations `{0,2}` when the
coefficient is dyadically odd and `{1,2}` when it contains the squarefree
dyadic atom.  Their totally-positive representatives are respectively
`{1,2}` and `{2+sqrt(2),2}`.

This reproduces the familiar `a=-1+2sqrt(2)` list

```text
1, 2, 3+sqrt(2), 6+2sqrt(2),
```

with gate `3+sqrt(2)` on the last two classes.  It also catches the more
fragile norm-49 case: the two primes above 7 vary independently, so there are
eight classes.  Their odd gates range over

```text
1, 3-sqrt(2), 3+sqrt(2), 7,
```

each combined with both dyadic alternatives.  Inert rows use the rational
gate `p`; products such as the norm-126 rows correctly combine that gate with
one selected prime above 7 and the dyadic alternatives.

## Exhaustion of the first necessary packet

For every local norm `q`, the replay enumerates every `x=m+n sqrt(2)` inside
the exact two-embedding rectangle

```text
x_0^2 <= q*701/a,       x_1^2 <= q/(-a_1).
```

It then imposes, exactly,

```text
x_0 > 0,
2*a*x/q in Z[sqrt(2)],
y^2 = q+a*x^2 >= 0 at the conjugate embedding,
h = a*x^2/q <= 701 at the distinguished embedding.
```

For an odd spatial gate `d`, the target is normalized as `y^2/d^2`, not by
the rational ideal norm of `d`.  Grouping all classes at equal height and
discarding only impossible binary Leech Gram data gives the same first packet
for all 57 directions.  The largest first height is
`335+237sqrt(2)` for `a=7+9sqrt(2)`, so the exact cutoff 701 exhausts every
earlier control level in the manifest.

Writing `y=u+sqrt(2)v`, a target `A+Bsqrt(2)` forces

```text
u^2+2v^2=A,       u.v=B/2.
```

The clean-room filter enumerates every even `(u^2,v^2)` allowed by the Leech
minimum and Cauchy.  It does not test a hand-picked list of short
combinations.  For a positive-definite binary Gram matrix of determinant
`Delta`, a combination of norm below four must satisfy

```text
i^2 <= 3*v^2/Delta,       j^2 <= 3*u^2/Delta,
```

so a finite rectangle exhausts **all** integer pairs `(i,j)`.  Rank-one data
are decided by `gcd(u^2,v^2,u.v)`.  Every surviving triple in all 57 packets
passes this exhaustive minimum-four test.

These are complete **necessary** packets.  A surviving Gram triple is not
claimed to be Leech-realized, Vinberg-acute, or an actual wall.

## Uniform fixed-direction replay

The direction norms `4,6,8,10` are nonempty without relying on a shell-count
flag: the audit archives explicit standard scaled-MOG vectors.  A type-`4^2`
minimal vector `P`, a disjoint one `R`, and the all-positive vector `Q` on the
standard Golay octad `{1,2,6,7,8,10,12,23}` give

```text
P, P+Q, P+R, P-Q
```

of respective Leech squares `4,6,8,10`.

For each necessary Gram type and fixed direction square `D`, the replay
enumerates every integral pairing allowed by the exact `3x3` positive
semidefinite Gram determinant and the minimum-four constraints on
`d+/-w,d+/-v`.  Dropping other Leech constraints only enlarges this set, so
the resulting `M_i` is a rigorous upper bound for every actual wall in the
component.

If a component has spatial divisor `d_i`, its correctly normalized slab bound
is

```text
B_i = a*x_i/d_i.
```

Thus `c_max=min_i(B_i/M_i)`.  For every row, at least one explicit direction
has `D*c_max^2-a>0`.  The archived rational shrink
`c=((n-1)/n)c_max` makes both the minimum wall margin and the hyperideal excess
strictly positive.  The tightest limiting clearance is independently the
norm-119 direction `a=3+8sqrt(2)`:

```text
D = 10,
c_max = 53/47 + (5/94)sqrt(2),
D*c_max^2-a = (21588-15022sqrt(2))/2209 > 0.
```

This supplies an open hyperideal patch for each finite first-packet prefix.
It cannot be promoted to a statement about the chamber after later mirrors.

## Exact comparison and tests

After independent generation, a structural comparison found equality of:

- all 32 atom records;
- all 57 ideal lineages and canonical coefficients;
- all local norm classes and spatial divisors;
- all 57 first heights, components, normalized targets, and surviving triples;
- every field of all 57 strict first-prefix witnesses.

The combined owner and clean-room suite passes:

```text
PYTHONPATH=.:src pytest -q \
  tests/test_silver_leech_manifest.py \
  tests/test_silver_leech_manifest_cleanroom.py

14 passed in 123.04s
```

The clean-room suite after adding the explicit MOG direction witnesses passes
separately with `7 passed`.

## Pinned SHA-256

```text
0f6efe6f454fb1133010d6f5694270ae54a9c7f96cb5f0814e3e090e3aff1069  src/coxeter_polytopes/search/silver_leech_manifest.py
7770b981852ab3ee0c58077e6501f59899ba3d0b40421562c7548b8fbf2f3043  src/coxeter_polytopes/search/silver_leech_manifest_cleanroom.py
98c72b6fcd7dae8eeb62e6ffb4306b3567bf844236c8bfecae3ff3499572bdf1  scripts/audit_silver_leech_low_norm_manifest.py
ceea5f3526e1fe64da53b443994a192f3253ed61cb7fc04c3085273e9a01b1a6  scripts/audit_silver_leech_low_norm_manifest_cleanroom.py
6025514976ae12735d94b9daa7088c29752c24baff2c051d3993d4597962ab71  tests/test_silver_leech_manifest.py
a6c777c6435132f3d5a87b74c5f1a2b651d14b98b46a6c11458d63593ec61ec9  tests/test_silver_leech_manifest_cleanroom.py
560f415844f20091fd8a7ef2e4db9667d9130148786ed58af201f0a948e13ffa  data/enumerations/silver_leech_low_norm_manifest.json
765c51828312befda80e76eca8cb454671eb0f1c1fe1f532554e210cf87f6065  data/enumerations/silver_leech_low_norm_manifest_cleanroom.json
0c7f2b7223b9bfa7c5b793dcbe22a5cc60e7510cafaa8a356077f3e84824214c  results/silver_leech_low_norm_manifest.md
```
