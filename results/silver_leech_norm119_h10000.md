# Silver Leech H24, norm-119 necessary prefix through H=10000

For

\[
L=\langle-(3+8\sqrt2)\rangle\oplus\Lambda_{24}
\]

over `Z[sqrt(2)]`, the complete arithmetic control prefix through `H=10000`
still leaves a strict hyperideal direction.  This is a bounded-prefix
obstruction, not a global nonreflectivity or reflectivity theorem.

Reproduce the archive and tests with

```console
PYTHONPATH=src python3 scripts/audit_silver_leech_norm119_h10000.py
PYTHONPATH=src pytest -q tests/test_silver_leech_norm119.py
```

## Local normalization, independently rechecked

Put `a=3+8sqrt(2)`.  Exactly,

\[
N(a)=-119,
\qquad
a=(3+\sqrt2)(5+2\sqrt2)(-1+\sqrt2).
\]

The first two factors are totally positive prime generators of norms 7 and
17; the last factor is a sign-changing unit.  At each odd prime dividing
`(a)`, the reflection-integrality valuation is either zero or one.  In the
second case the spatial vector lies in that prime times the Leech lattice.
The two odd choices are independent.  The even unimodular block supplies
dyadic factors 1 and 2.  This gives exactly eight local mirror classes:

| canonical `q` | dyadic factor | spatial divisor `d` |
|:--|--:|:--|
| `1` | 1 | `1` |
| `2` | 2 | `1` |
| `3+sqrt(2)` | 1 | `3+sqrt(2)` |
| `13-5sqrt(2)` | 1 | `19+11sqrt(2)` |
| `5+2sqrt(2)` | 1 | `5+2sqrt(2)` |
| `6+2sqrt(2)` | 2 | `3+sqrt(2)` |
| `26-10sqrt(2)` | 2 | `19+11sqrt(2)` |
| `10+4sqrt(2)` | 2 | `5+2sqrt(2)` |

The apparently mismatched product row is an associate identity, not a change
of spatial divisor:

\[
13-5\sqrt2=(19+11\sqrt2)(1+\sqrt2)^{-2}.
\]

For every class, writing the spatial vector as `y=d w` gives the simultaneous
normalizations

\[
w^2=\frac{q+a x^2}{d^2},
\qquad
z\mathbin{\cdot}w\le \frac{a x}{d}.
\]

The code rechecks these identities for all eight classes.  An independent
clean-room replay of the preceding norm-128 manifest also agrees on this
normalization and on the norm-119 first packet.

## Complete arithmetic census

The Vinberg height convention is

\[
H=\frac{a x^2}{q}.
\]

Every algebraic-integer control satisfying the two-embedding bound,
`2ax/q in Z[sqrt(2)]`, and nonnegative conjugate spatial target was enumerated
through `H<=10000`.  The result has:

- 358 arithmetic controls;
- 347 controls with no surviving necessary Leech Gram type;
- 7 nonempty heights and 11 nonempty mirror-class components;
- 176 raw binary Gram types;
- 12 types excluded by an explicit norm-1, norm-2, or norm-3 integral
  combination;
- 164 surviving necessary Gram types.

All controls with odd spatial divisor `d != 1` are empty through this cutoff.
That is a computed result, not an assumption: the census includes 16 controls
in the `q=3+sqrt(2)` class, 26 in `q=5+2sqrt(2)`, 33 in
`q=6+2sqrt(2)`, 53 in `q=10+4sqrt(2)`, 70 in
`q=13-5sqrt(2)`, and 140 in `q=26-10sqrt(2)`.

The complete nonempty component table is:

| packet | exact height | `q` | `x` | normalized target | raw | surviving |
|--:|:--|:--|:--|:--|--:|--:|
| 0 | `243+172sqrt(2)` | `1` | `3+2sqrt(2)` | `244+172sqrt(2)` | 4 | 3 |
| 0 | `243+172sqrt(2)` | `2` | `4+3sqrt(2)` | `488+344sqrt(2)` | 9 | 8 |
| 1 | `486+344sqrt(2)` | `2` | `6+4sqrt(2)` | `974+688sqrt(2)` | 11 | 9 |
| 2 | `972+688sqrt(2)` | `2` | `8+6sqrt(2)` | `1946+1376sqrt(2)` | 3 | 1 |
| 3 | `1417+1002sqrt(2)` | `1` | `7+5sqrt(2)` | `1418+1002sqrt(2)` | 13 | 12 |
| 3 | `1417+1002sqrt(2)` | `2` | `10+7sqrt(2)` | `2836+2004sqrt(2)` | 26 | 26 |
| 4 | `2834+2004sqrt(2)` | `2` | `14+10sqrt(2)` | `5670+4008sqrt(2)` | 36 | 36 |
| 5 | `3563+2520sqrt(2)` | `1` | `11+8sqrt(2)` | `3564+2520sqrt(2)` | 10 | 9 |
| 5 | `3563+2520sqrt(2)` | `2` | `16+11sqrt(2)` | `7128+5040sqrt(2)` | 19 | 18 |
| 6 | `4737+3350sqrt(2)` | `1` | `13+9sqrt(2)` | `4738+3350sqrt(2)` | 15 | 13 |
| 6 | `4737+3350sqrt(2)` | `2` | `18+13sqrt(2)` | `9476+6700sqrt(2)` | 30 | 29 |

The JSON records all 358 controls, including every empty control and every
raw/surviving Gram triple.  The last enumerated control is at
`4975+3518sqrt(2)` (approximately 9950.20); completeness up to the rational
cutoff 10000 follows from the exact two-embedding enumeration rather than
from stopping at that last value.

## Exact fixed-direction obstruction

For a fixed Leech vector `d0` of norm `D`, and every surviving binary Gram
type for `w=u+sqrt(2)v`, a positive-semidefinite `3 by 3` Gram determinant
together with the Leech minimum four gives an exact upper bound `M` for
`d0.w`.  The optimized feasible scaling is

\[
c_{\max}=\min_i\frac{B_i}{M_i},
\qquad B_i=\frac{a x_i}{d_i}.
\]

The oracle enumerates the exact integer interval cut out by the Gram
determinant.  Tests compare it to the older independent rectangular oracle on
a representative from every nonempty component and every
`D in {4,6,8,10}`.

At `H=10000`, directions of norms 4, 8, and 10 no longer certify a hyperideal
point.  Their limiting excesses are strictly negative.  Norm 6 survives.  Its
limiting component is packet 4,

```text
H = 2834+2004sqrt(2)
q = 2
x = 14+10sqrt(2)
B = 202+142sqrt(2)
M = 129+93sqrt(2),
```

and hence

\[
c_{\max}
=\frac{B}{M}
=\frac{118}{219}+\frac{52}{73}\sqrt2,
\]

\[
6c_{\max}^2-a
=\frac{77231}{15987}-\frac{18088}{5329}\sqrt2
>0
\quad(\text{approximately }0.0306571).
\]

To make every wall inequality strict, the archive chooses the smallest
integer shrink denominator that remains hyperideal, `n=936`, and sets

\[
c=\frac{935}{936}c_{\max}
=\frac{55165}{102492}+\frac{935}{1314}\sqrt2.
\]

For a fixed norm-six Leech vector `d0`, the point `z=c d0` has exact minimum
wall margin

\[
\frac{57515}{51246}-\frac{19541}{34164}\sqrt2>0
\quad(\text{approximately }0.313435)
\]

over all 11 nonempty components, and exact Lorentzian excess

\[
6c^2-a
=\frac{8428441993}{1750768344}
-\frac{38203717}{11222874}\sqrt2
>0
\quad(\text{approximately }2.31\cdot10^{-5}).
\]

Thus `z` lies strictly beyond the light cone and strictly inside every
necessary slab in the complete `H<=10000` prefix.  Since the necessary packet
is a superset of the actual realized, primitive, and Vinberg-accepted walls,
the same point survives every actual wall in that bounded prefix.

## Exact scope

Surviving Gram triples are necessary data only.  They have not all been
realized as Leech pairs, tested for primitive roots, or filtered by earlier
Vinberg walls.  The direction certificate is nevertheless valid for this
bounded prefix because it satisfies the larger necessary wall set.

No conclusion is made beyond `H=10000`; in particular, this is not a proof of
global nonreflectivity, termination, reflectivity, compactness, or novelty.

## Frozen artifacts

| artifact | SHA-256 |
|:--|:--|
| `src/coxeter_polytopes/search/silver_leech_norm119.py` | `128358c8bae27b3ad81b94d0aa9382de71641fff323b6dda820f1cbd0a6ec649` |
| `src/coxeter_polytopes/search/silver_leech_manifest.py` | `0f6efe6f454fb1133010d6f5694270ae54a9c7f96cb5f0814e3e090e3aff1069` |
| `scripts/audit_silver_leech_norm119_h10000.py` | `e05cc29884d7830906d3e21946a95be42217b9761279803b85bf9d8d4fbcbd6f` |
| `tests/test_silver_leech_norm119.py` | `d68a1b14308b365106a412fc40e71261ccc88dfb7bef7cb58ef6026145f0fb35` |
| `data/enumerations/silver_leech_norm119_h10000.json` | `3e845cada57b93a1aeecdaaf15310ca5fdaf1a97a05f4e9948e96bd9d0c52b13` |
