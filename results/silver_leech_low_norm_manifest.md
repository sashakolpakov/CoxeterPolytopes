# Silver Leech H24 coefficient-direction manifest through norm 128

This is a finite exact screen for

\[
L_a=\langle-a\rangle\oplus\Lambda_{24},\qquad
a=A+B\sqrt2,quad N(a)<0,quad |N(a)|\le128,
\]

For this prioritization screen, coefficients are identified under
multiplication by a field square; a final totally-positive unit square
`(1+sqrt(2))^(2k)` chooses the exact Pell-domain representative.  This is an
equivalence of diagonal quadratic spaces, not an asserted integral isometry
between lattices obtained by multiplying a coefficient by a nonunit square.
Since `Z[sqrt(2)]` is a PID and every totally-positive unit is a square, the
resulting squarefree coefficient-ideal directions are the squarefree products
of prime ideals.  There are **57** such directions of ideal norm at most 128.

The census includes the ramified ideal `(sqrt(2))`, inert ideals `(p)` of
norm `p^2`, and the two conjugate primes separately whenever `p` splits.
Thus rational-content examples missing from the earlier primitive slice are
present: the inert `p=3` class has canonical coefficient `-3+3pi`, while its
product with the dyadic prime has coefficient `3pi` and norm `-18`.

Reproduce with

```console
PYTHONPATH=src python3 scripts/audit_silver_leech_low_norm_manifest.py
PYTHONPATH=src pytest -q tests/test_silver_leech_manifest.py
```

The JSON artifact records every coefficient, every locally possible mirror
norm, its odd-prime spatial divisibility gate, and the complete first
necessary binary Leech Gram packet.  “Necessary” is deliberate: a Gram type
can still fail Leech realization or Vinberg acuteness, so this census does
not manufacture roots or turn an empty/loose screen into a theorem.

## Local mirror classes

Write `pi=sqrt(2)`.  If
`|N(a)|=2^e product p_i^e_i`, every odd prime-power ideal dividing `(a)` has
reflection valuation either zero or its full coefficient valuation.  In the
latter case reflection integrality forces `y in d Lambda`, and the exact
shell census uses `y=d z`, hence target `y^2/d^2`.  The dyadic alternatives
are represented by `1,2` for odd coefficient norm and by `pi,2` for even
coefficient norm.  Taking every product yields the complete finite local
class manifest.

As a check, `a=2pi-1` gives exactly

```text
1, 2, 3+pi, 6+2pi,
```

with the last two classes gated by `y in (3+pi)Lambda`, the same ideal as
`(2pi-1)`.

## First-packet ranking

Ten of the 57 first necessary packets consist entirely of rank-one scalar
Leech-shell components; seven consist of a single scalar shell:

| coefficient `a` | `|N(a)|` | first height | complete scalar shell |
|:--|--:|:--|--:|
| `-1+pi` | 1 | `14+10pi` | 10 |
| `pi` | 2 | `8+6pi` | 6 |
| `-1+2pi` | 7 | `5+4pi` | 4 |
| `-6+7pi` | 62 | `19+14pi` | 4 |
| `-4+7pi` | 82 | `29+21pi` | 6 |
| `3+6pi` | 63 | `33+24pi` | 4 |
| `7pi` | 98 | `31+22pi` | 4 |

The first packets of the first three already have strict hyperideal witnesses.
None of the four new single-shell rows closes at its first packet either.  For
`a=-6+7pi`, `3+6pi`, and `7pi`, the normalized
packet is the full minimal shell; a fixed norm-six direction pairs with it by
at most three.  For `a=-4+7pi`, the packet is the full norm-six shell and a
fixed minimal direction pairs by at most three.  Choosing `9/28` of the
normalized slab bound gives strict wall margin and strict positive
Lorentzian norm in every case.  In particular, the two closest earlier rows
have exact excesses

\[
\frac{17457}{196}-\frac{827}{14}\pi>0
\quad(N=62),\qquad
\frac{5009}{98}-\frac{211}{7}\pi>0
\quad(N=82).
\]

Thus the first packet of every scalar-shell row in this norm-128 manifest is
still infinite-volume.  This is a bounded-prefix statement: later mirrors
can cut the displayed open hyperideal patch.

The closest rank-two priority in the fixed-direction obstruction ranking is
the pair of norm-89 squareclasses.  In canonical coordinates their exact
first-packet certificates are:

- `a=-3+7pi`: first height `117+83pi`, `q=2`, `x=4+3pi`, target
  `236+166pi`, five surviving Gram types.  For any fixed minimal Leech vector
  `d`, the exact `3x3` Gram/minimum oracle gives
  `max d.y <= 21+14pi`.  Thus `z=(4/3)d` has wall margin
  `2+pi/3>0` and hyperideal excess `91/9-7pi>0`.
- `a=3+7pi`: first height `219+155pi`, eight surviving types and
  `max d.y <=29+21pi`.  The witness `z=(9/5)d` has margin
  `(9-4pi)/5>0` and excess `249/25-7pi>0`.

Consequently neither norm-89 first necessary packet closes the chamber, even
if every necessary Gram type is realized and accepted.  Later packets may
still close it.  The first class is the cheaper next-packet lane (five types
rather than eight).  All
coefficient and witness statements here use the canonical fundamental-domain
representatives; unit-square associates must transform the time and spatial
coordinates together.

## Uniform first-prefix result for all 57 directions

The same calculation was run over the complete manifest, using fixed Leech
directions of norms `4,6,8,10`.  For each same-height component `i`, let
`M_i` be the rigorous `3x3`-Gram/minimum-four upper bound on `d.y`, and let
`B_i` be its correctly normalized slab bound (including cancellation of any
odd-prime spatial divisor).  Put

\[
c_{\max}=\min_i B_i/M_i.
\]

For **every one of the 57 rows**, at least one direction norm `D` satisfies

\[
D c_{\max}^2-a>0.
\]

The artifact does not stop at this limiting inequality.  It records an
explicit integer `n` and the strict witness
`z=((n-1)/n)c_max d`; its minimum wall margin and Lorentzian excess are both
stored as positive exact elements of `Q(sqrt(2))`.  Thus every complete first
necessary packet in the norm-128 ideal-squareclass manifest leaves an open
hyperideal patch.  Equivalently, none of the 57 prefixes closes at its first
packet.  This does not eliminate any coefficient from a full Vinberg search:
higher packets may cut the patch.

The tightest obstruction is the norm-119 class `a=3+8pi`.  Its best direction
has `D=10`,

\[
c_{\max}=\frac{53}{47}+\frac{5}{94}\pi,
\qquad
10c_{\max}^2-a=
\frac{21588-15022\pi}{2209}>0.
\]

That makes `a=3+8pi` the ranked next-packet continuation if this Leech family
is pursued beyond the first packet; the next alternatives by exact clearance
are the already-studied unit class and the norm-113 class `7+9pi`.

## Exact scope

The manifest is exhaustive for squarefree coefficient-ideal directions of
ideal norm at most 128 and for the stated local norm/divisibility model.  For
all 57, the finite prefix consisting of the complete first necessary packet
has an open hyperideal patch.  No conclusion about the final chamber follows
without exhausting later mirrors.

It is not a census of integral-lattice isometry classes obtained by inserting
nonunit square factors, and it is not exhaustive over larger squarefree ideal
norm, different spatial Niemeier lattices, or the cubic field.  It supplies a
ranked, reproducible seed queue; it does not prove or disprove reflectivity,
termination, finite volume, or novelty for any unexpanded continuation.
