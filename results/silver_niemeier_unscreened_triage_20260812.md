# Silver-field triage of the unsearched rooted Niemeier hosts

## Outcome

No new compact or finite-volume Coxeter polytope was found in this screen.
The new exact host is

\[
  \langle-(1+\sqrt2)\rangle\oplus N(E_6^4)
  \quad\text{in }\mathbf H^{24}.
\]

Its full ternary-tetracode glue produces a genuinely coupled second Vinberg
packet, so this is not merely another disconnected-root-system proxy.
Nevertheless an explicit hyperideal ray survives every wall through the
complete height-15 prefix.

The triage covers the 19 rooted rank-24 Niemeier systems left after excluding
the already screened Leech, `A1^24`/Golay, `D24+`, `E8^3`, and `A24` hosts.
It identifies `A17 E7`, `A15 D9`, and `D12^2` as the next three exact glue
backends worth implementing.

## Exact common silver packet

For any even unimodular Niemeier block `N`, take the unit coefficient

\[
  a=1+\sqrt2.
\]

For a primitive root `r=(x,y)`, reflection integrality gives
`q(r) | 2 I(r)`.  Even unimodularity and the ramification
`(2)=(sqrt(2))^2` leave exactly the canonical norm representatives `1,2`.
The universal factorwise affine wall is

\[
  r_R=(2,(1+\sqrt2)\theta_R),\qquad r_R^2=2,
  \qquad H(r_R)=2+2\sqrt2,
\]

where `theta_R` is the highest root of a component `R`.

The other already-derived small coefficient classes transfer uniformly to
every even unimodular Niemeier block.  They define the natural coefficient
queue after the unit screen:

| field | timelike coefficient `a` | complete canonical mirror norms | extra local gate |
|:---|:---|:---|:---|
| `Q(sqrt(2))` | `1+sqrt(2)` | `1, 2` | none |
| `Q(sqrt(2))` | `sqrt(2)` | `2, 2+sqrt(2)` | dyadic valuation gate |
| `Q(sqrt(2))` | `2sqrt(2)` | `2, 4+2sqrt(2)` | top class has divisible spatial part |
| `Q(sqrt(2))` | `2sqrt(2)-1` | `1, 2, 3+sqrt(2), 6+2sqrt(2)` | the norm-seven classes have spatial part in `(3+sqrt(2))N` |
| cubic seven | `theta` (unit-square equivalent to `theta-1`) | `1, 2` | none |
| cubic seven | `2theta` | `2` | the whole form is even |
| cubic seven | `theta(2-theta)` | `1, 2, 2-theta, 2(2-theta)` | prime-seven gate |
| cubic seven | `2theta(2-theta)` | `2, 2(2-theta)` | dyadic and prime-seven gates |

The table follows from the same primitive pairing-ideal argument used in the
existing `D24+` and `E8^3` backends; no property specific to those root
systems enters beyond even unimodularity.  The present numerical triage uses
only the first silver row, which is the cheapest complete comparison.

If no earlier glue wall is inserted, the dominant section cut by these walls
is the product of affine Weyl alcoves.  At a finite node `i` its vertex energy
is exactly

\[
  \frac{4(C_R^{-1})_{ii}}{m_i^2},
\]

where `m_i` is the highest-root mark.  Convolving these rational energy lists
gives the table below.  Since the timelike coefficient `1+sqrt(2)` is
irrational, none of these factorwise vertices is ideal.

This is a **triage invariant**, not a claim that the displayed factorwise
packet is the first full Niemeier packet: glue vectors can give earlier or
simultaneous coupled walls.  That glue gate is resolved completely only for
`E6^4` in this report.

## Ranked host table

`Coupling support` is the support of a shortest explicitly certified word
which couples distinct root components.  A dash means that this report has
not decoded that host's glue code, not that coupling is absent.

| priority | host | components | roots | glue index | coupling support | affine rays | finite | hyperideal | glue status |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| 1 | `A17 E7` | 2 | 432 | 6 | 2 | 144 | 2 | 142 | forced cross-component discriminant glue |
| 2 | `A15 D9` | 2 | 384 | 8 | 2 | 160 | 2 | 158 | exact cyclic glue; screened through H=50 |
| 3 | `D12^2` | 2 | 528 | 4 | 2 | 169 | 3 | 166 | forced cross-component discriminant glue |
| 4 | `A12^2` | 2 | 312 | 13 | 2 | 169 | 1 | 168 | forced cross-component discriminant glue |
| 5 | `A11 D7 E6` | 3 | 288 | 12 | -- | 672 | 3 | 669 | three-factor mixed glue |
| 6 | `A9^2 D6` | 3 | 240 | 20 | -- | 700 | 2 | 698 | three-factor mixed glue |
| 7 | `D10 E7^2` | 3 | 432 | 4 | -- | 704 | 4 | 700 | three-factor binary glue |
| 8 | `D8^3` | 3 | 336 | 8 | -- | 729 | 4 | 725 | three-factor glue |
| 9 | `A8^3` | 3 | 216 | 27 | -- | 729 | 1 | 728 | three-factor glue |
| 10 | `E6^4` | 4 | 288 | 9 | 3 | 2,401 | 5 | 2,396 | explicit ternary tetracode; screened below |
| 11 | `A7^2 D5^2` | 4 | 192 | 32 | -- | 2,304 | 3 | 2,301 | four-factor mixed glue |
| 12 | `D6^4` | 4 | 240 | 16 | -- | 2,401 | 5 | 2,396 | four-factor glue |
| 13 | `A6^4` | 4 | 168 | 49 | -- | 2,401 | 1 | 2,400 | four-factor glue |
| 14 | `A5^4 D4` | 5 | 144 | 72 | -- | 6,480 | 2 | 6,478 | five-factor mixed glue |
| 15 | `D4^6` | 6 | 144 | 64 | -- | 15,625 | 7 | 15,618 | many-factor code |
| 16 | `A4^6` | 6 | 120 | 125 | -- | 15,625 | 1 | 15,624 | many-factor code |
| 17 | `A3^8` | 8 | 96 | 256 | -- | 65,536 | 1 | 65,535 | many-factor code |
| 18 | `A2^12` | 12 | 72 | 729 | -- | 531,441 | 1 | 531,440 | many-factor code |
| 19 | `D16 E8` | 2 | 720 | 2 | -- | 153 | 3 | 150 | decomposable `D16+` plus `E8`; glue does not couple the factors |

The top four positions use exact discriminant constraints, rather than a
guessed code distance:

- `A12^2` has prime discriminant 13 in each factor, so its index-13 glue is
  the graph of an anti-isometry and every nonzero coupling word has support
  two.
- `D12` has no proper even integral same-rank overlattice, so the index-four
  glue of `D12^2` must cross the two factors.
- In `A17 E7` and `A15 D9`, the odd ranks rule out an even unimodular factor;
  their index-six and index-eight glue therefore contain cross-factor words.
- `D16 E8` is deliberately demoted despite its small affine census: the
  index-two spinor extension is internal to `D16`, giving the decomposable
  spatial lattice `D16+` plus `E8`.

The priority within the forced two-factor tier is the exact number of
remaining hyperideal affine-product rays.  This is a search-cost heuristic,
not a theorem that a higher-ranked row cannot close earlier.

## Full `E6^4` tetracode screen

Write an E6 weight in Dynkin coordinates `z=C_E6 y`.  The discriminant
functional used here is

\[
  \rho(z)=(1,2,0,1,2,0)\cdot z\pmod 3.
\]

The four residues lie in the self-dual ternary tetracode generated by

\[
  (1,0,1,1),\qquad(0,1,1,-1).
\]

It has nine words: zero and eight words of support three.  Membership is
imposed independently on the rational and `sqrt(2)` coefficients, exactly
as required for `N(E6^4) tensor Z[sqrt(2)]`.  The height-zero census is
exactly four dominant norm-two orbits, hence precisely the `E6^4` root
system and 288 spatial roots.

There are eight complete reflective control levels through rational height
15.  The five simultaneous height packets are:

| height | raw candidates | accepted |
|:---|---:|---:|
| `(1+sqrt(2))/2` | 0 | 0 |
| `1+sqrt(2)` | 0 | 0 |
| `2+2sqrt(2)` | 4 | 4 |
| `(7+5sqrt(2))/2` | 0 | 0 |
| `7+5sqrt(2)` | 44 | 44 |

The four first walls are the independent affine E6 caps.  Their exact
product geometry has 2,401 rays: 5 finite, no cusps, and 2,396 hyperideal.

The second nonempty packet is more informative.  Its 44 roots are all
mutually Vinberg-acute and all survive every earlier wall.  They split as

- 20 roots in the root-lattice residue; and
- 24 roots with nonzero tetracode residue, genuinely coupling three E6
  factors.

Thus the glue acts early and substantially.  It still does not close the
chamber.  Let `z=(1,w)` have Dynkin weights

\[
  (0,0,0,0,0,3/5)
\]

in each of the four E6 blocks.  Since the relevant inverse-Cartan diagonal
is 2,

\[
  z^2=-\left(1+\sqrt2\right)+4\cdot2\left(\frac35\right)^2
     =\frac{47}{25}-\sqrt2>0.
\]

The exact audit finds negative pairing with all 48 positive-height walls;
the height-zero inequalities are nonpositive.  Hence this ray is hyperideal
in the complete height-15 chamber prefix.  The prefix is not finite-volume,
so it is not a compact example.  Later walls may still cut this ray; no
global nonreflectivity statement is made.

An independent cleanroom reconstruction, using its own integer box, Schur
enumerator, tetracode join, and pairing routine, returned PASS with no
discrepancy.  In particular it recovered the witness-pairing histogram

| exact pairing | walls |
|:---|---:|
| `-1-(6/5)sqrt(2)` | 24 |
| `-1-(3/5)sqrt(2)` | 12 |
| `-11/5-(9/5)sqrt(2)` | 8 |
| `-(4/5)-(4/5)sqrt(2)` | 4 |

so strictness does not depend on a numerical sign comparison.

## Next exact searches

The next bounded implementations should be:

1. `D12^2`, whose index-four glue is especially small;
2. `A12^2`, using its prime index-thirteen graph glue; and
3. the remaining many-factor glue hosts in the ranked order above.

The silver-unit `A17 E7` and `A15 D9` lanes have now been implemented.  The
`A15 D9` chronology is glue-complete through height 50 for the unit seed: its
three nonempty packets contain `2,15,191` walls, but a strict D9-axis
hyperideal ray survives all 208 positive-height walls.  The distinct
prime-seven coefficient has also been completed through height 63; its dense
second packet leaves 1,655 walls, but a strict D9-axis ray survives all 1,666
positive walls.  See `results/silver_a15_d9_h50.md` and
`results/silver_p7_a15_d9_h63.md`.

For each remaining host, the same dominant-weight engine can enumerate the
silver unit seed first.  The glue decoder should then be reused for the cubic-seven unit
seed `a=theta-1`, where even unimodularity again gives canonical mirror norms
`1,2`.  The cubic-seven control chronology must be recomputed rather than
inferred from the silver table.

## Reproduction

- E6 backend: `src/coxeter_polytopes/search/silver_e6_tetracode.py`
- Triage backend: `src/coxeter_polytopes/search/niemeier_host_triage.py`
- Audit: `scripts/audit_silver_e6_tetracode_h15.py`
- Tests: `tests/test_silver_e6_tetracode.py`,
  `tests/test_niemeier_host_triage.py`
- Archive: `data/enumerations/silver_e6_tetracode_h15.json`

```console
PYTHONPATH=src python3 scripts/audit_silver_e6_tetracode_h15.py
PYTHONPATH=src pytest -q tests/test_niemeier_host_triage.py tests/test_silver_e6_tetracode.py
```

The focused suite passes `7/7` tests.  The archive SHA-256 is
`875d49ffd81d5b8a56ab4de34dfb083cbbbbb4a69e4a5854a501452a077509bc`.
