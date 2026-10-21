# Prime-seven silver `N(A15 D9)` through height 63

## Outcome

The exact Vinberg prefix for

\[
 \langle-(2\sqrt2-1)\rangle\oplus N(A_{15}D_9)
 \quad\text{in }\mathbf H^{24}
\]

is not compact and does not have finite volume through rational cutoff 63.
The complete bounded chronology has 1,666 positive-height walls.  Its dense
second packet contains 1,655 accepted walls, including 1,412 with nontrivial
Niemeier glue and 1,252 which genuinely couple `A15` to `D9`.  Nevertheless,
the explicit `D9`-axis ray below is hyperideal and strictly satisfies every
wall inequality.

This is an exact bounded result, not a global nonreflectivity theorem; a wall
above the cutoff could cut the witness.

## Field and local mirror norms

Put

\[
 \pi=\sqrt2,\qquad \epsilon=1+\pi,\qquad
 a=2\pi-1,\qquad \delta=3+\pi.
\]

Then `epsilon*a=delta`, `N(a)=-7`, and `N(delta)=7`.  The spatial lattice is
the even unimodular Niemeier lattice with cyclic glue

\[
 \langle(2,1)\rangle
 \subset (\mathbf Z/16)\oplus(\mathbf Z/4).
\]

The exact index-eight glue and coefficientwise membership functional are the
ones proved in the companion unit-seed report
`results/silver_a15_d9_h50.md`.

For a primitive mirror `r=(x,y)`, reflection integrality gives

\[
 (q(r))\mid2I(r),\qquad I(r)=(ax,(y,N)).
\]

Even unimodularity excludes support away from `(pi)` and `(delta)`.  At
`(pi)`, the possible valuation is zero or two.  At `(delta)`, reflection
integrality either gives valuation zero or forces
`y in delta*N(A15 D9)` and valuation one.  Class number one and the fact
that all totally positive units are squares then give the complete canonical
norm list

\[
 1,\quad2,\quad\delta,\quad2\delta.
\]

For a `delta`-supported norm, the enumerator first writes `y=delta*y_0` and
uses the reduced shell.  The reduced target is divided by `delta^2`, whereas
the reduced reflection norm is divided by `delta`.  Reconstruction and the
glue membership test are performed in the full spatial lattice.  Norm-two
classes are additionally divided by `pi` in the full glue basis to remove
imprimitive mirrors.

## Complete control chronology

There are 21 controls at 14 simultaneous heights through cutoff 63.  The
first nonempty packet occurs at

\[
 H=5+4\sqrt2.
\]

Its norm-one control contributes two walls.  Its norm-two control contributes
nine primitive walls after two `pi`-divisible copies of the norm-one mirrors
are removed.  All eleven walls are mutually strictly acute.

Every intervening packet through `H=28+20sqrt(2)` is empty.  The next
nonempty packet is

\[
 H=31+22\sqrt2.
\]

The norm-one and norm-two controls give 15 and 1,655 primitive candidates,
respectively.  Exact comparison with the eleven earlier walls yields:

| status against earlier walls | count |
|:---|---:|
| candidates | 1,670 |
| rejected | 15 |
| accepted | 1,655 |

Across all `1,670*11=18,370` candidate/old-wall comparisons, the exact sign
census is 18,337 negative, 18 zero, and 15 positive.  The 15 candidates with
a positive pairing are precisely the rejected walls.

## Exact dense-packet acuteness

All inverse-Cartan denominators divide 16.  Multiplying the Lorentz pairing
matrix by 16 therefore gives integer pairs `A+B*sqrt(2)`.  The batched checker
constructs those coefficient matrices with signed 64-bit integer products;
the bounded entries are far below overflow.  Every sign is then decided from
`A`, `B`, and `A^2-2B^2`, without floating point.

Among the

\[
 \binom{1655}{2}=1,368,685
\]

same-height survivor pairings, exactly 1,367,834 are negative, 851 are zero,
and none are positive.  Hence the complete dense packet is simultaneous and
Vinberg-acute.  An independent implementation reproduced the candidate-old
and same-height integer-pair matrices and all sign counts.

The total bounded chamber prefix has

\[
 24+11+1655=1690
\]

walls, including the 24 height-zero `A15+D9` simple walls.

## Strict hyperideal witness

Let `z=(1,w)` have zero `A15` Dynkin weights and only the seventh `D9`
Dynkin weight nonzero:

\[
 w_{D,7}=\frac{13}{25}.
\]

In the numbering used here,
`(C_D9^-1)_(7,7)=7`.  Therefore

\[
 z^2=7\left(\frac{13}{25}\right)^2-(2\sqrt2-1)
    =\frac{1808}{625}-2\sqrt2>0.
\]

Dominance gives every height-zero inequality.  All 1,666 positive-height
pairings are strictly negative.  The closest is

\[
 -\frac{18}{25}+\frac15\sqrt2<0.
\]

Thus the `H=63` prefix has an explicit strict hyperideal ray and infinite
volume.  This prime-seven coefficient produces a much denser packet than the
unit coefficient on the same host, but it does not remove the geometric
escape.

## Reproduction

- Backend: `src/coxeter_polytopes/search/silver_p7_a15_d9.py`
- Audit: `scripts/audit_silver_p7_a15_d9_h63.py`
- Tests: `tests/test_silver_p7_a15_d9.py`
- Certificate: `data/enumerations/silver_p7_a15_d9_h63.json`

```console
PYTHONPATH=src python3 scripts/audit_silver_p7_a15_d9_h63.py
PYTHONPATH=src pytest -q tests/test_silver_p7_a15_d9.py
```

The focused suite passes `4/4`.  This is a new exact negative coefficient
screen over `Q(sqrt(2))`, not a new compact Coxeter polytope.

The frozen certificate SHA-256 is
`f4642be55151e1c656d1c838ed5fd30c5626a194d1b4aa3c243811dece80e9b7`.
