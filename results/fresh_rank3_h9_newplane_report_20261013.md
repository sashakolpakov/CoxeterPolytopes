# Two-weight H9 slice: bounded binary failure and verified free-ratio cuts

A bound-independent projection certificate now excludes finite-volume
completion by any additional original-compatible walls in the surviving
four-dimensional slice or its symmetry image. This closes this particular
pair of slices; it is not an ambient nonreflectivity theorem.

The new binary line selected by a 28-subspace screen has **no integral
norm-one or norm-two roots with `0<t<=1,000,000`**. An independent arithmetic
implementation repeats the entire bounded search and agrees. This is only
one line of spatial ratios inside the surviving four-dimensional slice.

Freeing those ratios is productive: the parent search subsequently found
335 compatible reflective cutters in a smaller coefficient box. The
independent audit verifies every cutter, including **127,970 exact
Coxeter-angle checks** against the 382 original walls. This does not say
that all 335 cuts can be used together; constructing mutually compatible
collections and testing their volume are separate tasks.

## Exact subspace selection

The starting cone is case 0 of the
[rank-three construction](fresh_rank3_h9_cones_20261013.json), with 382 walls
and escape v of norm `(128-79phi)/12>0`. The
[parent screen](fresh_rank4_h9_slice_screen_20261013.json) considers all 28
spans of time, the last scalar, and two standard E8 fundamental weights.

For 27 spans it gives exact nonnegative linear combinations of the
restricted original wall-pairing functionals equal to the restricted
functional of v. If r belongs to such a span and has nonpositive pairing
with every existing wall, the identity forces `(r,v)<=0`. Thus no compatible
real normal in any of those 27 spans can cut this particular escape. This
argument is independent of arithmetic and coefficient bounds.

The remaining pair is `(4,7)`, with zero-based fundamental-weight indices.
Their positive Gram is `[[12,9],[9,8]]`. The screen gives a geometric cutter
with time `4+3phi`, spatial part `2*omega4+2phi*omega7`, and scalar `1+phi`.
Its norm is `16-9phi>0`, its escape pairing is `9-11phi/2>0`, and all 382
wall pairings are nonpositive. It is primitive because its scalar coordinate
is a unit, but its reflection coefficients are not all integral. Hence it
is not an arithmetic reflective wall. The independent verifier replays all
28 entries and this primitivity check.

## Binary line search

Let

```
d = (12+10phi,24+20phi,36+30phi,30+24phi,
     24+18phi,16+12phi,8+6phi,18+16phi,1+phi).
```

These are the eight E8 coordinates followed by the scalar coordinate.
They satisfy `d^2=beta=82+107phi`. The line search used roots `(x,t*d)` and

`beta*t^2-2phi*x^2=q`, with `q in {1,2}`.

The scalar coordinate of d is a unit, so this binary module is saturated.
Its compatible ratio interval is

`4+3phi <= x/t < -13/4+15phi/2`.

The right inequality is exactly the condition for cutting v. The search
used `0<t<=1,000,000` and the necessary conjugate bound
`beta'*t'^2<=q`. Modulo two, q=1 forces `t=a+bphi` with a even and b odd;
q=2 forces both coefficients even. All remaining exact square equations
were tested, with these counts:

| q | Bounded t values, including zero | After parity | Integral x solutions |
| --- | ---: | ---: | ---: |
| 1 | 224,518 | 56,130 | 0 |
| 2 | 317,517 | 79,380 | 0 |

The producer took 23.317 seconds. The independent replay constructs its own
alphabet directly from `|t'|<1` and uses an integer norm-and-trace square
test instead of the project's algebraic square-root routine. It repeats the
complete stated bound. Since there are no integral norm solutions, there
are no cuts to submit to the angle or volume gates in this binary run.

An optional installed PARI/GP relative-norm probe returned rational
solutions for both norms. They are independently verified:

| q | t | x |
| --- | --- | --- |
| 1 | `50474/3+27223phi` | `(448166+725152phi)/3` |
| 2 | `(34+50phi)/9` | `(284+455phi)/9` |

These are nonintegral and give no reflective roots. They show that the
rational norm equation is solvable; they do not resolve the integral Pell
problem. No larger cutoff was run.

## Audit of the free-ratio search

The [parent search](../scripts/fresh_rank4_h9_free_search_20261013.py) instead
uses `(x,t*omega4+u*omega7,z)` with bounds `t<=300`, `u<=600`, `z<=600`,
`t,u>0`, and `z>=0`. Its spatial energy is

`E=12t^2+18tu+8u^2+z^2`.

The inverse two-weight Gram gives exact conjugate coordinate bounds
`t'^2<=8q/15`, `u'^2<=4q/5`, while `z'^2<=q`. The scalar parity follows
from the even E8 form. The producer's `la+lb*phi` equals twice the spatial
pairing with v:

`2L=12t+(5+3phi)u+(-1+phi)z`.

Its `ga+gb*phi` is exactly `(2L)^2-8phi*(E-q)`. Since L is positive for
the stated parameter signs, this is the correct prefilter for an escape
cut. The radicand is `(E-q)/(2phi)`, with coefficients
`((eb-ea)/2,ea/2)`. Its field norm must be a nonnegative integer square.
These identities and the integer-array sign filter agree with independent
field arithmetic on deterministic samples.

For these fixed bounds all parameter coefficients have absolute value below
602. Direct bounds keep the initial quadratic and pairing arrays inside
int64; the producer additionally checks magnitudes before its sign squares
and radicand norm squares. The producer now explicitly rejects command-line bounds above
300/600/600; the audit concerns this supported range.

The [independent accepted-root audit](fresh_rank3_h9_free_audit_20261013.json)
reconstructs all 335 normals from their parameters and confirms:

- 86 roots have norm 1 and 249 have norm 2, with no duplicates.
- Every root is integral and has integral reflection coefficients.
- Every root satisfies the parameter bounds, strict interior, positive
  escape pairing, and its recorded height.
- Every one of its 382 pairings is nonpositive and passes the exact golden
  Coxeter-angle gate. Values for angles such as `2pi/5` are excluded.

The independent audit does not repeat the entire free-ratio parameter-box
enumeration and does not assert mutual compatibility among the 335 cuts.
A subsequent bound-independent certificate closes finite-volume completion
using this slice and its symmetry image, as detailed below.

## Bound-independent obstruction after the 335 cuts

The finite dictionary contains the 335 direct cuts and their 335 block-swap
images. An independent integer-pair audit recomputes all 224,115 dictionary
pairs, verifies all 340 conflict edges, and checks the twelve saved selections.
Each selection is a maximal compatible subset of this finite dictionary,
though the twelve do not exhaust all maximal subsets. The resulting cones
have 854 through 934 walls and all have exact spacelike escapes.

More strongly, use case 1 (`decreasing_height`) of the
[completion archive](fresh_rank4_h9_completion_20261013.json). Its exact
escape has norm

`901480/358801-2767458phi/1794005 > 0`.

Let `S=span(time,omega4,omega7,scalar)`. For an ambient vector w, its
pairing functional restricted to S has coordinates

`ell(w)=(-2phi*w[0], w[5], w[8], w[9])`.

The [exact certificate](fresh_rank3_h9_rank4_obstruction_20261013.json)
gives identities

```
ell(v)  = sum(lambda_i * ell(w_i)), i in {4,28,91,277},
ell(Uv) = sum(mu_i * ell(w_i)),     i in {29,186,226,252}.
```

All eight coefficients are strictly positive and are stored literally in
the certificate. Here w_i refers to the original 382-wall cone and U is the
integral involution permuting that wall set. Exact Fraction arithmetic checks
every coefficient sign and every coordinate of both identities, as well as
the spacelike norm and all 1,052 finite-dictionary pairings.

For any r in S with `(r,w_i)<=0` for all original walls, the first identity
implies `(r,v)<=0`. For any r in U(S), apply the second identity to Ur and
use the permutation property. Thus this same spacelike v survives **every**
original-compatible real wall from S or U(S). No reflection-integrality,
norm, height, or coefficient bound is needed for this conclusion.

Any full-dimensional cone obtained by such additions therefore retains a
spacelike feasible direction and cannot have finite volume. This makes a
larger box or a complete search of the current box unnecessary for the
finite-volume goal in these two slices. Other ambient directions remain
unexcluded. The 27 earlier slice identities concern a different fixed
escape; they must not be combined into a claim that every two-weight span
is simultaneously excluded for this new escape.

The separate construction audit independently rebuilt the four basis vectors
from the inverse E8 Gram, checked the restriction coordinates against all 382
walls, and replayed both identities. Its [report](fresh_rank4_h9_selection_audit_20261013.md)
provides a second finite-dictionary escape certificate as well.

- [Whole-slice certificate and replay](../scripts/fresh_rank3_h9_rank4_obstruction.py).
- [Independent dictionary graph and twelve-cone replay](../scripts/fresh_rank3_h9_dictionary_audit.py), [audit](fresh_rank3_h9_dictionary_audit_20261013.json).

## Files and validation

- [Binary search](../scripts/fresh_rank3_h9_newplane_search.py), [result](fresh_rank3_h9_newplane_20261013.json).
- [Independent binary and screen replay](../scripts/fresh_rank3_h9_newplane_verify.py), [audit](fresh_rank3_h9_newplane_audit_20261013.json).
- [Rational probe](../scripts/fresh_rank3_h9_newplane_rational_norms.gp), [exact rational witnesses](fresh_rank3_h9_newplane_rational_norms_20261013.json).
- [Independent 335-cut audit](../scripts/fresh_rank3_h9_free_audit.py), [result](fresh_rank3_h9_free_audit_20261013.json).
- [Regression tests](../tests/test_fresh_rank3_h9_newplane.py): square tests, array signs, polynomial filters, and forged-cut/projection rejection; all six pass.

```sh
PYTHONPATH=src:. python3 scripts/fresh_rank3_h9_newplane_search.py --limit 1000000
PYTHONPATH=. python3 scripts/fresh_rank3_h9_newplane_verify.py \
  --output results/fresh_rank3_h9_newplane_audit_20261013.json
PYTHONPATH=. python3 scripts/fresh_rank3_h9_free_audit.py \
  --output results/fresh_rank3_h9_free_audit_20261013.json
PYTHONPATH=. python3 scripts/fresh_rank3_h9_rank4_obstruction.py --verify
PYTHONPATH=. python3 scripts/fresh_rank3_h9_dictionary_audit.py \
  --output results/fresh_rank3_h9_dictionary_audit_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_rank3_h9_newplane.py
```
