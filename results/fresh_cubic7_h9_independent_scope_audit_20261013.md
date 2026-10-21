# Independent audit of the fresh cubic H9 search

The three height-100 prefixes have respectively 24, 17, and 14 walls. Every
saved wall and every saved spacelike escape passes the independent exact
checker. These are nonclosed bounded prefixes, not nonreflectivity results.
No compact or new-dimensional example is established.

## Primitive norm alphabet

The transfer from `E8^3` to the odd positive block `E8+<1>` is valid: both
are unimodular over `O_K=Z[theta]`, where
`theta^3+theta^2-2theta-1=0`. For primitive `r=(x,y)` and
`q=(r,r)`, integral reflection is equivalent to
`q | 2 I(r)`, with `I(r)=(a*x,G*y)`. Outside primes dividing `2a`, local
unimodularity and primitivity force `I(r)` to be the unit ideal.

For `a=theta`, this leaves norm ideals `(1),(2)`. For `a=2theta`, if the
pairing ideal is divisible by `(2)`, spatial coordinates are even and `x`
is a local unit. Then `q=-2theta*x^2 mod 4`, giving valuation exactly one,
so no norm-four class is missing. For `a=theta*(2-theta)`, the same argument
at `delta=2-theta` gives valuation one at that prime, leaving
`(1),(2),(delta),(2delta)`. The parity change in the positive block does not
affect any of these arguments.

The cubic discriminant is 49. Minkowski's bound is `14/9<2`, proving class
number one. The units `-1`, `theta`, and `theta+1` have independent sign
vectors at the three real places; hence the unit-signature map has full
rank and every totally positive unit is a square. The displayed norm
generators are therefore complete up to rescaling roots by units.

The norm classes cannot share a positive control height: equality would
make the quotient of two distinct classes a field square, whereas its
rational field norm is a nonsquare quotient of `1,8,7,56`. Thus the script's
sequential level processing does not split an equal-height packet across
different norm classes.

## Height-zero selection and scope

The exact dominant census in the search finds the scalar norm-one orbit and
the E8 norm-two orbit, and no other orbit. Therefore the nine negative
simple coordinate roots select the expected `E8+A1` height-zero chamber.
The code retains exact reflection-integrality filters in both blocks.
The census implementation was reviewed; it was not independently
reimplemented in this audit.

The numerical linear program only proposes a point. A spacelike vector
satisfying all archived wall inequalities, together with a strict future
timelike interior point, certifies an open ideal-boundary region in the
bounded prefix and hence infinite volume of that prefix. It says nothing
about later walls or the full reflection group.

## Independent verification

[`audit_fresh_cubic7_h9_no_import_20261013.py`](../scripts/audit_fresh_cubic7_h9_no_import_20261013.py)
imports only the Python standard library. It uses generic polynomial
convolution and reduction modulo the cubic, independent rational Gaussian
inversion, rational interval isolation of all three real roots, and a
separately specified E8 edge list. It verifies:

- host coefficient, embedding signs, and expected norm alphabet;
- integrality, exact norms, integral reflection coefficients, and mutual
  acuteness of every saved root;
- the nine initial simple roots and the declared height bound;
- exact nondecreasing level heights, each identity `height=a*control^2/norm`,
  and every level's height bound, including levels accepting no wall;
- consecutive level indices and wall counts, with each row's `new_walls`
  consumed in order from `all_walls`, and each such wall's time coordinate,
  norm, and height matched exactly to its row;
- a strict future timelike interior;
- the saved positive escape norm and every wall inequality.

For all three height-100 files, an altered escape norm, a timelike vector
falsely presented as an escape, an altered level height, reordered levels,
and walls assigned to the wrong level are all rejected. The order-corruption
control moves each root packet and its pairing data together, updating
indices and cumulative counts; chronology is still rejected. Height-10
controls needing two levels or two nonempty packets are explicitly marked
inapplicable when those conditions are absent.
The check is independent of search and field code;
it does **not** independently certify exhaustive enumeration up to the
height bound.

The three height-10 controls passed in
[`fresh_cubic7_h9_h10_independent_audit_20261013.json`](fresh_cubic7_h9_h10_independent_audit_20261013.json).
All three final height-100 artifacts passed in
[`fresh_cubic7_h9_h100_independent_audit_20261013.json`](fresh_cubic7_h9_h100_independent_audit_20261013.json).

```sh
python3 scripts/audit_fresh_cubic7_h9_no_import_20261013.py \
  results/fresh_cubic7_h9_seed0_h100_20261013.json \
  results/fresh_cubic7_h9_seed1_h100_20261013.json \
  results/fresh_cubic7_h9_seed2_h100_20261013.json \
  --output results/fresh_cubic7_h9_h100_independent_audit_20261013.json
```
