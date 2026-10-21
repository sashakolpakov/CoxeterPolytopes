# Independent complete contact-domain review at the H20 null direction

**PASS:** all norm-two acute Coxeter walls through the specified null
direction have one of 125 finite contact profiles on the 13 spherical
incident roots, or belong to two continuous single-contact families.
This is a necessary-domain theorem, not an escape-cutting exclusion.

The [independent checker](../scripts/fresh_h20_null_link_profiles_20261016_review.py)
and [certificate](fresh_h20_null_link_profiles_20261016_review.json)
pass generation and fresh saved replay. They agree with the
[primary census](fresh_h20_null_link_profiles_20261016.md) on every accepted
profile and every exact cost. The
[independent physical cusp review](fresh_h20_affine_a17_completion_20261016_review.md)
supplies the old null incidence and its positive affine G2 relation.

If r.nu=0, the positive relation and acuteness force r to be orthogonal
to all three affine G2 roots. The remaining 13 old incident roots form
A2+A5+A4+2A1. In the positive Euclidean quotient nu-perp/R*nu, let
d_i=-r.n_i>=0 and J be the inverse of their Gram matrix. Projection gives
`d^t J d<=r^2=2`.

The reviewer constructs that Gram directly and inverts it independently.
All inverse entries are nonnegative. Its four smallest diagonal entries
are 1/2, 1/2, 2/3, 2/3. Nonzero Coxeter contacts have magnitude at least
one, so four nonzero contacts would cost at least 7/3>2. Dotted contacts
have magnitude greater than two and are impossible by the minimum
diagonal 1/2.

A nonisolated root has inverse diagonal at least 2/3, so its contact
magnitude is at most sqrt(3). An isolated root accompanied by another
nonzero contact satisfies d_i^2/2+1/2<=2 and has the same upper bound.
These cases have labels at most six. Only a single contact at isolated
root 17 or 31 can exceed sqrt(3); the interval [sqrt(3),2] covers every
remaining finite label and the parallel endpoint.

The independent enumeration checks all 19,604 assignments with one,
two or three nonzero contacts and labels three through six. It uses
recursive exact sign comparisons in the ordered tower
Q(sqrt(2))(sqrt(3))(sqrt(5)), rather than the producer's radical interval
comparison. Direct Gram inversion and symmetric quadratic-form summation
also differ from the producer's path-inverse formula and ordered sum.

| Number of nonzero contacts | Accepted profiles |
|---:|---:|
|0|1|
|1|29|
|2|73|
|3|22|

There are 22 exact equality cases with cost two; they are retained and
listed in the certificate. Both continuous tails deliberately include
non-Coxeter values and overlap the finite census at label six. No
finite-label cutoff is imposed on the actual candidate walls.

The theorem does not assert that any retained profile is geometrically
realizable, satisfies the other old wall contacts, or cuts an escape.
Those are separate obligations.

```sh
python3 scripts/fresh_h20_null_link_profiles_20261016_review.py --verify
```
