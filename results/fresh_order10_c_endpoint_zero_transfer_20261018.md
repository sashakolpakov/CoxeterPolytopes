# The order10 affine C endpoint unit is unnecessary

In a whole varying affine C9 factor of order10, assign weight zero to
the pair of doubled terminal leaves. Every quadrilateral with such an
actual corner receives a full retained unit at another corner. The C
family cost decreases from45 to44 without any new finite coupon. All
other weights remain unchanged.

This is a local improvement to the
[completed dimension58 rule](fresh_a15_d15_dimension58_20261018.md).
The order10 family totals become A45, B45, C44, D45. Their common upper
envelope remains45/10=9/2; neither the finite forest budget nor the
current dimension bound improves from this change alone. The
[ambient58 frontier](fresh_ambient58_payment_frontier_20261018.md)
therefore remains an unchanged diagnostic for a possible bound57.

The [producer](../scripts/fresh_order10_c_endpoint_zero_transfer_20261018.py)
and [certificate](fresh_order10_c_endpoint_zero_transfer_20261018.json)
reconstruct all physical quantities below. This is the shorter-resistance
case of the independently proved
[order11 C endpoint theorem](fresh_order11_c_endpoint_zero_transfer_20261018.md),
with a fresh exact check at resistance9. It does not use the conditional
order10 theorem that assumes absence of smaller affine corners.

## Current weight and complete physical profiles

The [original simultaneous rule](fresh_deep_finite8_dimension64_20261018.md)
assigns a unit to every pair in an ideal factor of order at most11.
The subsequent changed ideal roles in the completed58 chain have order
at least11; no order10 pair was changed. The frozen frontier also records
the current rank9 envelope9/2. The new source saves all45 old and new
pair weights and checks reflection invariance. Only the unique terminal
pair(0,9), at distance9, changes. All distance-at-most-eight weights are
unchanged.

Use an actual quadrilateral star W+a+b+c+d, with opposite pairs
{a,c},{b,d}. At the selected corner the whole varying affine factor F is
the path0,...,9, with doubled bonds0--1 and8--9 and marks a=0,b=9.
Every physical root has squared norm2. By the
[actual factor-fiber framework](fresh_cusp_factor_payment_framework_20261017.md),
the other supporting roots form a positive-definite diagram Z orthogonal
to F. The positive kernel vector of F is

    nu=(1/sqrt2,1,1,1,1,1,1,1,1/sqrt2).

The opposite heights h=−c·nu and k=−d·nu are strictly positive. If
c·b, c·d or d·a is nonzero, a distinct actual corner has marked distance1
and pays a retained full unit. If c,d share a contact coordinate in Z,
their actual corner has marked distance2 and also pays a full unit.
These alternatives already prove the claim. Otherwise the three
adjacent-black products vanish and the nonnegative spectator contact
vectors x,y have disjoint supports.

Deleting either selected endpoint from F leaves a B9 principal diagram.
Its norm-two inverse has diagonal1 at the simple endpoint nearest the
omitted mark, diagonal2 at the next node, and values strictly greater
than2 everywhere else. The source verifies both inverse diagonals and
positive leading determinants. A one-contact positive-semidefinite
extension therefore permits only a simple or doubled contact at that
endpoint, or a simple contact at its next node. A label at least5 has
squared native magnitude greater than2 and fails the endpoint bound.

Two nonzero old contacts are also impossible: both must occur at the
only individually admissible positions, whose inverse diagonals are1
and2. The inverse is entrywise nonnegative, so the combined projection
is at least3>2, already contradicting principal semidefiniteness. The
source checks the complete inverse sign and the exact two-contact
minimum. Thus the principal inequalities give exactly the following
four profiles, without any exception from another finite or affine
type. A principal obstruction cannot be repaired by adding spectators.
No bound on their number or contacts is imposed.

Let t denote the magnitude of the contact with the omitted old endpoint.
It can be zero or any native magnitude at least1, including parallel and
dotted products. The profile data are

|Contact on B9|Height h|Moment l|Constant c0|Native t range|Maximum g=2l/h+c0/h²|
|---|---|---|---:|---|---:|
|none|t/sqrt2|0|2|t≥1|4|
|endpoint3|t/sqrt2+1|1|1|t≥0|3|
|endpoint4|t/sqrt2+sqrt2|sqrt2|0|t≥0|2|
|next-node3|t/sqrt2+1|2|0|t≥0|4|

For the contact-free profile, positive height forces t>0 and hence t≥1.
Every other profile has positive height already at t=0. The producer
checks each upper bound on the whole corresponding continuous ray:
after t=t_min+u, the numerator of maximum−g has nonnegative coefficients
and the denominator is strictly positive for u≥0. In particular, the
endpoint4 and next-node3 affine alternatives remain in the calculation.
No alternative affine payment is spent or presumed positive.

## Strict determinant gap for arbitrary spectators

Eliminate the eight white roots of F and scale a,b by their positive
kernel coefficients. The reduced old-mark matrix is

    (1/9) [[1,−1],[−1,1]].

The source reconstructs this matrix, all four heights/moments/norm
constants, and all16 reflected cross identities. The exact determinant
identity after eliminating Z is

    0 = 1−(g_c+g_d)/9
          +(x/h−y/k)^T G_Z^(-1)(x/h−y/k)/9.          (1)

The actual-star factor-fiber theorem supplies the zero determinant:
after the supporting whites are eliminated, the four side roots lie
in the three-dimensional face quotient. The producer independently
expands the symbolic four-by-four determinant to verify(1).

The spectator quadratic form is nonnegative because G_Z is positive
definite. Every complete profile has g≤4, so the right side of(1) is at
least

    1−8/9 = 1/9 >0.

This is impossible. The16 exact ordered bounds have uniform minimum1/9.
Consequently an adjacent-black contact or a shared spectator coordinate
must occur, and a distinct actual corner supplies the full distance1
or2 unit identified above. This conclusion covers arbitrary native
opposite products, arbitrary supporting rank and all permitted spectator
attachments. It neither infers an empty spectator diagram nor relies
on a selected finite catalogue.

## Current composition and counting scope

For an arbitrary quadrilateral, first test whether any actual corner is
one of the newly zero C10 endpoint pairs. If so, the independent proof
above supplies a distinct full short unit. Those distance1/2 payments
are unaffected by every current ideal change, same-component finite
deletion and whole-component coupon transfer. If no changed C10 corner
occurs, every actual payment equals its completed58 value, so the full
baseline face theorem applies pointwise. This also handles any older
proof that happened to choose a C10 endpoint payment: the existence of
that actual corner activates the first branch. There is no circular
appeal to the removed unit.

The [triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
uses finite paths of length at most7, ideal paths of length at most5 or
protected ideal paths6/7, and whole affine factors of order at most9.
The changed order10 distance9 pair is outside these outputs. Triangles
therefore retain total at least2 and quadrilaterals total at least1.
All earlier finite6/finite7 and deep finite8 routes remain available;
the only new output needed here is a short unit.

No finite weight or coupon changes, so every previously certified finite
upper bound remains valid. The old rank9 all-family envelope9/2 is still
valid and attained by the unchanged A/B/D families. The current ideal
profile certificate and its low9 concentration argument may therefore be
retained as they stand. A future reduction of that all-family envelope
will require a new concentration justification; no such reduction is
claimed here. This packet does not regenerate a forest budget or an
ambient-dimension certificate.

## Verification

Final source-bound generation42bd6f and separate fresh saved replayf33543
both terminated exit0/PASS. The parent and a peer completed full
mathematical source/report reads. The peer identified an imprecise draft
cycle description; it was replaced before the final runs by the direct
two-contact inverse-projection exclusion above. The peer then read the
final source/report delta and passed the corrected complete profile
domain. The parent also read the corrected source/report delta and
passed the direct inverse exclusion. These reviews covered the native-ray bounds,
Schur signs, arbitrary spectator term, current actual-corner priority,
unchanged short units, triangle interface and unchanged all-family
envelope. Neither reviewer executed the producer or independently
replayed the inverses, profiles, determinants or historical census.

The producer checks the two B9
principal inverses, all four complete physical profiles, unbounded-ray
cost bounds,16 cross identities and strict pair bounds, symbolic Schur
identity,45 old/new pair records and exact family costs. The all-label
finite/affine classification and actual-star factor-fiber theorem remain
explicit inherited mathematical inputs. No external refereeing,
proof-assistant verification or new finite-volume example is claimed.
Direct metadata audit passes15 source bindings,8 report links and source
syntax. Source, certificate and report are final/frozen.

    python3 scripts/fresh_order10_c_endpoint_zero_transfer_20261018.py
    python3 scripts/fresh_order10_c_endpoint_zero_transfer_20261018.py --verify
