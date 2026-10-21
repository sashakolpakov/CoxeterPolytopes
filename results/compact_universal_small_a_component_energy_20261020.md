# Universal small-A component energy at a farthest vertex

20 October 2026. Complete paper candidate, awaiting independent review.
No computation, census, optimizer, build or formal verification was run.
There is no arithmetic hypothesis.

**Universal component bounds.** At a vertex farthest from an interior
point of a compact hyperbolic Coxeter polytope, every A3 component has
energy greater than 1/4 and every A4 component greater than 2/7.
The other finite components are arbitrary.

**Forest consequence.** In dimension n >= 14 the farthest vertex cannot
have whole finite diagram A1^a + A2^r + A3^s + A4^t. Other vertices are
unrestricted. This eliminates a pointwise family; it is not a uniform
dimension ceiling and does not include B, D, E, F or H components.

The farthest-point, reciprocal-shadow and no-common-opposite mechanisms
were co-developed by root, batch_guard and resource_cleanup. This source
uses the concise universal no-common-opposite proof at
`results/compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`,
SHA256 4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d.
It supersedes the scope, not the bytes, of the reviewed A3/A1/A2 note
`results/compact_a3_block_energy_and_low_rank_forest_exclusion_20261020.md`,
SHA256 70ab9f7e9fad80976b772e09388746db79c7bd161fc856bdc4ae0d38618b5068.
The universal small-height isolation below replaces its spectator lists.

## 1. Actual-wall tools and energy

Use outward unit normals n_i at the finite simple vertex p, positive
Gram B, inverse C, d_i = C_ii and u_i = omega_i/sqrt(d_i). Fix a strict
interior future unit point z and choose p maximizing E = -z.p over the
actual vertices. With z = E(p-x),

    ||x||^2 < 1,       x.u_i >= tau_i = tanh(L_i/2).     (1)

Here the actual adjacent endpoint is cosh(L_i)p - sinh(L_i)u_i.
The energies of the orthogonal finite components add to ||x||^2.

For an actual outside unit wall c write height A = -p.c > 0, h = A^2,
contacts beta_i = -c.n_i >= 0, and

    v_i = (C beta)_i/sqrt(d_i),       T_i = v_i/A,
    beta^T C beta = 1+h.

The actual-wall alternatives are

    T_i > 1: c closes edge i, tanh L_i = 1/T_i;
    T_i < 1: c does not close i, tanh L_i >= T_i;
    T_i = 1: impossible in a compact chamber.            (2)

To verify the nontrivial alternative, if a DISTINCT wall d with contacts
gamma and height D closes i, Schur complementation gives

    beta^T C gamma >= (C beta)_i (C gamma)_i/d_i.

Indeed C - C_.i C_i./d_i is zero on row and column i and equals the
entrywise nonnegative inverse of B_(I\{i}) on the other indices.
Acuteness gives AD >= beta^T C gamma; the actual endpoint formula
tanh L_i = D sqrt(d_i)/(C gamma)_i proves the shadow in (2).
This also proves the two-wall law

    tanh L_i tanh L_j >= K_ij := u_i.u_j
        if their actual opposite walls are distinct.    (3)

The no-common-opposite theorem in the pinned source proves that one
actual wall cannot close directions in two whole finite components.
Its direct proof is short: project the common opposite positively off
the retained old walls. For u_i.u_j = 0 its projected normal is
A p - (A/r_i)u_i - (A/r_j)u_j, r_i = tanh L_i. Its product with z
is strict negative, giving x_i/r_i + x_j/r_j < 1. But (1) gives a
lower bound tau_i/r_i + tau_j/r_j = 1+(tau_i^2+tau_j^2)/2 > 1.

We also need a universal cost fact. If c closes an edge in one component,
each OTHER whole component is retained in the finite white endpoint.
Within such a component c has at most one contact: two contacts and
the connecting old path form a cycle, whose unit Gram is nonpositive
on its all-ones vector. A nonzero contact contributes cost q = d beta^2
at least 1/4; in a nontrivial component q >= 1/3. For the latter, an
adjacent old bond has cosine at least 1/2, so its two-node inverse
diagonal is at least 4/3, and the full inverse diagonal is no smaller.

Consequently, if h < 1/4 and c closes any edge, c has NO contacts to
other whole components. Otherwise their cost q >= 1/4 > h makes
T_contact = sqrt(q/h) > 1; by (2) the same c closes that other contacted
edge, contrary to the no-common-opposite theorem. We call this the
small-height isolation consequence. It controls arbitrary spectators.

## 2. Small energy forces a single opposite in A3 and A4

Set H(t) = 2t/(1+t^2). It is increasing on (0,1), and log H is concave,
since its second derivative is -1/t^2-2(1-t^2)/(1+t^2)^2 < 0.
If a component has energy at most e, two of its directions of correlation
k satisfy tau_i+tau_j <= sqrt(2e(1+k)). Thus

    tanh L_i tanh L_j
       <= H(sqrt(e(1+k)/2))^2.                         (4)

For A3 adjacent k = 1/sqrt(3). At e = 1/4 the right side is
32(1+k)/(9+k)^2 < k: after clearing denominators this is 74k > 39.
It follows from 74^2 > 3*39^2. Hence (3) forces adjacent directions
to share their actual opposite, and ONE wall closes all three.

For A4 adjacent correlations are sqrt(3/8), 2/3, sqrt(3/8). At e = 2/7
the right side of (4) is 28(1+k)/(8+k)^2 < k. The cleared difference
is k^3+16k^2+36k-28, increasing for k > 0; at k = sqrt(3/8) positivity
is equivalent to 291sqrt(6) > 704. Squaring gives 508086 > 495616.
Again ONE wall closes all four directions.

For a common opposite, (1) bounds every actual tangent square by

    R_3^2 = 16/25,       R_4^2 = 56/81,                (5)

respectively. In particular h <= R_k^2 min_i v_i^2.

## 3. Universal A3 energy greater than one quarter

In path order,

    C = [[3/2,1,1/2],[1,2,1],[1/2,1,3/2]].

Assume energy <= 1/4, so c closes all three. Adjacent contacts are
impossible: at the third edge endpoint the retained white triangle is
nonpositive (affine at minimum contacts 3,3), contradicting its finite
Gram. If there are two contacts they therefore occur only at endpoints.
Both 3 give the forbidden affine four-cycle. In every other case their
values u,w have u,w >= 1/2 and one >= 1/sqrt(2). For the middle derivative,

    q - (16/25)v_2^2
      = (59/50)(u^2+w^2) + (9/25)uw
      >= 177/200 + 9sqrt(2)/100 > 1,

where q is the own-component cost. Since h >= q-1 this contradicts (5).

It remains to consider a single contact beta. If at an endpoint, its
own cost q = (3/2)beta^2 satisfies min_i v_i^2 = q/9. Thus

    q-1 <= h <= 16q/225,       q <= 225/209.

Therefore beta^2 <= 150/209 < 3/4: only native labels 3,4,5 remain,
and ALL have q < 1. Also h <= 16/209 < 1/4, so small-height isolation
removes every spectator. Then h = q-1 < 0, impossible.

For a middle contact q = 2beta^2 and min_i v_i^2 = q/3. Hence

    q <= 75/59,       beta^2 <= 75/118 < cos^2(pi/5).

Only labels 3,4 remain. They have q <= 1 and h <= 16/75 < 1/4;
again isolation gives h = q-1 <= 0. This proves the UNIVERSAL A3 bound.
No list or gap theorem for arbitrary finite spectator costs was used.

## 4. Universal A4 energy greater than two sevenths

The inverse in path order is

    C = (1/5)[[8,6,4,2],[6,12,8,4],
              [4,8,12,6],[2,4,6,8]].

Assume energy <= 2/7. A common wall c closes all four. Two contacts not
at the two old endpoints form a cycle in a retained white diagram after
deleting an old node outside their connecting path. This is nonpositive,
contrary to the finite endpoint. Thus multiple contacts can only be
at the two endpoints. Two 3 contacts give the forbidden affine five-
cycle. Otherwise u,w >= 1/2 and one >= 1/sqrt(2). For v_2,

    q - (56/81)v_2^2
      = (32/27)u^2 + (344/243)w^2 + (20/81)uw
      >= (230+15sqrt(2))/243 > 1.

All coefficients are positive; the displayed minimum allows the larger
contact at either endpoint. This contradicts h >= q-1 and (5).

For a single endpoint contact, q = (8/5)beta^2 and min_i v_i^2 = q/16.
Consequently q <= 162/155 and beta^2 <= 81/124 < cos^2(pi/5).
For the last comparison use 5*31^2 > 69^2. Only labels 3,4 remain,
both q < 1. We also have h <= 7/155 < 1/4. Small-height isolation
therefore gives the impossible equality h = q-1 < 0.

For a single inner contact, q = (12/5)beta^2 and min_i v_i^2 = q/6.
We obtain q <= 243/215 and beta^2 <= 81/172 < 1/2. Thus only label 3
remains, with q = 3/5 and h <= 28/405 < 1/4. Isolation gives the same
contradiction. This proves the UNIVERSAL A4 bound.

## 5. Universal short A1 classification

Call an A1 direction short if L < L_* = asinh(1/sqrt(2)). A closing wall
has h < 1/2. Write b for its released contact square and Q for its full
white cost. Then

    h = b+Q-1,       h < b/3,       Q < 1-2b/3.         (6)

Every positive spectator component has one contact of cost q >= 1/4.
If q > h, (2) makes c close that contacted direction, forbidden by the
no-common-opposite theorem. Equality q = h gives T = 1 and is also
impossible. Thus all spectator costs satisfy q < h.

There cannot be two spectators: then Q >= 1/2, while (6) implies
h < (1-Q)/2 <= 1/4, less than either spectator cost. If there is one,
its q = Q < h = b+q-1 forces b > 1, and (6) gives q < 1/3. The only
such nonzero cost is an A1 contact 3, with q = 1/4. With no spectator,
h = b-1 > 0 likewise forces b > 1. Therefore the ONLY short A1 walls
at ANY farthest vertex are

    I: isolated dotted A1 contact,       0 < h < 1/2;
    J: dotted A1 contact plus A1_3,       1/4 < h < 3/8. (7)

For J, in fact h > 1/3. Its heavy contacted edge contributes more than
A/2 to x.v_c, as follows from its reciprocal actual length. If h <= 1/3,
the shadow on the light A1 gives another contribution at least A/2:
its ratio is 1-sqrt(1-1/(4h)) >= 1/2. This contradicts the strict
interior inequality x.v_c < A. Its heavy short energy consequently
exceeds 9/100 (tanh^2 L > 4/13 > 3600/11881), while its light A1 shadow has energy
greater than 1/4 because h < 3/8 gives tanh^2 L > 2/3 > 16/25.

There is at most ONE short A1 direction. Any two distinct short walls
are orthogonal: their heights have product < 1/2 and spatial products
are nonnegative, so the Coxeter gap excludes any negative nonright
product. An I wall cannot be orthogonal to another I or J: disjoint
support gives a negative product, and any overlap contributes at least
1/2, greater than their height product. For two J walls, a shared heavy
contact contributes > 1, a heavy/light overlap > 1/2, and a sole light/
light overlap contributes 1/4, LESS than their height product since
both h > 1/3. All give nonzero products; disjoint support does too.
Each J closes only its heavy A1 direction. This proves the assertion.

Nonshort A1 energy is greater than 1/10, since
tanh^2(L_*/2) = (sqrt(3)-sqrt(2))^2 > 1/10. If J occurs, its two A1
components pay more than 9/100+1/4 > 2/10, so the full A1 energy is
at least a/10. Only I may create a debit; its height is treated below.

## 6. Universal A2 budget, with one named possible exception

Put phi = (1+sqrt(5))/2. Set sinh^2 L_0 = 31/50 and e_2 = 319/2000.
If both directions of an
A2 block have length at least L_0, their correlation 1/2 gives energy
at least (tau_1+tau_2)^2/3 > e_2. Here tau_0 > 173/500, and
4*173^2/(3*500^2) > 319/2000.

Consider a wall closing an A2 direction shorter than L_0. It cannot
have two contacts in that block: two 3 form affine A2, while the next
pair 3,4 has h >= sqrt(2)/3 and white residual <= 3/4, so
sinh^2 L >= 4sqrt(2)/9 > 31/50. The quotient increases with either
contact and any spectator cost. The same lower bound covers larger
pairs. If the closing wall contacts only the other node, its contacted
direction is shorter, so classify that direction first.

Let b = beta^2 be its single contacted square. The own cost is 4b/3,
and the complete spectator cost Q satisfies

    h = 4b/3+Q-1,       0 < h < (124/243)b,
    Q < 1-(200/243)b.                                  (8)

Every positive spectator cost must be LESS than h, by the universal
no-common-opposite theorem. For label 3, h < 31/243 < 1/4; hence no
spectator exists, but own cost 1/3 cannot give h > 0. For label 4,
h < 62/243, so every spectator must be an A1_3 of cost 1/4. Writing
Q = k/4, only k = 2 gives 0 < h < 62/243, but then h = 1/6 < 1/4,
again a contradiction. For label 5, Q < 1/2 allows at most one
spectator; its cost Q exceeds h because own cost phi^2/3 < 1. It is
therefore forbidden too. Label 6 contains the affine 3--6 path.

For native labels >= 7, b >= cos^2(pi/7) > 81/100, so (8) gives
Q < 1/3. Thus Q is zero or a single A1_3. The latter is short only
for b < 729/800, precisely the possible labels 7,8,9,10. Indeed
b_10 = (5+sqrt(5))/8 < 729/800, whereas
b_11 > 1-(2/7)^2 = 45/49 > 729/800, using pi < 22/7.
For a dotted release b > 1, Q < 43/243 < 1/4, so Q = 0.

We have therefore exhausted every short A2 profile, with arbitrary
finite spectators:

    isolated A2 contact, finite m >= 7 or dotted;
    A2_m + A1_3, finite m = 7,8,9,10.                   (9)

For isolated finite m the same wall closes the mate by (2), since its
mate T^2 = b/(4b-3) > 1. At m >= 8 its actual mate has
tanh^2 L >= 6sqrt(2)-8 > 400/841, hence half-tangent greater than 2/5;
the strict radical comparison follows from 2*841^2 > 1188^2. The block
energy exceeds 4/25 > e_2. For isolated dotted b, shortness gives
b < 243/200, and the distinct mate's shadow gives
tanh^2 L >= b/(4b-3) > 81/124 > 16/25, so its energy exceeds 1/4.
For the mixed profiles in (9), b <= b_10 < 11/12. Their mate white
cost b+1/4 > 1 excludes the same wall, and the shadow gives
tanh^2 L >= b/(4b-9/4) > 11/17 > 16/25. These too pay more than 1/4.

Thus EVERY A2 block pays more than e_2, except possibly an isolated 7
block. Its height h = 4b/3-1 lies in (0,1/12), where b = cos^2(pi/7).
The elementary bounds 623/1000 < cos(2pi/7) < 5/8 follow by substituting
in 8q^3+4q^2-4q-1, increasing on [1/2,1]. They justify all bounds on b
used above without a numerical approximation or census.

## 7. A uniform floor on the remaining directions

For a whole A_k block, k = 1,2,3,4, let K be the Gram of its normalized
dual directions. Its K^-1 1 is entrywise positive, and

    T_1 = 1,       T_2 = 4/3,
    T_3 = 5-2sqrt(3),       T_4 = (28-8sqrt(6))/5,
    T_k = 1^T K^-1 1 > 5k/12.                          (10)

These follow directly from K^-1 = diag(sqrt(d)) B diag(sqrt(d)).
For T_4 > 5/3 the strict comparison is 59 > 24sqrt(6), since
59^2-24^2*6 = 25. The smaller ranks satisfy their inequalities directly.
If all the block's actual half-tangent squares are at least theta,
Cauchy with the positive vector K^-1 1 gives energy >= theta T_k.
Thus a collection of these blocks of total rank R costs more than
(5R/12)theta. The coefficient is uniform over all four component types.

Suppose an isolated A2_7 block occurs in a whole A1..A4 forest, n >= 14.
Let c be its isolated wall, b = cos^2(pi/7), h = 4b/3-1. Consider the
actual opposite d of ANY direction in another component. It is distinct
from c. At d's finite endpoint the whole old A2 is white; it has at
most one contact, of cost q >= 1/3, or has no contact there.

If there is a contact, the spatial product with c is at least beta*gamma
times the off-diagonal inverse entry 2/3. Its square is at least bq/3.
Acuteness therefore gives D^2 >= (b/(3h))q. The squared derivative of
the chosen direction is at most its own component cost <= 1+D^2-q.
The ratio D^2/(1+D^2-q), increasing in D^2, gives

    tanh^2 L >= b/(b+6h) = b/(9b-6) > 5/9,

where monotonicity in q uses 1/3 <= q < 1, and b < 5/6.
If d has no contact to that A2, its product with c is -sqrt(h)D < 0,
so the Coxeter gap gives D^2 >= 1/(4h); hence
tanh^2 L >= 1/(1+4h) > 5/9. Every OTHER direction consequently has
half-tangent square greater than 1/5. Their total rank n-2 >= 12 and
(10) already require energy greater than (n-2)/12 >= 1. This is
impossible. The isolated 7 exception is therefore absent in this forest.

Next let an isolated dotted A1-short wall I have height h <= 1/5.
For any other opposite d, if its released-A1 contact is gamma = 0,
the same gap gives D^2 >= 1/(4h). If gamma > 0, then gamma^2 >= 1/4
and acuteness gives D^2 >= (1+h)gamma^2/h. In either case, since its
chosen direction derivative squared is <= 1+D^2-gamma^2,

    tanh^2 L >= 1/(1+4h) >= 5/9.

For the nonzero case, clear the inequality using
4hD^2 >= 4(1+h)gamma^2 >= 1-gamma^2. Thus all OTHER directions have
half-tangent square at least 1/5. Their rank n-1 >= 13 and (10) again
contradict energy less than 1. Such an I must therefore have h > 1/5.
Its actual isolated edge then has half-tangent square greater than
3/70: at h = 1/5 this is equivalent to 6*67^2 > 5*73^2.

## 8. The dimension-fourteen forest contradiction

Write n = a+2r+3s+4t >= 14. Section7 excludes isolated A2_7, so all
A2 energies exceed e_2 = 319/2000. The universal component bounds give
A3 energy > 1/4 and A4 energy > 2/7. If there is no isolated A1-short
wall I, Section5 gives A1 energy at least a/10, including the funded J
case. Since every per-rank rate is at least 1/14, their total is
strictly greater than n/14 >= 1, impossible.

If I occurs, it is the only short A1 direction and Section7 gives

    E_A1 > a/10-2/35.

Summing all components now gives the explicit lower bound

    ||x||^2 > n/14 + a/35 + (233/14000)r + s/28 - 2/35. (11)

For n >= 15 this exceeds 1 even after discarding the positive count
terms. At n = 14, a >= 2 makes the last A1 debit nonnegative. If a = 1,
then 2r+3s+4t = 13 forces s odd and hence s >= 1; its extra 1/28
exceeds the remaining debit 1/35. The case a = 0 cannot contain I.
Every possibility contradicts (1), proving the forest consequence.

The main new reusable mechanism is small-height isolation, obtained
from the actual cross-component opposite obstruction. It avoids
repeating finite spectator inventories. Larger A components and the
other irreducible finite types remain unclosed at the farthest vertex;
no bound for all compact Coxeter polytopes has been inferred.
