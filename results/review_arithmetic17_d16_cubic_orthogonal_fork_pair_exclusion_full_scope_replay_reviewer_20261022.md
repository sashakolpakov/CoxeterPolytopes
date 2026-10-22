# Full outside-author review: cubic D16 plus actual I

22 October 2026. **PASS at the four explicit written-input boundaries.**
This is a complete manual mathematical review of
`arithmetic17_d16_actual_i_cubic_orthogonal_fork_pair_exclusion_candidate_20261022.md`,
SHA256 `071de627b87d6a6be80d93af3310e1eef6a7dfc8b0109d11b31c688a14a78ed4`,
8,349 bytes, 184 lines. I read the entire proof and the relevant complete
input sources, checked the literal input identities, and independently
rederived the new argument. No mathematical worker, worker import,
enumeration, symbolic/numerical process, test, Lean or TeX ran. File reads
and streaming metadata hashes are not mathematical execution evidence.
The reviewed source and all bound inputs remain unchanged.

The author is the foundations agent, not this reviewer. I previously
reviewed the old-I contact theorem, authored the separate quartic
norm/unit exclusion, and contributed inherited actual-foot and field
normalization work. I also supplied a preliminary crosscheck of the new
cycle argument before receiving the frozen proof. These contributions
are disclosed: this report independently checks the author's complete
new cubic exclusion; it does not claim independent authorship or a fresh
referee proof of my own earlier foundations.

## 1. Exact premises and scope

All four literal hashes match the current bytes:

| Input under `results/` | SHA256 | Bytes | Interface checked |
|---|---|---:|---|
| `arithmetic17_d16_actual_i_cubic_old_i_contact_three_exclusion_candidate_20261022.md` | `45ee50551d069126849ab3fc4393e1ad99e5b31beaa834c71162110ba4b3fec0` | 13909 | In actual degree three, both actual forks are D-isolated, dotted, and orthogonal to the old I root |
| `arithmetic17_d16_actual_i_joint_head_height_and_shorter_forks_candidate_20261022.md` | `3a5c21e33622c62bb77a81d8cdbe394e1179001cc7a19a984cbbf6ebd6030da4` | 8316 | Both actual fork cosh lengths are strictly less than 37/4 |
| `arithmetic17_d16_actual_i_actual_field_degree_four_candidate_20261022.md` | `f0b9b9d8746719ac5a7c1c6af480bfe5fb6581be02376381329a3bbeb1af4ba4` | 6487 | An isolated dotted fork has integral T=4gamma², actual K=Q(T), and every other conjugate lies strictly between zero and 1/2 |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | 9509 | Actual K-valued root lines, invariant Cartan cyclic integrality, and positive definite nonphysical forms |

The first input's conclusion is exactly strong enough: it removes old-I
contacts for both forks in degree three, without assuming anything about
the opposite of I. I reread its full argument, including the unequal
parameter cases; its separate accepted scope is not enlarged here.
The second input's fork bound is independent of its optional final
nonisolated-role corollary. The third input's principal-Gram descent
identifies the actual field, rather than just a convenient subfield.
An extra embedding fixing Q(T) would retain the negative principal
Schur complement inside a positive definite conjugate form, so no hidden
larger actual field is allowed. The historical first-type wording of
the fourth source is not an additional subclass restriction on arithmetic
reflection lattices.

The theorem reviewed here excludes the entire **actual degree-three**
branch of the whole D16 plus actual isolated short-I forest in dimension
17. It does not exclude the whole forest or prove an arithmetic ceiling
16. The already accepted quartic exclusion is a separate input to the
current research status, not a premise used to prove this cubic theorem.

## 2. Actual normalizations and the two cycles

I checked the unit-to-raw conversions. Normalize the connected old D
roots over K to square two and each cap to raw own product minus one.
The cap has raw square 1/(2gamma²)=2/T. Its old-D projection has square
four; the two fork duals have mutual product 7/2. Orthogonality to the
old I root removes its spatial component regardless of the old-I raw
square. Because the entire spatial forest is D16 plus I, both remaining
cap components lie in the same one-dimensional time line. The resulting
actual Gram identity is therefore

    (P-7/2)² = 4(2T1-1)(2T2-1)/(T1*T2).

This is justified for actual cap roots, not for arbitrary completed Gram
matrices. The two walls are distinct because their nonzero old-D own
contacts occur at different fork roots. Their parameters need not agree.

The actual five-edge path through the fork roots and their common
branch root has Cartan factors

    -T1, -1, -1, -1, T2*P.

Their product is X=T1*T2*P. The actual cap two-cycle is
L=T1*T2*P². Both are integral invariant cyclic products. If P is zero,
both assertions hold immediately. No individual raw Cartan entry is
declared integral, and no relation between disconnected D/I root scales
is assumed.

Putting Q=T1*T2 and expanding the determinant identity gives

    L-7X+49Q/4 = 16Q-8T1-8T2+4,
    L-7X = 15Q/4-8T1-8T2+4.

Thus 15Q/4 is integral and Q/4=4Q-15Q/4 is integral. This cancellation
uses only integral closure and elementary ring arithmetic; it needs no
assumption that O_K modulo two is reduced and no ramification case list.

## 3. The norm bound and the complete coefficient strip

Old-I orthogonality gives C²=2T, so C<37/4 gives T<1369/32.
At the other two embeddings, 0<T_sigma<1/2. Consequently each norm
n_i is a positive integer at most ten. The cubic norm of Q/4 is
n1*n2/64, an integer. The product is positive and at most 100, hence
equals 64. The only factor pair with both factors at most ten is (8,8).
Both cap parameters therefore have norm eight.

For either parameter let its two small conjugates be x,y. Then Txy=8
and x,y<1/2 imply T>32. The function t+2sqrt(8/t) is increasing on
this range, so the integral trace a is strictly greater than 33 and
strictly less than 44. With b=T(x+y)+xy, the independent necessary
bounds are

    34<=a<=43,  b²>32(a-1),  2b<=a+31.

For the second inequality use T>a-1 and b>2sqrt(8T). For the third,
P(1/2)=b/2-a/4-63/8 is negative; converting this strict inequality
to its integer form gives exactly the displayed non-strict bound.
I checked all ten trace rows. The lower/upper b intervals are respectively

    (33,32), (33,33), (34,33), (34,34), (35,34),
    (35,35), (36,35), (36,36), (37,36), (37,37).

Thus the five pairs in the source are the complete necessary list, not
sampled polynomial choices. Direct substitution in the cubic discriminant
formula gives, in the source's order,

    -17131, -16124, -11651, -2704, 11821.

All values and signs check. A primitive totally real cubic generator has
positive discriminant. Hence both parameters have minimal polynomial
X³-43X²+37X-8. Each actual parameter has one physical conjugate above
four and two conjugates below one half. The polynomial therefore has a
unique large root, and both physical parameters must be that root.
Injectivity of the physical embedding gives T1=T2 in K. Equality is
derived only here; it has not been assumed in the cycle or norm argument.

## 4. The final integral-closure contradiction

With equality established, Q/4=T²/4 is integral. Since T/2 belongs to
K and satisfies a monic polynomial over O_K, integral closure gives
U=T/2 in O_K. I independently checked the reproduced no-even-T lemma.

Its conjugates satisfy U<1369/64 and 0<x,y<1/4. The positive integral
norm is less than 1369/1024, hence one. Therefore U>16. Setting
r=1/x and s=1/y gives r,s>4 and rs=U, so

    r+s < 4+U/4.

The trace of the integral unit U^(-1) is strictly between eight and
2409/256<10, so equals nine. The minimal polynomial of U consequently
has the form X³-aX²+9X-1. Its negative value at one quarter gives
a>81/4, whereas a<U+1/2<1401/64<22. Thus a=21. Direct evaluation of
the discriminant gives minus 864, contradicting total reality. This
exhausts the forced equal-pair case and completes the proof.

## 5. Verdict, status wording, and limits

**Full mathematical SOURCE PASS, with no correction required.** Every
new cycle identity, normalization, integer bound, coefficient-table row,
discriminant, equality step, and integral-closure implication checks.
All four literal premises match and their used scopes are sufficient.
There is no hidden common fork parameter, relative D/I squareclass,
contact restriction on the I opposite, auxiliary farthest vertex, or
addition of energies from different vectors.

The source's standalone statement that degrees two/four are unresolved
is read with its explicit qualifier **by this note**. It earns no claim
that the quartic branch remains globally open: the separately accepted
quartic theorem already excludes it. After root accepts this cubic
theorem, composing those results with the accepted degree bound leaves
only actual degree two for this forest. The forest itself remains open;
no global bound, example, reflectivity result, runtime or Lean credit
follows from this review. No frozen manuscript, guide, receipt, or source
was edited, and this report does not itself perform root acceptance.
