# Uniform golden first binary edge: outside mathematical review

20 October 2026. Reviewer: batch_guard. Complete independent paper audit.
No scientific execution, root/divisor search, Lean/TeX job or source edit.

**Verdict: PASS. No gap found in the uniform first-edge theorem.**

Bound source:
`results/golden_e8_binary_family_uniform_first_edge_20261020.md`,
SHA256 `71cfa11bcadf5fffae017d712d168eeec44eecf7d5b74f595bcebca7ad92af2c`.
The inherited full-base theorem remains at SHA256
`ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd`;
its outside review is
`review_golden_e8_binary_family_full_base_facet_batch_guard_20261020.md`,
SHA256 `8bf116807614a38d6215e5a3a3da8c49f4fae52165cce99d88fda171115b6580`.
I did not author the new unit-spacing argument and rederived it directly.

## All divisor classes and both embeddings

The unit (r,c) basis makes lambda,mu integral. For arbitrary totally
positive d|b, integrality against c and gcd(d,2)=1 imply mu=d*nu.
The resulting equation

    (b/d)lambda^2-b*lambda*nu+d*nu^2=1

is exact. The sign lambda*mu<0 follows by evaluating the oriented
mirror at p0 and p1. Hence T=-lambda*nu is a positive nonzero algebraic
integer. Nothing in this reduction restricts d to1 or b or uses
squarefreeness.

Physical AM--GM gives T<=1/(b+2sqrt(b)). The negative-product conjugate
case would give a positive nonzero integer N(T)<1, and is impossible.
The other sign yields |T'|<=1/(2sqrt(b')-b'), with a strictly positive
denominator because0<b'<4. The possibility of a zero conjugate factor
is excluded by lambda,nu being nonzero field elements.

Writing q=sqrt(b),q'=sqrt(b'), I independently expanded

    (b+2q)(2q'-b') = qq'/[(q-2)(q'+2)].

It follows exactly from (b-4)(4-b')=1. The stated comparison with1/2
has the correct direction: it reduces to q(2-q')<2(q'+2). The
decreasing q/(q^2-4) bound at b_min=3+phi proves this strictly,
including k=-1. Thus -2<N(T)<0 and integrality gives N(T)=-1.

For P=phi^(2k+1), U=P*T is a totally positive unit. Its physical
upper bound is strictly below1, so golden unit spacing gives
U<=phi^-2. Conversely N(T)=-1 and the conjugate bound give

    U>=P(2q'-b')=q'/(2+q')>phi^-2.

The equality uses P(4-b')=1. The strict final comparison follows from
b'>4/phi^2, since b'>=4-phi and3phi-4>0. These are uniform inequalities,
not a numerical scan, a fundamental-relative-unit assumption, or a
factorization-dependent norm pool.

## Actual edge and exact endpoint germ

The preceding full-base theorem supplies the actual E8+A1(r) corner.
The integral norm2 c mirror gives a finite endpoint on its outgoing
binary axis. Local finiteness supplies a first mirror. At an earlier
finite endpoint it would positively complete the retained actual E8
diagram; finite spherical classification forces orthogonality to E8,
so it is binary and covered by the scalar contradiction. At p1 every
binary zero normal is proportional to c. The same finite-extension
argument excludes a nonbinary additional actual zero wall. Hence the
endpoint germ is exactly E8+A1(c), with no omitted simultaneous wall.

This settles the one outgoing edge for every k>=-1 and every allowed
primitive divisor norm. It does not imply compactness, finite volume,
a second chamber automorphism, uniform nonreflectivity, or exclusion
of the infinite family. All earlier frozen sources remain unchanged.
