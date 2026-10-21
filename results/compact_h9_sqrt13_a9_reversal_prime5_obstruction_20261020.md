# The quadratic13 A9 reversal fails at prime5

20 October 2026. Own paper obstruction, no execution. The actual A9
corner is supplied by
`compact_h9_sqrt13_a8_edge_interval_20261020.md`, SHA256
`e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228`.
The complete primitive norm alphabet1,2 is inherited from the reviewed
B9 route b9226479. This note obstructs one proposed period mechanism;
it does not obstruct the lattice, compactness, or other chamber maps.

Put eta=(3+sqrt13)/2 and L=<-eta>+I9. The primitive future apex is

    Q=(9eta+1;Y,...,Y),    Y=5eta+2,    Q^2=-10eta.

It is primitive because N(9eta+1)=-53 and N(Y)=9 are coprime: a common
prime ideal dividing every coordinate would have rational norm prime
dividing both integers. The actual A9 roots are the chain
(a,alpha1,...,alpha8), all of norm2.

Let w0 be the longest element of this finite A9 reflection group. It is
an integral lattice isometry, since every defining reflection has an
integral norm2 normal. On the root span it sends each simple root to
the negative of its diagram reverse; it fixes Q.

Consequently the UNIQUE future-preserving real isometry extending the
A9 diagram reversal while fixing its vertex is

    T=J_Q*w0,
    J_Q(x)=-x+2(x,Q)Q/(Q,Q).

Indeed J_Q is -identity on the nine-dimensional root span and identity
on Q. A future isometry preserving the Q ray fixes Q itself, because
Q has nonzero norm. There is no free complementary sign or scalar
that can repair this extension.

J_Q is not integral. For any spatial coordinate vector e_i,

    J_Q(e_i)=-e_i - [Y/(5eta)]Q.

The coefficient has a genuine prime5 denominator: eta is a unit and
N(Y)=9 is coprime to5. Since Q is primitive, its coordinates cannot
cancel this denominator simultaneously. Equivalently, the integral
rank-one formula in a unimodular lattice would force (Q,Q) to divide2,
which -10eta does not. As w0 is integral with integral inverse, T is
not an integral automorphism either.

There is a stronger arrangement conclusion, rather than merely failure
of one lattice basis. J_Q(e_i) has norm1 and is not integral. If its
mirror were nevertheless an integral reflection mirror, the complete
norm alphabet would give an integral normal lambda*J_Q(e_i) of norm1
or2. Since both vectors are defined over K, lambda lies in K. It must
satisfy lambda^2=1 or2. The first only allows signs, which do not clear
the denominator. The second is impossible because sqrt2 is not in
Q(sqrt13). Thus this image is NOT a mirror of the full arrangement.
Because w0 permutes that full arrangement, T cannot permute it either.

The automatic E8 diagram-reversal mechanism used in earlier excluded
hosts therefore does not transfer to this A9 corner. The existing B9
chamber involution remains integral of order2. No claim is made about
the order of any different explicitly constructed chamber map; a new
period would require an actual map and a new full-germ proof.
