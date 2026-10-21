# Root review of the uniform first binary edge

20 October 2026. Complete independent written review **PASS** of
`golden_e8_binary_family_uniform_first_edge_20261020.md`, SHA256
`71cfa11bcadf5fffae017d712d168eeec44eecf7d5b74f595bcebca7ad92af2c`.
No scientific computation, factorization, worker import, Lean or TeX
was used in this review. The accepted full-base source is unchanged.

I checked the source in full, including all proper-divisor norm classes
and the actual first-endpoint argument. In the unit binary basis (r,c),
integrality against c gives d|2mu, hence mu=d*nu, since d|b and b is
coprime to2. The norm equation is exactly

    (b/d)*lambda^2-b*lambda*nu+d*nu^2=1.

At p0, c has product -phi*u<0 and r has zero product. At p1, r
has product -phi*u*b<0 and c has zero product. A distinct mirror
separating these two endpoints consequently has lambda*nu<0.
Neither factor is zero. Thus T=-lambda*nu is a nonzero algebraic
integer positive at the physical embedding. The choice of a totally
positive divisor d is legitimate because its square norm and b are
totally positive.

Physical AM--GM, with product of the two diagonal coefficients equal
to b, gives T<=1/(b+2sqrt(b)). If the conjugate product lambda'*nu'
were negative, the same bound at that embedding would make the
nonzero positive integer N(T) strictly less than1. Thus T'<0 and
the conjugate inequality is |T'|<=1/(2sqrt(b')-b'). The denominator
is positive since 0<b'<4.

I independently checked the delicate uniform norm bound. Put q=sqrt(b)
and q'=sqrt(b'). The identity (b-4)(4-b')=1 gives

    (b+2q)(2q'-b')=q*q'/((q-2)*(q'+2)).

For this product to exceed1/2 it suffices, and is equivalent, to
q*(2-q')<2*(q'+2). The displayed chain in the source proves this:
q/(q^2-4) is decreasing for q>2; its value at the lower endpoint
b=3+phi is sqrt(5phi+4)<4, whereas 2*(q'+2)>4. Every denominator
is positive. Hence |N(T)|<2, with N(T)<0, so N(T)=-1 exactly.

P=phi^(2k+1) has norm-1 and negative conjugate also at k=-1.
Therefore U=P*T is a totally positive unit. Since U<P/b<1, the
explicit golden unit group gives U<=phi^-2. Conversely T=1/|T'|
and the conjugate bound give

    U>=P*(2sqrt(b')-b')=sqrt(b')/(2+sqrt(b')).

This exceeds phi^-2 precisely when sqrt(b')>2/phi. The uniform
lower bound b'>=4-phi suffices, as
4-phi-4/phi^2=3phi-4>0. Thus the source's contradiction is strict
throughout the family. It neither enumerates divisors nor assumes
that b is squarefree or prime.

The geometric endpoint statement also has the required scope. The
accepted full Bugaenko base supplies an actual E8+A1(r) corner.
The norm2 vector c is an integral reflector, and p1=(2,-phi*u;0)
is timelike with norm -2phi*b. The segment to it is finite. Local
finiteness supplies an earliest cutting mirror. The retained eight
actual walls form E8 at that finite endpoint; the finite spherical
classification has no connected proper extension of E8. The actual
opposite is therefore orthogonal to E8 and binary, rather than merely
being chosen from a convenient binary subset. The scalar argument
excludes every earlier one. A binary wall through p1 is proportional
to c. Hence the full germ is E8+A1(c), including simultaneous walls.

The imported inputs remain the reviewed full-base theorem, golden
PID/unit facts, local finiteness and finite spherical classification.
The conclusion settles this one outgoing edge for all k>=-1, with
all allowed divisor norms retained. It does not settle other edges,
compactness, finite volume or nonreflectivity of the infinite family.

Before the source was written I suggested reducing an intervening
binary root to the integral norm equation. The author supplied the
small-norm product and unit-spacing argument. The full proof above
was reconstructed here; this contribution is disclosed rather than
described as unrelated independent invention.
