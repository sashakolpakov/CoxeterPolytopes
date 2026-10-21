# The first binary edge is uniform in the golden E8 family

20 October 2026. Complete paper argument, awaiting independent review.
No scientific execution, factorization, Pell search, Lean or TeX job.
The preceding full-base source and all frozen evidence are unchanged.

## Statement and inherited setting

Let O=Z[phi], phi^2=phi+1, and for every integer k>=-1 put

    u=phi^k,   P=phi^(2k+1),   b=4+P,
    L_k=E8(O) orthogonal-sum [[-2phi,phi*u],[phi*u,2]].   (1)

Take r=(u,2;0), c=(0,-1;0). Their Gram is
`[[2b,-b],[-b,2]]`, and they form an O-basis of the binary summand,
with determinant -u. The FULL base facet r is the compact Bugaenko
chamber by
`results/golden_e8_binary_family_full_base_facet_20261020.md`, SHA256
`ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd`.
In particular the intrinsic E8 vertex p0=(1,0;0) is an actual ambient
E8+A1(r) vertex.

**Theorem.** The outgoing edge releasing r and retaining the whole E8
has its first endpoint at c, for every k>=-1. There is no intervening
primitive binary reflection mirror of ANY allowed divisor norm. Its
endpoint is

    p1=(2,-phi*u;0),    p1.c=0,    p1^2=-2phi*b,

with actual whole germ E8+A1(c). Squarefreeness of b is not assumed.
This does not supply another chamber symmetry or prove uniform
nonreflectivity of the family.

## 1. A candidate mirror has one integral norm equation

Use the preceding note's necessary primitive norm result: every
primitive reflector has square2d with d in O totally positive and d|b.
Also b is coprime to2,3,5; hence d is coprime to2. A binary normal can
be written lambda*r+mu*c with lambda,mu in O. Its reflection integrality
against c gives

    d | (-b*lambda+2mu),     so mu=d*nu, nu in O.

Put e=b/d, also totally positive. Dividing its norm equation by2d gives

    e*lambda^2-b*lambda*nu+d*nu^2=1.                   (2)

A mirror cutting the OPEN segment from r to c must satisfy
lambda*mu<0 physically. Indeed orient it negative at p0; its product
there is mu*(p0.c), so mu>0. For it to separate p0 from p1 requires
lambda<0, since its product at p1 is lambda*(p1.r). Changing overall
normal orientation gives the same sign conclusion. Since d>0,

    T=-lambda*nu>0,        T in O, T!=0.               (3)

No assumption about the factorization or size of d has been used.

## 2. Both embeddings force N(T)=-1

Write q=sqrt(b)>2 and q'=sqrt(b') in(0,2). The prime denotes conjugation,
not differentiation. In fact

    b>=3+phi,     b'>=4-phi,     (b-4)(4-b')=1.        (4)

Physical arithmetic-geometric mean in(2), using ed=b, gives

    T<=1/(b+2sqrt(b))<1/b.                            (5)

If lambda'*nu'<0, the conjugate equation gives
`|T'|<=1/(b'+2sqrt(b'))`. Thus

    0<N(T)<=1/[(b+2sqrt(b))(b'+2sqrt(b'))]<1,

since N(b)>=11. This contradicts integrality of the nonzero norm.
Neither conjugate factor can be zero. Consequently lambda'*nu'>0,
T'<0, and conjugate arithmetic-geometric mean instead gives

    |T'|<=1/(2sqrt(b')-b').                            (6)

We claim the product of the denominators in(5),(6) is strictly greater
than1/2. The relation in(4) gives the exact identity

    (b+2q)(2q'-b') = q*q'/[(q-2)(q'+2)].              (7)

Its being greater than1/2 is equivalent to
`q*(2-q')<2(q'+2)`. But

    q*(2-q')=q/[(q^2-4)(q'+2)]
             <q/(q^2-4)
             <=sqrt(3+phi)/(phi-1)
             =sqrt(5phi+4)<4<2(q'+2).

Here q/(q^2-4) decreases for q>2, and phi<2 proves the last radical
bound. Therefore(5)--(7) imply

    0<|N(T)|<2,       N(T)<0,       N(T)=-1.           (8)

This is a norm bound for a single algebraic integer, not an unbounded
unit or divisor enumeration.

## 3. Golden-unit spacing contradicts the two bounds

By(8), T is a unit, positive physically and negative at the other
embedding. Since P=phi^(2k+1) has norm-1, U=P*T is a totally positive
unit. The physical strict bound in(5) gives

    0<U<P/b<1,       hence U<=phi^(-2).                (9)

The last step uses the explicit golden unit group: every totally
positive unit is phi^(2j), j in Z. On the other hand N(T)=-1 gives
T=1/|T'|. From(6) and(4),

    U>=P(2sqrt(b')-b')=q'/(2+q')>phi^(-2).            (10)

For the strict comparison, q'^2>=4-phi>4/phi^2, since
`4-phi-4/phi^2=3phi-4>0`. Thus q'>2/phi, exactly the inequality
equivalent to the final one in(10). Equations(9),(10) contradict.
No primitive binary reflector of any divisor norm crosses the segment.

## 4. Why this is the actual first edge

The actual E8+A1(r) corner supplies the outgoing binary axis. The norm2
mirror c is integral and meets it at the finite p1 displayed above.
Local finiteness gives a first mirror at or before p1. At a finite
point on this axis, the old eight actual walls retain their E8 simple
diagram. No connected finite Coxeter diagram properly extends E8;
thus the first actual opposite is orthogonal to all of E8 and is binary.
The scalar barrier excludes every earlier mirror. At p1 a binary
normal orthogonal to the same timelike point is proportional to c,
so there is no distinct simultaneous binary wall. The full finite germ
is E8+A1(c), not a merely formal positive definite candidate.

The firstness proof does not filter roots by an assumed finite prefix,
assume all primitive norms are2 or2b, or discard proper divisors when b
is composite. It uses neither an increasing-height search nor a Pell
fundamental unit. The inherited golden PID/unit facts and classical
finite E8 extension input remain explicit external premises, as in the
full-base note and its pinned binary audit
`bugaenko_facet_binary_extension_and_collar_audit_20261020.md`, SHA256
`774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`.

Uniform base calibration and this one finite edge do not imply that
the remaining full chamber is compact, finite volume, or nonreflective.
Transverse E7/other edge continuations and a second actual symmetry
remain separate obligations. Closed literal members remain closed;
no new host, class-wide exclusion, example or dimension ceiling is
promoted by this awaiting-review edition.
