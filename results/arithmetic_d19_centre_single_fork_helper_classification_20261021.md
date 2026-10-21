# A complete weak single-fork helper classification at the D19 centre

21 October 2026. Supporting written proof by root, pending full
outside-author review. No scientific, enumeration, optimization, TeX or
Lean process. This is a lemma for the existing arithmetic paper, not a
new manuscript or a dimension18 theorem.

## Statement and inputs

Let a farthest vertex of a compact arithmetic first-type hyperbolic
19-polytope have WHOLE finite diagram D19. Use unit old normals

    n_i=(e_i-e_(i+1))/sqrt2 (i1..18),
    n19=(e18+e19)/sqrt2,
    omega17=sqrt2 P17,
    omega18=(P18-e19)/sqrt2, omega19=(P18+e19)/sqrt2.

Suppose the actual opposite wall c at centre17 has exactly one white
contact: a simple3 contact with fork f18 or19. Write

    c=A p-v, v=gamma omega17+omega_f/2, A>0.

If h17=x.omega17<=sqrt34-21/40, then precisely the following necessary
possibilities remain:

    own5: gamma=phi/2 and actual K=Q(sqrt5);
    isolated from the other whites, dotted own contact:
          gamma=(1+sqrt2)/2 and actual K=Q(sqrt2).

This statement does not assert realization of either possibility. It
does not cover a left-white contact or a second fork contact. It does
not assert that either remaining case is impossible.

The geometric inputs are compact acute simplicity, the actual all-wall
nonpositive products, native Coxeter angles, the farthest-vector
normalization and the actual-edge Schur formula. Their internal proof
home is `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`
(3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84)
and `compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`
(e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c).
Field normalization, cyclic integrality and the all-conjugates criterion
are used exactly as proved in
`arithmetic_uniform_classical_weak_helper_obstruction_20261021.md`
(535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905).

## Eliminate the small own labels

If gamma=0, the actual tree contains an affine star with arms(1,2,5):
centre17, short fork other than f, the arm f--c, and five left old nodes.
Its positive null combination is feasible and future, contradicting
compactness. Own3 gives the affine3/3/3 triangle c,17,f. Own4 gives
the affine star centred at17 with leaves c,16 and the OTHER fork,
whose squared cosines sum to1. Own6 gives the affine6--3 path
c--17--otherfork. These principals retain all actual pair products.
Own parallel gives the feasible null sum c+n17 and is impossible.

The remaining possibilities are own5 or a native contact with
gamma>9/10, including every finite label>=7 and every dotted contact.
In the own5 case the native3 helper normalizes c to square2 over K.
The full old D19 plus c has entries in Q(phi), has rank20 and is
Lorentzian. A nonphysical embedding of K fixing Q(phi) would preserve
that Lorentzian Gram. Thus actual K=Q(phi), giving the first row.

## The smaller physical strip

Assume gamma>9/10 and put t=2gamma. With q=sqrt34 and c0=21/40,
the exact edge formula has

    cosh L/q=(gamma+1/4)/sqrt(3/4).

The sufficient ratio for h17>q-c0 is

    R(q)=1/c0-(q-c0)/(q(2q-c0)).

This function increases here and q<35/6. At35/6 it is
40/21-3822/46795<73/40: denominator clearing has positive margin75215.
Furthermore73/40<19/(6sqrt3), since575532<577600. Therefore gamma>=4/3
would make the actual edge ratio at least19/(6sqrt3)>R(q), contradicting
weakness. Consequently

    9/5<t<8/3.                                             (1)

The simple3 helper makes t an algebraic integer over the actual field,
without a square-root extension. For the norm-two old D19 inverse,

    (G^-1)_17,17=17, (G^-1)_17,f=17/2, (G^-1)_f,f=19/4.

At EVERY nonphysical embedding strict Schur positivity therefore gives

    17 t_sigma^2+17 t_sigma+11/4<0,
    -4/5<t_sigma<-1/5.                                    (2)

Indeed the roots are -1/2 plus or minus sqrt(3/34), strictly inside
the displayed rational interval. Every nonidentity embedding of Q(t)
lifts to K, so (2) holds for all other conjugates of t.

## A finite norm-and-trace argument

Let d=[Q(t):Q]. All nonphysical conjugates are negative; write them
as -a_j with1/5<a_j<4/5. Every norm below is a nonzero integer.

For d>=6, its magnitude is less than
(8/3)(4/5)^5=8192/9375<1, impossible.

For d5 the positive norm is less than2048/1875<2, so it is1.
The four a_j have product1/t. Hence t>625/256 and
t-sum a_j>625/256-16/5>-1. AM--GM gives sum a_j>=4/t^(1/4)>t,
because t^5<(8/3)^5<256. The integer trace is strictly between-1
and0, impossible.

For d4 the negative norm is-1. The trace lies in(-1,1): the lower
bound follows from t>9/5 and sum a_j<12/5; the upper bound follows
from their product1/t>3/8 and hence sum a_j>2. Its value is0, so
t=a+b+c and t^4>=27. The integral second coefficient is
-(a^2+b^2+c^2+ab+ac+bc), strictly greater than-96/25>-4 and
at most-(2/3)t^2<=-2sqrt3<-3. No integer is available.

For d3 the norm is1. The trace is positive, since t>9/5 and
a+b<8/5; AM--GM gives t-a-b<=t-2/sqrt(t)<2. It is1.
Thus a+b=t-1, ab=1/t, 9/4<t<13/5. The integral second coefficient

    S=1/t+t-t^2

is strictly decreasing and lies between-1227/325>-4 and
-341/144<-2. Its only possible integer is-3. But the resulting
polynomial X^3-X^2-3X-1=(X+1)(X^2-2X-1) is reducible, a
contradiction to minimal degree3.

For d2 the trace is strictly between1 and3, hence2. Its nonzero
integral norm t(2-t) lies in(-16/9,9/25), leaving only-1.
The positive root is t=1+sqrt2. For d1 the only integer in(1) is2,
the already forbidden actual parallel pair.

Thus gamma=(1+sqrt2)/2>1. The full rank20 Lorentzian Gram is now
defined over Q(sqrt2); a nonphysical embedding fixing this field
would contradict positive definiteness. Therefore actual K=Q(sqrt2).

## Exact time class for the remaining dotted branch

In square-two normalization the time coefficient H of c satisfies

    H^2=17t^2+17t+11/4=(283+204sqrt2)/4,
    Norm(H^2)=-3143/16=-7*449/16.

This is useful compatibility data, not an exclusion by itself. Another
actual K-valued wall must have a K-square time ratio; norm ratios are
then rational squares. A complete smaller-loss inventory and comparison
with its possible weak time classes are still required before claiming
that all other prefix indices are strong.

The proof has finite degree and coefficient cases, with explicit strict
endpoints. It uses no field census, optimizer, script result or assumed
whole-D20 theorem. WholeD19 and the global ceiling18 remain open.
