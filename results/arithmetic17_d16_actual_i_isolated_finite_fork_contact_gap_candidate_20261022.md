# D16 plus actual I: a short contact-gap proof for finite fork opposites

22 October 2026. Root-authored manual candidate, **pending complete
outside-author review**. No mathematical worker, numerical/symbolic process,
enumeration, scientific import, test, Lean or TeX ran. The earlier complete
two-finite-fork proof is preserved. The present argument proves a stronger
individual contact restriction and gives a shorter route to its joint
conclusion. It uses no finite-beta field census or squareclass enumeration.

## Statement and accepted inputs

At an actual farthest vertex of a compact arithmetic Coxeter17-polytope, suppose
the ENTIRE finite forest is D16 plus actual isolated short I. If an outgoing
fork opposite is isolated within D16 and has finite positive own contact,
then exactly the following necessary profiles remain:

| Own contact | Actual field | Contact with old I |
|---|---|---|
| cos(pi/5) | Q(phi) | zero or cos(pi/5) |
| cos(pi/12) | Q(sqrt3) | zero |

Here phi=(1+sqrt5)/2. The table is a necessary restriction, not an existence
claim. In particular two isolated finite-own fork opposites cannot occur
together, even with unequal contacts with I. If both fork opposites are
isolated within D16, at least one must have dotted own contact.

The argument uses these four source-bound interfaces under `results/`:

| Source | SHA-256 | Interface |
|---|---|---|
| `arithmetic17_d16_actual_i_strong_five_finite_head_and_height_reduction_candidate_20261022.md` | `71b49bf5272c31374702a91c271e6f203f2abc7a773b7f609b53f3a39a5c0412` | Actual D16/I coordinates and unconditional0<h<23/100, invariant J strip |
| `review_arithmetic17_d16_actual_i_height_and_head_reduction_root_acceptance_20261022.md` | `8e32ced5182f452b43302cf457d2d004bc69772e8800a4e6bb3b1a40846a5bec` | Accepted low-height scope and acyclic proof order |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | K-valued root-line normalization, integral invariant cycles, all nonphysical embeddings |
| `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md` | `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84` | Actual-wall acuteness and feasibility of positive actual-normal combinations |

Compact acute simplicity, the native wall-product alphabet, the admissible
arithmetic quadratic-form model and the integer rings of the two quadratic
fields are explicit foundations. That model is automatic for arithmetic
reflection lattices. No prescribed lattice, D/I relative root-length
squareclass, new centre floor or unreviewed fork theorem is assumed.

## The actual field and actual I height

Use unit normals, with old fork duals omega_f,omega_h satisfying

    omega_f²=omega_h²=8, omega_f.omega_h=7.

Let e be the old I unit root, orthogonal to D16. Its actual opposite is

    c_I=sqrt(h)p-sqrt(1+h)e,
    0<h<23/100,
    J=4(1+h) in O_K,
    4<J<123/25, 0<J_sigma<4 at every other embedding.       (1)

An isolated finite-own fork opposite has

    c=A p-gamma*omega_f-beta*e,
    s=gamma², A²=8s+beta²-1,
    gamma=cos(pi/m), m>=3,
    beta=0 or cos(pi/n), n>=3, 0<=beta<1.                  (2)

The beta pair is finite because c and old e belong to the actual finite
endpoint. In particular a nonzero beta is at least1/2. The wall c is distinct
from c_I: its nonzero old-fork product is -gamma, whereas c_I is orthogonal
to every old D root. We may therefore use their actual acuteness and native
product alphabet throughout.

Here is the finite-own field reduction in full. Normalize the connected old
D root lines over K to common square1 and rescale c's K root line to own
product-1. Its square is1/s. The old-D-plus-cap principal has rational old
entries and Schur complement1/s-8<0, independently of beta. Its Gram lies
over Q(s). If actual K properly extended Q(s), a nonphysical embedding fixing
s would preserve this indefinite principal inside a positive definite form.
Thus K=Q(s). Every other conjugate satisfies0<s_sigma<1/8.

The two-cycle makes t=4s integral. For own3,4,6, t is rational, so K=Q;
the integral J interval(4,123/25) is then empty. Own5 gives
s=(phi+1)/4 and K=Q(phi). For every finite own at least7, physically3<t<4
and every other0<t_sigma<1/2. Degree at least3 would give a positive nonzero
integral norm less than4*(1/2)^(d-1)<=1, impossible. Degree1 has no integer
in(3,4). Degree2 has norm1 and integer trace4, hence t=2+sqrt3 and own12.
This exhausts all positive finite own parameters.

The J strip now determines h without further tables. In O_Q(phi)=Z[phi],
write J=a+b phi. The embedding gap b sqrt5 lies in(0,123/25), so b is1
or2. For b1 the physical interval forces a3; for b2 it forces a1, whose
conjugate is negative. Thus J=3+phi and h=(phi-1)/4. In O_Q(sqrt3)=Z[sqrt3],
the gap2b sqrt3 forces b1 and then a3. Thus J=3+sqrt3 and
h=(sqrt3-1)/4. Both caps, when present, share the same ACTUAL field.

## A decreasing native product removes the foreign contacts

Put B²=8s-1>0. The magnitude of the product with the ACTUAL c_I is

    F(beta)=-c.c_I
      =sqrt(h*(B²+beta²))-beta*sqrt(1+h).                 (3)

Acuteness requires F>=0. Its derivative is strictly negative, since

    F'(beta)=sqrt(h)*beta/sqrt(B²+beta²)-sqrt(1+h)<0.

Its sign is the sign of hB²-beta², by comparing the two positive terms
in(3). If0<F<1, actual walls must have F=cos(pi/k), k>=3. In particular
there is no native value in(0,1/2) or in(1/2,1/sqrt2).

For own5, the exact height gives hB²=s. Thus F>=0 forces beta<=sqrt(s).
For beta<sqrt(s), F>0. Moreover F(1/2)<1/2: after squaring the two positive
sides, the difference needed is

    [(1+sqrt(1+h))/2]²-[s+h/4]
       =1/2+sqrt(1+h)/2-s>1-s>0.

Every nonzero beta with beta<sqrt(s) therefore gives0<F<1/2 and is
impossible. The only remaining contacts are beta0 and beta=sqrt(s)=cos(pi/5).
They give respectively F=cos(pi/5) and F0.

For own12, the exact height gives hB²=1+h. Consequently F>0 for all beta<1.
Every nonzero finite beta equals1/2 or is at least1/sqrt2. At beta1/2,

    1/2<F(1/2)<1/sqrt2.                                 (4)

For the lower inequality, the difference of the relevant positive squares
is(1+2h-sqrt(1+h))/2>0; squaring1+2h gives the positive margin3h+4h².
For the upper inequality use h<.23:
sqrt(1+5h/4)<1.135 and sqrt(1+h)/2>1/2, hence F(1/2)<.635<.7<1/sqrt2.
The seal1.135²=1.288225>1+5*.23/4=1.2875 is rational.

At beta1/sqrt2, F<1/2. Indeed

    [1/2+sqrt((1+h)/2)]²-[1+3h/2]
      =sqrt((1+h)/2)-1/4-h>.7-.48>0.

Monotonicity now puts every beta>=1/sqrt2 in0<F<1/2. Both possibilities
violate the native alphabet by(4). Therefore only beta0 survives for own12.
This proves the individual table without classifying beta's field degree
or compatible squareclasses.

## The two finite forks contradict actual acuteness

If both isolated fork opposites are finite-own, the field reduction makes
their own parameters the same gamma. The two walls are distinct, since each
has a nonzero contact with one old fork and zero with the other. Their
product is

    c1.c2=7s+beta1*beta2
             -sqrt((B²+beta1²)*(B²+beta2²)).              (5)

Equal betas give1-s>0, contradicting acuteness. In the own12 case the table
already forces beta1=beta2=0. In the own5 case the only unequal possibility
is beta1=0,beta2=gamma, or its interchange. Using phi²=phi+1,

    B²*(B²+s)-(7s)²=(phi-6)/16<0.

Both compared square roots are positive, so(5) is again strictly positive.
The two finite-own forks cannot coexist.

For completeness an isolated own0 would have time square beta²-1<0.
Parallel own yields the nonzero future feasible null c+n_f, a positive
combination of ACTUAL walls, contradicting compactness. Thus if both fork
opposites are isolated, at least one has dotted own contact. Nonisolated
endpoint3 roles and isolated dotted roles are not excluded. The whole
D16+actual-I family remains open, and no global bound changes.

The earlier longer proof, including its complete(4,12) compatibility case,
is retained at
`results/arithmetic17_d16_actual_i_two_isolated_finite_forks_exclusion_candidate_20261022.md`.
It is not an input here. The new individual contact restriction is stronger;
the final joint exclusion has a shorter proof. There is no new execution,
Lean verification, manuscript integration or example claim.
