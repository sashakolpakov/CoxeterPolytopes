# Arithmetic D16 plus actual I: the actual field has degree at most four

22 October2026. Root-authored complete manual candidate, **pending full
outside-author mathematical review**. No worker, mathematical import,
enumeration, numerical/symbolic process, test, Lean or TeX ran. The six
manuscripts and all global dimension ceilings remain unchanged.

**Theorem.** If a compact arithmetic Coxeter17-polytope has an actual
farthest vertex with ENTIRE finite forest D16 plus actual isolated short
I, its ACTUAL admissible defining field has degree at most4.

This is a field-degree reduction for the entire stated forest, including
two isolated dotted forks. It is not whole-forest nonexistence, an
arithmetic16 ceiling, a prescribed field list or a new example. The
admissible quadratic-form model is automatic for arithmetic reflection
lattices and introduces no extra first/simplest-type restriction.

## Exact premises and boundaries

All source paths in the table are under `results/`.

| Source | SHA-256 | Restricted interface |
|---|---|---|
| `arithmetic17_d16_actual_i_fork_cosh_and_three_white_role_reduction_candidate_20261022.md` | `cdd5b294c7ec4cb8f32c847bc50b0eb3a8fb5c3104acddb811409d9e2a010a26` | Both ACTUAL fork lengths C<11; exact complete retained-D white roles |
| `arithmetic17_d16_actual_i_all_forks_dotted_and_nonisolated_degree_four_candidate_20261022.md` | `6cb2a4662222dd7984f8e29f8ed532eee22456d2726034e4f115e49b65252387` | Both fork own contacts dotted; a nonisolated fork forces ACTUAL field degree at most4 |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | Actual K root lines, integral invariant cycles, nonphysical definiteness and principal-Gram descent |

The second source has a full outside-author review9c7c158f; root
acceptance is a separate obligation. The new isolated-fork argument
below does NOT use that source's nonisolated degree conclusion to
bound an isolated fork. No finite script output or relative root-length
squareclass is assumed. Elementary algebraic facts used explicitly are
that an algebraic integer has a monic integer minimal polynomial, its
norm is a nonzero integer when it is nonzero, and the arithmetic and
geometric mean inequality for positive real numbers.

## 1. A small algebraic lemma

Let T be an algebraic integer generating a totally real field K. Suppose
its physical embedding and EVERY other embedding satisfy

    4<T<121/2, 0<T_sigma<1/2 otherwise.              (1)

Then [K:Q]<=4.

Write d=[K:Q]. Its norm is a positive integer. If d>=7, then

    0<N(T)<(121/2)*(1/2)^(d-1)<=121/128<1,

impossible. It remains to exclude d5 and d6. Let P be T's monic
integer minimal polynomial. Since T>4, irreducibility implies P(1/2)
is nonzero, so

    |P(1/2)|>=2^(-d).                              (2)

For d5, call the four small conjugates x1,...,x4. Their product is
N(T)/T>1/61. Since

    61*(7/20)^4=146461/160000<1,

the arithmetic mean of the xi is greater than7/20. The positive
numbers yi=1/2-xi consequently have arithmetic mean less than3/20,
and hence product less than(3/20)^4. But then

    |P(1/2)|=(T-1/2)*product yi
       <61*(3/20)^4=4941/160000
       <5000/160000=1/32,

contradicting(2).

For d6, the five small conjugates have product N(T)/T>1/61. Since

    61*(2/5)^5=1952/3125<1,

their arithmetic mean is greater than2/5. The five positive deficits
1/2-xi have mean less than1/10. It follows that

    |P(1/2)|<61*(1/10)^5=61/100000<1/64,

again contradicting(2). This exhausts every degree above4 and proves
the lemma. No coefficient search or enumeration has been used.

## 2. An isolated dotted ACTUAL fork satisfies the lemma

Choose either original fork f. If its opposite is isolated within D16,
the accepted actual role and own reduction give, in UNIT notation,

    c=A p-gamma*omega_f-beta*e,
    gamma>1, 0<=beta<1, omega_f²=8,
    A²=8gamma²+beta²-1.

The retained-I contact beta is finite or zero. It is NOT assumed zero,
and no relation between its raw scale and D is imposed. The actual
outgoing edge has unit spatial direction-omega_f/sqrt8. Orthogonality
at its ACTUAL endpoint therefore gives

    C²=8gamma²/(1-beta²).

Set T=4gamma². Its two-cycle is an integral invariant, so T in O_K.
The dotted condition and inherited C<11 imply

    4<T=(C²/2)*(1-beta²)<121/2.                     (3)

For the field and nonphysical bounds, normalize the connected old D
roots OVER K to square2. The actual c root line can be normalized to
own product-1. Its raw square is1/(2gamma²)=2/T; its only old-D
contact is that own-1. Projection onto the raw D dual of square4
gives the restricted principal Schur complement

    2/T-4.                                        (4)

This uses a K-valued actual root line, not a physical unit normal
declared to belong to K. All principal entries belong to Q(T), and
(4) is physically negative. Any nonphysical embedding fixing Q(T)
would retain this indefinite principal in a positive definite form.
Thus the ACTUAL field is K=Q(T).

At every other embedding the raw diagonal2/T_sigma and the Schur
complement2/T_sigma-4 are strictly positive. Hence

    0<T_sigma<1/2                                 (5)

for EVERY other embedding. Equations(3),(5) give exactly the algebraic
lemma, so an individually isolated dotted fork forces [K:Q]<=4.
Both fork orientations are included because their dual squares are
the same; their own parameters and foreign-I contacts need not agree.

## 3. Composition covers the entire forest

Choose either actual fork. Its own contact is dotted by the explicit
source6cb2 interface. If it is isolated within D16, Section2 gives
actual degree at most4. If it is nonisolated, the independently proved
endpoint3-role degree bound of6cb2 gives actual degree at most4.
These are the complete retained-D alternatives from cdd5, including
all allowed foreign-I contacts. Thus no fork configuration in the
stated entire forest can have actual field degree above4.

The proof neither assumes the two forks share a wall nor requires
them to be distinct. It uses a single chosen actual edge and never
sums independent component energies. Field degrees2,3,4 and all
remaining geometric cases still require further work. There is no
global dimension improvement, example, execution or Lean claim.
Full outside-author review and root acceptance are required before
promoting this new whole-forest field-degree conclusion.
