# Arithmetic D16 plus actual I: quartic forks have zero old-I contact and norm-two invariants

22 October 2026. **Complete manual candidate, UNACCEPTED pending full root and outside-author review.** No worker, mathematical import, enumeration, numerical/symbolic process, test, Lean or TeX ran. Frozen papers, accepted sources, public guides and global ceilings are unchanged.

**Theorem at the exact source boundaries below.** Suppose an actual farthest vertex of a compact arithmetic Coxeter17-polytope has ENTIRE finite forest D16 plus actual isolated short I, and its ACTUAL defining field K has degree4. Both fork opposites are isolated within D16 and have dotted own contact. In addition:

1. BOTH are orthogonal to the **old I wall**;
2. for their independently varying own magnitudes gamma1,gamma2, the integral generators `T_i=4gamma_i^2` satisfy

       K=Q(T_i), 4<T_i<1369/32,
       0<(T_i)_sigma<1/2 at every other embedding,
       N_K(T_1)=N_K(T_2)=2.

This does not assert orthogonality to the I opposite, equality of fork parameters, whole-forest nonexistence, an arithmetic16 ceiling or an example. The quartic two-isolated-dotted-fork configuration with these necessary conditions remains open.

## 1. Exact inherited boundary and conventions

All paths below are under `results/`.

| Source | SHA-256 | Restricted use |
|---|---|---|
| `arithmetic17_d16_actual_i_all_forks_dotted_and_nonisolated_degree_four_candidate_20261022.md` | `6cb2a4662222dd7984f8e29f8ed532eee22456d2726034e4f115e49b65252387` | Every actual fork own is dotted; complete retained-D roles |
| `arithmetic17_d16_actual_i_quartic_nonisolated_forks_exclusion_candidate_20261022.md` | `651bb43b1062034eb7b132143f54e6e292af2ff0b1e7daa9b1f0aeb82210a698` | A nonisolated fork forces actual degree at most3; hence quartic forces both D-isolated |
| `arithmetic17_d16_actual_i_joint_head_height_and_shorter_forks_candidate_20261022.md` | `3a5c21e33622c62bb77a81d8cdbe394e1179001cc7a19a984cbbf6ebd6030da4` | BOTH actual fork cosh lengths are less than37/4 |
| `arithmetic17_d16_actual_i_actual_field_degree_four_candidate_20261022.md` | `f0b9b9d8746719ac5a7c1c6af480bfe5fb6581be02376381329a3bbeb1af4ba4` | Isolated-fork raw square2/T, actual-field descent and strict nonphysical bound; no field census |
| `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | Actual K root lines, invariant-cycle integrality, nonphysical definiteness and principal-Gram descent |
| `review_arithmetic17_d16_actual_i_dotted_forks_height_and_field_four_root_acceptance_20261022.md` | `50e65569883340da5a5269679d1a53e645a576d34a8bceca433f437ccd830ddd` | Recorded acyclic accepted scope of6cb2,3a5c andf0b9 |

The new651b input has a separate full outside-author report, SHA `71f46be4a8dd6ecfcff65f97c5752391a73746d02f7bbe8cf7882058d76d29a2`; any remaining acceptance obligation for it is retained explicitly. This note does not self-certify an upstream candidate. The contact argument below itself applies to any individually D-isolated dotted fork over a quartic ACTUAL field;651b is used only to apply it to BOTH forks in the stated entire forest.

I authored earlier retained-white and actual-foot supporting work, and reviewed several inherited D16 arguments. Root proposed pursuing the quartic isolated-fork lane. The contact/norm argument below is my new manual derivation and requires another author's review. Standard external foundations remain the admissible quadratic-form model automatic for arithmetic reflection lattices, actual K root lines, integral Cartan cyclic products and strict positive nonphysical forms. No additional first/simplest-type restriction, relative D/I length squareclass, script result, realized field list or parameter equality is assumed.

Use unit D16 normals and fork duals omega15,omega16 with

    omega_f^2=8, omega15.omega16=7.

The old I unit wall e is orthogonal to D16. At an individually isolated fork f its ACTUAL opposite has

    c=A p-gamma omega_f-beta e,
    gamma>1, beta=0 or cos(pi/m), integer m>=3,
    A^2=8gamma^2+beta^2-1.

The old-I contact is finite because e is retained at the compact adjacent endpoint. It is not assumed zero. Contact with the I opposite is unrestricted by this description and will not be discarded.

## 2. A full-principal conjugate inequality, with no disconnected scale assumption

Put `T=4gamma^2`. The own two-cycle makes T integral. The connected old D roots can be normalized OVER K to square2, and the actual c root line to own product−1. Its raw square is2/T and its only old-D contact is that own−1. The raw D dual has square4, so the restricted principal Schur complement is2/T−4. Physical negativity gives the inherited actual-field descent `K=Q(T)`, and strict nonphysical positivity gives `0<T_sigma<1/2` for EVERY other embedding.

Retain the ACTUAL old I root as well. Write its arbitrary totally positive raw square as ell in K. Its product with c has square `2beta^2*ell/T`. Thus its additional Schur cost is exactly `2beta^2/T`, independently of ell. The rank-two invariant makes `4beta^2` integral in K, so this calculation never declares disconnected physical unit normals to be K vectors.

The old D, old I and cap are independent: physically the cap has positive time A, and the old D+I is a Euclidean17-frame. Independence is preserved at every embedding. Strict nonphysical definiteness of their full principal consequently gives

    2(1-beta_sigma^2)/T_sigma-4>0,
    0<T_sigma<(1-beta_sigma^2)/2.                 (1)

For a finite native beta every conjugate of beta² lies strictly in(0,1). At beta0 the same inequality has its evident meaning. On the physical side the ACTUAL edge formula, already derived in f0b9, is

    C^2=8gamma^2/(1-beta^2)=2T/(1-beta^2).

Using the independently proved C<37/4 gives

    4<T<(1369/32)*(1-beta^2).                     (2)

Multiplying(1),(2) over the four embeddings yields the complete norm restriction

    1<=N_K(T)<(1369/256)*N_K(1-beta^2).           (3)

Every factor is positive. All old-I root scales have cancelled before taking the norm.

## 3. All finite old-I labels at least4 are impossible

Suppose beta=cos(pi/m)>0. The actual invariant beta² lies in K. For m>=3,

    E=Q(beta^2)=Q(zeta_m+zeta_m^(-1)),
    [E:Q]=phi(m)/2.

Therefore phi(m)/2 divides4, so phi(m) is2,4 or8. This gives exactly

    m in{3,4,5,6,8,10,12,15,16,20,24,30}.         (4)

Here is elementary completeness of this small integer list. Any prime p dividing m has p−1 dividing phi(m), hence p is2,3 or5. A power of3 or5 beyond the first would give an odd divisor of phi(m), impossible. Writing m=2^a3^b5^c, with b,c in{0,1}, the constraint phi(m)<=8 gives respectively m=2,4,8,16 when b=c0; m=3,6,12,24 when b1,c0; m=5,10,20 when b0,c1; and m=15,30 when b=c1. Remove m2, whose beta is zero. This proves(4) without a field census.

For completeness the norms can be evaluated without approximating any angle. In E,

    4sin^2(pi/m)=(1-zeta_m)(1-zeta_m^(-1)),
    N_E(4sin^2(pi/m))=Phi_m(1).

The second identity follows by taking the full cyclotomic norm on both sides and using positivity of all real factors. It also follows directly by pairing inverse roots in Phi_m. Passing to K gives

    N_K(1-beta^2)=Phi_m(1)^(4/[E:Q])/256.

The complete values needed in(3) are:

| m |3|4|5|6|8|10|12|15|16|20|24|30|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| phi(m)/2 |1|1|2|1|2|2|2|4|4|4|4|4|
| Phi_m(1) |3|2|5|1|2|1|1|1|2|1|1|1|
| 256 N_K(1-beta²) |81|16|25|1|4|1|1|1|2|1|1|1|

These cyclotomic evaluations are the elementary prime-power rule Phi_(p^a)(1)=p and Phi_m(1)=1 otherwise; they may equivalently be checked from these twelve explicit elementary cyclotomic polynomials. For EVERY m>=4, the last row is at most25. Equation(3) would then give

    1<=N_K(T)<34225/65536<1,

impossible. Thus the entire finite native old-I list leaves ONLY beta1/2, besides beta0.

## 4. Old-I contact3 is impossible too: an integral reciprocal shift

Assume beta1/2. Equations(1),(2) sharpen to

    4<T<4107/128<65/2,
    0<T_sigma<3/8 otherwise.

The positive integral norm satisfies

    0<N_K(T)<(4107/128)*(3/8)^3
                =110889/65536<2.

Therefore N_K(T)=1 and T is a unit. Let the three reciprocal nonphysical conjugates be Y1,Y2,Y3. Then

    Y_i>8/3, Y1*Y2*Y3=T,
    T>512/27, 1/T<27/512.

For three positive numbers greater than a, their sum is strictly less than `2a+product/a^2`: replacing two numbers y,z by a,yz/a increases their sum, since `(y-a)(z-a)/a>0`, and repeating once proves the inequality. At a=8/3 this gives

    Y1+Y2+Y3<16/3+9T/64.

The trace of the integral unit T^(-1) is therefore an integer strictly greater than8 and less than

    16/3+585/128+27/512=15293/1536<10.

It equals9. As the physical reciprocal is positive, each Y_i is consequently less than `9-16/3=11/3`. Thus the nonphysical conjugates of the integral element

    alpha=T^(-1)-3

lie in(-1/3,2/3), whereas its physical absolute value is less than3. It is nonzero, because T>4. Its nonzero integral norm would satisfy

    0<|N_K(alpha)|<3*(2/3)^3=8/9<1,

impossible. This excludes beta1/2 in the full quartic isolated-fork branch. Combining Sections3–4 proves old-I orthogonality for BOTH fork opposites under the stated quartic hypotheses.

## 5. A zero-contact isolated fork has norm invariant1 or2

Now beta0 and equations(1),(2) give

    4<T<1369/32<(7/2)^3,
    0<T_sigma<1/2 otherwise.

Write N=N_K(T), a positive integer, and x1,x2,x3 for the three small conjugates. Their product is N/T. The integral element1−2T is nonzero. Its absolute norm is

    (2T-1)*product(1-2x_i).

AMGM first on the x_i and then on their positive deficits gives

    product(1-2x_i)
       <=[1-2*(N/T)^(1/3)]^3.

The quantity in brackets is positive because each x_i<1/2. Hence

    0<|N_K(1-2T)|
       <2*[T^(1/3)-2N^(1/3)]^3.                  (5)

If N>=3, then `N^(1/3)>7/5`, since `(7/5)^3<3`. The right-hand side of(5) is strictly less than

    2*(7/2-14/5)^3=343/500<1,

contradicting the nonzero integral norm. Therefore EVERY remaining isolated quartic fork has N_K(T) in{1,2}. No coefficient enumeration or prescribed quartic field is used.

## 6. The two actual forks force BOTH norms to be2

The two D-isolated fork opposites c1,c2 are distinct: one has nonzero own product at old fork15 and zero at16, the other the reverse. Choose their ACTUAL K root lines normalized to raw own product−1, denoted r1,r2. Their squares are2/T1 and2/T2. Put `P=r1.r2 in K`.

The old raw D roots have square2 and fork inverse entries4 on the diagonal and7/2 off the diagonal. Since both old-I contacts are now zero, subtracting their D projections leaves two vectors on the SAME physical time line. The residual2x2 determinant is zero. Equivalently,

    (P-7/2)^2
       =(2/T1-4)*(2/T2-4)
       =4*(2T1-1)*(2T2-1)/(T1*T2).              (6)

This is an identity among K elements, obtained from the actual walls. It is not an assertion about an arbitrary formal rank-two Gram or an assumed equality of their time parameters.

The actual five-edge Cartan cycle

    r1 -> n15 -> n14 -> n16 -> r2 -> r1

has factors `-T1,-1,-1,-1,T2*P`. Its invariant cyclic product therefore gives

    X:=T1*T2*P in O_K.                            (7)

If P0, the same assertion is immediate. No assertion that P itself is integral is needed. Clearing(6) by(7) gives the exact identity

    (2X-7T1*T2)^2
       =16*T1*T2*(2T1-1)*(2T2-1).                (8)

Reducing(8) modulo2O_K shows `(T1*T2)^2 in2O_K`. Taking the quartic norm gives

    16 divides [N_K(T1)*N_K(T2)]^2.

By Section5 each positive norm is1 or2. Their product is therefore1,2 or4, and only4 has its square divisible by16. Consequently

    N_K(T1)=N_K(T2)=2.

This is a joint restriction on the TWO actual fork walls. It does not assume their parameters coincide, compare disconnected physical unit scales over K, or use a selected local model as an ambient example.

## 7. Exact conclusion and remaining obligation

Under the stated quartic actual-field hypothesis, both forks are D-isolated with dotted own contacts, their old-I contacts vanish, and their independent invariants T1,T2 both have norm2 with every nonphysical conjugate strictly between0 and1/2. The proof has excluded all positive old-I native labels and both norm-one alternatives, not the remaining norm-two quartet geometry.

No contact with the I opposite was assumed zero. No new farthest energy has been added and no incompatible head/I/fork payments have been summed. In particular these arithmetic constraints do not yet exclude a remaining pair of norm-two invariants or the full D16+actual-I forest. Closing that residual requires another complete actual geometric/arithmetic argument, not absence of a script output or a one-wall relaxation. The source is a supporting candidate pending full independent review and root acceptance, not an additional paper, a new bound or a new example.
