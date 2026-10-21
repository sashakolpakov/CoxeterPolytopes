# Whole A18 over the actual silver field: complete fresh half-loss inventory

21 October 2026. Complete manual supporting proof, **pending full root and
outside-author review**. No computation, enumeration job, optimizer, Lean or
TeX run, manuscript edit or status promotion. This proves a necessary model
inventory, not whole-A18 exclusion or a new dimension bound.

## Statement and exact input boundary

At an ACTUAL farthest vertex of a compact arithmetic first-type hyperbolic
18-polytope, assume the WHOLE finite diagram is A18 and the ACTUAL field is
K=Q(sqrt2). Write s=sqrt2, O_K=Z[s]. Use unit normals, their duals omega_i,
and the SAME spatial vector x with ||x||<1 and x.n_i>0. Put

    D_i=2i(19-i)/19, C_ij=2min(i,j)(19-max(i,j))/19,
    h_i=x.omega_i, h0=h19=0.

Call i in5..14 weak if h_i<=sqrt(D_i)-1/2. Every such actual opposite is
one of precisely the following44 NECESSARY models:

| Role | Parameters | Number |
|---|---|---:|
|One far white3 helper|i5..14, helper e1 or18, own t=2gamma=1+s|20|
|Adjacent white3, own4|(i,j)=(5,6),(6,7),(13,12),(14,13), t=s|4|
|Isolated dotted|i5..14, beta=4gamma²=3+2s or6+4s|20|

The list is the SAME44 as the silver portion of the old .52 necessary
inventory, but the proof below rechecks EVERY role at the new threshold.
It does not assert that every listed model is weak, realizable, or mutually
compatible; further actual-point and common-time tests can remove models.

| Source under `results/` | SHA-256 | Exact use |
|---|---|---|
| `arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md` | `b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707` | Loss-independent N19 finite-white types, actual own0/parallel obstructions, raw normalization and explicit conjugate Schur quadratics |
| `arithmetic_whole_a18_actual_silver_field_reduction_20261021.md` | `fb185ba3e154863520a62729348339ebde51443f5f94c77ff0abaa528c96485a` | Establishes the actual-field hypothesis at this exact whole-A18/rank18 scope; no D19 transfer |
| `arithmetic_whole_a18_quadratic_weak_cap_inventory_and_remaining_interface_20261021.md` | `4710025219a3d8a2d5a485de27d3a9f6962fa61289ef9e6120ee880ee4eec618` | Comparison with the old silver44 only; its .52 removals are NOT imported at .50 |

External premises are finite spherical classification, compact acute
simplicity, native Coxeter angles, first-type cyclic integrality and positive
conjugate forms. No preselected field list or literal-lattice reflectivity
theorem is used.

## 1. Full geometric white types and the new edge cutoff

Deleting i retains A_(i-1)+A_(18-i), both of rank>=4. Its actual opposite
with these roots forms a finite rank18 endpoint. A finite diagram is a tree,
so it has at most one contact to each white path. The complete types are:
isolated; one white3 at an end, D-position or E-position; one white4 at an
end; or two white3 contacts at the two path ends. The single-path augmented
rank>=5 excludes H/F exceptions. A two-path join has connected rank18>8;
its two long arms preclude D and an internal4 precludes B, leaving an A
path with both contacts3. There is NO foreign component at this whole vertex.

Own0 remains impossible independently of any loss. One far-end3, next-to-
end3 or endpoint4 would give positive A19,D19 orB19 in signature(18,1).
Other internal3 positions contain an affine E principal; internal4 has the
affine squared-cosine star(1/2,1/4,1/4). Two white3 contacts give an affine
cycle. With no white contact the new spacelike normal would lie in the
timelike orthogonal complement of all18 old roots. An own parallel contact
also remains impossible: the nonzero future null sum of the actual unit
pair is feasible against every actual wall, contradicting compactness.

At half loss a sufficient strong-edge target, for q=sqrt(D_i), is

    T(q)=2-(q-1/2)/[q(2q-1/2)].

It increases for sqrt(140/19)<=q<31/10. Consequently a weak index has

    R=cosh(L_i)/q <=T(q)<3274/1767<1853/1000.        (1)

The last cleared difference is251. This is only a necessary weak-edge test,
not an equivalence asserting that a retained row really has weak height.

## 2. Every nonisolated white role at the changed threshold

For one white3 on a side with m-1 nodes, m=i or19-i>=5, number its position
j from that side's far end. Its exact endpoint Schur quantities are

    delta=1-j(m-j)/(2m), sigma=gamma+j/(2m), R=sigma/sqrt(delta).

Normalize the old roots to square2. A native3 helper permits a K-valued new
normal also of square2, helper product-1 and own product-t, t=2gamma. The
old3 path gives cyclic integrality t in O_K. The full old18-plus-new Gram
has entries in Q(t) and is Lorentzian; an embedding of K fixing a proper
Q(t) would repeat it in a positive definite conjugate space. Hence the
ACTUAL field equals Q(t), with no square-root extension silently adjoined.

For finite own contacts t=2cos(pi/l), the degree formula for the real
cyclotomic field gives l3,4,5,6 as the COMPLETE possibilities of degree<=2.
Label5 generates Q(sqrt5), label6 Q(sqrt3), so neither is this actual field.
Own3 with a white3 gives an all3 affine cycle. The only remaining finite
own label is4, t=s. If its helper is nonadjacent to i, the principal star
at i with its two old neighbors has squared bonds1/2,1/4,1/4 and is affine.
Thus only ADJACENT own4/white3 is retained; its exact conjugate filter is in
Section3. This disposes of EVERY finite near/D/E label, including7/8/9/12,
by the actual-field normalization rather than a borrowed height cutoff.

For dotted own, t>2. The following fresh bounds cover all white roles.

**Far endpoint3.** Here j1, delta=(m+1)/(2m), and

    R²-(4gamma-2)=(gamma-1-1/(2m))²/delta>=0.

Equation(1) therefore gives t<68/25. The loss-independent FULL-Gram Schur
test is D_i*(t')²+2C_ie*t'+D_e<4. At N19 its exact endpoint boundary tests
give -.772<t'<.48, i.e. -193/250<t'<12/25. These are the explicit inverse
quadratics of b57 Section4, not an embedding of an assumed unit normal.
Since t>=phi>1.6 and t'<.48, write t=a+b*s with b a POSITIVE integer.
Now 2b*s=t-t'<68/25+193/250=873/250<4s, so b1. The physical interval
1.6<t<2.72 then forces a1. Therefore t=1+s, and no other dotted or finite
far parameter is missed. The necessary twenty index/helper models are
retained, with no height-dependent deletion.

**Near endpoint3.** For a dotted contact, delta<=3/5 and
sigma>=gamma+2/5. At gamma>=21/20,

    R²>=(29/20)²/(3/5)=841/240>(1853/1000)².

Thus a weak dotted parameter would have2<t<21/10. The independent conjugate
quadratic has t' in(-3/2,0): Q(-3/2)=(148+15i-i²)/38>=81/19>4,
Q(0)>=120/19>4, and Q(-1)=36/19<4 after reflection to j=i-1.
Then 0<2b*s=t-t'<3.6 forces b1, whereas
2-s<a<2.1-s contains no integer. Hence NO near dotted model remains.

**D-position3.** In either orientation delta=2/m. A dotted contact has
R>=(m+1)/sqrt(2m)>=6/sqrt10>1.853. Finite own3/4/5/6 and every higher
finite label were already covered above; therefore no D-position is weak.

**E-position3.** The complete finite white possibilities are(m,j)=(6,3),
(7,3/4),(8,3/5). For own>=4 their lower ratios are respectively
2gamma+1/2, sqrt7*(gamma+3/14), and4*(gamma+3/16).
At gamma>=1/sqrt2 the first is>=sqrt2+1/2>1.9 and the other two are
larger. Own3 is the affine cycle and own0 was handled above. Thus ALL E
roles are strong, finite or dotted, at the new threshold.

**Endpoint4.** Far own3 contains affine F4; near own3 has R>2. For every
own>=4, finite or dotted, R>=(m+1)/sqrt(2m)>=6/sqrt10>1.853.

**Two endpoint3 contacts.** Own0/3 contains an all3 affine cycle. Every
own>=4 has R>=q/sqrt2+1/q>2>1.853, using the exact two-white Schur formula.

These exhaust the finite white TYPES and every native/dotted own contact.
No new role is discarded merely because it was strong at loss.52.

## 3. The four and only four adjacent own4 models

For j=i-1 and t=s, the square-two time is F=(t*Omega_i+Omega_j)²-2,
where Omega_i=omega_i/sqrt2. At the nonphysical embedding its numerator is

    19F'=f_i=59i-3i²-58-2(i-1)(19-i)*s.

Admissibility requires F'<0. This quadratic in i is concave. At5 and12
its values are162-112s>0 and218-154s>0; the respective squared margins
are1156 and92. Therefore EVERY i5..12 is inadmissible. At13,14 the values
are202-144s<0 and180-130s<0. Reflection supplies the other orientation.
Exactly the four rows in the theorem remain; their physical time squares
are positive, with

    19F=180+130s for(5,6),(14,13),
    19F=202+144s for(6,7),(13,12).

For completeness every far t=1+s model has, with k=i(19-i) and a=19-i
for helper1 or a=i for helper18,

    19F=3k+2a-20+(2k+2a)*s.                         (2)

Its conjugate is k(3-2s)+2a(1-s)-20<90*.172-20<0 and its physical
value is positive. Thus all20 are legitimately retained NECESSARY models;
none is asserted to extend to a compact chamber.

## 4. Isolated finite, own0/parallel, and the full expanded dotted domain

For an isolated finite contact put beta=4cos²(pi/l). The K-valued raw
normal has square2beta and marked product-beta. The full Lorentzian Gram
locks K=Q(beta); beta is integral and all nonphysical roots obey

    0<beta'<4/D_i<=19/35.                           (3)

The real cyclotomic degree condition leaves only l3,4,5,6,8,10,12.
Labels3/4/6 have rational beta (and affine geometric obstructions); labels5
and10 generate Q(sqrt5), label12 generates Q(sqrt3). The ONLY possible
silver isolated finite label is8, beta=2+s. Its beta'=2-s violates the
per-index condition beta'<38/k because the raw diagonal k/19>=70/19 is
greater than2+s=2/(2-s). Thus NO isolated finite cap is half-loss weak
or even admissible in this actual-field/rank window.

For an isolated dotted contact beta=4gamma²>4 is integral, has actual field
Q(beta)=K, and obeys(3). Set u=4/beta. Its actual single-foot estimate bounds
height loss by .191u+sqrt(u(1-u)). For u<=1/4 this is at most
.04775+sqrt3/4<.5. Therefore any HALF-loss weak isolated dot satisfies

    4<beta<16, 0<beta'<19/35, beta=a+b*s in Z[s].

This is the EXPANDED physical bound at.50, not the old40/3 bound. Since
beta=2a-beta', the integer a is3..8 and b is positive. Its complete
integer interval is (a-19/35)/s < b < a/s. Direct endpoint comparisons give

| a |3|4|5|6|7|8|
|---|---:|---:|---:|---:|---:|---:|
| integers b in that open interval |2|none|none|4|none|none|

Both beta=3+2s and6+4s have physical value in(4,16). Their conjugates are
less than.2 and.4, respectively, so BOTH obey the stronger per-index
condition beta'<38/k for every k<=90 (since38/90=19/45>.4).
They therefore give exactly twenty necessary index models, with raw times

    F_i=beta*(k*beta-38)/19, c²=2beta,
    c.old_i=-beta, projection beta*Omega_i.         (4)

This covers the complete continuous isolated-dotted class and all its
quadratic integral parameters. No finite domain was enumerated by a job.

## Resulting interface

Equations(2),(4) and the four near4 times specify every raw time and root
normalization of the44 necessary models. The positive physical time must
be used in any subsequent shared-time/Cartan comparison. Neither a square
norm ratio alone nor this model list implies geometric coexistence.
Conversely every actual .50-weak cap has been covered, so the old44-model
necessary compatibility graph applies to the NEW weak set once its own
exact field/time/Cartan completeness is supplied. Any energy exclusion
must still treat isolated weak5/14 and the genuine silver far-helper pair
honestly, using one original vector. No such whole-A18 conclusion or
dimension improvement is asserted by this inventory note alone.
