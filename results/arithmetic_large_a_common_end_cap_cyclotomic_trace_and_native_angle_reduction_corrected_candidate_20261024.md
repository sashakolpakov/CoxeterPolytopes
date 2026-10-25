# Preserved normalization clarification: the native-six affine null

24 October2026. FULL corrected edition of the frozen304-line candidate
`arithmetic_large_a_common_end_cap_cyclotomic_trace_and_native_angle_reduction_candidate_20261024.md`,
SHA256 `1bc5b56b2239576f6d665db6886c689972486aae24dfc5598460357e7631f0dc`.
Composition identified that the Section6 native-six display needed its
root normalization stated: Section1 gives3c+2n1+n2, while fresh square2
normalization gives sqrt3*c+2n1+n2. Both are now explicit, with all-wall
checks. No angle/field conclusion, trace argument or other input changed.
The original ambiguous display is preserved. This full corrected source
still requires ENTIRE different-author review and root acceptance.
No execution, geometric Lean or dimension improvement is claimed.

---

# Arithmetic large A: a finite common-end-cap alphabet from a trace polynomial

24 October2026. **Complete handwritten supporting CANDIDATE, pending
different-author full review and root acceptance.** This is an ARITHMETIC
local reduction. No assertion below is imported into the general bounds.
Root derived the field equality, trace-polynomial degree bound, complete
angle exhaustion and both adjacent-end continuations. Historical novelty
is not asserted. There is no new
dimension bound, field census, enumeration or realized example.

## Statement and precise external boundary

Let an actual vertex of a compact arithmetic hyperbolic Coxeter polytope
have a whole A_r component, r>=11. Let c be an ACTUAL different facet
whose only nonzero contact within this component is with an END root,
at finite native angle pi/m, m>=7. Every contact with other old components
is allowed. No farthest-point or small-energy assumption is needed.

**Theorem.** The only possible labels are

    m in{7,8,10,12,18}.

Moreover the ACTUAL ground field is exactly

    K=Q(2+2cos(2pi/m)).

Its degree is consequently at most three. Sections4--5 also prove:

* for adjacent end contacts(m,3), with m>=4, only m=4,5,6,9 is possible;
* for adjacent end contacts(m,4), with m>=3, only m=3 is possible.

Here the second contact is at old node2, the first at old node1, and
there are no other contacts inside A_r. All foreign components remain.
These are necessary possibilities;
no label or lattice is claimed realized. In the common-end-pair inventory
the theorem applies to the isolated native-seven-plus L branch, but its
statement does not require that c close two edges.

The external arithmetic model is the automatic admissible quadratic-form
model for arithmetic reflection lattices: original K root lines, a totally
real field K, and positive-definite conjugate forms at every nonphysical
embedding. Thus first type is not an extra subclass. The precise original
root normalization and cycle-integrality boundary is recorded in
`arithmetic_uniform_classical_weak_helper_obstruction_20261021.md`, SHA256
`535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905`;
only that model/normalization interface is used, not its helper theorem.

For the published boundary, see Vinberg's
[original1967 paper](https://doi.org/10.1070/SM1967v001n03ABEH001992).
The automatic quadratic-form reduction is restated as Lemma2.1, and
the rational-field compactness exception is stated immediately before
it, in Belolipetsky's [Section2, pages5--6](https://arxiv.org/pdf/1506.03111#page=5).
Those references supply only the external model and compactness facts;
the trace polynomial and native-angle exhaustion below are derived here.

Other external elementary facts are cyclotomic irreducibility and degree,
Euler's totient formula, the elementary bound3<pi<4, and ordinary real
trigonometric/calculus identities. The trace identities and the finite
totient exhaustion needed here are proved below. There is no arbitrary
face heredity, root-scale identification between disconnected components,
published ground-field list or use of the general dimension ceilings.

## 1. Original K normalization forces the actual field

Put N=r+1. Old3 bonds allow K-valued rescaling along the connected A_r
chain to one common norm b0. It is totally positive: every nonphysical
ambient form is positive definite. Scale the form by2/b0, so the whole
old Gram is the rational Cartan matrix A_r, with diagonal2 and bonds-1.
This global scalar is totally positive. No square-root field extension
or relation with a disconnected old component has been assumed.

For an original c0 on the actual K root line put a=(c0,n1)!=0.
Replace c0 by-c0/a, whose product with n1 is-1. This is a K rescaling
and is positive at the physical place. The finite-angle two-cycle is

    M=4cos²(pi/m)=2+2cos(2pi/m) in K,
    (c,c)=2/M, (c,n1)=-1, (c,ni)=0 for i>1.

Consequently the ENTIRE A_r+c principal Gram is defined over F=Q(M).
This statement retains every foreign contact: none enters this principal.
Its Schur complement against the old chain is

    S=2/M-r/(r+1).

Since m>=7, M>3; since r>=11, 2(r+1)/r<=24/11<3. Thus S<0:
this is an actual physically Lorentzian principal, not just a finite
angle which happens to generate a subfield. In particular the selected
vectors are independent. If K/F were proper, a nonidentity embedding
of totally real K fixing the physical F would leave this same indefinite
principal unchanged, contradicting positive-definite ambient conjugacy.
Hence ACTUAL K=F.

At each other embedding the diagonal2/M and the Schur complement are
strictly positive. Therefore EVERY other conjugate satisfies

    0<M_sigma<2(r+1)/r<=B, B=24/11.                 (1)

This is the exact source of the field restriction. No condition(1)
is available for a general nonarithmetic Coxeter polytope.

## 2. A quadratic trace polynomial bounds the cyclotomic degree

Let zeta be a primitive m-th root, M=2+zeta+zeta^-1 and
d=[F:Q]=phi(m)/2. The element M is an algebraic integer. Summing
over the primitive residues, paired by k and-k, gives

    Tr(M)=2d+mu(m),
    Tr(M²)=6d+4mu(m)+c_m(2),                       (2)

where c_m(2)=sum_(gcd(k,m)=1) zeta^(2k). For completeness the
primitive-root sum is mu(m) by inclusion-exclusion over divisors.
If m is odd, multiplication by2 permutes its primitive residues,
so c_m(2)=mu(m). If m=2s with s odd, it equals mu(s). If4 divides
m, reduction to residues modulo m/2 has exactly two lifts and gives
c_m(2)=2mu(m/2). In every case

    mu(m)>=-1, c_m(2)>=-2.

Use the quadratic test polynomial P(X)=X²-BX. By(1) it is NEGATIVE
at every nonphysical conjugate. Physically3<M<4, so
P(M)<4(4-B)=80/11. Thus its full trace satisfies

    Tr(P(M))<80/11.

On the other hand(2) gives

    Tr(P(M))=(18/11)d+(20/11)mu(m)+c_m(2)
              >=(18d-42)/11.

Therefore18d<122 and d<=6. This is a complete analytic degree bound,
not a scan of fields or an assumed finite coefficient cutoff.

## 3. The complete totient exhaustion is elementary and finite

We now have phi(m)<=12. Every odd prime p dividing m has p-1<=12,
so p is3,5,7,11 or13. For11/13 the complete possibilities are
11,22 and13,26. For7 they are7,14,21,28,42. For5 they are
5,10,15,20,30. The remaining2^a3^b possibilities are

    1,2,4,8,16; 3,6,12,24; 9,18,36.

These lists follow directly by multiplying p^(a-1)(p-1): another
odd prime or higher exponent would exceed12. They exhaust all
prime-factor combinations; this is not a program's output. Restricting
to m>=7 leaves exactly

    7,8,9,10,11,12,13,14,15,16,18,20,21,22,24,26,28,30,36,42.

For every entry outside the claimed five, choose this explicit primitive
residue k:

|m|9|11|13|14|15|16|20|21|22|24|26|28|30|36|42|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|k|2|2|2|3|2|3|3|2|3|5|3|3|7|5|5|

Each has gcd(k,m)=1 and1<k<=7m/30<m/2, so it defines a
NONPHYSICAL embedding of F. Monotonicity of cosine gives

    cos(2pi k/m)>=cos(7pi/15)=sin(pi/30)>1/11.

The last strict bound has an exact elementary seal. For0<x<1,
sin x>x-x³/6, as follows by differentiating three times from0.
Since1/10<pi/30<2/15,

    sin(pi/30)>1/10-(2/15)³/6=2017/20250>1/11.

Hence this conjugate M_sigma=2+2cos(2pi k/m)>24/11, contradicting(1).
Only7,8,10,12,18 survive. Their real cyclotomic degrees are respectively
3,2,2,2,3, and their other primitive angles all exceed pi/2, giving
0<M_sigma<2. No further removal or realization claim is made.

## The abstract trace lemma

The preceding trace argument has a useful abstract form: if q>=7,
X=2+2cos(2pi/q), every nonphysical conjugate lies in(0,24/11), and
the physical X lies between24/11 and4, then q is7,8,10,12 or18.
Sections2--3 proved precisely this implication without a geometric
extra hypothesis. We now use it for the adjacent profiles as well.

## 4. Adjacent contacts(m,3): only4,5,6,9

Suppose c has only contacts m at n1 and3 at n2 inside A_r, m>=4.
Normalize c to square2 over K using its native3 contact, so its raw
products are-t,-1, where t=2cos(pi/m). The actual triangle cycle
makes t integral; the old A_r+c Gram is defined over Q(t). Its
Schur complement is

    -[(t+2)(rt-2)]/(r+1).

It is negative physically. The same totally-positive-scale/proper-
extension argument therefore gives ACTUAL K=Q(t), and at every other
embedding

    -2<t_sigma<2/r.

For X=t+2=2+2cos(2pi/(2m)), all nonphysical conjugates are consequently
in(0,24/11). Its physical value is between3 and4. Apply the abstract
trace lemma at conductor q=2m>=8: only q8,10,12,18 remain. Thus
m=4,5,6,9. Foreign components have never entered this principal or
been assumed zero. The all-three m3 case, excluded from this statement,
has its positive actual affine-A2 null and cannot occur compactly.

## 5. Adjacent contacts(m,4): only3

Now the native4 at n2 allows K rescaling to c²=4 and raw products
(c,n1)=-t,(c,n2)=-2, where t=2sqrt2*cos(pi/m). In detail, if c0
has n2-product a, its four-bond identity is a²=(c0,c0); multiplying
by-2/a gives this normalization OVER K. The triangle cycle makes t
integral and the full principal is defined over Q(t). Its Schur
complement is-P_r(t)/(r+1), with

    P_r(t)=rt²+4(r-1)t+4r-12=r(t+2)²-4(t+3).

Physically t>0 and P_r(t)>0, so ACTUAL K=Q(t). At every other
embedding P_r(t_sigma)<0. In particular

    t_sigma<-1,                                    (3)

because P_r(-1)=r-8>0 and P_r is strictly increasing on[-1,infinity).
Put M=t²/2=2+2cos(2pi/m), F0=Q(M) subset K. If[K:F0]=2 and
F0!=Q, a nonidentity embedding of F0 extends with BOTH signs of
sqrt(2M_sigma); its positive extension contradicts(3). Thus a proper
quadratic extension is possible only when F0=Q. In that case m is
3,4 or6 by the elementary totient formula.

The following explicit cyclotomic automorphisms cover EVERY other m,
without a cutoff. Write S=Q(zeta_(2m)). If4 does not divide m, choose
a primitive residue k modulo lcm(8,2m) with

    k=1 modulo2m, k=5 modulo8.

The Chinese remainder conditions agree modulo gcd(8,2m), which divides4.
This automorphism fixes S but sends sqrt2 to-sqrt2. Hence sqrt2 is
not in S. Since2cos(pi/m) is a nonzero element of S, t is not in S,
and in particular t is not in F0. The quadratic-extension observation
then leaves only m3 or6 in this case.

If8 divides m, the primitive residue k=m-1 modulo2m fixes M and
sqrt2 but negates cos(pi/m), so t is not in F0. Here F0!=Q, already
a contradiction. If m=4s with s odd and s>=3, choose instead the
primitive residue k=2s-1 modulo8s. It sends cos(pi/(4s)) to
sin(pi/(4s)), while sqrt2 is sent to one of its two signs. Therefore

    |t_sigma|=2sqrt2*sin(pi/(4s))
              <=2sqrt2*sin(pi/12)=sqrt3-1<1.

This is a nonphysical embedding of K: its value differs from the
physical t>=sqrt2. It contradicts(3). The only unhandled multiple
of four is m4, with t2 and K=Q.

Use here the explicit standard external compactness premise: a compact
arithmetic hyperbolic n-polytope with n>=4 cannot have ground field Q.
This is Meyer's rational isotropy theorem in rank at least five plus
the isotropy/nonuniformity criterion for the quadratic-form lattice.
It was also the user's stated arithmetic boundary; it is not proved
or formalized here. Since n>=r>=11, it excludes m4.

For m6, K=Q(sqrt6) and the other t=-sqrt6 gives

    P_r(-sqrt6)=10r-12-4(r-1)sqrt6
                 >=98-40sqrt6>0,

since r>=11, 10-4sqrt6>0 and98²-1600*6=4>0. This contradicts
nonphysical positivity. Only m3 remains. No realization or validity
at arbitrarily large r is asserted for that remaining possibility.

## 6. The resulting literal arithmetic cap alphabet

An isolated end-five cap has M=(3+sqrt5)/2>24/11, so the Section1
principal argument gives ACTUAL K=Q(sqrt5), with arbitrary foreign
contacts. For isolated end-six, use the Section1 normalization c²=2/3,
(c,n1)=-1. Its actual positive affine3/6 null is3c+2n1+n2.
Equivalently after rescaling c to square2 it is sqrt3*c+2n1+n2.
All selected products vanish; its positive c coefficient makes it future,
and every unselected actual wall has nonpositive product. Compactness
therefore excludes this nonzero feasible null. All isolated native7+
labels are exactly covered by Sections1--3. The H foot profile4@1,
3@3 normalizes to square2 along its simple3 helper; its actual cycle
puts sqrt2 in K, and its Lorentzian A_r+c principal gives K=Q(sqrt2).
The C profile4/3 at3/4 has the same field from its actual334 triangle.
These assertions use no geometry of an auxiliary arithmetic face.

Thus AFTER a general actual-wall inventory has reduced a common end
or3/4 pair to its literal L/H/F/C profiles, their arithmetic continuation
has only the following ACTUAL ground fields:

    Q(2cos(2pi/7)), Q(sqrt2), Q(sqrt5), Q(sqrt3),
    Q(2cos(pi/9)).

For L, isolated m is5,7,8,10,12,18; adjacent(m,3) has m4,5,6,9;
adjacent(m,4) has m3. H/C use Q(sqrt2), and F uses Q(sqrt5).
The fields have degree two or three. This is a local consequence of
proved geometric support, not a global five-field theorem for all
compact arithmetic polytopes or all A components. It supplies a
finite arithmetic alphabet before any whole-diagram search and retains
every foreign contact until its separate geometric/arithmetic proof.

## Final scope and verification status

The proof applies with every foreign component and every actual degree
retained, and uses no numerical energy lemma. It can prune a literal
arithmetic cap catalogue AFTER its geometric support is proved.
It does not classify all common caps, remove a whole arithmetic15 group,
prove arithmetic14 or change any general bound. The general bounds
continue to assume NO arithmeticity. No scientific import, worker,
enumeration, optimization, numerical/symbolic process, test, Lean or TeX
ran. Frozen sources, manuscripts and successful verification records
remain unchanged. Full different-author review is still required.
