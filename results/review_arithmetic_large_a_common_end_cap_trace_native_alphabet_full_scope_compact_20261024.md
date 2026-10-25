# Arithmetic large A common-end caps: full different-author source review

24 October2026. **FULL mathematical SOURCE PASS for the corrected source
identified below, with its explicit external arithmetic premises.** The
predecessor's native-six normalization display is preserved and is not
credited as an unchanged correct identity. No execution, general-case
field restriction, global five-field theorem or dimension improvement
is established by this review.

## 1. Exact target, complete read and independence boundary

I read the ENTIRE original304-line proof, then the ENTIRE corrected323-line
proof, and compared their complete textual difference. The final target is:

| Item | Exact source | SHA256 | Bytes/lines |
|---|---|---|---|
| Corrected target | `arithmetic_large_a_common_end_cap_cyclotomic_trace_and_native_angle_reduction_corrected_candidate_20261024.md` | `178dacabfdc05da8e2ee4a93a4695ee33215e66f21241fcd24f8fb327c4f70b9` |14728/323|
| Preserved predecessor | `arithmetic_large_a_common_end_cap_cyclotomic_trace_and_native_angle_reduction_candidate_20261024.md` | `1bc5b56b2239576f6d665db6886c689972486aae24dfc5598460357e7631f0dc` |13633/304|
| Sole inherited mathematical pin | `arithmetic_uniform_classical_weak_helper_obstruction_20261021.md` | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` |Section1 model/normalization only|

These identities were freshly authenticated. I read the pinned Section1
in full, together with the arithmetic quadratic-form boundary in
`docs/manuscript_assumptions.md`. The helper obstruction theorem of that
older note is not a premise of this target and was not substituted for
any new trace, degree, conductor or angle argument.

Root authored the new proof. I did not supply its field-equality argument,
trace polynomial, totient table, rejection residues or cyclotomic CRT
arguments before its freeze. My prior work on actual geometric component
and farthest-vector interfaces elsewhere is inherited background, not an
independent review of my own older conclusions. During this audit I found
the native-six normalization mismatch; Root wrote the corrected full copy.
That correction and its mathematical meaning are independently checked
below. There is no new coauthorship claim for the other arguments.

This is a source/math review, not an independent runtime implementation or
literature reproving exercise. The automatic arithmetic quadratic-form
model, its original K root lines and nonphysical positive definiteness,
cyclic integrality, and the stated Meyer/nonuniformity compactness input
remain explicit external premises. They are not supplied by a general
nonarithmetic Coxeter argument or by a Lean axiom audit.

## 2. Original-field normalization and exact field equality

All old A_r roots lie on actual K root lines. Across an old native3 bond,
the identity `(n_i,n_j)^2=n_i^2*n_j^2/4` permits a K-valued rescaling
to equal norm; its physical sign can be chosen positive. Propagation
along the CONNECTED chain gives common norm b0. This norm is positive
physically and at every other embedding. Thus the global form scaling
by2/b0 is totally positive and gives the rational square-two Cartan
Gram. It imposes no square-root extension and no scale identification
with a disconnected spectator.

For the one-end contact, the actual K rescaling `c=-c0/(c0,n1)` gives
`c.n1=-1` and `c^2=2/M`, where `M=4cos^2(pi/m)` belongs to K by
the two-cycle. The raw A_r inverse entry at1 is r/(r+1), so the exact
Schur complement is `2/M-r/(r+1)`. For r>=11 and m>=7 it is negative:
M>3 whereas `2(r+1)/r<=24/11<3`.

Consequently this actual principal is nonsingular and Lorentzian, with
all its Gram entries in F=Q(M). If K/F were proper, another real
embedding of totally real K fixing the physical F would keep this
indefinite Gram unchanged inside a positive-definite conjugate ambient
space. This is impossible. The conclusion is ACTUAL K=F, rather than
merely F contained in K. Every nonphysical conjugate then has
`0<M_sigma<2(r+1)/r<=24/11`: positivity of its new-root diagonal gives
the lower bound, and strict Schur positivity gives the upper bound.

Foreign contacts never enter this selected principal. Their arbitrary
values therefore do not invalidate the argument, nor are they silently
assumed zero. Independence is supplied by the nonzero Schur complement,
not by a guessed ambient diagram or arithmeticity of an auxiliary face.

## 3. Complete trace and degree argument

For M=2+zeta+zeta^-1 and d=phi(m)/2, pairing primitive residues gives

    Tr(M)=2d+mu(m),
    Tr(M^2)=6d+4mu(m)+c_m(2).

I checked the real-field factor-of-two normalization. Inclusion-exclusion
gives the first primitive-root sum. The second satisfies exactly the three
cases in the source: mu(m) for odd m, mu(s) for m=2s with odd s, and
2mu(m/2) when4 divides m. Thus mu(m)>=-1 and c_m(2)>=-2.

Set B=24/11 and P(X)=X^2-BX. All nonphysical values of P are strictly
negative. Physically 3<M<4 and P is increasing there, giving
`P(M)<80/11`. The trace identity is

    Tr(P(M))=(18/11)d+(20/11)mu(m)+c_m(2)
              >=(18d-42)/11.

Combining these strict inequalities gives18d<122 and therefore d<=6.
No coefficient scan, field census or unproved finite cutoff is used.
The abstract version also applies when the physical value is only in
(24/11,4), since P is increasing on that entire interval.

## 4. Full totient exhaustion and every rejected label

I independently reconstructed phi(m)<=12 from the prime-power formula.
An odd prime can only be3,5,7,11 or13. The11/13 possibilities are
11/22 and13/26; the7 possibilities are7/14/21/28/42; the5 possibilities
are5/10/15/20/30. The remaining2^a3^b cases are exactly the three rows
listed in the source. A further odd prime, power, or two-power violates
the product bound. Restricting to m>=7 gives exactly the displayed20
labels, without an omitted conductor.

For each of the15 rejected labels I checked the specified residue is
coprime to m, differs from ±1 modulo m, and satisfies
`1<k<=7m/30<m/2`. It therefore gives a genuinely nonphysical embedding
of Q(M). The source's uniform comparison is correct:

    cos(2pi*k/m)>=cos(7pi/15)=sin(pi/30)>1/11.

The elementary sine bound has the exact rational seal
`1/10-(2/15)^3/6=2017/20250>1/11`, using3<pi<4. It follows that
the selected conjugate is greater than24/11, contradicting the full
conjugate bound. This checks ALL15 rejections, including30 at the
endpoint k=7.

The five surviving labels are exactly7,8,10,12,18. Their degrees are
3,2,2,2,3; the other primitive angles are in(pi/2,pi), so their other
M-conjugates lie strictly between0 and2. This records necessary
survivors, not existence or validity at every larger r.

## 5. Adjacent native3: normalization, conductor and all labels

The native3 contact at node2 allows square-two K normalization with
contacts -t,-1, t=2cos(pi/m). The actual triangle Cartan cycle is ±t,
so t is integral. Using the raw inverse entries at1/2 gives exactly

    S=-(t+2)(rt-2)/(r+1).

It is physically negative for m>=4. The same proper-extension argument
gives ACTUAL K=Q(t), and nonphysical positivity gives precisely
`-2<t_sigma<2/r`. For X=t+2 the correct conductor is2m, not m.
Every nonphysical X lies in(0,24/11); physically3<X<4. Applying the
already proved trace lemma at2m>=8 leaves conductors8,10,12,18, hence
m=4,5,6,9. The omitted m3 is the actual affine all-three triangle:
its positive null combination is future and feasible against every
other actual wall, contradicting compactness. Foreign components are
retained throughout this principal argument.

## 6. Adjacent native4: full unbounded conductor coverage

For a native4 at node2, its original square relation is `a^2=c0^2`
after the old norm-two normalization. Rescaling by-2/a is OVER K and
gives c^2=4 and contacts -t,-2, with t=2sqrt2*cos(pi/m). The unequal-
norm triangle cycle still equals ±t, rather than ±t/2. The exact
Schur numerator is

    P_r(t)=rt^2+4(r-1)t+4r-12=r(t+2)^2-4(t+3).

It is physically positive, so ACTUAL K=Q(t). At a nonphysical embedding
P_r(t_sigma)<0. Since P_r(-1)=r-8>0 and its derivative is positive on
[-1,infinity), every nonphysical t_sigma is strictly below-1.

For F0=Q(t^2/2), K/F0 has degree at most2. If it is quadratic and
F0 is not Q, a nonidentity F0 embedding extends with BOTH signs of
sqrt(2M_sigma); the positive choice contradicts t_sigma<-1. This
uses total reality and all conjugates, not one selected conjugate.

I checked every case in the source's no-cutoff conductor split:

* If4 does not divide m, the CRT congruences k=1 mod2m and k=5 mod8
  agree modulo their gcd. They produce a primitive cyclotomic residue
  fixing Q(zeta_(2m)) and negating sqrt2. Hence t is outside that field,
  in particular outside F0, and only rational F0, m3/6, remains.
* If8 divides m, k=m-1 is primitive modulo2m, fixes M and sqrt2, and
  negates cos(pi/m). Thus K/F0 is proper, while F0 is not Q, a
  contradiction by the preceding extension argument.
* For m=4s, odd s>=3, k=2s-1 is primitive modulo8s. It sends the
  cosine to sine and sends sqrt2 to a sign of itself. The resulting
  nonphysical value has absolute value at most sqrt3-1<1, contrary
  to the full nonphysical restriction.

These three cases exhaust all m other than the remaining m4, together
with m3/6 already retained in the first case. For m4, t=2 forces K=Q;
the stated EXTERNAL Meyer/isotropy-nonuniformity premise excludes a
compact arithmetic chamber in ambient dimension>=4. Here n>=r>=11.
It is not an internal general-Coxeter compactness theorem. For m6,
the conjugate -sqrt6 instead gives

    P_r(-sqrt6)>=98-40sqrt6>0,

where the increasing-in-r coefficient is positive and
`98^2-1600*6=4`. This contradicts conjugate positivity. Only m3 is
left as a necessary possibility; the source explicitly avoids claiming
it works for arbitrarily large r.

## 7. Corrected six-null and remaining literal profile fields

The sole original finding was substantive as a normalization identity:
Section1's end-six root has c^2=2/3 and c.n1=-1. Its positive null is
`3c+2n1+n2`, not `sqrt3*c+2n1+n2` in that normalization. The corrected
source states both forms: after c is freshly rescaled to square2 and
contact -sqrt3, the latter coefficients are correct.

I checked zero products with the THREE kernel roots c,n1,n2, null
square, and nonzero future orientation from its positive c coefficient
and p.c<0. The next old root n3 has negative product, and every other
actual wall has nonpositive product by acuteness. Thus foreign contacts
do not obstruct feasibility; they need not vanish. This is an actual
feasible null, so compactness excludes end-six without arithmeticity.

The corrected copy adds a15-line historical clarification and replaces
only the native-six paragraph with the two explicit normalizations and
all-wall/future check. The complete textual diff has no other change.
The original304-line source bytes remain at their literal path and hash;
no silent edit or pass for that old display is asserted.

For end-five, M=(3+sqrt5)/2 exceeds24/11, so the same actual principal
gives K=Q(sqrt5). For H=4@1/3@3, raw square-two normalization gives
principal Schur numerator `8-3r-2sqrt2*(r-2)`, divided by r+1; it is
negative for r>=11. Its actual cycle puts sqrt2 in K, so the same
field-equality reasoning gives K=Q(sqrt2). For C=4@3/3@4 the numerator
is `26-8r-6sqrt2*(r-3)`, again negative. Its actual334 triangle puts
sqrt2 in K. These checks retain every foreign contact.

## 8. Verdict and exact permitted consumer

The full corrected source is mathematically sound on its explicit
arithmetic model and external compactness premises. It proves the one-
end m>=7 theorem, both adjacent-end continuation statements, and their
conditional literal-cap alphabet. All degrees, finite labels, foreign
contacts and proper extensions have been covered in the stated domains.
No unresolved branch or further correction was found.

The five listed fields are a conclusion ONLY AFTER the geometric
inventory has reduced the actual cap to the specified L/H/F/C profiles.
This review supplies no such inventory for arbitrary caps, no global
five-field theorem, no realized example, no whole arithmetic15 group
deletion and no arithmetic14 or general dimension change. General
ceilings remain free of arithmeticity. No worker was imported or run;
no enumeration, scientific calculation process, test, parser, Lean or
TeX compilation was launched. Root acceptance remains separate.
