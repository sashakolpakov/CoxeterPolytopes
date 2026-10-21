# Independent review: isolated opposite side with no long-arm contact

19 October 2026. **PASS for the stated restricted local theorem.** Full
paper/source review only: no checker, matrix calculation program, search,
optimizer, Lean or build was executed, and no pinned source was changed.

Reviewed source:

* `results/arithmetic_rank20_isolated_side_no_long_contact_20261019.md`,
  SHA256 `af2ff3c9c4667b037ea93ed1afa8f1dcaf364c0ac41e894e3541b2bcde1df879`.

I also checked its named native eight-field interface in
`fresh_arithmetic_high_label_20261016.md` and compared the earlier
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`
for the normalization and finite-contact formulas. The older rank19/20
threshold is not being imported as a smaller-rank exclusion.

## Actual geometric and arithmetic scope

The claim concerns a compact arithmetic simplest-type **20-dimensional**
polytope over one of the eight displayed fields. The complete old vertex
is `D_r + T_s`, with `r+s=20`, `10<=r<=16`, and the stated D/E spectator.
The opposite wall `d` contacts only the short leaf `b`; `c` is orthogonal
to `b,d`, the entire spectator, and the long arm. Its only allowed old
contacts are the fork `a` and the other short leaf `e`.

These hypotheses are essential. There is no arbitrary extra old factor
or supporting white omitted from the calculation. This review does not
exclude long-arm or spectator contacts, and does not prove a full fork
refund, simultaneous refunds, or a new dimension ceiling.

The inherited arithmetic inputs remain explicit: cyclic algebraic
integrality, the determinant-ratio interpretation of unit-normal
expressions, admissible strict positive definiteness at nonphysical
embeddings, and the native eight-field scope. Unit-normal contacts are
not silently treated as elements of the ground field.

## 1. Schur reduction, strictness and field generation

Write `N=r-2`, `L=N-1`. Eliminating the actual positive white support
`A1(e)+A_L+T_s` gives

    s0 = r/(4N),  t = 1-gamma^2,  X = x+gamma/2,
    (X^2/t-s0)(beta^2-1)=1/4.

The coefficient `s0` follows from the A_L endpoint inverse and the
isolated e root. The off-diagonal correction is `-gamma/2`. Since b
and d are orthogonal to the support, their side normals in the actual
two-face are unchanged. Acute actual-face projection and compactness
give `beta>1`; this is not an unsupported assertion that both ambient
opposite pairs must be dotted.

There are 22 actual roots in Lorentz dimension 21. The white principal
is positive definite, so their rank defect gives exactly the displayed
four-root determinant equation. The inverse of the old D_r block gives

    qc = 2N*x^2+2N*x*gamma+(r/2)*gamma^2,
    qd = r*z/8,   z=4*beta^2,
    delta=z-4=1/(R-s0),  R=X^2/t.

Each added wall is strict at the simple old vertex. Its component in
the timelike line perpendicular to the entire old vertex is nonzero;
therefore its old projection cost is strictly greater than one. At
every nonphysical embedding the corresponding independent principal
is positive definite, giving `qc'<1`, `qd'<1` and

    0<z'<8/r,  -4<delta'<-4+8/r<0.

An embedding fixing z would preserve the physical negative determinant
of the old vertex plus d. Thus `K=Q(z)=Q(delta)`. The same argument
justifies each later field-equality step for a scalar determining qc.
In particular, in the `j=3` branch the triangle cycle puts gamma in K
and qc is a rational polynomial in gamma; its physical strictness
also gives `K=Q(gamma)`. This makes the trace degree there the actual
degree of the ground field, as required.

For positive gamma the scalar cycles are correctly normalized:
`w=4x^2`, `v=4gamma^2`, and `q=4x*gamma` are integral. In the norm-two
Gram, the triangle cycle is `-q`; `q^2=wv`. This permits all embedding
inequalities without asserting that x or gamma separately belongs to K.

## 2. Zero and finite contact cases

With gamma zero, `qc=Nw/2` and `K=Q(w)`. For finite nonzero x,
`1<=w<4` and `0<w'<2/N`. At degree at least two,

    0<Norm(w)<4(2/N)^(d-1)<=1.

The strict upper bound includes the endpoint `N=8,d=2`. Degree one
forces Q and is outside the stated field list. The case x zero is
geometrically impossible because c would be orthogonal to the entire
positive old vertex.

For positive gamma and x zero, the cosine-square trace forces degree
at most three with threshold `2/r<=1/5`. The complete small-totient
label list in the source is correct. Labels 7,9,14,18 and 10 have a
nonphysical square above 1/5; the rational labels and the square-root
2/5 fields are excluded. Only label12 over Q(sqrt3) remains.

For finite `j>=4`, completing the gamma-square gives
`qc'>=(4N/r)*(x^2)'`. An embedding fixing the physical square would
violate admissibility, so `K=Q(x^2)`, with all other squares below
`r/(4N)<=5/16`. The exact trace excludes degree five and above,
including equality at the degree-five threshold. The complete
totient-at-most-eight list is correct; each of the ten explicitly
excluded labels has the stated coprime exponent and conjugate above
5/16. The surviving square fields leave only j12 in the eight-field
scope. Its nontrivial sqrt2 square class permits e-label4 or12;
equal12 violates `qc'<1` already at r10.

For j3, the exact factorization

    qc-1=(gamma+1)(r*gamma+r-4)/2

forces every nonphysical gamma into `(-1,-1+4/r)`. The ordinary cosine
trace then reduces to labels 3,4,5,6,7,9. The stated affine, field and
conjugate tests leave only label6.

I rederived the common value

    R0=(7+4sqrt3)/4

for all three residual profiles `(0,12),(1/2,6),(cos(pi/12),4)`.
Writing `A=3r-7`, one has
`delta=2(r-2)/(A+2(r-2)sqrt3)`, and its trace gives exactly

    Tr(z)=4+4(7r-15)/(3r^2-6r-1).

This is strictly between four and five for every r at issue:
`3r^2-34r+59` is positive at r10 and increasing thereafter. Thus the
integrality contradiction closes all these finite and zero profiles.

## 3. Divisibility and the dotted-contact norm bound

For gamma zero and dotted x, the source's integer Bezout step is
valid over any ring of integers. For
`A=(N/g)w-r/g`, coprimality with `N/g` and
`A*delta=4N/g` imply `(A)|(4)`. The physical absolute value exceeds
four, and each nonphysical absolute value is strictly above
`N/g>=4`. This contradicts the norm divisor bound. No principal ideal
domain or integrality of x itself is used.

For gamma positive, put `u=v-4`, and use the source's integral B,C.
I independently expanded the determinant identity to obtain

    B=N*C+2u,  B*delta=-4N*u.

At a prime ideal, with `a=val(2)`, `n=val(N)`, `t0=val(2u)`:
if `n>t0`, then `val(B)=t0`. Otherwise the integral delta gives

    val(B)<=a+n+t0<=a+2t0=val(8u^2).

All valuations are nonnegative, so the first case also satisfies the
bound. Thus **(B) divides (8u^2)** with the correct direction. Taking
absolute norms and using the other delta conjugates, all strictly
between -4 and zero, gives

    |Norm(delta)| >= (N/2)^d/U,
    delta > 4(N/8)^d/U,   U=|Norm(u)|.

The physical bound `R>=25/12`, together with `s0<=5/16`, gives
`delta<=48/85`. The strictness needed in the later comparisons is
retained; none of the boundary inequalities are used in the wrong
direction.

## 4. Complete native label/norm table and the last prime obstruction

The finite e label is not subject to an arbitrary cutoff.
`F_ell=Q(cos(2pi/ell))` lies in K, so `phi(ell)<=10`. The displayed
totient list, after removing 7/14, is complete. The field containments
used for its prime-power labels agree with the eight-field scope:

|ell|possible K|U|
|---:|---|---:|
|5|K15,K20|25|
|8|K16,Q(sqrt2,sqrt3)|4|
|9|K9|3|
|11|K11|11|
|16|K16|2|

The norm formula `|Norm_F(u)|=Phi_ell(1)` follows by pairing the
cyclotomic conjugates of `1-zeta`; extension to K raises it to the
stated relative degree. Non-prime-power labels have U1. The rational
prime-power labels 3/4 are correctly reserved for the sharper final
argument rather than mistakenly assigned U1.

The lower bounds exclude U1 and labels8,9,16. For label11,
`gamma>2/3` gives `delta<80/231<4/11`. For label5 the exact comparison
`gamma>21/26` gives `delta<3/16`; if N is at least nine the norm lower
bound exceeds one quarter. At N8, the absolute integral norm lies
strictly between ten and twelve, hence equals eleven. Its sign is
negative, since precisely the other three conjugates are negative.

The final residue-degree obstruction is correct. In each real field
K15 or K20, 11 is unramified. Its residue degree is its order modulo
the real quotient by plus/minus one: 11 is not plus/minus one and
its square is one in both moduli. Therefore all primes above11 have
residue degree two and norm121. No integral principal ideal can have
norm11. This contradicts `Norm(delta)=-11` without any class-number
assumption.

## 5. Rational e-labels and conclusion

For ell4, `Z=q+1` is integral, `R=Z^2/4`, and
`qc=N*Z^2/4+1/2`. This gives the stated nonphysical bound and the
same Bezout divisor4 contradiction.

For ell3, `Z=2q+1` is integral, `R=Z^2/12`, and
`qc=N*Z^2/8+1/4`. Taking `g=gcd(N,6)` gives the divisor12. The
five ranks 10,12,13,15,16 have nonphysical absolute values strictly
greater than or equal to the required threshold12, with strictness
coming from `qc'<1`; the physical values also exceed12. At r11,
`B=3Z^2-11` is coprime to3, improving the divisor to4. At r14,
`B=2Z^2-7` is coprime to2, improving it to3. Their conjugate bounds
contradict these smaller divisors. This covers every r10 through16.

No missed native field, finite/zero contact case, valuation case,
conjugation-sign issue, or norm-degree mismatch was found. The theorem
passes as a restricted actual-quadrilateral exclusion with the stated
arithmetic and geometric premises. The source correctly leaves all
long-arm/spectator branches and the full refund composition open.
