# The old-seven defect forces every other irreducible rank-two block above one fifth

21 October2026. New complete manual supporting proof, pending FULL root
and outside-author review. No scientific, symbolic, numerical, enumeration,
Lean or TeX job ran. No frozen source, manuscript, receipt or ceiling
changed. This is one geometric lemma in the existing compact family.

**Proposed theorem.** At a vertex farthest from a strict interior point
of a compact hyperbolic Coxeter polytope, suppose a whole I2(7) component
has an ACTUAL opposite wall supported only there, with one native3
contact and zero contact at its other node. Then EVERY OTHER whole
irreducible rank-two component has actual projection energy greater
than1/5. Arbitrary finite spectators are allowed; no arithmeticity is
assumed.

Consequently this typed OLD I2(7)/NEW isolated3 defect cannot coexist
with both a whole A_r component,r>=14, and another irreducible T2.
In particular it excludes the ENTIRE retained general dimension20
family A16+that defect+ordinary T2. Other retained families and the
ambient compact ceiling20 remain separate.

## Inputs and same-vector normalization

| Input | SHA256 | Interface used |
| --- | --- | --- |
| results/compact_universal_small_a_component_energy_20261020.md | 801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6 | Section1's actual two-wall law, finite white endpoints, inverse positivity and universal positive foreign cost>=1/4 |
| results/compact_farthest_rank_four_forest_exclusion_20261020.md | 68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6 | The precise OLD I2(7)/NEW3 defect: isolated support, squared height in(0,1/3), actual energy>9/100 |
| results/compact_large_an_rank14_five_sevenths_uniform_transfer_20261021.md | 2af6d0032c431610a7ece98c3fc8d590e79b1bd75ade7d3c77486c2d51bb9f27 | Actual A_r energy>5/7 for every r>=14, arbitrary finite spectators |
| results/compact20_arithmetic18_remaining_farthest_forests_20261021.md | a49d2c52ce5b5cb0a77b655cb0552229d9c415aab2df7cf4e0612c88cd796ad8 | The retained general A16/old-seven/T2 family, not an arithmetic premise |

The author has shared earlier actual-point/framework contributions;
these inherited foundations are not claimed as independently authored
here. The new rank-two cutoff and full native/foreign case proof follow.
External premises remain compact acute simplicity, native finite angles
and the product gap for DISTINCT actual outward unit wall normals:
a strictly negative product must be at most-1/2.

Write the interior point as E(p-x), ||x||²<1, and e_T=||x_T||² for
orthogonal whole finite components T. For an actual edge i let
tau_i=tanh(L_i/2)>0 and u_i its normalized dual. Then

    x_T.u_i>=tau_i,
    tanh(L_i)tanh(L_j)>=u_i.u_j
       whenever the ACTUAL opposite walls are distinct.

All edges/endpoints below are actual. No formal Gram completion,
field criterion, or changed A-component threshold is introduced.

## 1. The typed wall forces a high opposite in every other block

Let D be the typed I2(7), and write its isolated wall as
c_D=sqrt(h_D)*p-v_D, supported spatially only in D. With
d7=csc²(pi/7), its single3 contact gives

    h_D=d7/4-1 in(0,1/3).

Let c be an actual opposite for any direction in another component.
Its finite white endpoint retains the WHOLE old D block. A positive
contact there contributes at least d7/4>1 to the white Schur cost,
contradicting positive definiteness. Thus c is orthogonal to every
old root in D. This covers all finite/dotted magnitudes and all
additional contacts; it is not an assumed isolated-opposite profile.

Write c=sqrt(h)*p-v. Its D projection is zero, so c is distinct from
c_D and

    c_D.c=-sqrt(h_D*h)<0.

The actual Coxeter product gap implies h_D*h>=1/4. Consequently

    h>=1/(4h_D)>3/4.                                (1)

Every actual opposite in the other component obeys this bound.

## 2. Energy at most one fifth forces one common rank-two opposite

Suppose contrary e_T<=1/5 for a whole T=I2(m),m>=3. Put
k=cos(pi/m), d=1/(1-k²). Its two normalized duals have correlation k.
For s_i=x_T.u_i,0<s_i<1, put F(s)=2s/(1+s²), increasing on(0,1).

If m>=4, distinct opposites would give

    k<=F(s1)F(s2)<=F(sqrt(1/5))²=5/9<1/sqrt2<=k,

impossible. For m3 use log F concavity and
(s1+s2)²<=e_T*||u1+u2||²<=3/5. Distinct opposites would imply

    1/2<=F(s1)F(s2)<=F(sqrt(3/20))²=240/529<1/2.

Here(log F)''=-1/s²-2(1-s²)/(1+s²)²<0. Thus ONE actual wall c
closes both directions for EVERY old label m>=3. Since each direction
has tau_i²<=e_T<=1/5, both actual lengths satisfy

    cosh L_i<=3/2.                                  (2)

## 3. Complete common-wall contact and spectator formulas

Let b,g>=0 be c's two contacts inside T. Let Q>=0 be its ENTIRE
foreign Schur cost. The two actual finite endpoint Grams give

    Q+b²<1, Q+g²<1,
    h=d(b²+2kbg+g²)+Q-1>3/4,
    C_b=cosh L_b=sqrt(d)*(b+kg)/sqrt(1-Q-g²),
    C_g=cosh L_g=sqrt(d)*(g+kb)/sqrt(1-Q-b²).       (3)

The strict endpoint inequalities show that every positive b,g is a
FINITE native cosine below1; parallel/dotted own contacts have thereby
been excluded from this LOW-ENERGY common-wall branch, not omitted.
At least one contact is positive, since otherwise h=Q-1<0. Every
positive foreign cost is>=1/4 by the retained whole-component inverse
bound. Thus Q is zero or at least1/4; all its components, native
contacts and multiplicities remain included in(3).

### One positive contact

By reversal take b=beta>0,g0. The common mate in(3) satisfies exactly

    C_g²=1+h/(1-Q-beta²).

Equations(1)–(2) imply1-Q-beta²>3/5, hence Q+beta²<2/5.
As beta²>=1/4, a positive Q>=1/4 is impossible. Therefore Q0,
and the native gap between3 and4 forces beta=1/2. Now

    h=d/4-1>3/4, so d>7.

But d8=4+2sqrt2<7. Thus the old label m>=9. Using sin t<t and
the elementary pi²<10,

    d>=csc²(pi/9)>81/pi²>81/10>31/4,
    C_g²=(d-1)/3>9/4,

contradicting(2). This covers every single-contact native label and
every possible foreign forest.

### Two positive contacts with Q>0

Both contacts are at least1/2. If their maximum is at least1/sqrt2,
take g to be that maximum. Since k>=1/2,d>=4/3 and Q>=1/4,

    C_b² >= (4/3)(1+1/sqrt2)²
           =2+4sqrt2/3>9/4.

This contradicts(2). Otherwise both contacts are exactly1/2. If
old m3, (3) gives h=Q<3/4, contradicting(1). If m>=4, then
k>=1/sqrt2>2/3 and

    C_b² >= (1+k)/(2(1-k))>5/2>9/4.

Thus no two-contact branch with positive foreign cost survives.

### Two positive contacts with Q0

If the larger contact g>=cos(pi/5)>4/5, then b>=1/2 and(3) gives
C_b²>3>9/4, already using only k>=1/2,d>=4/3. If both contacts
are4, then C_b²=(1+k)/(1-k)>=3>9/4. The remaining finite list
is therefore just(3,3) or(3,4), up to reversal.

For(3,4) and old m>=4, (3) gives C_b>=2>3/2. For old m3 its
squared height is h=sqrt2/3<3/4, contradicting(1).

For(3,3) and old m>=5, k>=cos(pi/5)>4/5 gives
C_b²=(1+k)/(3(1-k))>3>9/4. Old m3 has h0, while old m4 has
h=1/sqrt2<3/4. Both contradict(1).

These are ALL two-contact finite labels: labels>=5, both4, mixed3/4
and both3. Hence every common-wall possibility contradicts(1) or(2).
This proves e_T>1/5 with arbitrary finite spectators.

## 4. The complete A16 mixed-family consequence

The inherited general theorem gives e_A>5/7 for any whole A_r,r>=14,
and the ACTUAL old-seven defect itself gives e_D>9/100. Adding the
new e_T>1/5 for the SAME original x yields

    ||x||²>=e_A+e_D+e_T>5/7+9/100+1/5=703/700>1.

Thus A_r plus this typed defect plus another irreducible rank-two block
is impossible for every r>=14, even with further orthogonal spectators.
In particular this excludes all native labels of the retained ordinary
T2 in the general A16 mixed dimension20 family, without changing the
A16 loss threshold or imposing arithmeticity.

The theorem concerns an ACTUAL isolated3 wall on OLD I2(7), not a bare
old label7 or the different OLD A2/NEW7 defect. The A18-plus-defect
family has no additional rank-two block and is not excluded here.
Other general dimension20 and arithmetic dimension18 families, full
composition and any improved ambient ceiling remain separate. The six
papers and all frozen evidence are unchanged.
