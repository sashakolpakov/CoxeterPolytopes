# Native5/D19: a strict relaxation left by the selected centre and fork pressures

22 October2026. Manual obstruction supplement by sharp_ra_foundations_audit, at root's request. **This is not a Coxeter-polytope realization or a counterexample to any proposed dimension theorem.** It records precisely why the listed prefix/head/centre/fork constraints alone do not close the remaining centre-helper branch. No enumeration, mathematical import, test, scientific, Lean or TeX job was used. All frozen inputs and papers are unchanged.

Exact geometric notation and the native5 head are those of
`results/compact_whole_d19_native_five_head_fork_constraints_foundations_20261022.md`, SHA256 `4c6247e36ae06ab6ee9ad00d39ecd8db3ffd4bf11f59a2d83a0b1ce3eadbbfa1`. The fresh complete-centre reduction is the pending supporting candidate
`results/compact_whole_d19_native_five_head_complete_centre_reduction_foundations_20261022.md`, SHA256 `f71a0edacabfa1e4c21e48cd58cc763515738b8527b9fc9e17f8d007342bb445`; its theorem is not assumed in constructing this relaxation. Scalar increment upper bounds below are the elementary bounds already used in Section5 of the first input. No whole-B or larger-D geometric conclusion is transferred.

## 1. Three residual caps have legal selected products

Use unit old D19 normals, centre17 and forks18/19, with the same omega notation. Put a=phi,

    c_H=A_H p-(a/2)omega1, A_H²=(a-1)/2,
    c=A p-(6/5)omega17-omega1/2, A²=2543/50=50.86,
    b18=sqrt37 p-2omega18, b19=sqrt37 p-2omega19.

These are unit spacelike normals with positive time. Their old contacts have precisely the stated native5 head, pure far3/dotted-own centre, and two isolated dotted forks. The selected mutual magnitudes satisfy

    -b18.b19=3,
    -c.b_f=A sqrt37-209/5>1.57,
    -c_H.c=A_H A-(17/10)a>1.20,
    -c_H.b_f=A_H sqrt37-a>1.75.

For example A>7.131, sqrt37>6.0827, A_H>.555 and a<1.619 prove the strict lower bounds by rational multiplication. Every cap/cap product is therefore an allowed dotted value; none lies in a forbidden finite native gap. Both fork cosh lengths are sqrt38<10. The selected four-cap Gram data do not assert that these are facets of a compact polytope.

## 2. An explicit norm<.996 vector

Set f_i=sqrt(i)-.58/sqrt2 and d_i=sqrt(i+1)-sqrt(i). Define h=.3931, z=.0325 and

    x1=h+z=.4256,
    x2=x3=x4=h-z/3,
    x5=f5-4h,
    x_(i+1)=d_i for i5..12,
    x14=x15=x16=.1414,
    x17=2369/22500, x18=13/360, x19=0.            (1)

This vector is weakly chamber-ordered. The only equalities of adjacent coordinates are inside the first four block and the14..16 block. Since d12>.14144>.1414, every other ordering gap is positive. Arbitrarily small ordered perturbations make all old-wall products positive; the strict margins below and the norm upper bound survive. One can also raise the first13 coordinates slightly so that every prefix floor used below becomes strict.

The first four squares sum to4h²+4z²/3<.619518774. Moreover x5<.2536, the increment squares from5 through12 sum to<.2394, 3(.1414)²=.05998188, x17²<.1053² and x18²<.0362². Thus

    ||x||²<.619518774+.06431296+.2394
              +.05998188+.01108809+.00131044
           <.99561215<.996.                      (2)

For the increment bound, use d_i²<10/(40i+19); upper numerators at denominator10^6 are(45663,38611,33445,29499,26386,23867,21787,20041), summing239299. These are eight displayed rational inequalities, not an executed search.

Write S_i=x.P_i. We have S_i=f_i exactly for5<=i<=13. At14 the .1414 increment exceeds d13, and the next two .1414 increments exceed d14,d15, so S_i>f_i for14<=i<=16. Also

    S17=f13+.4242+2369/22500>3.72488>f17,
    S18=S17+13/360>3.7609.

Thus every .58 prefix floor at indices5..17 holds; this statement does NOT include indices1..4. The exact native5 H4 head requirement also holds strictly:

    S4=1.5724>2sqrt(a-1).

The last inequality follows from a<1.6181 and1.5724²>2.4724. These constraints remain compatible with(2).

The future unit point proportional to p-x is strictly inside all four selected cap halfspaces as well: its spatial products are respectively<.488<A_H, <6.627<A and<5.324<sqrt37. Thus the displayed relaxation does not obtain its small norm by putting the proposed interior point beyond a selected cap.

## 3. The listed actual-style edge and core inequalities all hold

The cap c at centre17 has actual-edge candidate cosh

    C=(34*(6/5)+1)/sqrt18=209/(15sqrt2)<9.853.

Its original-edge farthest pressure is

    S17>=sqrt17 sqrt((C-1)/(C+1)).

The right side is<3.724, while(1) gives S17>3.72488. This tests the edge inequality in addition to the corrected-core inequality; neither is silently omitted.

For the centre's future core c+(5/6)n17, put a0=5/12, b0=47/72 and c0=25/72. The normalized displacement and demand are

    V=P16+a0P1+b0e17+c0e18,
    G=sqrt(2543)/12-5sqrt22/72<3.877.

Direct substitution gives

    x.V>3.6196+.17733333+.06873024+.01253858
       >3.8782>G.

Here S16=f13+.4242>3.6196. The pressure uses the actual first-coordinate term, not a hypothetical lower bound obtained by subtracting unrelated pressures. The centre is still .50-weak: S17<3.725<sqrt17-.50/sqrt2.

Each isolated fork has the original-edge pressure

    S18 +/- x19 >=sqrt19 sqrt((sqrt38-1)/(sqrt38+1)),

whose right side is<3.71; (1) gives S18>3.7609. Its individual future core b_f+n_f/2 has normalized displacement and demand

    V_f=P17+(3/4)(e18 +/- e19),
    G_f=sqrt(37/2)-sqrt(3/8)<3.69.

The left side is>3.75. Finally b18+b19 is itself a future feasible negative core of square -4. Its normalized pressure is

    S18>=(sqrt37-1)/sqrt2<3.595,

also satisfied. Both signed fork orientations have therefore been checked. This is a finite stated set of constraints, not a claim that every feasible point arising from these normals, or every unknown actual wall, has been checked.

## 4. What an additional argument must supply

The relaxation simultaneously has the correct old chamber signs after a small perturbation, all .58 prefix floors5..17, the exact native5 head H4 pressure, legal selected cap products, both fork-cosh bounds, the centre actual-edge/own-core pressures and both fork actual-edge/own-core/joint-core pressures. Its norm remains strictly below1. Therefore those constraints alone cannot imply a norm contradiction, even for an unequal-parameter theorem that includes the symmetric delta18=delta19=2 instance. A new actual-wall or feasible-point inequality is needed. This does not invalidate any restriction proved from a larger set of actual data.

A possible additional route is the ACTUAL adjacent centre vertex q, whose finite forest is A17+2A1. It is not asserted to be farthest. If the original vector is z=E(p-x) and q=Cp-sP17/sqrt17, then

    tau=(-z.q)/E=C-(s/sqrt17)S17<=1.

For an actual edge from q in either retained A1 direction n_f, with endpoint r=Dq-sqrt(D²-1)n_f, the ORIGINAL farthest condition gives the exact inequality

    tau*D-sqrt(D²-1)*x.n_f<=1.                   (3)

This identity is legitimate for the same original z without reclassifying q as farthest. Any quantitative use of the A17 component theorem, an I/A1 theorem, a no-common-opposite statement or a nonshort A1 bound at q must first establish its hypotheses with this slack. None is transferred to q in this supplement. No polytope, finite-volume completion, new dimension ceiling, actual-field restriction, runtime or Lean result is claimed.
