# B18 short anchors and weak prefixes: exact necessary comparison

21 October 2026. Mathematical specification for the new worker
`scripts/check_b18_short_anchor_prefix_signature_census_20261021.py`.
**Pre-execution specification; no output or B18 exclusion is asserted.**
This supporting computation belongs to the existing compact-arithmetic
research family, not a new manuscript. Root authored this specification
and worker; an outside-author source review is required before execution.

## 1. Geometric inputs and exact conclusion sought

The accepted short-edge theorem
`arithmetic_whole_b18_short_cap_inventory_and_degree_three_reduction_20261021.md`
(SHA256 `0d9cffcc58700e38985e31865380457bbbf885d0c329cb6c53ea69c089f675e1`)
and the complete .48 weak inventory
`arithmetic_whole_b18_point48_weak_prefix_inventory_and_cap_census_domain_20261021.md`
(`04dfdc05a18b6aa2d7d7ba71a0e6b3824601dce14a797e113a038f209ff003cb`)
are the geometric premises. They concern the ACTUAL first-type field of
a hypothetical compact polytope with whole B18 at a farthest vertex.
They supply degree two or three, the complete lists below and at least
one weak prefix index. They are not re-proved by this worker.

The proposed computed proposition is a COMPLETE NECESSARY comparison of
these lists: every actual short-anchor/weak-prefix pair occurs among
the retained output pairs. Quadratic pairs pass exact same-field and
time-square tests. Cubic pairs pass only discriminant/norm squareclass
tests: every such match is retained UNRESOLVED. The computation supplies
neither cubic field isomorphisms, mutual native Gram admissibility,
actual realization, a whole-B18 exclusion, nor dimension seventeen.

In fact the finite silver anchor and the silver far-helper prefix at
index five give a genuine time-square match. This fixture must survive.
An empty quadratic output would be an error, not a desired conclusion.

## 2. Complete finite domains

Write a monic quadratic as X^2-TX+N, and a monic cubic as
X^3-TX^2+SX-M. Every range includes its integer endpoints. The anchor
parameters B, t and z are integral generators of their ACTUAL fields;
z is the near-helper parameter t+2, not the unshifted t.

| Anchor parameter | Quadratic outer coefficient box | Cubic outer coefficient box |
|---|---|---|
| Isolated B | T=5..37; N=1..8 | T=5..37; S=1..16; M=1 |
| Far helper t | T=3..6; N=-3,-2,-1,1,2 | T=2..6; S=-7..4; M=-2,-1,1,2 |
| Near helper z | T=5..7; N=-2,-1,1,2,3 | T=5..7; S=-4..8; M=-2,-1,1,2 |

These give 299 quadratic and 924 cubic raw entries, before root tests.
Section 7 of the complete weak inventory proves these coefficient
bounds from trace, symmetric sums and norm, using all conjugates.
No field database, sample or bound on a chosen defining extension is used.

For B require one physical root in (4,338/9) and every other root in
(0,2/9). Put a=(sqrt(19)-1)/9 and b=(sqrt(19)+1)/9. For far t require

    2sqrt(2)<t<(13sqrt(19)-1)/9; every other root in (-b,a).

For near z require

    2+2sqrt(2)<z<(13sqrt(19)+1)/9; every other root in (-a,b).

The single finite anchor has theta=sqrt(2), raw norm 2, projection
Omega17+(theta/2)Omega18, and squared time 24+17theta.

An isolated weak prefix at i has parameter alpha and raw norm 2alpha.
In degree two use i=5..16, T=5..16, and integers 1<=N<32/i.
The twelve N upper bounds are 6,5,4,3,3,3,2,2,2,2,2,1: 420 entries.
In degree three use T=5..16 and

    i=5: S=1..12,M=1..2; i=6: S=1..10,M=1;
    i=7: S=1..9,M=1.

There are 516 entries (288+120+108). In both degrees retain one physical
root in (4,16), all other roots in (0,2/i). Add exactly the twelve
quadratic typed rows in Section 6 of the inventory: own4/right-long3
at i5; far3 at i5..10; golden far-D at i5,6; golden near-D at i5;
golden far3 at i5; and native12/SHORT18 helper at i16. The last row is
essential and is not removed by any older .61 threshold argument.

The outer boxes may contain false positives. Keeping one is harmless;
omitting a genuine member is not. No enumeration count beyond the raw
loop sizes is asserted until the reviewed program actually runs.

## 3. Exact root tests, with termination and completeness

All polynomial arithmetic is over the rational numbers. A quadratic
is irreducible precisely when its positive integral discriminant is
not a square. For the monic cubic the rational-root theorem reduces
irreducibility to testing BOTH signs of every positive divisor of
abs(M). Negative constants are allowed; none of the domains has M=0.
An irreducible polynomial of degree two or three has no rational root
and is squarefree in characteristic zero.

Use the derivative and negative Euclidean remainders to construct its
Sturm sequence. Euclidean degrees decrease strictly. Variation counts
at rational endpoints give the number of distinct roots in each open
interval; endpoints are never roots. Require exactly one physical and
d-1 other roots in the stated disjoint intervals. This also proves
total reality and identifies a unique physical root.

The irrational anchor bounds are tested exactly as follows. For far t,
first count roots in the rational supersets (2,25/4) and (-3/5,3/8).
At the physical root require

    t^2-8>0, 81t^2+18t-3210<0;

at every other root require (9/2)t^2+t-1<0. The last inequality has
roots -b,a; the physical inequalities, together with t>2, give exactly
the physical interval above. For near z use rational supersets
(4,13/2) and (-3/8,3/5), physical tests

    z^2-4z-4>0, 81z^2-18z-3210<0,

and nonphysical test (9/2)z^2-z-1<0, whose roots are -a,b.
Here 3211=169*19, so squaring the positive expressions 9t+1 and 9z-1
introduces no extra root branch. The outer intervals are genuine
supersets, as sqrt(19)<35/8, a<3/8 and b<3/5.

To evaluate a test polynomial at an algebraic root, reduce it modulo
the irreducible minimal polynomial. A zero remainder is exactly zero.
A nonzero remainder of smaller degree cannot vanish at that root.
Bisect rational intervals, keeping the unique root with Sturm counts,
until rational Horner interval evaluation excludes zero. This halts:
the isolated simple root has positive separation from the other roots,
the test is continuous and nonzero there, and the interval widths tend
to zero. Initial isolation of d-1 roots likewise terminates by positive
separation. Thus this is an exact decision, not a floating tolerance or
a resource-limit pruning rule. Every domain member passes precisely
its stated open-interval tests.

## 4. Raw times, independent norms and signed squareclasses

With old B18 raw duals Omega_i=P_i, their products are min(i,j).
If a raw cap has projection v=sum c_i Omega_i and raw norm N0,
its squared time is U=v^2-N0. The worker recomputes this identity
modulo the minimal polynomial for EVERY retained model. The isolated
prefix time is alpha(i alpha-2). Anchor times are

    B(9B-2)/2, (9/2)t^2+t-1, (9/2)z^2-z-1,

and the finite silver time given above. The shifted near expression
equals (9/2)(z-2)^2+17(z-2)+15.

For any time polynomial f, form the multiplication matrix of f in
Q[X]/(p), in the basis 1,X,...,X^(d-1). Its determinant is the field
norm. Gaussian elimination and the full permutation determinant
(at most six terms) must agree. For isolated and dotted-anchor rows
also check an independently expanded closed norm formula. If p is
quadratic, use its T,N; if cubic, use T,S,M. The formulas are

    prefix quadratic: N(i^2 N-2iT+4);
    prefix cubic:     M(i^3 M-2i^2 S+4iT-8);
    B quadratic:      N(81N-18T+4)/4;
    B cubic:          M(729M-162S+36T-8)/8.

For a far/near helper put e=+1/-1, respectively. Then

    quadratic: [81N^2+18eTN+40N-18T^2-4eT+4]/4;
    cubic: [729M^2+162eMS+360TM+116eM-162S^2
            -36eTS-80S+36T^2+8eT-8]/8.

These formulas also allow negative N and M. The norm is nonzero,
negative in degree two and positive in degree three, since the time
is physical-positive and negative at every other real embedding.

For a reduced rational norm n/d, d>0, its signed rational squareclass
is sign(n) times squarefree(abs(n)*d), NOT the squarefree numerator
alone. Indeed n/d=sign(n)*squarefree(abs(n)*d)*(k/d)^2. Thus denominators
such as 8 in cubic anchor norms cannot be silently dropped. Trial
division yields the unique squarefree part and verifies the displayed
factorization exactly; its increasing integer loop terminates.

## 5. Necessary matching and exact quadratic square decisions

Every actual pair has the SAME ACTUAL field and physical embedding.
Since its two time coordinates lie in one one-dimensional K-line,
U_prefix/U_anchor is a square in K. Necessarily their degrees agree,
their polynomial discriminants have the same rational squareclass,
and their time norms have the same signed rational squareclass. The
discriminant statement follows from discriminant(order)=
discriminant(field)*index^2 for an integral generator; the time
statement follows by taking the field norm of a K-square.

Group prefix models by these three signatures. Every compatible pair
must lie in the same group. For CUBICS retain EVERY signature match
with an explicit UNRESOLVED field-map/K-square status: neither
condition follows from these signatures. Their failure is not guessed.

For quadratics the common positive squarefree discriminant D identifies
K=Q(sqrt(D)). An integral generator with polynomial X^2-TX+N maps to
(T+k sqrt(D))/2 where discriminant=k^2 D and k>0. This is its unique
larger, hence physical, root in every domain above (also for the typed
sqrt(2),phi,sqrt(3) rows). Evaluate both times exactly in this basis.

For x=a+b sqrt(D)>0, a square y=u+v sqrt(D) must satisfy

    u^2=(a+r)/2, r=+sqrt(Norm(x)) or -sqrt(Norm(x)),
    2uv=b.

First require the rational norm to have a rational square root.
Try BOTH signs r and exact rational square roots of (a+r)/2.
If u!=0 put v=b/(2u). If u=0, necessarily b=0; test v^2=a/D.
Verify y^2=x exactly and choose its positive physical sign. Conversely
every K-square has one of these cases; choosing the positive rational
u loses only the irrelevant overall sign. Thus the decision is complete,
including pure-radical roots. No square root is silently adjoined to K.
The sign of a+b sqrt(D) follows from rational signs and a^2-Db^2,
not a floating approximation. Quadratic nonsquares can be rejected;
all squares remain for later geometric analysis.

As an independent necessary survivor, the finite silver time U0=24+17s
and far3/index5 time U5=16+12s, s=sqrt(2), satisfy

    U5/U0=12-8s=[2(s-1)]^2.

Their output row must have positive square root -2+2s. This fixture
prevents misreporting the finite silver branch as an empty comparison.
It asserts no actual completion or whole-B18 realization.

## 6. Execution, output and logical boundary

The worker visits the four finite coefficient domains, all twelve typed
prefix rows and the finite anchor; applies the exact tests above; saves
every retained model with minimal polynomial, physical-first rational
root boxes, raw projection/norm/time, discriminant and exact norm; then
saves every anchor's signature targets, rejected quadratic nonsquares
and retained pairs (with a positive root or unresolved cubic status).
All IDs are local row identifiers, not geometric equivalence classes.

All loops are finite or have the proved decreasing/isolation termination
arguments. No symmetry quotient, heuristic pruning, numeric root
solver, optimizer, native-angle guess or unreported survivor deletion
occurs. An actual geometric pair maps into the proved coefficient box,
passes the exact root tests and necessary signatures, and, in degree
two, its actual time square passes the complete quadratic decision.
This proves the precise necessary completeness conclusion in Section 1.
The finite loop/decision arguments prove termination independently of
the process cap. A resource refusal or timeout is never an empty result.

Run once under `run_bounded_research_check_20261018.py --parent-lock`,
only after full outside-author source review, with one thread, 59 seconds,
1 GiB process-tree cap, the shared lock, 40 GiB disk/6 GiB RAM start gates
and the unchanged 8 GiB/3 GiB abort gates. Imported rational/Sturm helpers
are pinned at `29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932`;
their CLI and old A18 scientific routines are not called. Source, this
specification and both geometric premises are hashed before and after.
Write one NEW compact result JSON and one NEW execution ledger; never
overwrite a historical result or launch concurrent proof work.

The exact retained counts and pair types must be reported from that
actual terminal output, independently read and reviewed afterwards.
A later proof must handle every surviving weak pattern using actual
geometry of the SAME original farthest vector. No full Lean, new bound
or new example claim follows from this computation alone.
