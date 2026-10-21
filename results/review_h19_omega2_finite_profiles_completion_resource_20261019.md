# Independent review of the smaller twenty-seventh H19 mirror

19 October 2026. **Written/source PASS; no flaw found in the stated
bounded construction and infinite-volume conclusion.** Reviewed report:
`results/h19_omega2_finite_profiles_completion_20261019.md`, SHA256
`5568426ae65b447a9aef5ea97b943499e6ac908de0bcdcc037efb7ffccd8043a`.
I read that report, the complete prepared search source
`scripts/search_h19_omega2_finite_profiles_20261019.py`, SHA256
`02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e`,
and its preparation report. All were preserved unchanged.

This review used paper algebra and source reading only. It did not
rerun the search, checker or vertex census. The independently supplied
generation/replay timings and resource ledgers remain execution
evidence supplied by the coordinator, not measurements made here.
The report's then-pending fixed27-wall checker is a historical status;
any subsequent execution is recorded separately.

## 1. The bounded search does what its scope says

The finite-profile and congruence inputs are the separately pinned
52/51 ->40 ->34 proofs. The component-choice enumeration keeps actual
node labels, applies the Schur q<2 test, and checks the exact identity
K=30(2-q). Its modulo3 filter correctly uses the original norm
equation, not the vacuous K=R² reduction. The prime-D test is correct.

The coefficient box covers every x=a+bphi with0<x<=600 and |x'|<=2:
b=(x-x')/sqrt(5) is between-1 and301, and the six integers around
floor(b/phi) contain every possible a. The floor formula uses integer
square roots and is exact. The subsequent field inequalities remove
the conservative excess; no approximate height is used.

The field-square test exhausts c+e*sqrt(5) by both signs of the
rational square root of its rational norm and the two candidates
c²=(a plus or minus sqrt(N))/2. It handles c=0 and verifies each
returned square. Both field-root signs are then tested. T/y native
integrality is checked after the scalar Pell equation, rather than
assumed from it. Every integral candidate receives all26 signed
native pair and angle tests, plus the actual endpoint checks.

The angle menu is complete for these norm classes by the pinned
field-angle proof. An unreachable extra m=4/6 equality in the code
does not admit a forbidden pair: a K-valued product with those
normalized squares cannot exist in the relevant norm classes.
Old-center positivity is not used as a rejection filter; the
p+epsilon*c construction instead supplies a strict center when
necessary. No full-height or other-star completeness is claimed.

## 2. Direct check of the smaller mirror

For the displayed root q the native products give exactly

    Gq=(-8-11phi,1;0,1,phi,0,0,0,0,0),
    q²=2,  (p,q)=-16-23phi,
    (c,q)=-(68+95phi)/10.

Its coordinates are integral. They are primitive, since its two
binary coordinates sum to the unit phi. Norm2 therefore gives a
full-lattice reflector. Its profile beta=1,eta_1=1,T=phi has
qW=7/6, R=20 and K=25, agreeing with A=20+30phi.

I independently substituted its sparse covector into EVERY retained
old and later root. The resulting full26-wall ledger is:

|Wall|Product with q|
|---|---:|
|old0,old2|-1|
|old1,old4,old5,old6,old7,old8|0|
|old3,old13,old19|-phi|
|old9|-2-phi|
|old11,old15|-1-2phi|
|old12,old16,old17|0|
|old14|-2phi|
|old18|-1-phi|
|r,z|-b*phi³|
|k|-7-11phi|
|w|-b*phi^4|
|a|-3-4phi|
|d|-12-20phi|
|s|-3-7phi|

There are exactly two angle3, three angle5, nine right, and twelve
strictly ultraparallel pairs, with no positive or parallel entry.
The preserved common strict center and acute projection give actual
facets for all27 inequalities, not just a compatible normal set.
The q mirror crosses the old chamber interior, so it is genuinely
additional to the previously generated reflection subgroup.

## 3. The new finite endpoint and 85-vertex lower bound

The endpoint parameter and square check exactly:

    tau=(16+23phi)/(20+30phi)=(11-2phi)/10,
    (p+tau*omega2)²=71/2-22phi<0.

For the sign, N(71-44phi)=-19 while its conjugate is positive.
Future orientation follows from its product p² with p. The old
first-hit parameter was(13+8phi)/33; subtracting tau gives
(-233+146phi)/330>0, for instance using phi>8/5. Thus the new
endpoint lies strictly before every old positive-wall first hit.
The complete saved old first-hit proof makes every old inequality
feasible there, with exactly the eight white walls active.
The new q equality is exact. The white-plus-q Gram is2A4+A1,
positive definite of rank9, so this is an actual finite vertex.

The products with r,w are strictly ultraparallel, so the two whole
compact H8 hyperplanes and their41+41 vertices are retained on the
correct sides. The projected strict center fixes that side; disjoint
hyperplanes cannot change sign on the connected hyperbolic sheet.
The old p is strict on q. The earlier vertex p_s=p+D*u*s is also
strict since

    (p_s,q)=(p,q)+D*u*(-3-7phi)<0.

The new endpoint is strict on old3, whereas p,p_s lie on it, so it
is distinct from both. All three lie strictly off the two caps.
The justified count is therefore41+41+1+1+1=85, conditional only on
the inherited H8 completeness for those two41-element sets. It is
not a complete ambient vertex count.

## 4. The two search outcomes cannot be appended together

The four coordinates/covector entries used in the displayed product
of q_L and q give4+5phi>0 by direct expansion. Thus their selected
inward orientations cannot coexist in an acute extension. The report
correctly chooses only the smaller q and does not count both results
as simultaneous walls. Existence and compatibility of the other
individual bounded survivor are inherited from its saved search
output, not a second full manual pair ledger in this review.

## 5. The retained escape really gives infinite volume

For E=p+((8+21phi)/44)*omega3, direct expansion gives

    E²=(1557-955phi)/484>0,
    (E,q)=(5phi-18)/11<0.

The positive square follows already from phi<13/8. All26 preceding
wall inequalities are inherited from the fixed26-wall escape ledger;
the last displayed product supplies the missing twenty-seventh.
The future orientation is fixed by (p,E)=p²<0. Perturbing toward
the strict center leaves a strict feasible spacelike point; the
segment crosses the future null cone inside every wall halfspace.
Hence an open ideal patch, and infinite hyperbolic volume, follows.

The larger finite-vertex count and closure of this initial edge do
not prove global boundedness. No new finite-volume H9 example,
complete mixed-norm chamber, nonreflectivity theorem or dimension
record is claimed. The source's limited execution and mathematical
status are consistent with these conclusions.
