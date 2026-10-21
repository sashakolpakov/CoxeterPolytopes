# Spectator contacts in the isolated-d quadrilateral row

19 October 2026. Root-agent arguments independently checked and
written out by batch_guard. Paper only: no scientific job,
enumeration, calibration, Lean, TeX or build was run. Earlier
source-bound files are preserved unchanged.

This note excludes the remaining e-plus-spectator and
path-plus-spectator continuations of c. It then audits the complete
isolated-d row over the eight target fields. It does not install
a simultaneous refund or a dimension bound.

## 1. Complete actual-quad hypotheses and beta>1

The entire old compact finite vertex is D_r+T_s, r+s=20,
10<=r<=16. The spectator T_s is D_s,4<=s<=10, or E6,E7,E8.
Let a be the D_r fork, b,e its short leaves, and l_1,...,l_L
the long arm. Write N=r-2>=8,L=N-1 and let

    W=A1(e)+A_L+T_s

be the COMPLETE white support of the actual quadrilateral a,b,c,d.
Assume d contacts ONLY b among the old vertex roots, and

    c.b=c.d=0.

In particular b,d are orthogonal to W and a.d=0. Let
x=-a.c>=0 and beta=-b.d. The native actual-wall condition is

    x=0, or x>=1/2.

More precisely 0<x<1 is a native cosine cos(pi/j),j>=3;
x>=1 is also allowed. This uses the ordinary acute Coxeter-wall
fact, not an assumption that ambient opposite walls are dotted.
For two actual acute walls with Gram product in(-1,0), the positive
two-normal projection of an interior point gives an actual ridge,
so their finite angle is a Coxeter angle. This supplies the gap.

Here beta>1 is a CONSEQUENCE of the stated actual geometry.
Project the finite white span off to the hyperbolic plane of the
actual quad. The b,d normals are unchanged, unit and orthogonal
to the projected c normal. Thus its two consecutive corners bc,cd
are right angles. The c side has positive finite length ell, so
b,d are the two supporting geodesics perpendicular to that common
geodesic at distinct points. With the outward orientations their
product is -cosh(ell). Hence

    beta=cosh(ell)>1.

This reasoning needs the ACTUAL compact face and its positive side
length. A formal four-root pattern alone would not justify this
geometric assertion. It does not apply the same dotted conclusion
to a,c, whose projections can change.

Actual cyclic integrality gives

    delta=4beta^2-4 in O_K,       delta_0>0.

The full actual ground field K is totally real. In the new cases
below A=2x belongs to O_K from the indicated actual cycle. We do
NOT assume that a unit-normal product x or beta belongs to K in
every other isolated-d branch.

## 2. Complete finite inventory when c touches T_s

The finite corner W+b+c contains the full positive-definite W+c.
A new c can contact at most one node in each connected old factor:
two contacts would create a cycle. Since T_s has a fork, every
finite connected extension is D/E and every new bond is simple.
The complete one-root table, with t_T the square remaining after
projection off the entire T_s, is

|spectator|contact|t_T|
|---|---|---|
|D_s|long endpoint; any equivalent leaf if s=4|1/2|
|D5|short leaf|3/8|
|D6|short leaf|1/4|
|D7|short leaf|1/8|
|E6|an endpoint of a length2 arm|1/3|
|E7|the length3 arm endpoint|1/4|
|E8|none|none|

This is the finite D/E continuation table, with residuals given
by the unit-normal determinant ratios. The full proof and entire
spectator qualification are in
`arithmetic_rank20_fork_spectator_continuation_research_20261019.md`,
SHA256 `95b0102731bb12a43416197260d65ff2998cc0a92db189722bee8ad43bf2d754`.

There are exactly three possibilities for the remainder of c's
white contacts:

1. no e contact and no long-arm contact;
2. an e3 contact, with no long-arm contact;
3. a simple contact to a long-path endpoint, with no e contact.

For case2, the connected factor T_s+c+e contains a fork and has
rank at least6, so an e4 or larger bond is not finite. Subtracting
the independent e projection cost1/4 gives exactly

|spectator contact|complete factor|t=t_T-1/4|N|
|---|---|---|---|
|D_s long endpoint, including D4 symmetry|D_(s+2)|1/4|8,...,14|
|D5 short leaf|E7|1/8|13|
|E6 length2 endpoint|E8|1/12|12|

The other t_T values give nonpositive t; they are not finite.

For case3 the joined path-plus-spectator component has rank18,
so it must be D18: the spectator must be D_s, c must attach to
its long endpoint (any leaf for D4), and to one endpoint of A_L,
both simply. An internal path attachment would make a second fork;
an E spectator has two long fork arms and cannot fit in D18.
The COMPLETE W+c diagram is D18+A1(e), of rank19, not D19+A1.
Its residual is

    t=1-1/2-(N-1)/(2N)=1/(2N).

No e contact can be added: c would make a second fork in a
connected rank19 diagram. Likewise there is no finite simultaneous
e+path+T continuation. This exhausts c's contacts with the complete
white support without suppressing outside spectator contacts.

## 3. Common Schur identity

Let t be c's square after projection off W, and write the residual
a.c product as -(x+kappa). The inherited unit-normal inversion and
actual corank identity give

    s=(N+2)/(4N),
    R=(x+kappa)^2/t,
    (R-s)*delta=1.                                 (1)

The old-full-vertex-plus-c projection cost is

    qc=1-t+2N(x+kappa)^2,

and its nonphysical conjugates satisfy qc_i<1. These are the
actual arithmetic/geometry interfaces of the previous sources,
not a formal finite-diagram realizability claim. The variables
used below are integral from actual cycles; no explicit integral
basis or false field-generation premise is required.

## 4. e3 plus spectator: all three residuals

The actual a,e,c triangle gives A=2x in O_K; at x=0 it is zero
and integral as well. Set V=2A+1=4x+1. Here kappa=1/4, so

    R=V^2/(16t).

### Residual t=1/4

Put

    C=NV^2-(N+2)=2D,
    D=2N*A(A+1)-1.

Equation(1) gives D*delta=2N. Since D=-1 modulo2N, its
principal ideal is coprime to(2N). As D divides2N, D is a unit.
If x=0 then D=-1 and delta=-2N, contrary to delta_0>0.
Otherwise the native gap gives A_0>=1 and D_0>=4N-1>1.

The nonphysical inequality qc_i<1 gives V_i^2<2/N. Hence

    C_i<-N,       D_i<-N/2<-1.

Every conjugate has absolute value greater than1, contrary to
abs(Norm(D))=1. This works over every totally real K.

### Residual t=1/8, N=13

Set

    C=2NV^2-(N+2)=8N*A(A+1)+N-2.

Equation(1) gives C*delta=4N=52. Since C=104A(A+1)+11,
it is11 modulo52. The rational integers11,52 are coprime,
so C is coprime to52 and must be a unit. Physically A_0>=0
gives C_0>=11>1, even at x=0. Nonphysically qc_i<1 gives
V_i^2<1/N, so C_i<-N=-13. This again contradicts its unit norm.

### Residual t=1/12, N=12

Set

    3NV^2-(N+2)=2D,
    D=18V^2-7=72A(A+1)+11.

Equation(1) gives D*delta=24. The residue D=11 modulo24
is coprime to24, so D is a unit. Physically D_0>=11>1.
Nonphysically qc_i<1 gives V_i^2<1/18, so D_i<-6.
This is the same contradiction. Thus ALL e3-plus-spectator
cases are excluded, including zero, native finite and dotted x,
without a ground-field list.

## 5. Long path plus spectator bridge: both endpoints

The old-path cycle gives A=2x in O_K (also at x=0). Write
h=1 or N-1 according to the contacted endpoint, and put

    Z=NA+h,
    kappa=h/(2N),       t=1/(2N),
    R=Z^2/(2N),
    C=2Z^2-(N+2),       C*delta=4N.

The proper-block inequality is Z_i^2<1 at each nonphysical
embedding. In the far role h=1,

    C=ND,       D=2N*A^2+4A-1,       D*delta=4.

D is congruent to-1 modulo2; since it divides4, it is a unit.
At x=0 it equals-1 and gives delta=-4, impossible. At x>0
the native gap gives A_0>=1 and D_0>=2N+3>1.

In the near role set U=A+1, so Z=NU-1 and

    C=ND,       D=2N*U^2-4U-1,       D*delta=4.

Again D is odd and therefore a unit. Now U_0>=1 gives
D_0>=2N-5>1, including at x=0. At either endpoint,

    C_i<-N,       D_i<-1

by Z_i^2<1. The product of the absolute conjugates exceeds1,
contradicting that D is a unit. Both bridge roles are excluded
over every totally real K, without a field-generation assumption.

## 6. Audit of the complete isolated-d row over the eight fields

For clarity, the composition uses all values of x, all white
contacts of c, and the entire field set

    Q(sqrt3), Q(sqrt6), K9, K11, K15, K16, K20, K24.

The complete branch table is:

|c's white contacts|x range|source settling the branch|
|---|---|---|
|no spectator, no long arm; e contact absent or native|zero, finite native, dotted|no-long source below|
|no spectator, a long contact (with all allowable e contacts)|zero, finite native|finite-long source below|
|same long-contact inventory|dotted x>=1|dotted-long unit closure below|
|spectator only; no e or path|zero, finite native, dotted|continuation source Section3|
|spectator plus e3; no path|all x|Section4 here|
|spectator plus path endpoint; no e|all x|Section5 here|

The first three named sources are

* `arithmetic_rank20_isolated_side_no_long_contact_20261019.md`, SHA256
  `af2ff3c9c4667b037ea93ed1afa8f1dcaf364c0ac41e894e3541b2bcde1df879`;
* `arithmetic_rank20_isolated_long_finite_product_20261019.md`, SHA256
  `9644c76993b26d97803686787c5a87c284051184aa08cd8d2d7d6250a06e1ece`;
* `arithmetic_rank20_isolated_d_dotted_long_unit_closure_20261019.md`, SHA256
  `5c583ce01db154908d5d0e139c61f4986cf1491cd4f3f991f0c7e26c8a5a92c6`.

The continuation source's Section3, bound in Section2 above, was
reread in full for this audit. It excludes x=0 by the physically
independent old-plus-c block, finite native x by the small positive
norm of w=4x^2 and the exclusion of K=Q, and dotted x by the
coupled w,z norm argument. Its separate degree2,N=8 case is kept.
That branch does NOT incorrectly assume2x in K: it uses the
actual integral squares w=4x^2,z=4beta^2.

For every branch beta>1 follows from the actual common-perpendicular
c side as proved in Section1. No unproved ambient-opposites
dottedness premise remains for b,d. What remains conditional is the
explicit actual-face framework, arithmetic cyclic integrality and
strict nonphysical positivity, and acceptance of the named complete
prior branch proofs; these are not replaced by numerical evidence.

The finite inventory in Section2 excludes all remaining mixtures:
multiple contacts into an old connected factor create a cycle,
and e+path+spectator creates a second fork. For c orthogonal to
the spectator the complete long inventory retains the two N=8
E8 cases. Thus the table contains every c white-support possibility
in the isolated-d row, without omitting zero or finite x.

Subject to those stated interfaces and composed branch proofs,
the isolated-d row is now excluded over ALL eight fields. This
does not exclude the distinct first-d rows F/B/T, does not establish
all triangle or simultaneous modified-angle consumers, and does
not install a global half refund or a lower dimension ceiling.
