# Two bounded-profile roots and the smaller27th H19 wall

19 October 2026. The coordinator executed the prepared finite-profile
search once, then replayed it separately. Both passed the shared-lock
resource wrapper. The source and preparation report remain unchanged.
The smaller compatible root gives a27-wall Coxeter chamber, still of
infinite volume by the exact spacelike witness below. A dedicated fixed
27-wall checker is being prepared separately and has not yet run.

## Executed finite search

Source `scripts/search_h19_omega2_finite_profiles_20261019.py`, SHA256
`02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e`,
checked405 component choices,51 nonzero Schur profiles,40 after the
original-equation modulo3 filter, and34 after the prime-D filter.
Its1,806 conservative coefficient pairs contained1,073 x in the
requested domain0<x<=600, |x'|<=2. There were36,482 profile/x exact
field-square tests,24 field-square solutions and48 signed branches.
Forty-four failed T/y integrality; four were integral norm-two roots.
Two passed every retained wall and the actual old3 endpoint test.

Generation/replay ran in0.650/0.654 seconds with process-tree peaks
57,065,472/56,836,096 bytes. Start resources were186.65GiB disk and
7.51/7.81GiB available RAM. The wrapper owned the single lock,59-second
limit and1GiB RSS cap; the child added its45-second alarm. This lane
did not run a second search or any checker body.

## The smaller root and its actual endpoint

The profile44 survivor has x=1+3phi. Write it q; its coordinates are

    q=(1+3phi,-1-2phi;
       7+10phi,14+20phi,20+30phi,16+24phi,
       12+18phi,8+12phi,4+6phi,10+15phi).

Its sparse covector is

    Gq=(-8-11phi,1;0,1,phi,0,0,0,0,0),

and q^2=2. Integral coordinates and this integral covector prove full
lattice reflection integrality. The contact profile has beta=1 and
eta_1=1, all other prescribed eta zero. Thus q_W=7/6 and K=25.
The free coefficient is T=phi, A=20+30phi and

    (p,q)=-16-23phi,
    (c,q)=-(68+95phi)/10<0.

Here c is the already fixed strict centre p+(0,1/10;C^(-1)1/10).
The executed certificate checks all26 pair products: two angle3,
three angle5, nine right and twelve strictly ultraparallel. In
particular q is strictly ultraparallel to both full compact H8 facet
walls r,w and has negative product with the old finite vertices p,p_s.

The actual new old3 endpoint is

    p_q=p+tau*omega2,
    tau=(16+23phi)/(20+30phi)=(11-2phi)/10,
    p_q^2=71/2-22phi<0.

The search certificate checks the future time cone, all26 old
inequalities, the q equality, and exactly the eight white old walls
active. The enlarged white Gram is positive definite, of type
2A4+A1. Thus p_q is an actual finite vertex, distinct from p,p_s.

The two inherited full H8 facet vertex sets have41 vertices each and
are disjoint. Appending an acute wall strictly ultraparallel to a
compact facet preserves its entire hyperplane: the projected strict
centre gives one negative point on that hyperplane, and a disjoint
hyperplane cannot change sign on its future hyperbolic sheet. Hence
both sets survive q. Together with p,p_s,p_q this gives at least85
actual finite vertices, using the inherited H8 completeness only for
the two41-vertex counts. This is not a full27-wall vertex census.

## The other survivor cannot be appended together with q

The profile22 survivor has x=44+72phi and covector

    Gq_L=(-202-326phi,0;0,1,16+28phi,0,0,0,0,1).

Its binary coordinates are44+72phi,-36-58phi; its E8 coordinates
at nodes1,2 are344+560phi and515+840phi. Taking their product
against the displayed sparse Gq gives, entirely by phi^2=phi+1,

    (q_L,q)
     =(-8-11phi)(44+72phi)+(-36-58phi)
       +(344+560phi)+phi(515+840phi)
     =4+5phi>0.

Therefore the two inward wall orientations cannot coexist in an
acute Coxeter extension of the current chamber. The27-wall choice
here appends only q; q_L remains a distinct bounded-search outcome.

## An unchanged spacelike escape still proves infinite volume

The fixed26-wall old4 first-hit point is

    E=p+((8+21phi)/44)*omega3,
    E^2=(1557-955phi)/484>0.

All its26 old products were checked in the saved26-wall certificate.
Since the q spatial coordinate at node3 is16+24phi,

    (E,q)=-16-23phi+((8+21phi)/44)(16+24phi)
         =(5phi-18)/11<0.

Thus E is also feasible for all27 walls. Its product with the future
timelike p is p^2<0. Perturbing E slightly toward the common strict
centre c makes every inequality strict while retaining positive
square and this time orientation. The segment from c to that point
crosses the future null cone in strict wall inequalities; consequently
the chamber contains an open patch of ideal boundary and has infinite
hyperbolic volume. This is a written geometric consequence of the
explicit feasible vector, not a final-cone census. The local checker
will independently reproduce its complete27-product ledger.

No finite-volume example, complete root enumeration, maximal reflector
system, whole-lattice nonreflectivity or dimension-bound result follows.

## Evidence

The following are immutable executed outputs; hashes are bound here.

|Evidence|SHA256|
|---|---|
|`results/h19_omega2_finite_profiles_20261019.json`|`46ab5c5d83415941362c8ab22df5a42d7c96ed595069f01a8f5e3597ae445813`|
|`results/h19_omega2_finite_profiles_generation_20261019.json`|`c6a104f51c2a1e2acf6603bfac813d821c9c8a7b4e934a565083ff20fe29f24d`|
|`results/h19_omega2_finite_profiles_replay_20261019.json`|`29657d7fbb6a485460d9a5cfcca98a3b0bf7418b2358e1f67648c319cd707c46`|

The complete contact-profile proof, exact arithmetic, coefficient-bound
coverage and inherited norm-alphabet scope are in
`results/h19_omega2_finite_profiles_preparation_20261019.md` and its
two pinned profile reports. Only the stated bounded root search has
been executed; the all-height pure-pencil exclusions are separate
paper congruence theorems.
