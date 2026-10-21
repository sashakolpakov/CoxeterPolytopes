# D20 centre own3/helper4: the other fork contact must be dotted

21 October 2026. New manual geometric lemma, pending full independent review. No arithmeticity, field restriction, census, scientific/import, Lean or TeX job. No frozen manuscript, guide or earlier proof was edited. This is a restricted necessary condition, not a whole-D20 exclusion or a general dimension improvement.

## Inputs and statement

Inherited source pins, under `results/`:

|Source|SHA-256|
|---|---|
|`compact_whole_d20_complete_continuation_20261021.md`|`9cec3f696bd5e67fba56de98db45d433188d5b13ee88eb1be5be788a8d0fe3f3`|
|`compact_whole_d20_long_end_reduction_20261021.md`|`58a158998d06f0384aae4e0bbf4a0069d05900ce6c8ccaa6937dc3164708fc76`|
|`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`|`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`|

The inherited sources supply the actual corner/normal convention, ordinary finite spherical classification, and the isolated native5 head interface. I authored the first restricted continuation and used root's observation of the remaining centre wall; the new contribution below is the complete finite/parallel own-contact exclusion for the other fork. It uses no numerical floor or conclusion from the inherited tail theorem.

Let p be a vertex of a compact Coxeter20-polytope, with whole old corner diagram D20. The farthest-vertex hypotheses of the inherited continuation are allowed but are not needed for this lemma. Work with unit old normals

    n_i=(e_i-e_(i+1))/sqrt2 (1<=i<=19),
    n20=(e19+e20)/sqrt2,
    omega_i=sqrt2*P_i (i<=18),
    omega19=(P19-e20)/sqrt2, omega20=(P19+e20)/sqrt2.

Here P_i is the sum of the first i orthonormal spatial coordinate vectors at p; p has Lorentz square -1 and is orthogonal to their space. Suppose two ACTUAL opposite walls are:

1. the isolated native5 head `h=A0*p-g*omega1`, where `g=cos(pi/5)` and `A0^2=2g^2-1`;
2. the centre18 wall `c=A*p-omega18/2-omega19/sqrt2`, with own3 at n18 and helper4 at n19.

**Lemma.** The actual outgoing-fork20 opposite wall d has dotted own contact with n20: `gamma=-(d,n20)>1`. Its retained A19 white contacts are exactly one of the seven rows in the table below. In particular every finite own angle, a right own angle, and a parallel own contact are impossible in every row.

The wall c is allowed to be the fork19 opposite wall. No argument assumes those two walls distinct. The other fork20 opposite cannot be c: its white diagram retains both n18 and n19, while c,n18,n19 form an indefinite 3–3–4 triangle. Equivalently the one-helper finite white list below cannot contain c's two nonzero retained contacts.

## 1. Complete white inventory and common algebra

Deleting n20 leaves the long A19 path n1,...,n19. Compactness makes its replacement by the actual d a finite spherical rank20 diagram. Two retained contacts would give a cycle. A finite extension of an A19 path is therefore: an isolated A1; a simple3 attachment at an endpoint or next-to-endpoint (A20 or D20); or a double4 attachment at an endpoint (B20). Exceptional E/F/H diagrams cannot contain this rank19 path. This proves the complete white list, without a search or an arithmetic assumption:

|White contacts of d|Spatial part v_d|`T^2`, where `d=T*p-v_d`|
|---|---|---|
|none|`gamma*omega20`|`10gamma^2-1`|
|simple3 at n1|`gamma*omega20+omega1/2`|`10gamma^2+gamma-1/2`|
|simple3 at n2|`gamma*omega20+omega2/2`|`10gamma^2+2gamma`|
|simple3 at n18|`gamma*omega20+omega18/2`|`10gamma^2+18gamma+8`|
|simple3 at n19|`gamma*omega20+omega19/2`|`10gamma^2+9gamma+3/2`|
|double4 at n1|`gamma*omega20+omega1/sqrt2`|`10gamma^2+sqrt2*gamma`|
|double4 at n19|`gamma*omega20+omega19/sqrt2`|`10gamma^2+9sqrt2*gamma+4`|

These are actual unit normals, so T>0. The own parameter is either gamma=0, gamma=cos(pi/m) for an integer m>=3, gamma=1, or gamma>1. All outward products of distinct actual walls are nonpositive; every strictly negative product of magnitude below 1/2 is forbidden for Coxeter walls.

Put

    s=sqrt2, a=(1+s)/2, t=1+s/2=s*a,
    v_c=a*P18+(e19-e20)/2,
    A^2=13+9s, B=v_c.omega20=9t,
    D=10A^2-B^2=(17+18s)/2>0.

For a helper of coefficient beta at j, define `K=v_c.v_d=B*gamma+beta*b_j`, where

    (b1,b2,b18,b19)=(t,2t,18t,9+5s).

The actual wall product is `(c,d)=K-A*T`. Since K,A,T are positive in every case to which the square comparison is applied, its sign is controlled without extraneous roots by

    P(gamma)=A^2*T^2-K^2.

If P<0, the product is positive, contradicting acuteness. If P>0 and `Q=P-K-1/4<0`, then `K<A*T<K+1/2`, so `-1/2<(c,d)<0`, contradicting the Coxeter angle spectrum.

## 2. Isolated and near-fork rows: acuteness already excludes gamma<=1

For the isolated row, `P=D*gamma^2-A^2`. Hence for 0<=gamma<=1,

    P<=D-A^2=-9/2<0.

If gamma=0 its time square is negative anyway. Thus no isolated finite or parallel contact is possible.

For the three near-fork helper rows the exact polynomials are

|Helper|P(gamma)|P(1)|
|---|---|---|
|simple3 at n18|`D*gamma^2-9gamma-A^2-9/2`|`-18`|
|simple3 at n19|`D*gamma^2-B*gamma-A^2-1/4`|`-B-19/4`|
|double4 at n19|`D*gamma^2-s*B*gamma-A^2-1/2`|`-14-9s`|

Each has P(0)<0 and is convex. Thus P<0 throughout [0,1]. Here K remains positive even at gamma=0. All finite, right and parallel own contacts in these rows contradict centre-wall acuteness.

## 3. Simple3 at n2: a strict forbidden contact interval

Here

    K=B*gamma+t,
    P=D*gamma^2-gamma-t^2.

The zero-own case has T=0 and cannot be an opposite wall. At every finite positive own contact or parallel contact, 1/2<=gamma<=1. We have

    P(1/2)=(1+10s)/8>0,
    P'(gamma)=2D*gamma-1>0 on[1/2,1].

Also

    Q=D*gamma^2-(B+1)*gamma-t^2-t-1/4,
    Q(0)<0,
    Q(1)=3s-17/4<0,

where the last strict inequality follows from `sqrt2<17/12` (288<289). Convexity gives Q<0 on[0,1]. Thus every 1/2<=gamma<=1 produces the forbidden strict interval `-1/2<(c,d)<0`. No limiting native3 or parallel case survives the strict endpoint tests.

## 4. Simple3 at n1: the native5 head and the finite-label gap

Here `T^2=10gamma^2+gamma-1/2`. Set `g=cos(pi/5)`, so

    4g^2=2g+1, 0<g<1.

The spatial product with the native5 head is `g*(gamma+1)`. Its squared difference from the time product is

    F(gamma)=g^2*(gamma+1)^2-A0^2*T^2
            =(10-19g^2)*gamma^2+gamma+g.

It is concave, with

    F(0)=g>0, F(g)=(1-g)/8>0.

Therefore F>0 on[0,g]. Whenever T>0 both sides compared are positive; the actual product `(h,d)=g*(gamma+1)-A0*T` is positive, impossible. A nonpositive time square is itself impossible. This excludes every native own label m<=5.

Any remaining finite native label has m>=6, and hence `3/4<=gamma^2<=1`; the parallel value gamma=1 is included in the following test. With `C=(55+38s)/8`,

    K=B*gamma+t/2,
    P=D*gamma^2-gamma/2-C
      >=3D/4-1/2-C=2s-1>0.

Meanwhile

    Q=D*gamma^2-(B+1/2)*gamma-C-t/2-1/4,
    Q(0)<0, Q(1)=-(69+4s)/8<0.

Convexity again gives Q<0 for 0<=gamma<=1. Thus all native m>=6 and the parallel endpoint give `-1/2<(c,d)<0`, forbidden. The interval between g and cos(pi/6) is not silently covered by a continuous inequality: it contains no finite Coxeter own parameter. This is exactly where the integer native-label spectrum is used.

## 5. Double4 at n1: head acuteness excludes the whole interval

Here `T^2=10gamma^2+s*gamma` and the head's spatial product is `g*(gamma+s)`. The exact squared difference is

    F4(gamma)=g^2*(gamma+s)^2-A0^2*T^2
             =(10-19g^2)*gamma^2+s*gamma+2g^2.

It is concave. Its first endpoint is positive, and

    F4(1)=10-17g^2+s > s-4/3>0,

using `g^2=(3+sqrt5)/8<2/3` and `sqrt2>4/3`. Therefore the actual head product is positive for every gamma in[0,1] with T>0. At gamma=0, T=0 instead. This excludes every finite and parallel own contact in the last far-helper row.

## 6. Exact remaining scope

The seven finite-white rows are exhaustive, and all gamma<=1 Coxeter possibilities have now been excluded. Consequently the actual other fork20 wall has gamma>1 in one of those same seven rows. The proof uses the SAME original corner p and actual native5 head and centre walls throughout. It neither adds norm lower bounds from different farthest vectors nor pays a prefix coefficient twice; no prefix reference is used at all.

**All dotted rows remain open in this note**, including the isolated one. For the isolated dotted row alone, centre acuteness gives the necessary bound

    gamma^2 >= A^2/D=(206+162sqrt2)/359>1.

Equality is orthogonality of c and d and is not excluded here. The other dotted cases retain their actual centre/head products and the complete white inventory above; no assertion that these local data extend to a compact polytope is made. This is a reduction of the particular centre18 own3/helper4 residual in the native5-head branch, not its full exclusion, not a complete centre18 floor, and not a whole-D20 or general19 theorem.
