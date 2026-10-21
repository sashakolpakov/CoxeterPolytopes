# The second initial edge: 29 finite profiles reduce to 18

19 October 2026. Exact paper analysis, independent of a root search.
No scientific computation, enumeration, optimizer, Lean job or build
ran. The earlier pinned pencil and omega2-profile reports are unchanged.

**Result.** For a canonical short root closing the initial omega3 edge,
the complete Schur-positive finite white-contact list has 29 labeled
profiles, including the orthogonal profile. The prime-D test leaves
22, and a native divisibility condition at sqrt(5) leaves 18. All 18
have only simple finite bonds. This does not construct their roots
or rule out their simultaneous native and chamber constraints.

Long roots cannot close the edge: the absence of finite mixed bonds
forces white orthogonality, and its entire long pencil is obstructed
modulo sqrt(5). The entire orthogonal short pencil is also obstructed
modulo D. These two inherited arguments are proved in the unchanged
`compact_h9_discriminant19_edge_pencil_obstructions_20261019.md` and
`compact_h9_discriminant19_edge_profile_congruences_20261019.md`.

## 1. White support and the complete finite-profile list

Use O=Z[phi], phi²=phi+1, D=1+5phi, and the native form

    G=[[-2phi,phi],[phi,2]] + C_E8,
    E8 edges01,12,23,34,45,56,27.

The omega3 edge frees old4=-e3. Its eight incident white walls are

    W=A4(E8 chain0--1--2--7)+A3(E8 chain4--5--6)
      +A1(binary old0).

For a canonical short closing root, every finite absolute white
contact is0,1,phi. The preceding pencil report proves this full
contact alphabet, rather than imposing a label cutoff. Let beta be
the binary A1 contact. Within A4 write e4 for an endpoint contact
(node0 or7) and i4 for an internal contact (node1 or2). Within A3
write e3 for an endpoint contact (node4 or6) and m3 for its middle
contact (node5). All absent contacts are0.

The enlarged white Gram is positive definite exactly when

    qW=eta_W^T*W^{-1}*eta_W < 2.

At most one node per component can be contacted. Indeed two unit
contacts in A3 contribute at least2 (the two endpoints attain2),
and the same holds in A4; its inverse has nonnegative entries, so
larger contacts cannot help. Single unit-contact contributions are

    A4 endpoint:4/5, internal:6/5;
    A3 endpoint:3/4, middle:1;
    A1:1/2.

Thus the profiles with only contacts1 are:

* no contact:1;
* one contacted component:4+3+1=8;
* two contacted components: A4+A3 gives10, A4+A1 gives4,
  and A3+A1 gives3, totaling17.

In A4+A3 the only rejected pairs are the two A4 internal choices
with the A3 middle node, whose sum is11/5. Three contacted components
are impossible: their minimum sum is4/5+3/4+1/2=41/20>2.

A contact phi multiplies its contribution by phi². None on A4 is
allowed, since (4/5)phi²>2. An A3 endpoint contact phi is allowed
alone, in two labeled positions; its middle contact is not, and
even adding A1 to an endpoint makes the sum exceed2. A binary A1
contact phi is allowed alone, but adding the least expensive A3
or A4 contact makes the sum exceed2. Two contacts phi are also
excluded. These assertions follow from3/2<phi<5/3.

The exact complete count is therefore26 zero/1 profiles plus3
profiles involving phi, for29. No finite Coxeter diagram
classification or numerical census is used in this count.

## 2. Exact native norm and completed square

For q=(x,y,z), let eta_j be the absolute product with -e_j, j!=3.
Then

    phi*x+2y=beta,
    z=T*omega3+sum_{j!=3} eta_j*omega_j,
    omega_j=C^{-1}e_j,
    omega3=(8,16,24,20,15,10,5,12), omega3²=20.

All x,T,eta_j are in O. In addition to any scalar equation, native
integrality requires beta-phi*x in2O. Define

    R=sum_{j!=3}(C^{-1})_3j*eta_j,
    Q0=sum_{i,j!=3}eta_i*(C^{-1})_ij*eta_j.

The exact norm-two equation is

    beta²+2Q0+4TR+40T²-Dx²=4.

Multiply by10 and complete the square:

    A=20T+R,
    A²-10D*x²=K,
    K=40-10beta²-20Q0+R².                          (1)

The block inverse identity at node3 gives

    K=20*(2-qW)
     =40-16e4²-24i4²-15e3²-20m3²-10beta².         (2)

The scalar equation is necessary, with the additional native
congruence A=R modulo20O retained. Neither that congruence nor
the binary parity condition is replaced by (1).

## 3. Complete prime-D test: 29 profiles reduce to22

Modulo the prime ideal(D), the residue field is F19 and phi=15.
Equation(1) requires K to be a square. Its squares including0 are

    {0,1,4,5,6,7,9,11,16,17}.

The following table contains every profile. A label3 denotes
contact1, and5 denotes contactphi. All endpoint and internal
positions and the component labels remain distinguished.

|Contacted components|Multiplicity|Exact K|K modulo19|Prime-D result|
|---|---:|---:|---:|---|
|None|1|40|2|Excluded|
|A4 endpoint(3) alone|2|24|5|Retained|
|A4 internal(3) alone|2|16|16|Retained|
|A3 endpoint(3) alone|2|25|6|Retained|
|A3 middle(3) alone|1|20|1|Retained|
|Binary A1(3) alone|1|30|11|Retained|
|A4 endpoint(3)+A3 endpoint(3)|4|9|9|Retained|
|A4 internal(3)+A3 endpoint(3)|4|1|1|Retained|
|A4 endpoint(3)+A3 middle(3)|2|4|4|Retained|
|A4 endpoint(3)+binary A1(3)|2|14|14|Excluded|
|A4 internal(3)+binary A1(3)|2|6|6|Retained|
|A3 endpoint(3)+binary A1(3)|2|15|15|Excluded|
|A3 middle(3)+binary A1(3)|1|10|10|Excluded|
|A3 endpoint(5) alone|2|25-15phi|9|Retained|
|Binary A1(5) alone|1|30-10phi|13|Excluded|

The exclusions have multiplicities1+2+2+1+1=7, leaving22.
In particular the orthogonal case is obstructed at every height,
consistent with its separately proved pure-pencil argument.

## 4. Native ramified-prime condition: 22 reduce to18

Let pi=2phi-1=sqrt(5). The residue field O/(pi) is F5, where
phi=3 and D=1. Unlike the bare requirement that K be a square,
we now use the native congruence A=20T+R.

In the four groups involving only A3 or only the binary A1, the
omega3-row entries are15,10,5 on A3 and0 on the binary root.
Hence R belongs to5O, and A also belongs to5O. Write A=5B.
Every corresponding K is in5O; dividing (1) by5 gives

    5B²-2D*x²=K/5.

Reduction modulo pi therefore requires

    -2x²=K/5 in F5.                               (3)

This is stronger than reducing (1) modulo pi before division.

* A3 middle(3) alone has K/5=4, so x²=3, impossible.
* Binary A1(3) alone has K/5=6=1, so x²=2, impossible.
* A3 endpoint(5) alone has K/5=5-3phi=1 modulo pi,
  so x²=2, impossible at either labeled endpoint.
* A3 endpoint(3) alone has K/5=5=0 and requires x=0
  modulo pi; this is not a contradiction and remains open.

Thus exactly four more profiles are excluded by these arguments,
leaving18. For the other retained rows R is nonzero modulo pi;
the above divisibility-by5 argument does not apply to them.

The complete retained list is

|Contacted components|Multiplicity|Exact K|
|---|---:|---:|
|A4 endpoint(3) alone|2|24|
|A4 internal(3) alone|2|16|
|A3 endpoint(3) alone|2|25|
|A4 endpoint(3)+A3 endpoint(3)|4|9|
|A4 internal(3)+A3 endpoint(3)|4|1|
|A4 endpoint(3)+A3 middle(3)|2|4|
|A4 internal(3)+binary A1(3)|2|6|

The multiplicities sum to18. All remaining contacts are1;
no label5 contact survives these two local filters.

## 5. An exact factorization for two of the retained profiles

For an A3 endpoint(3) alone, beta=0 implies x=2u with u in O
and y=-phi*u. At node6 the endpoint has R=5 and Q0=2, so the
uncompleted norm equation reduces to

    5T(2T+1)=D*u².

At node4 it instead has R=15 and Q0=12, giving

    5(2T+1)(T+1)=D*u².

These are exact all-variable identities, not evidence that their
positive solutions exist or pass the chamber. They offer a smaller
factor or relative-unit route if a later bounded search is authorized.

## Scope

The reduction29->22->18 is a necessary finite-contact and arithmetic
analysis at the actual initial edge. No candidate root was enumerated
or constructed here. Remaining work includes A=R modulo20O,
beta-phi*x in2O, nonnegative native E8 covector coordinates, all
retained wall products and their Coxeter-angle constraints, and
actual finite endpoint position. No all-height exclusion of the18,
finite-volume completion, new dimension bound or nonreflectivity
claim follows from this packet.
