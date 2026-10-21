# Norm239: the next actual wall and its E7+A2 germ

22 October 2026. **Complete supporting proof candidate, pending root and
outside-author review of this composition.** The exact finite algorithm
already has separate root and outside-author source reviews, and all20
current units and their combiner now have actual guarded terminal PASS.
The author prepared the folding contract and worker; this note is not an
independent review of those foundations. Root owned every mathematical
execution. Preparation of this note used source reads, manual algebra and
small receipt/hash checks only; no worker was imported or rerun.

## 1. Literal module, inherited inputs and theorem

Let O=Z[phi], phi²=phi+1, and use the ORIGINAL module

    L=E8(O) perpendicular-sum [[-2phi,1],[1,8]].

The binary basis is X,Y. Native roots e0,...,e7 have norm2 and products
-1 on02,13,23,34,45,56,67, with all other off-diagonal products0. Put

    D=1+16phi, b=phi D, A=phi², E=7-4phi, F=3-phi=A E,
    q=(21+33phi)X-(17+25phi)Y,
    m=(3+7phi)X-(3+5phi)Y,
    H=Omega7=(2,3,4,6,5,4,3,2),
    J=Omega6=(4,6,8,12,10,8,6,3),
    R=phi(q-m)/2+phi³H,
    Q=2Aq+DH, W=m+(phi/2)R.

Thus q²=-2D, m²=R²=2, m.R=-phi, H²=2, J²=6 and H.J=3.
The inherited ACTUAL full germ at future Q is E7+I2(5), with simple
walls -e0,...,-e6,m,R. Release m and retain E7 and R. The outgoing
line is P(z)=Q-zW, where Q²=-2DE, W²=F/2 and Q.W=0.

The exact inputs are:

| File | SHA-256 | Supplied boundary |
|---|---|---|
| `results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md` | `879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02` | Original module, dyadic glue, actual base and binary edge |
| `results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md` | `2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15` | Actual full starting germ at Q |
| `results/compact_h9_prime239_second_h2_edge_finite_upper_contract_candidate_20261021.md` | `1d11ee290abb97a9da0aaab1d86496071de46d2cf74a02489457fb7f369f89b4` | Original-integral finite upper and quarter-integral folding isometry V |
| `results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md` | `b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce` | Primitive norm alphabet2/2b and full arrangement/local-finiteness input |
| `results/compact_h9_prime239_folded_second_h2_edge_checker_spec_20261022.md` | `4e070627629fd0e56aa3dff34de6cdd1f424c25a3677271dd9a6f2f7b09a5a29` | Complete finite domain and exact algorithm |
| `scripts/check_compact_h9_prime239_second_h2_edge_batches_20261022.py` | `6569c6083ce781400fa2052f5a863af526605e4d27e367b4f172f45a6d28c9e2` | Executed implementation of that domain |

Finite spherical Coxeter classification and recognition of the local
simple cone remain explicit geometric inputs. No compactness, finite
covolume or chamber-stabilizing property of V is assumed.

**Theorem at these inputs and the exact computed proposition below.**
The first new simple wall of this actual outgoing edge is

    S=(55+90phi)X-(43+69phi)Y+(6+10phi)H.             (1)

Its collision is z*=(664+108phi)/131. The endpoint is finite and has
the full ACTUAL simple germ E7+A2, with simple walls
`-e0,...,-e6,R,S`. This is one next-edge theorem for the literal module;
it is not a chamber completion, volume result, example or host exclusion.

## 2. The exact computed proposition and why its domain is complete

The folded contract proves an original-integral finite mirror upper
V³m at distance3l0, where V fixes E7,R and has active trace
t=563+914phi=2cosh(l0). The matrices V^j and V^-j, j0,...,3, have
entries in(1/4)O on the ORIGINAL basis; V^+-3 are integral. No power
is assumed to preserve the chamber. Every actual first simple wall
before the upper has a nearest phase j0,...,3 and signed folded distance
of absolute value at most l0/2.

Finite spherical extension of the retained E7+R and the primitive norm
alphabet give exactly five necessary profiles. A folded seed has

    s=(Tq-Km)/2+hH+epsilon J, s.R=-kappa,
    norm2: (epsilon,kappa)=(0,0),(0,1),(0,phi),(1,0),
    norm2b: (epsilon,kappa)=(0,0).

Both signs of h must be kept. For d1 or b define

    Lh=2A(2h+3epsilon)+2kappa/phi,
    M0=4h²+12epsilon h+12epsilon²-4d,
    Z=D[Lh²-(D-1)M0].

The contact and norm equations give K=DT-Lh and
`T=(D Lh +/-sqrt(Z))/(D(D-1))`. These are at most two EXACT field
roots per scalar row. In Q(sqrt5), writing Z=a+c sqrt5 reduces square
membership to both signs of sqrt(a²-5c²), then rational-square tests;
the c0, zero and pure-sqrt5 cases are separately retained and every
returned root is squared back. Thus a square norm is not mistaken for
field-square membership.

The contract's Cauchy and definite-conjugate estimates give |h|<342,
|h'|<=2 in norm2. For norm2b they give h/b in(1/4)O, |h/b|<53
and |(h/b)'|<43/50. Writing4h=a+c phi, or4h/b=a+c phi, gives
respectively |c|<=625 with at most17 a per c, or |c|<=97 with at
most7 a per c. Hence there are at most345732 scalar rows across all
twenty phase/profile units. These finite loops and at most two field
roots per row prove termination; the recorded runtime is separate.

Each seed is checked for the complete quarter-module condition and
the signed centered window. After reconstruction r=V^j s, the worker
checks ORIGINAL integral coordinates, exact norm, every reflection
matrix coefficient, positive collision at or before the finite upper,
and acuteness to all known actual walls -e0,...,-e7,r0,m,R, where
r0=X+2phi Y. These tests overlist mirrors but cannot discard an actual
first simple wall. All comparisons use exact ordered-field arithmetic.
Each unit retains its least collision and every canonically distinct
minimum tie. The bounded tie limit causes an error, never truncation.
The combiner requires every distinct phase/profile under identical
current bindings, rechecks each recorded original root, and compares all
unit minima exactly. An incomplete unit set earns no edge theorem.

**Actual finite output proposition.** Every one of those twenty units
ran to terminal PASS under the unchanged guard. The actual scalar total
is325636: each of the sixteen norm2 units has20017 rows and each of
the four norm2b units has1341 rows. The numbers of admitted rows are:

| Phase | Orthogonal | Simple | Golden | E8 | Norm2b |
|---:|---:|---:|---:|---:|---:|
|0|0|1|0|0|0|
|1|0|0|0|0|1|
|2|0|1|1|0|0|
|3|1|0|1|1|0|

The global exact minimum has ONE original root, namely(1), at z*.
The aggregate cone check reports one candidate and Gram E7+A2.
Later admitted mirrors in other units do not make this minimum tied.

The complete aggregate is
`results/prime239_folded_second_edge_complete_20261022.json`,
SHA-256 `6559027bf88d0ac040fac1d860f85fc408bbe401a5674f777da31efd048dea96`.
The actual execution conjunction is bound by
`results/prime239_folded_second_edge_runtime_manifest_20261022.json`,
SHA-256 `e8a64ac748aa6b85681bab9bea6cc24fbcf52e8fd6372190dab2e06e7468655f`.
It records all20 outputs, all20 terminal runner ledgers, their preserved
stdout, and the separate bounded combiner. Metadata inspection for this
note checked all71 structured binding occurrences/70 distinct paths,
all actual commands, output digests and all21 preserved stdout records.
It did not rerun the mathematical loops.

Units used74.703s in total, with maximum6.579s per unit and maximum
57950208B process-tree RSS. The separate combine used0.298s. Each
had59s/1GiB/one-thread limits; minimum recorded unit admission was
134.65GiB disk and6.23GiB available RAM. These are actual recorded
runs, not a runtime guarantee derived from the domain bound.

The algorithm's full outside-author source review is
`results/review_prime239_folded_second_h2_worker_and_contract_full_scope_compact_reviewer_20261022.md`,
SHA-256 `6657f95985c59311f8807191191218849da7f386b20e80467c6aa6e063f8da58`.
The root source review is
`results/review_prime239_folded_second_edge_contract_and_worker_root_20261022.md`,
SHA-256 `dc11c1585be150c25438e97918d68e2c7a30eeed6023a69d8145afceeb9c3887`.
They are pre-execution source reviews; the runtime manifest supplies the
subsequent actual execution evidence. Neither is a review of this new
mathematical composition.

## 3. Direct original-module verification of S

Set x=55+90phi, y=-43-69phi and h=6+10phi. Its ORIGINAL Gram
covector is

    GS=(-223-359phi,-289-462phi,0,0,0,0,0,0,0,h).

The binary contribution to S² is -270-440phi, whereas the E8
contribution is2h²=272+440phi. Hence S²=2. All its coordinates
and covector entries lie in O. Thus `v -> v-(v.S)S` preserves the
ORIGINAL lattice. Also S.R=-1 is a unit, so S is primitive. This
verifies an actual arrangement mirror without replacing L by the
quarter-module or a rationally similar lattice.

Direct pairings give

    S.(-e_j)=0 (j0,...,6), S.(-e7)=-6-10phi,
    S.R=-1, S.m=-5-10phi=-5phi³,
    S.r0=-1147-1861phi,
    S.q=-67-116phi,
    S.Q=-34-66phi=-2D A,
    S.W=-5-21phi/2.

These display every known-wall acuteness check, the correct orientation
S.Q<0, and the positive outgoing derivative -S.W>0. In particular

    -S.Q/(-S.W)=(68+132phi)/(10+21phi)
               =(664+108phi)/131=z*.

## 4. Firstness in the full chamber and the full finite germ

Let z_first be the actual first boundary on this outgoing chamber edge.
Local finiteness and the original-integral finite upper guarantee that
it exists before that upper, at a finite point. A first new SIMPLE facet
is in one of the complete phase/profile domains. It passes every
original-integrality, orientation and known-wall test. Therefore the
computed minimum is at most z_first. Conversely every admitted original
reflector is a mirror of the full arrangement. Such a positive mirror
cannot cross the interior of the chamber segment before its first
boundary. Thus every admitted collision is at least z_first. The two
inequalities prove z_first=z*, and uniqueness of the computed root
identifies S as the actual first new simple wall.

An explicit endpoint ray is

    P*=1/131 [(2206+2368phi)X-(1370+2030phi)Y
                    +(-695+830phi)H].              (2)

This equals Q-z*W. It is orthogonal to -e0,...,-e6,R,S.
Their Gram is the direct sum of the positive definite E7 Gram and
`[[2,-1],[-1,2]]`; its determinant is2*3=6. Thus these nine roots
are independent and their one-dimensional orthogonal complement in
the original Lorentz space is negative. Since P*.Q=Q²<0, (2) is
the future timelike ray of that complement. As a scalar cross-check,

    P*²=30(85255-52712phi)/17161<0.

At (2) the released walls remain strict:

    P*.m=-z*F/2<0, P*.(-e7)=(695-830phi)/131<0,
    P*.r0=D*(-1370-2030phi)/131<0.

Firstness proves feasibility against EVERY actual chamber wall by
closedness of the chamber; it is not inferred merely from these known
wall checks. Finite Coxeter chamber germs are simplicial. The eight
retained actual walls and the unique actual new simple wall therefore
give the FULL germ E7+A2 at this vertex. The independent exact cone
calculation agrees: one positive-definite nine-wall basis and one
minimum ray, with no unresolved tie or row-order selection.

## 5. Reproduction and remaining scope

The contract Section6 gives the guarded per-unit command. The manifest
binds every actual argv and output; reproduction requires all twenty
phase/profile commands through the SAME bounded runner, using NEW
destinations and fresh unchanged40GiB/6GiB admission, followed by a
separate bounded combine of all twenty NEW part files. The59s/1GiB,
one-thread watchdog and8GiB/3GiB abort floors remain mandatory. Nothing
in this proof authorizes rerunning already successful units or lowering
those gates.

The next edges leaving this E7+A2 vertex remain open. Neither the
integrality of a folding power nor the finite local diagram supplies a
full-chamber period. This note proves no finite volume, compact H9
example, nonreflectivity, form-class exclusion or dimension improvement.
The six delivered papers and all frozen sources/evidence are unchanged.
