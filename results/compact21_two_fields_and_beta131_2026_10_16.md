# Two remaining compact dimension-21 fields and a distinct H9 host

**New internally reviewed proof candidate:** compact arithmetic Coxeter
polytopes over eight additional fields have dimension at most 20. Together
with the earlier field-specific arguments, this leaves only
**Q(sqrt5) and Q(sqrt2,sqrt5)** unresolved in compact dimension 21.
The overall arithmetic bound-21 candidate is unchanged, and no new
finite-volume example has been found. External mathematical validation
and priority remain open.

The [preceding checkpoint](arithmetic_hub_counts_and_e8_obstruction_2026_10_16.md)
and its exact source/report hashes remain frozen.

## A shorter field obstruction

The [new theorem](fresh_arithmetic_high_label_20261016.md) applies whenever
the ground field K is not a subfield of F=Q(sqrt2,sqrt5) and does not
contain K7=Q(cos(2pi/7)). It is independent of the large Q(sqrt5) attachment
enumerations. The field hypothesis forces every minimal indefinite diagram
to have order two or three. Every triangle has either a bond of label at
least nine or a bond six/eight with nonright contacts at both endpoints.

These bonds are isolated in finite corners. If a quadrilateral has no
marked path of length at most two, its two opposite minimal indefinite
diagrams force a principal six-node configuration. Eliminating the two
opposite rank-two blocks leaves a negative-definite two-plane, impossible
in Lorentz signature. Distance-two angle weights, with an extra unit on
isolated high-label dihedral components, then give vertex total at most
2n-2. At dimension 21, the h-vector requires a strict average greater than
40, contradicting the vertex maximum 40.

The [independent review](fresh_arithmetic_high_label_20261016_review.md)
uses separate cyclotomic automorphism witnesses for both field
noncontainments, a complete 35,904-graph enumeration, exact scalar bounds
and all-rank counting identities. Its nine surviving graph supports agree
with the primary's differently organized 43,836-case enumeration. Both
implementations pass generation and fresh saved replay. The argument needs
no finite angle cutoff, full-star rank reduction or spectator catalogue.

The eight new fields are Q(sqrt3), Q(sqrt6), Q(sqrt2,sqrt3), and
K9, K11, K15, K16, K20. The [complete branch composition](fresh_compact21_field_reduction_20261016.md)
accounts for thirteen published dimension-at-least-14 branches: Q is
excluded separately by rational isotropy; eight fields use the new theorem;
Q(sqrt2) and K7 use their earlier complete reviews; two fields remain open.
Generation and saved replay pass, binding 96 recursive manifests. A
separate read-only scope review confirms the old Q(sqrt2) coverage/weight
pairing and that no dimension-22 field restriction is applied at 21.

A [fresh hostile audit](fresh_arithmetic_dimension21_hostile_20261016.md)
also finds no gap in the existing bound-21 geometric/coverage argument.
It makes the corank-one bridge explicit: a k-node quadrilateral core has
a positive-definite (k-2)-node corner and a negative direction, so its
nullity is at most one. This supports the dimension reuse already made by
the complete queues. Internal review is not external validation.

The remaining F branch has no arithmetic Lanner diagrams of order at least
four, by the published arithmetic Lanner classification. Nevertheless,
triangles (2,5,8) and (4,5,5) defeat the new bond criterion. The Q(sqrt5)
rule currently excludes dimension 22 upward, not 21. These are concrete
next targets; compact dimension 21 has not been excluded.

## A distinct dimension-9 construction with a calibrated compact facet

The new lattice over Z[phi] is

\[
E_8\perp\begin{pmatrix}-2\phi&-1\\-1&6\end{pmatrix},
\qquad\phi=(1+\sqrt5)/2.
\]

Its determinant has norm -131. This separates it from the previously
excluded golden determinant-norm-11/55 hosts even up to similarity. The
complete primitive reflection norm alphabet is 2 and 2(12+13phi).

The [height-20 pilot](fresh_h9_beta131_pilot_20261016.md) has 52 controls
and twelve actual walls. An [independent prefix review](fresh_h9_beta131_pilot_20261016_review.md)
reconstructs all roots, while a [different complete cone enumeration](fresh_h9_beta131_geometry_20261016_review.md)
checks all 220 nine-wall supports and confirms 22 finite vertices and
28 hyperideal rays.

The [height-100 continuation](fresh_h9_beta131_h100_20261016.md) has
**114 controls, 41 actual walls and 45 finite vertices**. The primary
complete cone has no ideal rays and 1,891 hyperideal rays. An exact strict
escape proves the current chamber still has infinite volume.

The [independent height-100 checks](fresh_h9_beta131_h100_20261016_review.md)
reconstruct all controls and roots, verify the 45 finite-vertex witnesses,
and confirm that wall 8 carries **the literal compact Bugaenko H8 example**.
Its eleven source roots and form agree exactly with the archived calibration;
a separate exhaustive 165-support calculation gives all 41 finite vertices.
Every other one of the 29 walls has a strictly future timelike projection
and is redundant throughout that facet. Four current ambient finite
vertices lie outside the facet.

All 66 saved integral corner maps fix the same future timelike point, so
the group they generate is finite. This rules out an infinite-order
symmetry from those particular maps, not from the full lattice. The
independent checks verify each map but do not independently repeat the
complete ambient cone or map-search domain. Primary and independent
generation and saved replays pass within these stated scopes.

## What prevented an unrestricted improvement

The unrestricted finite-volume bound candidate remains **163**, and the
published unrestricted compact bound remains **29**. The
[new exact rank-162 configurations](fresh_nikulin_rank162_local_obstruction_20261016.md)
show why simply adding total rank and common heights to the existing local
estimates cannot prove 162: two full Lorentz Coxeter Grams satisfy the
endpoint budgets, common height one and all present pairings below 14,
with the critical pairing exactly 335/24=14-1/24. One has strictly elliptic
endpoint frames. Both have certified strict null escapes and infinite
volume; neither is a finite-volume counterexample.

A [changed-omission diagnostic](fresh_nikulin_center_omission_scout_20261016.md)
also exposes cheap contact rows missed by the proposed improvement. It
does not change the 163 theorem candidate. Additional geometric incidence
or a complete improvement valid for all newly allowed contacts is needed.

The [archive checker](../scripts/audit_compact21_beta131_checkpoint_20261016.py)
and [integrity certificate](compact21_two_fields_and_beta131_2026_10_16_audit.json)
preserve the preceding archive and bind this checkpoint. Archive integrity,
independent computation, mathematical review, and external validation
remain separate claims.
