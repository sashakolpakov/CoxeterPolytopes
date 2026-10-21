# Full-cusp deletion rigidity and curvature-financed ideal priorities

21 October2026. **New source/manual mathematical proof, pending independent
review.** No scientific, enumeration, optimizer, Lean or TeX job ran; no
manuscript, accepted input, rule or certificate was edited. The dimension
ceiling remains57. The results below use the inherited actual-face and
product-of-simplices cusp framework, not an asserted ambient completion
of a locally realizable Gram star.

## 1. The spectator forest is a constrained deletion

Let an actual ideal corner of a quadrilateral vary in a whole irreducible
affine factor F. Write its supporting white forest as

    W = (F minus its two marked roots) orthogonal-sum Z.

The actual factor-fiber theorem says more than Z being positive definite:
for EVERY other whole affine cusp factor G, Z contains exactly G minus
ONE node. There are no additional supporting components outside that full
cusp diagram. This follows directly from the simplex-factor description
of the actual link edge: each inert simplex factor contributes one of
its vertices, hence exactly one omitted facet/root.

The following elementary deletion classification therefore applies to
the whole components of Z. Write B1=A1, D2=2A1 and D3=A3; an empty term
is omitted. For a large classical affine factor, deleting a node gives:

| Full affine family | Possible finite components after one deletion |
| --- | --- |
| A cycle | the connected A path of the same affine rank |
| C path with two double ends | one B path, or two B paths |
| B single-fork/double-end tree | one B or D component, or D_i+B_j |
| D double-fork tree | one D component, or D_i+D_j |

This table is a graph proof, without a rank census. Deleting an interior
spine node separates the two decorated ends. Deleting a fork center
leaves its two isolated leaves, which are the D2 convention above.
Deleting a fork leaf removes that fork and leaves the other decorated
end. An end of a B path may be isolated, giving B1. In each case retain
every native double bond. In particular a double-bond component is not
mistaken for a simple A path.

**Deletion rigidity.** If every other full cusp factor is classical of
rank at least9, a pure A component of Z has rank1, rank3, or rank at
least9. The ranks at least9 arise only by deleting a node of a whole
affine A factor of that same rank. No I2(5) component occurs; B2=I2(4)
can occur. More generally, without any rank restriction, a WHOLE A_s
component of Z with s>=9 must come from a WHOLE affine A_s factor of
the actual cusp: exceptional affine factors have rank at most8, while
the B/C/D deletion list has no such component.

The smaller statement that Z has no A2 already holds when all other
full factors are classical and none is affine A2. Small A1/A3 components
from the D2/D3 degeneracies must be retained; excluding them would be
false.

## 2. Exact small-spectrum and coupon consequences

Let x,y be the native nonnegative, disjoint contact rows of the two
opposite sides on Z, with both separate augmented edge Grams strictly
positive definite. The inherited complete strict spectator theorem gives,
for D=(x-y)^t G_Z^-1(x-y)<1,

    D in {0} union {s/(s+1):s>=1}
                       union {2-sqrt2, 1-sqrt5/5}.

A rational nonzero value forces ONE whole A_s component and unit leaf
or adjacent contacts; all other components are uncontacted. The two
irrational cases force whole I2(4)/I2(5), respectively. Deletion rigidity
therefore sharpens this at a cusp with the large-classical complement to

    D in {0,1/2,3/4,2-sqrt2}
                          union {s/(s+1):s>=9}.             (1)

For a ONE-SIDED contact row of energy less than9/10, only 0,1/2,3/4
remain. The nonzero row is a unit endpoint contact on a whole A1 or A3.
This is a restriction on actual completed cusps, not on arbitrary
positive-definite spectator forests in a local star.

In the accepted opposite-fork proof for affine D9, the one-extreme
unmixed side has energy<=3/4; with two extremes, select an unmixed side
of energy<=7/8. Under the same full-cusp condition it is empty or whole
A1/A3. Its separate finite recipient is therefore only A1, A2 or A4,
instead of the unrestricted A1 through A8 list. The analogous affine
D10 unmixed energy3/2 argument selects a side of energy<=3/4 and gives
the same recipient list, instead of A1 through A4. The joined finite-D
alternatives and all mixed-component branches are unchanged. These
statements localize those particular coupon branches; they do not erase
the same finite roles if another proof needs them.

There is also an unconditional large-component consequence. The affine
D9 proof's mixed whole-A15 spectator branch MUST originate at a cusp
containing another whole affine A15 factor. This holds with arbitrary
other cusp factors, including small and exceptional factors. Its chosen
finite outputs are the restored A17 endpoint payment or, in the other
old attachment choices, the protected D25 short-leaf distance8 unit.

## 3. A fully financed priority at every affine-A15 cusp in ambient57

Use affine ranks summing to56, V=product(r_i+1), the actual product-link
curvature kappa=(5G-56V)/2, and the CURRENT MANUSCRIPT envelope c(r).
The complete ambient57 diagnostic and the accepted low8 concentration
argument give, for EVERY rank56 profile,

    m = kappa/V - sum c(r_i) >= -261/2662 > -1/10.       (2)

This is a normalized assertion. Concentration needs no positivity to
deduce (2); no multiplication of a negative terminal margin is used.
The terminal domain and exact minimum are inherited checked propositions.

Suppose the actual cusp has a whole affine A15 factor. Its order16 cycle
has exactly80 retained distance-at-most-five pairs, hence normalized
charge5, whereas c(15)=411/64. Keeping this A15 rule saves91/64 against
the all-family envelope. A single such factor suffices for what follows.

At this cusp assign weight ONE to EVERY marked pair in EVERY varying
factor of rank at most11. Leave all higher-rank factors, including that
A15 factor, unchanged. This is an intrinsic condition on the FULL actual
cusp, not a condition on a displayed proper affine subdiagram. The only
possible increases relative to c(r) are

| Varying rank r | Full pair cost divided by order | Increase over c(r) |
| --- | ---: | ---: |
| 1 through8 | r/2 | 0 |
| 9 | 9/2 | 1/16 |
| 10 | 5 | 2/11 |
| 11 | 11/2 | 1/3 |

Each positive increase is at most r/33. These varying factors consume
at most56-15=41 ranks, so their total normalized increase is at most
41/33. If C_new is the complete ideal charge after this contextual
restoration, factor-fiber multiplicities give exactly

    C_new/V <= sum c(r_i) - 91/64 + 41/33.

Combining with (2),

    (kappa-C_new)/V > 91/64 - 41/33 - 1/10
                    = 839/10560 > 0.                    (3)

Moreover V>=16*42=672: the other ranks sum to41, and
product(s_j+1)>=1+sum s_j. Thus the ABSOLUTE surplus in (3) is

    kappa-C_new > 5873/110 > 43/128.                    (4)

Multiple A15 factors only improve the saving. This estimate allows all
other factors and all their multiplicities, not only the terminal
profiles or the41 negative diagnostic rows.

For any simultaneous rule whose other factor charges are bounded by
the same manuscript c(r), this proves that the contextual restoration
retains a positive absolute ideal surplus at these cusps. In particular
the reviewed order11 composition already has uniform actual surplus
43/128 elsewhere, and its factor envelopes are no larger than c(r).
Consequently that uniform43/128 bound survives the contextual restoration.
This is a new financing argument; it does not rerun its2,270 convolutions.

## 4. What the financing removes, and the remaining finite obstruction

Give an actual corner with the new weight1 its own priority: any
quadrilateral containing it is paid directly. Triangles containing the
relevant changed opposite-fork pair are already excluded by the original
triangle argument; all other triangle payments only increase.

The WHOLE-A15 spectator branch in Section2 necessarily has such a
priority corner: its selected affine D9 has rank9, and its cusp has
the required other affine A15. Thus that entire branch needs neither
its finite A17 endpoint restoration nor its finite D25 leaf unit.
This is an actual source-financed replacement with arbitrary other
supporting roots. It is not a claim that a finite coupon can be deleted
merely because an isomorphic unused corner occurs in a saved packet.

To remove an ENTIRE finite A17 or D25 entry, the complete current
chosen-output interface must still show that no other producer requires
it, or provide a replacement for each remaining demand. The theorem
discharges exactly this D9 producer. No current finite weight was edited,
and no new finite maximum is asserted. In particular it does not remove
the pure finite D16/D17/D23/D24 leading allowances or prove364 as a
global finite budget. The new low-coupon restrictions also cannot be
substituted into low8 concentration without a proof for the new typed
rule. The order25 local star remains realizable locally and unchanged.

The next mathematical target is now precise: use the FULL source cusp
to finance identified finite producers, rather than add their coupons
at every recipient of the same abstract type. The affine-D9 mixed-A15
producer is completely resolved above. Any further producer needs its
own forced full factor, quantitative curvature saving and complete
current-output retention. This is the missing composition boundary
before a ceiling56 consequence, rather than another local-rule census.

## Source bindings and status

All listed files are unchanged; hashes identify inherited inputs, not
new executions. The present arguments are subject to independent review.

| Input under results/ | SHA256 |
| --- | --- |
| fresh_cusp_factor_payment_framework_20261017.md | `c65cfc5a4409a8cc83fcf0a083a6c1a7a29bb7bc7d649c11435773877e70ccab` |
| fresh_spectator_energy_below_one_20261017.md | `90ceea392d5201c531fd9952455bd706bda5f80eb5a80a8fb8d4f840364f1f6e` |
| fresh_order10_d_opposite_fork_fractional_transfer_20261018.md | `3f2154ec10e3782918c36f62157986e1831983dbbfb36efeec3a0ff6069564cb` |
| order11_fractional_cusp_closure_research_20261018.md | `3cdd8e76fb472d45a36ff4ad4ff75af17ffaf4cd0e7eaf965759566a3f46fb0c` |
| fresh_order11_fractional_cusp_closure_finite_budget_20261018.md | `37da62033e2f06b419475df17bae2fa476a6c1680640ee4030448ad3417af89b` |
| fresh_ambient57_payment_frontier_20261018.md | `6b6a9582260ac280f50e299250d226eec1d2dbca0e444b114d21da5a1eb1d474` |
| finite57_structural_reentry_targets_20261021.md | `822953841f6b551e590496b6d617afb67c5f5eceba37ad1fd9566c26d9ef53f7` |
| finite57_uniform_large_cycle_four_rule_and_boundary_profile_20261021.md | `e2a637b6a96ec650fe161afa5216a131a9493a5653cb160f38fb5b6473e297b8` |

The manuscript read is `manuscripts/finite_volume57.tex`, SHA256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The old negative-cusp diagnostic is not the final held ideal frontier:
the reviewed order11 composition already closes that ideal side. Its
finite new-positive family upper2965/8 and predecessor formal maximum
11861/32 exceed364. Later physical refunds remain separate preserved
inputs; no combined updated finite maximum is inferred here. No new
dimension bound, polytope, scientific replay, external refereeing or
end-to-end geometric formalization is claimed.
