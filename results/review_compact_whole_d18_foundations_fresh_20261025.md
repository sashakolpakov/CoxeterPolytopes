# Independent review of the complete earlier whole-D18 appendix

25 October 2026. Reviewer: `/root/audit_sharp_ra`.

**Verdict:** the complete appendix has a coherent mathematical route at its
declared shared-lemma and finite-classification boundaries. I found one
incorrectly enlarged scalar assertion, F1 below, which must be repaired in
a preserved successor. I also recommend the explicit local height step in
F2. Neither finding requires a new bound, a new geometric premise or a new
search. The frozen source is not accepted unchanged by this report. No
whole-manuscript, numerical replay, Lean or TeX PASS is claimed.

## 1. Exact objects, independence and method

I read the entire 3,428-line, 127,358-byte fragment
`manuscripts/revisions/parts/compact_whole_d18_foundations_working.tex`,
SHA256 `9b0ba1b88efedb3a664e335a301d7d73a47341af3a23b3d0bde187e07bae0320`,
including every proof body, scalar table, composition paragraph and trailing
blank line. I read the entire focused ledger
`results/compact_whole_d18_foundations_transcription_ledger_20261025.md`,
SHA256 `69961018442080d98faaccccb2b34265a086d0fc78880b9b9ff1abb0180ebac7`,
the freeze receipt and the complete 55-file source-pin register
`results/compact_whole_d18_foundations_source_pins_20261025.json`, SHA256
`a086bd3fdba774ff1b31d9f96e835708ba55e0992dcd6125bf433eadb84b83ad`.
All 55 observed pins matched during this review. That metadata check did not
import a worker or evaluate a mathematical certificate.

The authoritative archive is
`results/source_archive/compact_whole_d18_foundations_transcription_20261025/`.
The 1,861-line pre-pair draft remains a separate preserved predecessor. I
did not edit either draft, the main manuscript, a historical proof source,
or any of the six canonical TeX/PDF pairs.

I previously reviewed the compact common foundation, the local A18 endpoint
fragment, and the B18 rewrite. Those contributions include identifying the
future-point perpendicularity/payment-strictness repairs and the two B18
exposition repairs. They are recorded in
`audit_compact17_foundation_transcription_fresh_20261024.md`,
`review_compact_whole_a18_endpoint_transcription_fresh_20261024.md`,
`review_compact_whole_b18_terminal_rewrite_fresh_20261025.md`, and
`review_compact_whole_b18_r1_r2_local_repairs_20261025.md` under `results/`.
I am a different author from this D18 transcription, but those already
reviewed common arguments are not being credited as a second fresh
independent derivation here. Root must retain its separate review role.

The present work was source-only. I launched no enumeration, symbolic or
scientific calculation, test, worker import, Lean, TeX or Sphinx job.
The mathematical checks below are direct reading and hand reconstruction
of the stated formulas, signs, domains and proof dependencies. No old
`PASS` or acceptance label was used as the proof of a new statement.

## 2. Required repair F1: the finite isolated subcritical cutoff

At lines 818--820, within `app:D18purecenter`, the text says that finite-g
monotonicity gives `delta < 8/5`. In its surrounding paragraph the positive
subcritical mutuals are **both** lambda = 1/2 and lambda = 1/sqrt(2).
The assertion, read over both values, is false as a scalar claim in the
displayed local parameter relaxation.

Set g = phi/2, an allowed native-five center own, and

```
a = 32*g^2 + 2*g - 1/2,
A = sqrt(a), k = 16*g + 1/2, a = 9*phi + 15/2 at this g,
F(delta,g) = A*sqrt(9*delta^2-1) - k*delta.
```

The elementary enclosures `.809 < g < .81`, `4.697 < A < 4.703`
and `13.444 < k < 13.46` give

```
F(8/5,g) < 4.703*4.695 - 13.444*(8/5)
          = .570185 < 1/sqrt(2),
F(3,g) > 4.697*8.944 - 13.46*3
        = 1.629968 > 1/sqrt(2).
```

The function is continuous and strictly increasing in delta > 1. Thus
its lambda = 1/sqrt(2) solution lies strictly between 8/5 and 3, still
within the permitted dotted-own and fork-cosine-less-than-twelve range.
This is a counterexample to the enlarged scalar cutoff, not a claim that
the entire forbidden geometric configuration is realized by a polytope.

The exact source F, Section 4, is more careful:
`compact18_d18_finite_own_weak_pure_centre_helper_complete_exclusion_candidate_20261023.md`.
It proves `delta_0 < 8/5` **only for lambda = 1/2**. For fixed g,

```
lambda/delta = A*sqrt(9-delta^(-2)) - k
```

increases with delta, and the solution delta increases with lambda.
Consequently BOTH subcritical values satisfy

```
lambda/delta > (1/2)/(8/5) = 5/16,
j-zeta = (k/a)*(lambda/delta) > (5/16)*(k/a).
```

The four strip endpoint bounds on k/a then yield exactly the printed
`.180, .170, .161, .153` column. The S16 lower table and its consumer
need this ratio bound, not a shared upper bound on delta. Replace the
unqualified cutoff by this source argument. The dotted-g `delta < 8/5`
claim later in the paragraph has its own stronger test and is not the
same error. Root and the author were notified; frozen bytes were retained.

## 3. Clarification F2: which local edge supplies the height

At lines 1940--1956, `app:D18localtwoA1` deals with a near-own-four,
Q=0 central opposite using its actual **seven- or nine-edge mate** of
length greater than nine. The following common display writes the bound
as a height in the central direction omega8. That does not follow from
the mate length alone without one more step.

The intended component-energy conclusion is already immediate: the mate
is in the SAME A16 component, and its normalized height exceeds .8,
which is more than required for e_A > .46. One clean repair is to state
the central omega8 display only for the C8 > 5 rows and state this mate
energy separately.

Alternatively strict-interior height concavity supplies the missing
central-coordinate step explicitly. With N=17 and D7=140/17,
D8=D9=144/17,

```
h8 >= (9/10)*h7, or h8 >= (8/9)*h9,
h8/sqrt(D8) > .8*min((9/10)*sqrt(140/144), 8/9) > .68.
```

Thus the displayed central conclusion is repairable at an already stated
geometric premise. The source V repeats the abbreviated central display;
this is not solely a new transcription discrepancy. No elimination fails.

## 4. Proof architecture and scope of the source comparison

The directed dependency order is valid after F1 is stated correctly:

1. Actual-normal positive projection, future orientation, original all-point
   maximum, positive-dual interpolation and the corrected common square.
2. The generic old-.61 prefix only in its legal rank range; the new actual
   D4/D3 tail inventory; strong old-five and genuine head isolation.
3. The fresh .58 inventory and its two newly admitted tail cases. These use
   early K11/K12 references, not later center/fork information.
4. The complete center-sixteen inventory and all five legal ending-fifteen
   reference contracts. Pure, finite nonpure, dotted nonpure and zero-own
   high-dihedral weak centers are all excluded before center strength.
5. Append actual strong sixteen, omitting a weak fifteen BEFORE pressure.
   Derive both fork bounds, the all-head envelope and the independent local
   neighbor lemmas.
6. Prove the standalone two-positive-dual A17/A18 and two-A1 comparisons
   before any II/IN/IF/FF endpoint consumer.
7. Close the all-head pair families; then close native-five and finite-seven/
   eight heads. The final dotted-head near/far/far-next argument remains a
   separate main-text consumer.

In particular the all-head FF proof does not import a finite-head reference,
the finite7/8 reference is extracted independently of source W's opening
near/far hypothesis, and the local two-dual proofs do not assume an isolated
pair exclusion that later invokes them. There is no general17/arithmetic15
cycle. The local A18 single-endpoint lemma is used as a local comparison
statement with its exact factor, not as a whole-original-farthest theorem.

I compared the substantive mathematical sections in all the ledger's short-ID
source families H/D/C/F/E/Q/R/N/P3/Z/L/T/B/N9/Env/I/II0/II-dot/IN-f/IN-d/V/
IF-f/J/FF/F2/G/W. This was comparison with the selected bodies used by the
fragment, not a fresh audit of every historical or deliberately pruned branch
in those files. The complete TeX bodies were read even where a source body
was reused. Important direct ranges include H78--454; D67--355; C68--329;
F65--444; E93--335; Q77--315; R83--370; N98--260; P3's actual-neighbor and
central proof; Z74--398; L65--413; T's terminal/local-neighbor/native-five
sections; B/N9's complete four-initial-case proofs; and all the selected pair
and finite7/8 mathematical bodies. The 55 pin matches include metadata-only
observations of canonical PDFs and earlier archives; those objects were not
mistaken for freshly read mathematical evidence.

In the observed main source `ac3f245a...`, I additionally reread shared
`app:Acentralwindow` (3241--3448), the relevant generic-prefix body beginning
4491, the reference/head-foot body beginning 6222, the full local-A16 regime
statement and foreign/ordinary inventory beginning 7724, and the precise
`lem:localA18` hypotheses beginning 1091. Their earlier scoped reviews remain
separate evidence; this report does not claim a fresh full-main-source audit.

## 5. Complete subsection review

| Frozen fragment | Evidence and mathematical assessment |
|---|---|
| 1--43, normalization | Unit D18 normals and all dual squares/crosses agree; x18 remains signed. All original reference uses concern one x. |
| 44--302, old head | M12 is legal at rank18; no M13/M14 geometric prefix is silently imported. Actual D4/D3 extra rows, both forks, the H4 point and both finite-side cores are present. Strong-five's dotted foot and each possible weak-fifteen case retain separate anchor prices and full crosses. Head isolation follows only afterward. |
| 303--452, fresh window | Old strength is not renamed .58 strength. Own-four uses the actual affine B principal. The new native-seven/fork row at14 and dotted/end row at15 are both rechecked on their whole domains. Early J6 weight .19 fits; late .21 is not used there. Later reference replacement does not add both payments together. |
| 453--572, center inventory | Full A15+2A1 white list, H3(3,5), both fork contacts and all own-zero rows are retained. Reciprocal B18/D18 endpoint exclusions are actual positive projections. The necessary weak list has exactly pure helper, single-fork-three positive own, and single-fork own-zero labels7--10. |
| 573--662, references/ellipse | Increment identities, nonadjacent omissions and five ending15 contracts are stated independently of the selected center. The ellipse is one paid-reference consequence. Its z0<=10P condition is checked at consumers, not inferred by optimizing an unrelated norm. |
| 663--1026, pure weak centers | Correct common cleared Gram and both signed fork tables; singular, subcritical and supercritical boundaries are separate. All far profiles and zero mutuals are included. F1 is a genuine statement-scope repair; the correct ratio argument retains the existing table. |
| 1027--1151, nonpure centers | Finite positive-own branch first proves fifteen strong with the FULL fork contact. Its later .13 payment uses P15, not an absent P14 price, and keeps the negative reference cross. Dotted and own-zero branches use actual A15+A2+A1 / A15+I2+A1 neighbors, without a farthest assertion there. |
| 1152--1314, local A15 | Actual central finite-white/null inventory includes own zero and both helper orientations. The full A2+A1 foreign alphabet is checked. High I2(m>=7) cannot be contacted because its cost exceeds one. Factor3/2 is used only in the high-dihedral local energy result. |
| 1315--1356, strong sixteen | Every weak-center branch has already failed. Reassembly omits a weak15 before pressure; the conservative old debit may remain. Net append>.016119 and terminal price>.127 survive the zero-sum early corrections. Both actual fork bounds<10 use this SAME reference. |
| 1357--1451, fork inventory | All seven original white profiles and own-zero/native/parallel/dotted branches are dispatched before the four-role remainder. Signed T<3.823 follows from terminal surplus, not a second optimized energy. Actual A18 endpoints use factor<11/8 and endpoint>.78; near/near acuteness is uniform over every head. |
| 1452--1570, large finite heads | Three bare references and the fourth weak-helper reference are not conflated with previously paid lower bounds. Entire ten-plus interval and five u-strips include the full head cross. Native-nine isolated-five rows pay one combined anchor debit and subtract their positive mutual cross. The helper row uses a native gap. |
| 1571--1711, envelope/cutoff | One unspent P16 coefficient gives every asymmetric direction. The all-point distance cap is derived from original maximality. Isolated-own15/8 cutoff uses an actual A17+A1 neighbor at factor8/5, both end orientations and a full-correlated H4 exceptional square. |
| 1712--1857, two duals | Both coefficients of the reverse unit remain positive and present. A17 uses sites8/9 to keep all side orders>=8. A18's right mates use three directions. Both strong-anchor capacities are individual; the exact negative residual correction is retained. These are standalone local statements before pair consumers. |
| 1858--2032, two A1s | Both foreign coordinates and the entire native cost alphabet are retained. Actual zero-own/null cases, pure helper ends and dotted intervals are covered. The native-five point retains its foreign component and uses one correlated square. F2 supplies the omitted central/mate link or a cleaner split. |
| 2033--2233, finite/orthogonal II | Positive mutuals give actual full dihedral neighbors with explicit factor; the A2 extension is supplied rather than applying B18's m3 lemma outside its factor. Orthogonal unequal owns give the actual A15+I2+A1 neighbor and separate factor1.7 local proof. No symmetry or new farthest choice is used. |
| 2234--2350, dotted II | Own-product and individual cutoffs force both owns dotted. The max-foot bound is the exact441225/6208, sufficient for beta>.99294; merely71 would not justify that pressure. The actual forward edge supplies an endpoint. Full B17/D17 contacts and both end orientations are checked before the two standalone A consumers. |
| 2351--2473, IN | Finite, parallel, dotted and orthogonal mutuals are distinct complete domains. Dotted proof uses the continuous rational quartic/tangent certificate, not an own-label sample. Orthogonal case retains the whole two-A1 reverse direction and invokes precisely the local factor33/20 result. |
| 2474--2654, nondotted IF | The product inequality permits equality. Zero mutual, B18 projection, H3 projection with its FULL P2 direction, affine boundaries and continuous finite7+ interval all remain separate. Native-three gives actual A18 with both reverse duals, using the earlier independent local lemma. |
| 2655--2914, FF/dotted IF | Zero FF is an actual positive D18 projection; native-three and parallel are actual nulls. Both FF owns are derived dotted; IF keeps finite far owns .92<h<1. Feet are feasible and unequal-own formulas retained. Released-wall sinh ordering proves the next endpoint is farther. All seven germs, both shadows, full A18 inverse cross and the final finite-own incidence are present. |
| 2915--3041, native-five head | Finite-cosine conjugate tests are elementary restrictions on an ACTUAL mutual, not arithmeticity assumptions on P. The entire near/far rectangle has positive derivatives and falls in the native gap. The six unordered pair classes invoke only earlier all-head statements. |
| 3042--3158, finite7/8 reference | Extracted without source W's NF premise. Three bare rows and fourth helper row use one vector each; combined head/U5 prices and the full positive cross are explicit. Every new debit is before seven, leaving P16 untouched. The old-A15 all-point cap keeps both signed final coordinates. |
| 3159--3377, finite7/8 forks | Entire finite far-next list5/7/8/9/10/11 and infinite far-three finite-own tails are covered by monotonicity and consecutive native gaps. Remaining NF mutuals0 or1/2 are fully dispatched using actual A18, actual finite dihedral neighbor, or actual dotted core. The all-point cap is used for a core that need not be a vertex. |
| 3378--3428, composition | Finite heads are exhausted. Only the stated earlier five all-head pair exclusions pass to the separate final dotted-head argument. No whole-D18 conclusion is assumed by a foundation used to prove it. |

## 6. Sensitive algebra, original walls and coefficient currency

The selected-point checks use one of three legitimate constructions:
positive finite projection `p + sum(lambda_i*m_i)`; a positive actual-normal
core with explicitly nonpositive selected products; or a positive future
null contradicting compactness. Corrected spacelike roots are used through
their finite projected point, not called points or new facets. Unknown wall
signs hold termwise because the summands are actual distinct walls. Selected
caps and released roots require their separate product checks; these are
present at the important two-cap feet, central cores and neighbor projections.

The actual original-shadow inequality used at the FF exits is legitimate
for the new wall, which is distinct from both original opposites. For a far
opposite the positive old inverse Schur complement gives
`v.Omega1 >= (v.Omega_f)/9`; together with
`F(g) < 3g+1/6`, acuteness gives `(v.Omega_f)/(3T)<1`.
For an isolated opposite it follows directly from `A<3g`. Reflection gives
the other shadow. This uses the full original contact vector; it does not
give an arbitrary new wall the finite-white restrictions of an original
opposite. In the finite-far-own case, distinctness uses the complete
nonpositive cycle, not an unavailable dotted-own argument.

The price checks of greatest consequence are retained:

- Old weak-five .07 foot: BOTH P1/P6 prices; its .088 P6 debit is below
  .08806. Large finite/dotted weak15 tails share P1 but cost no P6 and
  have zero residual cross with that foot. Native5/7/8 tails do not also
  receive the weak-five foot.
- Fresh early isolated-five J6 .19 debit, not .21; the later .21 version
  uses the larger K14/15 capacity and replaces the early reference.
- All ending15 head corrections have zero sum through seven. Appending
  sixteen does not change their prices or full crosses.
- Nonpure finite center: terminal .13, full residual square<1.174 and
  reference cross `-a*w15`; an absent P14 capacity is never spent.
- Native-nine and finite7/8 head+U5: the combined debits .85/6 and .9/7
  fit individually. Their cross `1-5/J` is subtracted, not discarded.
- A18 two-dual reference: each of its two possible terminal anchors has
  its own capacity; the exact residual square
  `1-2*K*L^2/(19*m)` is retained. No independently optimized energy or
  single-endpoint surrogate replaces the two-dual reverse.

## 7. Remaining proof and verification boundaries

The external inputs remain compact acute simplicity, actual-wall native
products and the finite/affine Coxeter classification at stated ranks,
together with elementary real algebra/cyclotomic facts for finite cosines.
This review does not freshly authenticate a literature classification.
No field of definition, integrality, arithmetic reflection hypothesis or
inheritance of Coxeter/arithmetic structure by an arbitrary face is used.

The tables here are explicit internal small inequalities and inverse-Gram
identities, not an imported large enumeration output. I checked their
formulas, parameter ranges, monotonic directions, signs and their uses,
and reconstructed the sensitive rational comparisons discussed above.
There has been no independent automated replay of every printed decimal
or all transitive historical source arithmetic. Such selected finite
verification remains a separate root-owned step, not an external theorem
to be silently assumed in a formalization.

After F1 and the F2 wording are repaired, a source-bound minimal-diff review
should verify the successor before integration. That can preserve the rest
of this full-body assessment; it does not require restarting the earlier
accepted A18/B18 proofs or rerunning successful modules without reason.
Root's own full-source review and acceptance remain required.

This appendix supplies earlier whole-D18 foundations only. The main final
dotted-head near/far/far-next proof, other exact18 consumers, higher-dimensional
general-only covers/consumers, their algorithm contracts and any conditional
Lean/TeX work remain at their separate declared scopes. No global ceiling,
canonical edition, current PDF, reproducibility catalogue or verification
status was changed by this review.
