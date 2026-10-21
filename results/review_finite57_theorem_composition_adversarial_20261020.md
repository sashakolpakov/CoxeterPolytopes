# Independent adversarial composition audit of finite-volume 57

Date: 2026-10-20. Scope: source-only mathematical composition, with no
scientific checker, Lean, TeX, or certificate generator executed. The
canonical source was not edited. The source's own claims were read and the
calculations below derived before consulting the companion lemma map.

Audited source: `manuscripts/finite_volume57.tex`, 156120 bytes, 3322 lines,
SHA-256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.

**Result: no mathematical composition blocker identified at the source's
declared exact-finite-output boundary.** This result concerns the implications
between the local rule, finite and ideal estimates, incidence formula, and
dimension-tail cases. It does not independently establish the saved exact
outputs. In particular, no new full-manuscript runtime PASS or end-to-end Lean
claim is made.

## 1. Dimension 58: the same charge reaches both sides of the contradiction

* Lines 90–116 define a charge on each actual corner, not on occurrences in
  proof diagrams. Nonnegativity handles every face with at least five sides;
  triangle and quadrilateral lower payments imply
  `D5 <= sum_v C(v)`. No strict face inequality is required.
* Lines 1253–1259 specify one finite/ideal function, and lines 1456–1469
  identify its finite component, cross, and context sums. At an ideal vertex
  only the varying affine factor contributes. Thus the finite maximization
  and cusp envelope bound the same corner rule used in the face theorem.
* The final negative guard has rank `13+3*15=58`. It leaves no extra component
  or repeated witness, and it reduces precisely two units to two halves.
  The positive guards also exhaust rank 58: `28+2*15=58` and
  `12+16+2*15=58` (lines 423–438, 634–643). Consequently they cannot overlap
  an uncounted extra context.
* The eight preceding exceptions are two choices for the rank-13 type times
  the four multisets of three rank-15 types. All have `78+3*100=378`.
  The B13 endpoint deletion lowers its four cases by one; the guarded D13
  reduction lowers the other four by one. The two increasing families cost
  `170+200+1/2=370.5` and `66+108.5+200+1=375.5`.
  The latter includes both A12 endpoint coupons. Every other forest changes
  pointwise downward (lines 624–645).
* The predecessor count is an explicit computed proposition, not an inference
  from example maximizing forests. The written cover fixes context-sensitive
  counts and refines all costly actual fibers. The residual proof requires
  every last-item upper inequality and a tight predecessor at every positive
  stored finite rank; positive item ranks rule out fabricated finite states.
  Rank zero is explicitly exempt from the no-predecessor/unreachable
  equivalence (lines 322–357, 374–408, 1564–1835).

The final local theorem's inherited consumers are not simply assumed to
retain their old units. The affine-D12 branch explicitly uses a distinct
retained A1/D12 payment or a guarded second half at a different corner
(lines 536–545, 581–589). The order of complete face inequalities is stated
at lines 1189–1249. At this composition level I found no call from a new
decreased role back to its own old unit. Exact chosen-role and domain
coverage assertions remain required finite inputs; their listed counts
alone would not suffice.

## 2. Nonsimple cusps, normalized descent, and the strict global correction

The rank used in a cusp link is 57, not 58. A product of simplices with
factor dimensions `r_i` has `sum r_i=57` and
`V=prod(r_i+1)` link vertices. An actual two-face corner corresponds to an
edge in this link. Such an edge varies in one factor: fixing its pair of
omitted facet labels leaves an independent choice of a vertex in every
other factor. Thus its multiplicity is `V/(r_i+1)`, exactly as used at
lines 272–287. This explains why summing affine charges without fiber
multiplicity would be wrong; the source does include that multiplicity.
The ideal rule has no cross-factor context payment (lines 1554–1559).

For the low8 move, `1 <= a <= b <= 7` ensures that both resulting ranks
remain at most 8. The envelope sum is unchanged. Convex order increases
`E|2S-57|`, hence decreases `G/V`, normalized curvature, and normalized
margin. Meanwhile
`(a+1)(b+1)-a(b+2)=b-a+1>0` decreases `V`. The source correctly uses the
positive terminal margin before multiplying: `V>=V'` and `m>=m'>0` imply
`Vm>=V'm'>=36`. No invalid inference from a normalized inequality to an
absolute surplus is made (lines 732–790). The integer square-rank potential
proves termination, and a move never reaches the non-linear rank-9 entry.

I independently recomputed the incidence constants from the definitions:

* A finite vertex has a 57-dimensional simplex link with 58 vertices and
  `E=57*58/2=1653` edges.
* `G=sum_{j=0}^{57} min(j,57-j)=2*(1+...+28)=812`.
* `(5G-57*58)/2=(4060-3306)/2=377`.
* At either global extremum the full link is uncut. Its excess over the
  minimum-index estimate is
  `E-G=sum_{j=29}^{57}h_j(2j-57)>=1+3+...+57=29^2=841`.
  All coefficients of a simplex/product-of-simplices h-polynomial are
  positive. Two distinct extrema therefore contribute
  `(5/2)*2*841=4205`, counted once (lines 801–827).

The induced-orientation assignment is consistent: below the level section
an endpoint receives its incoming-edge count `j`, and above it receives its
outgoing-edge count `57-j`. Hence every vertex contributes at least their
minimum. The published index distribution and incidence identities are
the external input; the inequalities and strict correction are internal.
Combining the final charge with this lower bound gives exactly
`0 >= 4205 + sum_ideal(kappa-C_actual) >= 4205+36*v_i > 0`
(lines 830–841). This also excludes an all-finite chamber; there is no
hidden assumption that a cusp exists.

## 3. Higher dimensions 59–163: coverage and concentration hypotheses

Writing `d=n-1`, a simplex link gives
`beta_n=n(n-6)/8` for even n and `(n-1)(n-5)/8` for odd n. These are the
coefficients printed at lines 1846–1849; the parity shift is correct.
The two-extrema argument remains strictly positive in every dimension,
so `U_n<=beta_n` is sufficient: equality causes no endpoint leak.

The intervals in the source are contiguous and exhaustive:
`59–65`, `66`, `67–68`, `69`, `70`, `71`, `72–74`, `75`, `76–77`,
`78–80`, `81–86`, `87–91`, `92–97`, `98–108`, `109–116`, `117–126`,
`127–136`, `137–146`, `147–156`, `157–163`.
They are not inferred by taking a Coxeter facet or by induction from 58.

For the early layers the stated finite maxima satisfy the common comparison.
In particular the delicate endpoints give
`U59=783/2=beta59`, `U60=404<405=beta60`,
`U61=420=beta61`, `U62=434=beta62`,
`U63=449<449.5=beta63`, `U64=464=beta64`,
`U65=479<480=beta65`, `U66=495=beta66`,
`U67=504<511.5=beta67`, and `U68=512<527=beta68`.
For n=69, `544=beta69`. These equality cases still retain the strict
global-extrema correction. The remaining listed finite comparisons have
the claimed direction (lines 1963–1974, 2231–2246, 2499–2510).

Every concentration move additionally requires nondecreasing envelope
cost, not only convex order. The source states this at lines 1879–1889.
For a linear/cap envelope, concentrating two parts below the cap or two
above it preserves cost and terminates with arbitrary cap parts and at
most one part on each side. The piecewise envelopes give their additional
positive duplicate-rank gains explicitly. For example, at n=66 the
rank-13 gain is `75/13+32/5-2*6=11/65>0`; the rank-12 gain would be
`16/3+6-2*(75/13)=-8/39<0` and is correctly not used. At 67–69,
`65/12+43/7-2*(75/13)=23/1092>0`. Thus the printed terminal domains do
not silently discard repeated rank 12 in the n=66 case
(lines 2206–2226). The 59–62 reduction is low-nine-only, leaving the
nonconstant high-rank envelope unaltered (lines 2491–2496).

The layer-specific simultaneous routing and exact positive minima remain
essential inputs. A proof of the final dimension-58 local rule alone is
not a proof of every earlier marked table, and a family total alone is
not a marked rule. The source supplies the layer selectors, complete
array/patch distinction and own-current-value routing requirement at
lines 2095–2166 and 2248–2520. This audit checked their composition and
types, not every distributed layer certificate.

## 4. All n>=164: original budgets precede residual budgets

The cusp branch of the narrow-facet argument uses the same ambient
`N=n-1`. One omitted node per affine factor gives an actual outgoing edge
in the chosen facet. The omission in its one distinguished factor avoids
the facet node. Rank-one avoidance is explicitly handled by repeating
the other omission and using target zero (lines 2541–2551).

At a finite endpoint, a terminating wall has a positive-definite pair with
the facet wall. At an ideal endpoint, the N-root support omits exactly
one node of each new affine factor. A missing node in the same factor
works except for affine A1; then a node of another factor works and is
orthogonal. Another factor must exist when an A1 factor occurs in the
rank range in question. The two-wall interior-point projection then
proves an actual ridge, so these terminating walls really are in the
framing set (lines 2553–2574). Merely being an arbitrary outgoing wall
would not have justified applying narrowness.

Both endpoint principal Grams are positive semidefinite, giving two
separate budgets `sum C0a<=2` and `sum C1a<=2`. Positive proportional
affine radicals give a common, strictly positive height ratio. In the
Lorentz identity the factor scaling cancels:
`sigma_a=k_a*b`, `tau_a=k_a*b'` imply
`sigma_a*tau_a*E_a = b*b' * ||u_a/b-u'_a/b'||^2`.
No common integral normalization is required (lines 2576–2626).

The full-budget classical unit-charge lemma explicitly starts with each
row budget at most two, not 3/2 (lines 2775–2793). Its B/C/D proof uses
positive inverse entries and the contact gap rather than a finite label
cutoff. Its A proof has an exhaustive finite KKT domain for orders 2–17,
checked at the actual wider ratio `21/20`, followed by positive all-order
polynomials beginning at 18 and 19 (lines 2882–2988). For example the
even endpoint margin is `128/[36*(18-2)]>0`; the odd endpoint margin is
`360/[36*19*(19-2)]>0`. Thus the residual cap is not used circularly to
justify the first failure.

The exceptional first-failure check has original-budget scope
`C0+C1<T`, not the residual 3/2 scope (lines 3056–3067). Its possibilities
leave residual budgets `(10/9,4/3)` or `(1,3/2)`/`(3/2,1)`, each coordinate
at most 3/2. Only then is the residual model enumeration used. Its ordered
packing bounds give `16/9+20/9=4` or `14/9+22/9=4`. Keeping a modified
residual record after a modified selected failure only enlarges the
domain. If no factor fails, the original two budgets themselves give
loss at most four (lines 3069–3077). This addresses both branches of
the failure selection.

For the other height range the independently proved coefficient-4/3
estimate sums to loss at most `(4/3)*(2+2)=16/3`. It is not used in place
of the stronger coefficient-one lemma needed above (lines 3079–3120).
There is no gap at height ratio `21/20`: both ranges include that value.

Finally, at `n=164`, `N=163`, the near-range product lower bound is
`M*(163-37)/9>=14`. The far-range bound is
`M*(163-43)/9 >= (21/20)*(120/9)=14`. The excess-rank terms are respectively
`t/9` and `7t/60`, positive or zero for `N=163+t`. Since framing narrowness
is strictly `<14`, equality at the threshold already contradicts it
(lines 2634–2644). No strictness was lost by using weak budget and
packing estimates.

## 5. A chosen narrow facet may have no cusp

This branch does not assume the facet is Coxeter. Its compact two-faces
are actual faces of the ambient chamber, and its vertices are finite
ambient vertices. The minimal-core argument gives two distinct paths in
a compact triangle, each of length at most eight, and one adjacent-mark
path of length at most 17 in a compact quadrilateral. The nonorthogonal
opposite cores have at most eight white nodes each, giving the bound
`8+1+8=17` (lines 3124–3156).

Ambient spherical distance-17 charge is at most `17n` per finite ambient
vertex: all pairs in rank at most 35 cost at most `17r`; larger finite
irreducibles are paths or D forks, with at most r pairs at any fixed
distance. Facet corner pairs form a subset of ambient corner pairs.
Thus the facet bound is `17(d+1)` for `d=n-1`, not `17d`
(lines 3157–3167).

The symmetric h-vector inequality is strict because the nonzero endpoint
coefficients lie above the middle minimum. Its required-payment lower
bound exceeds the finite charge by
`t^2/2+36t+3/2` at `d=143+2t` and
`t^2/2+(73/2)t+19` at `d=144+2t`. At t=0 these are 1.5 and 19, so
the branch excludes `d>=143`, equivalently `n>=144`
(lines 3167–3187). This covers every cusp-free narrow facet in the
unbounded range. Together with the cusp branch, the bounded ranges and
dimension 58, it excludes every n>=58.

## Remaining typed premises and what this audit does not certify

External geometric inputs are exactly the finite/affine classifications,
minimal-indefinite order ceilings, finite-volume cusp/simple-edge facts,
incidence/index formulas, and Nikulin's actual framing-wall narrowness
statement as declared at lines 857–864. Their historical citations were
not independently checked in this source-only task.

Internal finite outputs must retain all of the following types:

1. Complete local-domain generation, exact extension/seal predicates, and
   sums of **distinct current** corner contributions, for the dimension-58
   rule and every earlier ambient layer. Output counts or old corner units
   cannot replace these assertions.
2. The complete representative domination and two exceptional-family
   cover, followed by actual-type fibers and every residual transition
   inequality/tight predecessor used for the rank-58 maximum.
3. All exact terminal cusp profiles, convolutions and positive margins
   for 58 and 59–163, under their stated envelope and concentration rules.
4. At n>=164, original-budget small-A KKT outputs, the separate original-
   budget exceptional first-failure classification, residual support
   coverage and domination/packing bounds, and the independent arbitrary-
   ratio four-thirds finite outputs. These interfaces are not
   interchangeable even when they use the same grounded matrices.

The source gives domains and checking predicates for these internal
propositions. This audit found no additional internal geometric theorem
silently imported as a published axiom and no transitive hypothesis lost
between them. The appropriate next action is to retain the declared
distributed-data boundary and the existing scoped replay evidence; this
audit identifies no reason to rerun successful scientific batches or to
change the dimension claim.
