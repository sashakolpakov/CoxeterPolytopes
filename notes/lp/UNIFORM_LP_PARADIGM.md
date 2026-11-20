# Uniform LP recovery of the existing dimension bounds

Updated 18 November 2026. This is the main explanatory note for the recovered
bounds **57, 17 and 15**. It describes the common method, the actual LPs,
their retained mathematical inputs and the precise sense of uniformity.
The [execution results](existing_proof_lp_recovery_results_20261118.md)
give model sizes, timings, artifact paths and acceptance records.

The present task was to express the **existing proofs** through LP, keeping
their proved local inequalities and complete applicability domains. That
task is complete at those inputs. Deriving the local inequalities anew from
primitive Gram constraints, removing imported classification results and
improving the bounds are separate, deferred tasks.

## 1. The common idea

Each proof has a globally constrained quantity and local estimates for it.
The LP chooses a compatible linear combination of those estimates.

| Case | Global quantity | Local estimates and compatibility | Contradiction |
| --- | --- | --- | --- |
| Finite volume | Sum of charges on vertex–two-face incidences | One shared corner charge pays every small face and stays within every vertex allowance | Incidence requires a strictly larger total than the vertex allowances permit |
| Compact | Squared norm of one original farthest vector | Component and joint-component energy bounds, sharing the original orthogonal energy variables | Local bounds force total energy at least one, whereas geometry gives total energy less than one |
| Compact arithmetic | The same original-vector energy | Compact estimates plus arithmetic local bounds, field restrictions and matched arithmetic exclusions | The same energy contradiction, with complete arithmetic domains retained |

The finite-volume LP and the compact energy LP are different physical
models. They share the certificate principle: preserve the meanings of
local quantities, identify quantities that must agree, then find a valid
linear combination contradicting the global budget. The two compact cases
also share the **same implemented assembler**.

An abstract finite LP makes the common logic explicit. After adding slacks
and splitting free variables when necessary, write a necessary feasibility
problem as

\[
             M_N z=d,\qquad z\ge0.
\]

Here \(N\) is the dimension threshold, \(z\) contains local masses or
scaled quantities, and the rows of \(M_N\) express the retained conditions
and compatibility identities. Every hypothetical object in the covered
dimensions must give a feasible \(z\). A certificate is a vector \(y\)
with

\[
             M_N^Ty\ge0,\qquad d^Ty<0.
\]

Multiplying the feasibility equations by \(y\) would give
\(d^Ty=(M_N^Ty)^Tz\ge0\), a contradiction. This is the finite-dimensional
Farkas certificate principle. It is not a claim of strong duality for an
arbitrary infinite geometric domain. Analytic tail theorems and exact
extended formulations must justify any infinite family represented by
finitely many certificate conditions.

Here the necessary geometric feasibility problem is the **geometric
primal**. A feasible certificate on its dual side excludes an actual
object. The finite-volume minimax search for charge coefficients is also
an LP in its own right, with its own feasible and dual optimality witnesses.
Those two uses of “dual” refer to different explicitly identified problems.

### Four meanings of uniformity

1. **Uniform interpretation.** A shared variable has the same geometric
   meaning wherever it occurs. A face-side corner feature is the same
   feature used at its vertex. An energy component is always a projection
   of the same original vector.
2. **Uniform generation within a family.** A fixed theorem or rule is
   instantiated over its stated ranks, types, marks and parameter ranges.
   Different geometric and arithmetic families may use different rules.
   An unbounded dihedral family is retained symbolically only when its
   theorem covers that entire family.
3. **Uniform assembly.** The assembler introduces the local rows, identifies
   their shared quantities and enforces the global budget. It does not
   insert the desired global dimension bound as a row. Dimension-specific
   local theorems remain explicit input data.
4. **Uniform thresholds.** The underlying domains and rules stay fixed;
   increasing \(N\) only removes admitted dimensions. The exclusion therefore
   applies to every larger dimension, subject to the explicitly retained
   tail theorem. The current catalogues have a stated lower range; lowering
   the threshold below it requires additional justified domains and inputs.

Uniformity does not mean that linear programming has reproved the local
geometry, field theory or classification. It means that their precisely
stated outputs enter a common, reviewable linear argument.

## 2. Three different kinds of conditions

The distinction below prevents a restriction on the *search for a proof*
from being mistaken for a necessary condition on a polytope.

| Kind | Meaning | Examples |
| --- | --- | --- |
| Necessary geometric or combinatorial conditions | Every actual configuration in the stated domain satisfies them | Corner conservation; spherical vertex components; nonnegative orthogonal energies and their strict total budget; arithmeticity under its original-field hypotheses |
| Proved local consequences used as inputs | A retained theorem supplies an inequality or excludes a precisely specified local domain | A component energy floor; an estimate involving one defect and another component; a local null-vector contradiction; a complete necessary cover |
| Certificate conditions | They ensure that the chosen linear combination is a proof | Face-payment restrictions on searched charge coefficients; nonnegative multipliers; energy-slot capacities; sufficient frontier domination conditions |

The second kind is itself necessary for actual configurations satisfying
the theorem's hypotheses, but its geometric proof is retained rather than
rederived by the global LP. The third kind can be a sufficient restriction
on a certificate family. It need not characterize every possible valid
certificate.

Every imported statement needs its complete domain: actual rank and type,
marks, parameter interval, spectator assumptions, original-vector identity
and any field or lattice hypotheses. A source label alone does not prove
that it applies. The source reviews check those applicability matches.
The LP checks the final linear combination.

## 3. Finite volume: the charge and incidence LP

### 3.1 Corners, features and the searched charge

Fix a dimension \(n\ge58\). A corner \(\kappa=(v,F)\) is an incidence of
a vertex \(v\) and a two-dimensional face \(F\). Let \(s(F)\) be the
number of sides of \(F\). The fixed feature vector
\(\phi(\kappa)\in\mathbb Q^{26}\) records the constant feature
\(\phi_0=1\), finite/ideal baseline weights and specified conditions on
marked diagram components. Its [complete coordinate definitions](../finite_volume_parametric_corner_features_20261107.md)
are part of the input. Feature values are prescribed by the configuration;
their charge coefficients are searched by the LP.

Write \(a^{58},a^{\mathrm{tail}}\in\mathbb R^{26}\) for the coefficient
vectors at \(n=58\) and \(n\ge59\). They are called \(e,c\), respectively,
in the original implementation. This note uses \(a\) to distinguish them
from the energy variables in the compact cases. Put \(\epsilon_n=1\) if
\(n=58\), and \(\epsilon_n=0\) otherwise. The charge is

\[
w(\kappa)=
\bigl(\epsilon_n a^{58}+(1-\epsilon_n)a^{\mathrm{tail}}\bigr)
 \cdot\phi(\kappa)+\epsilon_n\lambda\chi_B(\kappa).
\]

Here \(\lambda\ge0\) is a searched coefficient and \(\chi_B\) is a fixed
indicator of a particular ideal-corner configuration. Its complete varying
affine factor is native \(B[13]\), the marks are \(\{10,11\}\) or
\(\{10,12\}\), and its complete corner diagram contains at least three
distinct connected components of types \(A_{15}\) or \(D_{15}\) disjoint
from both marks. Components created by deleting the marked nodes do not
count for this predicate. The [accepted composition](../finite_volume_common_law_dispatch_and_incidence_composition_20261109.md)
gives the exact native conventions and residual-node choices.

The 53 conserved feature coordinates are

\[
       (\epsilon_n\phi,(1-\epsilon_n)\phi,\epsilon_n\chi_B).
\]

They are counted on the same actual corner incidences on both sides of
the face/vertex sum. A nonnegative margin \(\rho\) is a separate searched
variable, giving 54 shared certificate variables in total. It is not a
54th corner feature. Complete native numberings and marks are coherent
across corners belonging to the same vertex.

### 3.2 Local payments and vertex allowances

The retained local proofs turn linear coefficient restrictions into:

\[
w(\kappa)\ge0,\qquad
\sum_{\kappa\text{ on }F}w(\kappa)\ge
\begin{cases}2,&s(F)=3,\\1,&s(F)=4.\end{cases}
\]

For \(s(F)\ge5\), nonnegative corner charges already pay the nonpositive
demand \(5-s(F)\). These statements retain the complete face-domain
coverage, exact provider associations and supporting-wall arguments of
the existing proof. The LP does not infer that coverage from a generic
continuous Gram assignment.

Set \(C(v)=\sum_{\kappa\text{ at }v}w(\kappa)\). A finite vertex has
allowance

\[
\beta_n=\begin{cases}
n(n-6)/8,&n\text{ even},\\
(n-1)(n-5)/8,&n\text{ odd}.
\end{cases}
\]

For an ideal vertex put \(r=n-1\). Its link is a product of simplex
factors of dimensions \(r_i\) with \(\sum_i r_i=r\). Define its
\(h\)-polynomial by

\[
\prod_i(1+t+\cdots+t^{r_i})=\sum_{j=0}^r h_jt^j,
\qquad V_L=\sum_jh_j=\prod_i(r_i+1).
\]

The ideal allowance is

\[
\kappa_n(v)=
\frac{5\sum_{j=0}^rh_j\min(j,r-j)-rV_L}{2}.
\]

The certificate conditions ensure
\(C(v)+\rho\le\beta_n\) at finite vertices and
\(C(v)+V_L\rho\le\kappa_n(v)\) at ideal vertices. The symbol
\(\kappa_n(v)\), an allowance, is distinct from a corner \(\kappa\).

Let \(B_n\) be the sum of these allowances over all vertices, and put

\[
                    D_5=\sum_{\dim F=2}(5-s(F)).
\]

Summing face payments, then counting the same corners by vertices, gives

\[
                 D_5\le\sum_v C(v)\le B_n
\]

even at \(\rho=0\). The retained incidence theorem gives

\[
                  D_5\ge B_n+\Delta_n,
                  \qquad \Delta_n\ge5(n-1)>0.
\]

The contradiction therefore requires a feasible charge certificate, not
a strictly positive value of the certificate margin \(\rho\).

### 3.3 The incidence feasibility problem behind the charge

For fixed dimension \(n\), denote the 53 conserved features above by
\(\Psi(\kappa)\). Let \(V\) range over allowed vertex/link types and
\(Q\) over allowed triangles and quadrilaterals. Define

\[
a(V)=\sum_{\kappa\text{ at }V}\Psi(\kappa),\qquad
a(Q)=\sum_{\kappa\text{ on }Q}\Psi(\kappa),\qquad d(Q)=5-s(Q).
\]

Write \(b_n(V)\) for the relevant finite or ideal allowance. Introduce
nonnegative masses \(x_V\) for vertices, \(y_Q\) for small faces, and
\(z_\kappa\) for the residual corners on faces with at least five sides.
The feasibility equations are

\[
\sum_Va(V)x_V-\sum_Qa(Q)y_Q-\sum_\kappa\Psi(\kappa)z_\kappa=0,
\qquad
\sum_Qd(Q)y_Q-\sum_Vb_n(V)x_V=1.
\]

An actual polytope gives these masses by scaling its corresponding counts.
If \(f_3,f_4\) count triangular and quadrilateral faces, the positive
normalizing factor is

\[
           L=2f_3+f_4-\sum_vb_n(v)>0.
\]

Indeed \(2f_3+f_4\ge D_5\ge\sum_vb_n(v)+\Delta_n\). Divide every count
by \(L\). The residual-corner variables conserve all features of the
remaining faces; their contributions must not simply be dropped.

The charge coefficient vector is
\(\theta=(a^{58},a^{\mathrm{tail}},\lambda)\). Its three certificate
families are

\[
b_n(V)-\theta\cdot a(V)\ge0,\qquad
\theta\cdot a(Q)-d(Q)\ge0,\qquad
\theta\cdot\Psi(\kappa)\ge0.
\]

Multiplying them by \(x_V,y_Q,z_\kappa\) and summing gives \(-1\ge0\),
using feature conservation and the second normalized equation. This is the
direct weak-duality explanation of why a charge certificate excludes the
polytope. Infinite supports may be read as finite-support nonnegative
measures; an actual polytope uses finitely many incidences. No infinite
strong-duality assertion is needed.

These equations explain the underlying necessary feasibility problem.
They are not a claim that this full geometric incidence primal was newly
exported or solved. The executed LP searches sufficient charge conditions
through the retained local proof and network blocks.

### 3.4 Why a network LP handles all finite vertex diagrams

A formal forest is a disjoint union of admitted finite Coxeter diagrams
whose ranks sum to 58. Every actual finite vertex gives such a forest;
realizability of every formal forest is unnecessary.

The retained exact construction represents these forests by paths in
finite directed acyclic networks. For an arc \(u\to v\), let \(q_{uv}\)
be its fixed feature contribution. Here \(u,v\) are network states, not
polytope vertices. Within each fragment \(S\), introduce private potentials
\(p_{S,u},p_{S,v}\) and impose

\[
 p_{S,v}-p_{S,u}\ge a^{58}\cdot q_{uv},\qquad p_{S,\rm start}=0.
\]

Terminal-entry rows add the remaining fixed features and bound the final
charge by a fragment allowance \(U_S\). All network fragments share the
same \(a^{58}\) and \(\rho\). Path inequalities telescope, while
longest-path potentials give the converse. This is an exact extended
formulation of the retained path budgets, without integer optimization
variables. The graph coverage and feature interpretation remain proved
inputs; a graph index alone is not a geometric justification.

### 3.5 The two roles of 377

First, incidence gives the finite-vertex allowance
\(\beta_{58}=58(58-6)/8=377\).

Second, the specified law family has an independently certified minimax
optimum of 377. Replace the fixed network allowance by a free common
upper variable \(U\), keep the declared local coefficient conditions and
\(\rho=0\), and minimize \(U\). Equivalently,

\[
       \min_{a\in\mathcal A}\ \max_F a^{58}\cdot\Phi(F).
\]

Here \(a\) denotes the full shared law data, \(\mathcal A\) is the
specified local-certificate family, \(F\) ranges over all retained formal
rank-58 forests, and \(\Phi(F)\) is their aggregate corner-feature vector.
The family includes the two stated A13 receiver restrictions as well as
the other retained local obligations. These are sufficient restrictions
on this certificate search, not a characterization of all possible laws.

To see the lower bound, write \(a_0,a_1\) for coordinates 0 and 1 of
\(a^{58}\). The retained corner/triangle rows give
\(a_0\ge0\) and \(3a_0+2a_1\ge2\). The formal forest
\(F_*=2D_{17}+D_{24}\) has feature vector
\(1653u_0+377u_1\), where \(u_i\) is a standard coordinate vector and
\(1653=\binom{58}{2}\). Thus the exact nonnegative-combination identity is

\[
U-377=(U-1653a_0-377a_1)
 +\frac{377}{2}(3a_0+2a_1-2)+\frac{2175}{2}a_0\ge0.
\]

A complete compatible rational feasible assignment supplies the reverse
inequality \(U\le377\). Matching certificates prove the optimum. No
numerical optimizer was needed, and this is not optimality over arbitrary
geometric charging rules. The capped system forces \(\rho=0\); the strict
incidence excess still gives the dimension contradiction.

### 3.6 All dimensions and actual evidence

The endpoint network is not the whole higher-dimensional argument. The
accepted composition also retains finite-vertex estimates for \(n\ge59\),
ideal-vertex estimates at 58 including the guarded cases, and ideal estimates
at 59 and every \(n\ge60\). The face coverage is valid at every \(n\ge58\).
Together they exclude every integer dimension at least 58. No reduction
to a 58-dimensional Coxeter face is assumed.

The actual coupled certificate consists of the exported local block and
verified network blocks; a single flattened global matrix was not exported.
The local block has 5,810 variables and 103,606 rows. Complete finite58
coverage and the matching eight-row minimax lower certificate are bound by
the [10 November acceptance](../../audit/joint58_complete_lp_execution_20261110.json).
The [15 November rerun](../../audit/lp_rerun_report_execution_20261115.json)
rechecked the declared certificate/assembly tasks; it did not reexecute all
1,228 original network leaf jobs. These records establish an exact LP
solution at the retained proof inputs, not a new geometric formalization.

## 4. The compact energy LP used for both remaining cases

### 4.1 One vector, complete branches and local rows

The retained compact geometry chooses one original farthest vertex and
one original spatial vector \(x\), with \(\|x\|^2<1\). Orthogonal incident
components give energies \(e_j=\|\operatorname{proj}_j x\|^2\ge0\), with
\(\sum_j e_j=\|x\|^2\). An aggregate slot may group several of those
orthogonal components. No energy of a new auxiliary or neighboring vector
is added to this account unless a retained theorem has already transferred
its conclusion back to the original vector.

A branch \(b\) is a complete admitted configuration domain at an actual
dimension. It records component types/ranks, actual defect ownership and
every hypothesis needed by its local rows. A complete necessary cover
means that every actual configuration matches at least one branch. The
branches need not be disjoint or realizable; one matching branch suffices.

Write the retained rational energy consequences as

\[
                      A_be_b\ge c_b,\qquad e_b\ge0.
\]

The column index of \(A_b\) is an original energy slot. A usual row has
entry 1 on the support of a local theorem, 0 elsewhere, and a lower bound
in \(c_b\). Strict lower bounds may safely be weakened to closed bounds:
the global upper budget remains strict. An irrational bound needs a proved
rational lower estimate; arbitrary numerical rounding is not a proof.

Alternative estimates on the same component act through their maximum.
Joint estimates remain simultaneous constraints. For example,
\(e_I+e_A\ge a\) and \(e_I+e_B\ge b\) do not imply
\(e_I+e_A+e_B\ge a+b\): that would count \(e_I\) twice.

### 4.2 The normalized feasibility problem

For each branch introduce mass \(\tau_b\ge0\), scaled energies
\(E_{bj}\ge0\) and deficit \(\delta_b\ge0\):

\[
 A_bE_b\ge c_b\tau_b,\qquad
 \delta_b=\tau_b-\sum_jE_{bj},\qquad
 \sum_b\delta_b=1.
\]

A covered energy point with \(s=1-\sum_j e_{bj}>0\) lifts by
\(\tau_b=1/s\), \(E_b=e_b/s\), and zero other branch masses. Conversely,
a feasible normalized LP has some \(\delta_b>0\); division by that
branch's positive mass gives an energy point below one. Thus strictness
is encoded without a chosen epsilon.

An imported local contradiction applies only to its complete matched
domain and gives \(\tau_b=0\). It may be a null-vector, native-contact,
acuteness or arithmetic contradiction, with no asserted energy floor.
Such a branch cannot represent an actual configuration satisfying that
theorem. It is not credited as a new geometric deduction by the LP.

### 4.3 The capacity certificate actually checked

For every nonexcluded branch the solver finds weights \(\alpha_b\ge0\)
and a common \(\eta\ge1\) such that

\[
                 A_b^T\alpha_b\le\mathbf1,
                 \qquad c_b^T\alpha_b\ge\eta.
\]

The vector \(\mathbf1\) has one entry per energy slot. Its capacity
inequality prevents any slot, including a shared defect, from being used
more than once. The certificate gives

\[
 \sum_jE_{bj}\ge\alpha_b^TA_bE_b
             \ge(c_b^T\alpha_b)\tau_b\ge\tau_b.
\]

Every deficit is consequently zero, contradicting their normalized sum.
The exact checker verifies all capacities, payments and variable bounds
as rational equalities or inequalities. `EXACT_LP_FEASIBLE` in the saved
results describes this **dual certificate problem**; the normalized
geometric-budget primal is certified infeasible.

The implementation's numerical objective maximizes the sum of local
payments subject to these constraints. It does not change the feasible
set or require a claim of numerical optimality. The accepted result is
the exactly checked rational feasible dual point.

## 5. Compact 17: uniformly generating the domains

The fixed catalogue covers every actual dimension 18–29. The independent
published compact ceiling 29 excludes dimensions at least 30. Call an
\(A_r\) component with \(r\ge7\), or a \(B_r,D_r\) component with
\(r\ge9\), dominant. The retained component theorems give energy at least
\(1/2\) for each dominant.

Every ordinary nondominant component costs at least its rank divided by
14. The retained actual-wall compatibility theorem permits at most one
of the following exceptional profiles, with no coexistence:

| Profile | Original ranks it owns | Retained cost or deficit |
| --- | --- | --- |
| \(I\) | One old \(A_1\) with its actual isolated dotted short opposite wall | Its nonnegative energy may lose at most \(1/14\) against the rank rate |
| \(J\) | Two old \(A_1\) factors in the actual common short-wall profile | Joint energy at least \(17/50\), so it pays both ranks |
| \(N_7\) | One old \(A_2\), with its actual isolated native-7 wall | Energy at least \(9/100\) |
| \(O_7\) | One old \(I_2(7)\), with its actual isolated native-3 wall | Energy at least \(9/100\) |

The typed profiles' rank-rate deficit is \(2/14-9/100=37/700<1/14\).
Consequently nondominant total rank \(S\) has the uniform aggregate floor
\((S-1)/14\). A bare seven bond is not an actual \(O_7\) profile. A marked
support never also receives the ordinary bound it replaces.

The generator covers all possibilities as follows:

1. **No dominant:** one aggregate slot has floor \((n-1)/14>1\).
2. **At least two dominants:** two actual dominant slots each have floor
   \(1/2\); an orthogonal remainder retains all other energies.
3. **Exactly one dominant of rank \(r\):** put \(S=n-r\). If \(S\ge8\),
   the dominant floor and spectator floor already total at least one.
   Stronger applicable component floors can similarly represent an entire
   context by two explicit aggregate rows.
4. **Contexts needing refinement:** generate every nondominant spherical
   spectator partition and every allowed actual defect profile. In this
   source data only remainders of rank at most four need nontrivial
   refinement. Retain all applicable component, joint and propagated rows,
   and only exactly matched local contradiction theorems.

The coarse contexts remain explicit branches in the LP; they are not
silently discarded by a separate proof of the answer. The LP certifies
their simultaneous energy contradiction. The historical eleven-case
rank-18 forest list is not assumed by this generator.

The ordinary rank-two types include \(A_2=I_2(3)\) and the full symbolic
family \(I_2(m)\), every integer \(m\ge4\). No finite native-label sample
replaces that family. Exceptional spherical types and rank-eight components
remain covered by the ordinary theorem even when they cannot fit a refined
small remainder.

The [complete compact family table](compact17_complete_energy_lp_recovery_20261118.md)
lists the source bounds, prefix rank conditions and local exits. The actual
catalogue has 829 branches, including 151 imported local contradiction
leaves. The capacity dual has 3,157 variables and 2,224 rows. It passes
the solve, exact check and fresh whole-model regeneration. This recovers
\(n\le17\) with the explicitly retained external ceiling.

## 6. Compact arithmetic 15: the same assembler with arithmetic inputs

At actual dimensions 16 and 17, the initial recovery retains the existing
**proved pruned covers** as inputs. Their 19 and 14 source groups expand
to 49 and 39 complete parameter-family branches, respectively. The LP
does not yet derive that arithmetic pruning from all spherical forests.

It does not search over number fields or infer integrality from a real
Gram matrix. Arithmeticity enters through the admitted domains and the
validity hypotheses of the retained arithmetic theorems.

The additional retained families are:

* Arithmetic local component and joint estimates, with their original
  field, root-line, integral-cycle, field-square and nonphysical-embedding
  hypotheses.
* Four named computed-output premises: quadratic \(A_{17}\) emptiness;
  the two cubic \(A_{17}\) anchor domains; the shared quadratic
  \(B/D_{17}\) output with two survivors; and the whole-\(B_{16}\) output
  with one survivor. The last two outputs are not emptiness statements:
  their distinct actual-domain inclusions and survivor consumers remain.
* The separate whole-\(D_{16}\) dyadic argument and every matched local
  arithmetic contradiction. A theorem for whole \(D_{16}\) is not a
  theorem for \(D_{16}+I\).
* The literal invariant-lattice and full reflection-group containment
  bridge wherever a local proof uses the golden-lattice theorem. A common
  field alone does not supply that bridge.

A simple exact-budget example is the covered domain
\(A_{11}+O_7+A_3\). Its retained local rows are

\[
 e_{A_{11}}\ge13/20,\qquad e_{O_7}\ge9/100,
 \qquad e_{A_3}\ge13/50.
\]

The last bound has the actual \(O_7+A_3\) hypotheses. Their disjoint
supports give total energy at least one exactly. Equality is sufficient
against the strict upper budget; no positive uniform margin is imposed.
Other branches use joint supports or explicit local exclusions, as listed
in the [complete arithmetic table](compact_arithmetic15_global_energy_lp_recovery_20261118.md).

Append the **same compact rank-18–29 catalogue**, using arithmetic implies
compact. No old global \(n\le17\) row is inserted into the arithmetic LP.
The combined catalogue has 917 branches and 167 local contradiction leaves
(the compact 151 plus 16 arithmetic leaves). Its dual has 3,324 variables
and 2,488 rows. Exact checking and fresh regeneration pass. The joined
rank-16–29 exclusion and the published compact ceiling recover \(n\le15\).

## 7. Which statements are inputs and which follow from the LP?

| Statement | Present role |
| --- | --- |
| Geometry of the original quantities and their global budget | Input theorem |
| Spherical/native classification, complete source domains and defect compatibility | Input theorems |
| Local face-payment implications, component/joint energy floors, local contradictions | Input theorems with full applicability domains |
| Finite58 network coverage and feature interpretation | Retained input construction and coverage theorem |
| Arithmetic rank-16/rank-17 pruning and field/data/lattice arguments | Explicit retained inputs |
| Compact cover from component-count classes and complete small remainders | Uniformly generated and source-reviewed cover; old eleven-case list not assumed |
| Compatible charge coefficients, network potentials or energy multipliers | Certificate outputs, exactly checked |
| Exclusion of all compact aggregate and remaining energy branches | Global LP consequence, using the explicitly retained local contradiction leaves |
| Optimum 377 within the declared finite-volume law family | Matching exact feasible and dual certificates |
| Full bounds 57/17/15 | LP certificate plus the retained incidence or energy theorem, complete coverage and stated tails |

A zero multiplier in one certificate only means that certificate did not
need the corresponding row. It does not prove that the row is logically
redundant. A stronger bound on the **same support and domain** does imply
redundancy of its weaker predecessor. For example, a valid prefix-plus-gain
floor dominates its corresponding prefix floor. No complete minimal-row
or independence classification has been computed.

## 8. Dimension thresholds and monotonicity

For compact recovery, freeze the complete catalogue and all rows. At
threshold \(N\), retain branches labelled by actual dimensions \(n\ge N\).
Increasing \(N\) removes branches and changes no theorem or cost. A feasible
normalized primal point at the higher threshold extends by zero masses to
the lower threshold. Thus infeasibility at \(N\) implies infeasibility at
every larger threshold within the catalogue. The published ceiling handles
all larger dimensions outside it.

The existing compact catalogue starts at 18 and the joined arithmetic
catalogue at 16. Neither claims to describe the missing lower dimensions.
Changing only a threshold number cannot justify a bound below those ranges.

For finite volume, the accepted certificate already covers \(n=58\) and
every \(n\ge59\) through its endpoint/tail proof. Restricting that covered
set to \(n\ge N\), \(N\ge58\), preserves the certificate. This is not a
claim that a new generic finite-volume model has been separately exported
and solved at every threshold. Lower thresholds need newly justified
domains and estimates. None of these arguments uses geometric inheritance
of the Coxeter or arithmetic conditions by arbitrary faces.

## 9. Evidence, navigation and deferred strengthening

The [results note](existing_proof_lp_recovery_results_20261118.md) is the
execution ledger for the three cases. The final
[root acceptance](../../audit/existing_proof_lp_recovery_completion_root_acceptance_20261118.json)
binds the old finite-volume completion, four new monitored jobs and the
independent source/execution reviews. Saved compact models are complete
primal/dual exports, catalogues and rational candidates. Fresh verification
regenerates their complete bytes and repeats exact candidate checking.
Numerical status alone is not the proof.

The detailed master design, compact table and arithmetic table were frozen
as reviewed inputs **before** the new jobs ran. Their preparation-era
sentences saying no execution had yet occurred are historical statements
at those source hashes. Preserve those bytes; use this overview, the
[public reading index](README.md) and actual execution records for status.

The exact new commands are saved verbatim in the four monitor records
linked from the root acceptance. Reproduction uses
`scripts/run_bounded_research_check_20261018.py --parent-lock`, fresh output
directories and ledgers, the recorded catalogue source hashes, and stages
`export-solve` then `verify`. The two compact cases use
`scripts/run_existing_proof_energy_lp_20261118.py`; arithmetic supplies both
catalogue modules. Successful unchanged jobs are not rerun for documentation.

For future SDP work, first identify whether a proposed change refines the
necessary geometric domain, expands the available certificate family, or
adds an obligation to a certificate search. These are different operations.
An additional valid necessary condition can shrink a counterexample set;
a richer multiplier cone can supply more certificates; another required
face-payment row can instead restrict the searched charges. Any hierarchy
must preserve the actual-object lift and its threshold nesting.

That future task remains deferred. In particular the incomplete continuous
finite-volume/H1 search is not established equivalent to the working
proof-based LP, and the compact partial-star relaxation is known to admit
abstract feasible points at every threshold. Neither is the completed
proof-recovery model described here. No further strengthening or numerical
experiment accompanies this documentation update.

The completed aggregate models are LPs: no new positive-semidefinite matrix
variables are being optimized in them. Positivity and lattice arguments
remain in their named local theorem inputs. A future SDP reformulation
must make explicit which of those inputs it replaces, strengthens or proves.
