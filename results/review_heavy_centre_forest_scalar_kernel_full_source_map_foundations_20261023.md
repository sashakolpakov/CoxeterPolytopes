# Full outside-author static review: heavy-centre scalar and forest logic

23 October 2026. **Complete mathematical/source/map PASS at the explicit
conditional boundary. No correction requested.** This is a PRECOMPILE review,
not an elaboration, object check or actual axiom audit. No mathematical worker,
import, enumeration, symbolic/numerical process, test, compiler, Lean, TeX,
Sphinx, plan or Git operation ran. Only source reads and bounded byte/hash
metadata were used. Neither reviewed target nor any inherited source changed.

## 1. Exact objects, complete read and independence

| Object | Bytes/lines | SHA256 |
|---|---|---|
| formal/CoxeterProofs/HeavyCentreForestScalarKernel.lean |12329/263|a839befdbd0454d70fd67cc458641469524dce1e69bdb1b4ba89e777c4fa2060|
| docs/heavy_centre_forest_scalar_kernel_logic_20261023.md |10117/172|db2f524ae6cb3beb230c15021d1b8466fd418fdde5cb4cb5ced3c1d2dc28e3c2|

I read the ENTIRE module and map, all182 lines of the new written059e
theorem, its entire2510 different-author report and root acceptance, and
the entire284-line inherited finite-centre source. I also inspected the
finite57 manuscript's actual-corner, all-wall projection and outgoing-edge
interfaces. This does not claim transitive regeneration of all3322 manuscript
lines, its scripts, or the inherited allocation tables.

Root authored the new forest mathematics and specified the two-local-bridge
interface. Replay implemented this module and map and authored the inherited
finite-centre classification. Composition supplied the separate full written
forest review. I previously reviewed the inherited finite-centre source and
contributed related geometric/ordered-ring foundations; I did NOT author the
new forest proof, this implementation, or its strong-induction proof. This is
a different-author review of the NEW formal source and its assumption map,
with those inherited contributions disclosed.

All SIX literal written bindings in the map freshly match current bytes:

| File | SHA256 |
|---|---|
| results/finite57_heavy_centre_forest_and_linear_exit_expansion_candidate_20261023.md |059e3fcf6c80700a0b11778638cbb47e758ca245f6cd9573c747e1ff8caa158f|
| results/review_finite57_heavy_centre_forest_linear_exit_expansion_full_scope_compact_reviewer_20261023.md |2510181660697bbc329f3c34762891d46f8fcc96f4a089aa818f96b809358b06|
| results/review_finite57_heavy_centre_forest_and_linear_exit_expansion_root_acceptance_20261023.md |97ff5c89037fc7e70e6b1393568750e52c95d50a6900d70b136c287ef43cefb5|
| manuscripts/finite_volume57.tex |c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0|
| results/finite57_complete_finite_centre_recipients_and_heavy_component_escape_candidate_20261023.md |994d3b9ed36a32e8698610140ee656b40a88507c7f0a7b6b1978faef21b00186|
| results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md |1c2e887918d4277efffd9c568aed5d9e203a61dce2ca7c98eca96665831792f9|

The preparation-preservation manifest is
results/source_archive/heavy_centre_forest_before_local_chain_induction_20261023/manifest.json,
SHAff4b711b610797450027d14cca79378325e160ef1b3122f8d396ebd64b311a3a.
Both stored gzip hashes match: b988761f69cb5977cc60f9fa7c592c53ba3fedac1da781ef08fa863a1a69ea3a
for the source and cf4b3de1e71cb5bd196ed4aa866742cd06decd3173ab0c46e1cceebb2aeeb348
for the map. Their decompressed hashes independently match828a17b988917fb380abeb952586f9d4c0293e0aea2c9ab9cb2bf31ed9e3bcd4
and e61b5f5c5a16d1d689e3a35c691e8d379b0b423256342d05ff518b1a88f4fce7.
These remain uncompiled predecessors, not inherited formal execution credit.

## 2. Primitive ordered-ring bounds, source lines31--90

The scalar domain is one generic linearly ordered commutative ring R,
with the same explicit Lean.Grind.OrderedRing interface as the already used
D15/A18 scalar pattern. It is not Nat or Int standing in for real lengths.
All scalar variables occur explicitly in theorem signatures. No real square
root, hyperbolic function or geometric object is constructed.

I rederived every difference-of-squares step in cone_primitive_bounds:

* c>=16 gives (c-16)(c+16)>=0 and hence c²>=256.
* If4t>=3, the nonnegative product(4t-3)(4t+3) would equal
  16t²-9=-1, using2t²=1. Thus4t<3.
* If s>=c, (s-c)(s+c)>=0 would equal s²-c²=-1. Thus s<c.
* If16s<=15c, its difference-of-squares product would give
  0<=225c²-256s²=256-31c², inconsistent with c²>=256.

All product-factor signs follow from the stated primitive data, including
c>0 and s>=0. No radical estimate or final cone gap is assumed. For
long_edge_cone_gaps, the two explicit nonnegative/positive products give

    4(c-st)-c = c(3-4t)+4(c-s)t > 0,
    16(s-ct)-3c = (16s-15c)+4c(3-4t) > 0.

Together with c>=16 these imply c-st>4 and s-ct>3, including strictness
at the endpoint c16. This is genuine derived scalar content.

## 3. Cosine and cleared cotangent rules, source lines93--161

cosine_rule_growth does not multiply by an uncontrolled-sign cosine.
With0<=sR<=cR, t>=0 and cosB<=t, the proof derives

    cR*t-sR*cosB=(cR-sR)t+sR(t-cosB)>=0.

Its exact decomposition is

    cD-4cR=cR(c-st-4)+s(cR*t-sR*cosB)>0.

The first term is strictly positive because cR>0 and the gap is derived;
the second is nonnegative. This remains valid for negative cosB. The
proof uses the supplied exact cosine equation and no hidden lower bound
on cosB, division by it, or unjustified reversal of an inequality.

For cotangent_rule_growth, the numerator decomposes as

    cothR*s-c*cosB=(s-ct)+(cothR-1)s+c(t-cosB)>3.

If cotC<=3, then(3-cotC)sinB>=0 and sinB<=1 imply
sinB*cotC<=3sinB<=3, contradicting the cleared cotangent equation.
Only0<sinB<=1 is needed; no inverse or division operation is used.
The theorem's positive sine input belongs to the geometric triangle
boundary, not to an unstated analytic instance. triangle_cone_step calls
these two theorems with the SAME primitive c,s,t and angle data, so it
does not combine incompatible long edges or triangles.

## 4. Two-local-bridge induction, source lines166--206

The two exact external interfaces are, for every1<=i<n,

    cone i -> path(i+1) != path0,
    cone i -> path(i+1) != path0 -> cone(i+1).

They have the correct order: tangent-angle control first excludes the
degenerate return; only then does the nondegenerate triangle propagate
the strict cone. Neither premise asserts no return for every path index
unconditionally, nor is a global noCycle assumption supplied.

I checked the complete strong-induction proof. At k1 it uses the two
base facts. Otherwise lower1<=k and k!=1 imply k>=2. Hence k-1 is
positive and strictly smaller, and k<=n implies k-1<n. The induction
hypothesis supplies cone(k-1); the first local bridge supplies the next
nonreturn, and the second supplies the next cone. The arithmetic identity
(k-1)+1=k then gives the stated conjunction. No predecessor no-return
conclusion is silently substituted for the new one. The final declaration
specializes to k=n using its explicit n>=1 hypothesis.

The invariant's _length_positive argument is logically redundant for a
possibly empty index interval, but harmless. Its final no-return corollary
does require and use n>=1. V and path are arbitrary, and the Nat values
are ONLY path indices. Quantitative4^(m-1) iteration and the geometric
realization of cone are not formalized by these declarations; the map
explicitly retains them at the written boundary.

## 5. Actual-centre normalization and incidence, source lines209--249

projected_centre_length_lower starts from delta>0,
(rank+1)delta1, shadow>=1+delta, wall<=0 and
delta*coshL=shadow-wall. Its identity

    delta(rank+2-coshL)=1+delta-shadow+wall<=0

forces coshL>=rank+2>=16. The rank here is an ordered-ring variable,
not a count used as a fake metric. In the written actual model it is
s14 or20, delta1/(s+1), shadow h^tB^-1h', wall the independently
norm-two actual a.a'. Those interpretations, the full support and
actual edge orientation are supplied outside Lean. The normalized
equation has the correct sign and factor: delta*coshL=shadow-wall.

All three Nat declarations are genuine linear incidence consequences:
3H=2E+boundary and E+components=H imply
boundary=H+2components; components1 gives H+2; the ACTUAL bound
boundary<=3receivers then gives H+2components<=3receivers. These
are counts only. The final inequality is correctly denominator-cleared,
not an implementation of ceiling notation, a receiver-capacity estimate
or a Hall theorem. It does not assert that all exiting targets are finite.

The forest Euler identity is explicitly an EXTERNAL input to these
counting declarations. The no-return declarations do not themselves
construct a graph, prove that a graph cycle produces an eligible path,
or establish its Euler identity. The map states this separation correctly.
There is consequently no claim of a completed geometric forest instance
merely because no global noCycle parameter occurs in the chain theorem.

## 6. Code/elaboration-risk inspection and exact formal scope

The module imports only Std and declares exactly ELEVEN named theorems,
with exactly the same eleven fully qualified footer requests. There is
no axiom declaration, sorry, admit, native_decide, external computation
oracle, old-kernel import or hidden geometric instance. The last three
declarations use omega only on Nat linear arithmetic. Scalar products
are introduced with OrderedRing.mul_nonneg/mul_pos, and the remaining
goals are polynomial normalization plus ordered linear consequences.

I inspected the pinnedv4.32.2 source-only APIs: the mul_nonneg/mul_pos
signatures in Init/Grind/Ordered/Ring.lean match these calls, and
Nat.strongRecOn in Init/WF.lean has the exact ind-case recursion shape
used here. The already verified ResidualBellmanCertificate uses the
same induction syntax. All grind-only calls have no redundant explicit
local theorem parameters. Names, binder order, predecessor use and
the simpa-only successor rewrite reveal no concrete static blocker.
This is not a prediction that substitutes for ACTUAL elaboration:
solver elaboration/typeclass/tactic performance remain untested here.

The map's geometric boundary is faithful to059e/2510 and the finite
centre input: original finite-volume57 polytope, whole2D17+D23 heavy
vertices, actual centre edges, complete heavy completion, Schur/shadow
identities, actual orthogonal tangents, hyperbolic trigonometry and angle
interpretation, two local bridges, graph cycle-to-path realization,
forest Euler identity and finite-target reverse ports all remain explicit.
Hash authentication of those sources does not create Lean instances.

For a subset of heavy vertices, its boundary may reach another heavy
vertex outside the subset. Such incidences do not furnish fresh debt-free
capacity. No ideal multiplicity, eight-unit financing, matching, universal
364 cap, dimension56, compactness/example or global ceiling is concluded.
All accepted allocations and the six frozen manuscripts stay unchanged.

## 7. Verdict and execution gate

FULL precompile source/math/assumption-map PASS at the exact hashes in
Section1; no edit requested. All eleven declarations and the entire map
were checked. This review certifies no compiler output, object or axiom
list: requested11, ACTUAL ACCEPTED OUTPUTS ZERO at this review boundary.
Only root may admit one capped module compilation and a SEPARATE actual
axiom audit after its resource/preclear checks. A refusal or failure earns
no formal credit. Even a future successful run will verify these conditional
ordered-ring, local-induction and incidence statements, not the external
real geometry, classification, actual graph/port instances or dimension56.
