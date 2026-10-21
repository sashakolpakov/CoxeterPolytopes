# Whole A18: exact necessary silver weak-cap pair graph

21 October2026. New root algorithm specification, pending outside-author
review and execution. No numerical output is asserted in this version.
This is supporting mathematics in the existing arithmetic paper family.

## Domain and geometric necessity

Use the accepted actual-field reduction fb185ba3 and complete quadratic
.52-weak inventory47100252, at WHOLE A18 in dimension18 only. The actual
field is K=Q(sqrt2). Write s=sqrt2, and raw old roots of square2 have
inverse Gram Ginv(i,j)=min(i,j)(19-max(i,j))/19. All coordinates below
are in their raw dual basis Omega_i. The COMPLETE necessary weak domain
is the following44 models; no model is asserted to be actually weak.

*20 far-helper3 models: i5..14, j1 or18, projection (1+s)Omega_i+Omega_j,
 raw norm2.
*4 near-own4/helper3 models: (i,j)=(5,6),(6,7),(14,13),(13,12),
 projection s Omega_i+Omega_j, raw norm2.
*20 isolated dotted models: beta=3+2s or6+4s, i5..14,
 projection beta Omega_i, raw norm2beta.

For each cap let P be its displayed projection, N its raw norm and
F=P.P-N its positive physical time square. The other embedding has
F'<0, while N and N' are positive. These are checked for every model.
These44 include every .52-weak actual opposite. An actual opposite is
unique for each outgoing edge, so two different models at the SAME
index cannot occur as two weak opposites. Different indices are tested
as follows; failure of a necessary test removes that pair.

All cap residuals occupy the SAME one-dimensional K-valued orthogonal
complement of the complete old basis. The physical time signs agree:
the old finite vertex satisfies each nonincident actual wall strictly.
Thus, if caps a,b occur together, q=sqrt(F_b/F_a) is in K, with its
positive physical sign. Their time product is F_a q>0, and their raw
mutual Gram product is R=P_a.P_b-F_a q. No arbitrary square root is
adjoined, and square norm ratios alone do not establish this condition.

All-wall acuteness requires R<=0. First-type Cartan cyclic integrality
also requires R in O_K=Z[s]. To see that it is a cycle, select one old
native helper contact for each nonisolated cap and its marked contact
for each isolated cap, and join the selected old nodes along the A18
chain. From the first cap to its old contact the Cartan entry is -1;
the last old node to the second cap has entry -1 for a norm2 cap or
-beta_b for a norm2beta_b cap. The returning entry is respectively R
or R/beta_b. All intermediate native3 entries are -1. Their product
is therefore +R or -R, and is integral. This argument covers mixed
norms as well as both equal-norm cases.

Independently, the TWO-cycle requires B=4R^2/(N_a N_b) in Z[s].
R-integrality alone does not guarantee this for two field-valued norms.
The worker imposes both necessary tests, including in the dotted case.
If R=0 the pair is orthogonal. Otherwise B>0. If B>4 it passes
the necessary dotted test. B=4 is parallel and impossible for a compact
chamber: the sum of the two outward UNIT normals is a nonzero null
vector feasible against every wall, so the chamber has an ideal point.
It is future since its product with an interior point is negative.
If B<4, the
walls intersect and must have an actual native angle. Since B is in K,
the possible finite angles are exactly m3,4,6,8, with B=1,2,3,2+s.
This is the usual degree<=2 real cyclotomic restriction: m5,10,12 have
different quadratic fields; m2 is the already separated zero case.
The positive conjugate2-s is not a native angle pi/m. This test is
necessary only. Passing does not establish an ambient chamber.

## Exact algorithm and proof

The worker uses pairs of Fractions for a+bs. Addition/multiplication,
conjugation and inverse via a^2-2b^2 are literal field operations.
For opposite-sign a,b, sign(a+bs)=sign(a)*sign(a^2-2b^2);
same-sign and zero cases are immediate. No floating-point comparison
or algebraic numerical tolerance occurs.

For a positive physical x=a+bs, a K-square root u+vs must satisfy
(u^2-2v^2)^2=a^2-2b^2. A necessary rational nonnegative square root r
of this norm is found by integer-square tests on numerator/denominator.
For BOTH signs of r, test u^2=(a+r)/2 for a rational square. If u!=0,
v=b/(2u). If u=0, a root is possible only when b=0 and a/2 is a
rational square, giving v. Check the candidate square exactly, then
choose its positive physical sign. Necessity of the equations and the
two possible norm signs proves completeness; direct squaring proves
soundness. The tests include pure rational and pure-s multiples.

Generate all44 models in the displayed finite loops. Enumerate ALL
946 unordered pairs; each gets one exhaustive first-failure reason or
a compatible result with the exact square-root, time product, Gram
and angle witness. No omitted pair, heuristic cutoff, optimizer or
height relaxation is used. Finite loops and terminating rational
integer-square arithmetic prove termination. Every actual weak pair
passes by the necessity arguments above; this is the precise coverage
conclusion, and does not assert every surviving pair occurs.

Connected components are computed by finite graph traversal, including
isolated vertices. No clique enumeration or new mathematical inference
from connectedness is included in this batch. At least one weak model
exists by the accepted .52 reference. Actual weak models must form a
clique with distinct indices; these observations are left for the
subsequent geometric argument.

The independent sanity witness far(5,18),far(7,1),t=1+s has
F5=(200+150s)/19, F7=(256+192s)/19, ratio32/25=(4s/5)^2,
time product(240+160s)/19 and R=-2-s. It MUST survive as dotted.
The worker asserts this witness rather than misreporting an empty graph.

## Execution boundary

Run only through the existing one-thread59s/1GiB parent-lock watchdog,
with40GiB disk/6GiB RAM start and8GiB/3GiB abort gates. The output
contains the complete model and pair list, compact graph, input/source
pins and explicit necessary-only scope. Existing successful workers,
certificates, six manuscript sources/PDFs and Lean objects are untouched.
No global ceiling changes without the remaining actual geometric proof.
