# Root static review: D18 centre scalar kernel

21 October 2026. **Source and mathematical-scope PASS; not compilation.**
No scientific worker, import, test, compiler, Lean or TeX process ran.
Root proposed the integer-square formulation; this review is not a
claim of independent authorship from that mathematical design. A
separate outside-author source review and actual compile/audit records
are required before execution credit.

Reviewed module `formal/CoxeterProofs/D18CentreScalarKernel.lean`, SHA-256
`998ae907c17886fd0c72a7bce74ce36126d128b6e51453966982025410cc9eb6`,
and map `docs/d18_centre_scalar_kernel_logic_20261021.md`, SHA-256
`dc257720eb80b1fef4f6562698a09524dfc621deea13600920fe2b5e3f6cf846`.
There are five requested declarations and five requested axiom prints,
with zero actual outputs at this preparation stage. Only Std is imported.

I expanded the exact identity by hand. Four times the leading
coefficient is387499200, twice it is193749600, and the completed
constant is28961120000. All three coefficients agree with
q(z)=96874800*z squared-1528000*z+6100. The square is nonnegative
by the ordered-ring laws; the strictly positive residual gives q>0
without a supplied application positivity premise.

I inspected the pinned Lean4.32.2 source for the two named APIs:
`OrderedRing.sq_nonneg` has an implicit element argument and
`mul_nonpos_of_nonneg_of_nonpos` takes the displayed two signs.
This static lookup supports the source formulation but is not a
kernel check or a claim that the tactics have executed successfully.

The supplied primitive analytic bound minus34190400 is exactly q.
Its strict positivity therefore derives10000*F>341904. In K10,
the primitive base and height floors clear to4901074<4900000*T;
in K11 they clear to4901534<4900000*T. Both imply T>1.
Those integer arithmetic identities and their strict margins check.

The map correctly assigns the radical chord, tangent inequality,
real parameter substitution, actual negative core, all-wall feasibility,
full cap classification, positive-reference payments and the original
norm budget to written application inputs. The final scalar lower
bound and T>1 are derived, not assumed. The weak5 references satisfy
the same base floors without reusing I energy. No geometric exclusion,
number-field theorem, whole dimension bound or complete D18 manuscript
is claimed as formalized by these five conditional scalar statements.

The source appears suitable for the next selected Lean batch AFTER
the promised pending manuscript replay slot and only when the unchanged
resource gate admits it. It has no actual compilation or axiom-audit
credit now. No prior successful module or dependency needs rerunning.
