# Full static review: Lorentz obstruction kernel

21 October2026. Reviewer: sharp_ra_foundations_audit. **PASS for the
mathematical statements, source proof structure and corrected assumption
map. Compilation and actual axiom outputs remain pending.** No Lean,
compiler, test, scientific import or mathematical job was run.

Reviewed in full:

| File | SHA256 | Bytes |
|---|---|---:|
| formal/CoxeterProofs/LorentzObstructionKernel.lean | `8713ee5909934f8b9263ff8cb77c2a0fe7752de603d36edd0168090524964690` |4523|
| docs/lorentz_obstruction_kernel_logic_20261021.md | `9a51f5b4af9cb82b558962aeae79d66ebdbeb00615a8d64ad86adffcc387ffaa` |3200|

I did not author the Lean module. I authored the earlier restricted
marked-D geometric note linked by the map and contributed to inherited
Lorentz/actual-foot foundations. This is an independent review of the new
algebraic kernel, not a claim to independently author-review that earlier
geometric classification or the pending full norm79 application.

## All seven declarations

1. **second_orbit_gap** (lines28–43). The assumptions are precisely c>0,
   q>=0, W=q+c, and `2q²<=W*(W+second+q)`. They imply W>0. If
   second<=-3c, the last bracket is<=2q-2c. Multiplication by W>=0
   and `(q+c)*(2q-2c)=2q²-2c²` contradict the Cauchy assumption,
   because c²>0. The source supplies the positive-square and multiplied
   inequalities explicitly; it does not assume a strict cutoff.
2. **norm_two_second_orbit_gap** (lines45–52). Substitution c2 in the
   preceding theorem gives exactly second>-6. Positivity of2 is valid
   in the generic ordered-ring context.
3. **golden_pair_beyond_cutoff** (lines54–57). The hypotheses phi>1 and
   second=-80-132phi imply second<-212<-6. The proof does not need or
   pretend to formalize the golden polynomial or an actual root product.
4. **no_finite_orbit_data** (lines60–69). A proposed finiteOrder supplies
   q,W and the primitive Cauchy data through the explicit bridge. The
   previous two strict inequalities contradict each other. The bridge
   does not contain the desired no-finite-order conclusion or the
   second-orbit cutoff. It is a substantive external geometric theorem
   parameter; no custom axiom is declared.
5. **schur_numerator_identity** (lines72–75). Direct ring expansion gives
   the stated polynomial identity. It is the determinant of the cleared
   two-by-two sum block: multiplying that block's determinant by
   (m+1)² yields (m+1) times the printed numerator. There is no division
   or unmentioned nonzero denominator inside this Lean statement.
6. **large_d_schur_numerator_nonzero** (lines77–94). At g0 the numerator
   is -2h-4<0. At g>=1, m>=3 gives mg>=m and
   `(m+1)g-2>=m-1>=2`. Thus
   `h*((m+1)g-2)+2mg-4>=2m-4>=2>0`. The proof explicitly constructs
   both required multiplication inequalities and the rearrangement.
   Only h>=0 is needed; a stronger native gap for h is unnecessary.
7. **full_white_rank_contradiction** (lines96–98). For natural numbers,
   rank<=n+1 and rank=n+2 are inconsistent. This checks the final rank
   arithmetic, while actual ambient/Schur rank interpretation is an
   input rather than a formalized Gram-matrix theorem.

All seven are true over the stated generic linearly ordered commutative
ring, apart from the last purely natural-number theorem. There is no
integer/finite-field proxy for real geometry. The hypotheses contain no
assumed dimension bound, chamber stabilizer or lattice exclusion.

## Mathematical mapping and external boundaries

I reconstructed the finite-order Cauchy application rather than relying
on its label. If a future Lorentz isometry F has finite order, averaging
a future timelike orbit supplies a fixed future unit vector y. For a
root z of square c, write z=a*y+w, q=a² and W=w²=q+c. If z.Fz=0,
then w.Fw=q and w.F²w=second+q. Euclidean Cauchy on w and
Fw+F^(-1)w gives

    4q² <= W*(2W+2(second+q)),

equivalently the input inequality. Averaging, the existence of this
Euclidean orthogonal decomposition and this Cauchy realization are
not proved by the module. The map states that boundary explicitly.
The norm79 written source's proposed two-step product fits the kernel
only after its separate actual root/isometry and full-germ premises
are established. Original-lattice integrality, arrangement preservation,
full chamber stabilization and the index/covolume consequence are all
outside this kernel, as the map correctly says.

For the marked-D application I compared the kernel with the linked
restricted geometric theorem. Its difference eigenvalues2+g,2+h,
white-block inverse, exact carried white support, actual long-end/center
roles, denominator positivity and the resulting full Gram rank are
external. The module checks the central polynomial/nonvanishing and
the final natural-number contradiction. It does not formalize those
geometric interfaces, a universal finite57 transfer, or any improvement
of the finite charge372.

## Authorized precision correction and preservation

The original map a8242e69 called W=q+c the Euclidean "norm". Since W
is w², I reported this precision issue and root authorized the sole
insertion "squared" at line12. The Lean source and all statements
remain unchanged. The original3192-byte map was preserved before the
edit, gzip-decoded and compared byte-for-byte with the original:

* archive: `results/source_archive/lorentz_obstruction_logic_before_squared_norm_precision_20261021.md.gz`,
  SHA256 `5e883cd0cff4925d2f613d8429643f638eafaa61cb4e0bff61c4010f3c159d5c`,1653 bytes;
* decoded original SHA256:
  `a8242e69f873c98536b91e00a26d063bf34535a2a2325b479fddabe80775495c`;
* preservation receipt:
  `results/lorentz_obstruction_squared_norm_precision_preservation_20261021.json`.

The exact unified diff has only that insertion. It corrects terminology,
not a Lean or mathematical premise.

## Static proof and execution status

The only import is Std. The OrderedRing multiplication API, local cast
instances and grind/omega usage match the existing source patterns read
in EndpointTriangleNormKernel and the tail/payment kernels. The two
`grind only` identities are universal polynomial normalizations. The
order proofs expose the needed positive products before their final
grind steps; the rw at lines41 and93 matches the displayed expressions.
I found no concrete tactic/API or elaboration defect. Tactic success
still requires the separate bounded compile; static inspection gives
no execution credit or axiom list.

The source has exactly seven theorem declarations and seven fully
qualified #print axioms commands, one for each theorem above. It
contains no sorry, admitted axiom or complete geometric theorem. Both
source and corrected map truthfully say compilation and the separate
actual axiom audit are pending. That pending status must be replaced
by actual terminal evidence before any verified-output claim.
