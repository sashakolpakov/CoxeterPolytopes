# New-dimension search: a new field, noncentral sections, and coupled polar caps

This round has not produced a new compact or finite-volume Coxeter dimension. It changes the cubic arithmetic field, the positive lattice blocks, the coefficients of H20 section normals, and the truncation equations of H9 lifts. Exact volume and actual-ridge checks determine the negative results; numerical failures to find an escape are left undecided.

## Arithmetic H9: conductor nine and a changed conductor-seven lattice

The [conductor-nine search](fresh_cubic9_dihedral_report_20261013.md) adds exact arithmetic in `Q(2cos(π/9))`, a complete primitive root-norm alphabet, and the full zero-height chamber `I₂(18)^4×A1`. Four H9 hosts were searched through height 300. Their accepted prefixes have **472, 204, 410 and 204 walls**. Every prefix has an independently checked strict spacelike escape, proving infinite volume. Two H8 variants, with 182 and 148 walls, fail the same test.

Independent unpruned enumeration repeats all spatial candidates and wall decisions through height 30, and all time/norm controls through height 300. A separate arithmetic implementation checks the saved larger wall sets and geometry. The same code exactly certifies the known compact `(2,4,18)` triangle as a calibration. The five arithmetic/certificate regression tests pass. No statement of full-host nonreflectivity follows from these bounded searches.

The [changed conductor-seven host](fresh_cubic7_mixed_blocks_report_20261013.md) replaces one binary block and a scalar by A3: `L=<-2τ>+B³+A3`. Determinant square classes distinguish this quadratic space from six previous saved hosts. Its zero chamber is `I₂(14)³×B3`, and complete enumeration through height 100 gives **67 walls with an exact infinite-volume escape**. A standalone Fraction verifier independently repeats the full search through that height. The lattice's integral overclosure `M=<-2τ>+B³+I3` has the same orthogonal and reflection groups; those two presentations are therefore counted once.

## H20: arbitrary real coefficients on elliptic three-wall supports

The [new section theorem](fresh_h20_bivariate_sections_20261013.md) excludes every interior H20 section of the frozen Borcherds H21 chamber whose normal lies in an **elliptic span of at most three original walls**, for arbitrary real coefficients. Unlike the previous central-section theorem, the section need not contain the old center. The proof covers 266,000 unordered triples in 15 symmetry orbits and 45 sign families. Independent replay verifies 58 actual-ridge certificates and the continuous polynomial sign arguments.

A constructed quartic coefficient repairs two actual angles to π/2 and π/8, but a third actual ridge has an acute angle strictly between π/3 and π/2. Thus it remains non-Coxeter. A [separate bounded search on nonelliptic supports](fresh_h20_bivariate_nonelliptic_20261013.md) finds two finite-volume quadratic sections with simultaneous π/2 and π/3 repairs. Both retain exactly verified obtuse actual ridges. This latter batch explores twelve declared parameter lines; it does not exclude all affine or Lorentzian supports.

## H9: a genuine π/10 polar-cap seed and compatible completions

Three alternative H8 lifts release different inherited walls before solving their cap equations. One admits a valid polar truncation meeting an old wall at π/10. A further cap over `Q(φ,√(49−24φ))` meets that polar at π/10. The [independent seed audit](fresh_h8_alternate_first_cap_independent_audit_20261013.json) checks the resulting twelve-facet H9 polyhedron: its Gram matrix has rank ten and inertia `(9,1,2)`, every pair is Coxeter or disjoint, and every declared wall has a strict facet witness. It is a valid Coxeter seed within the finite bank ruled out below.

Adaptive cap equations yield 144 coupled caps in addition to 348 base caps and the polar. Numerical escape search on their intersection was inconclusive. Pairwise compatibility resolves the issue: **no selection from this 493-cap bank, retaining the ten seed walls, can be a finite-volume Coxeter polyhedron**.

Five exact spacelike witnesses have complete cutter lists whose union contains only 54 caps. A separate [cyclotomic audit](fresh_h8_alternate_core_angles_independent_20261013.json) classifies all 1,431 pairs among them against every possible finite Coxeter label in their degree-at-most-eight fields. There are no omitted compatible finite-angle or dotted pairs. Allowing both parallel pairs as well, a [standalone cover replay](fresh_h8_alternate_core_cover_independent_20261013.json) proves that no compatible selection can cut all five witnesses. It checks the cutter lists against all 493 normals and records a 17-state Boolean proof. This closes the saved finite bank even for finite volume; it does not exclude cap equations outside that bank. The [full construction report](fresh_h8_alternate_20261013.md) and [independent obstruction review](fresh_h8_alternate_core_independent_report_20261013.md) record the seeds, equations, scopes and replays.

The independent source and certificate audits are saved beside their respective constructions. Later mirrors, different arithmetic hosts, larger section supports and cap equations outside the declared finite banks remain possible routes. None of the bounded failures is promoted to a universal dimension obstruction.
