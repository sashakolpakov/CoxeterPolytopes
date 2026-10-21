# Independent review of the symmetric-arc common escape

**PASS:** every coherent cap retaining the 133-wall H20 chamber and having
symmetric finite or parallel contacts on a consecutive three-wall chain
arc misses the saved simple escape. All 18 support contacts are arbitrary
Coxeter contacts. Any finite mutually compatible collection of such caps
therefore leaves infinite volume. Asymmetric or sparser nondotted contacts
are outside this theorem.

The [independent checker](../scripts/fresh_h20_symmetric_arc_escape_gate_20261016_review.py)
and [certificate](fresh_h20_symmetric_arc_escape_gate_20261016_review.json)
import no producer arithmetic. They reconstruct the rational inverse support
Gram with SymPy, symbolically derive the three-contact norm identity, and
use symbolic quadratic arithmetic for the complete finite support-profile
census. Both the [all-support theorem](fresh_h20_symmetric_arc_escape_gate_20261016.md)
and the earlier [single-support theorem](fresh_h20_dotted_arc_family_gate_20261016.md)
are covered. The independently checked chain nonabelianity and the
[134-wall cap's strict escape](fresh_h20_dotted_arc_cap_20261016_review.md)
are explicit dependencies.

For side magnitude `B` and middle magnitude `A`, the symbolic reconstruction
gives `y^2=A^2/2-B^2/(2C)`, with `C=a^2/(a^2+2)`. If `B<=A`, every chain
pairing would be at least minus two. Acuteness and the Coxeter condition
would make the entire chain nondotted, contrary to the reviewed chain
obstruction. Hence `B>A`, `B>=1`, and `y` is future timelike.

If every support contact were finite or parallel, its rational-Gram cost
would belong to the maximal abelian extension of the rationals, forcing
`C` into that field. Thus some support contact is dotted. The norm equation
and `a^2>280` give the all-label support budget `d^T J d<281/70`.
These are algebraic reductions, not bounded label searches.

The independently reconstructed A10 inverse is
`Jij=min(i,j)*(11-max(i,j))/11`. A dotted A10 contact can only be at one of
the two endpoints; its cost leaves less than one half for all other
coordinates, forcing them to vanish. The two continuous endpoint inequalities
are checked directly. The review also verifies all eight margins of the
earlier single-support result, including its isolated-root cases.

Otherwise a dotted contact outside A10 costs more than two, leaving an
A10 cost below `141/70`. Exact inverse-diagonal bounds exclude all interior
roots except the second endpoints, all endpoint labels at least five,
and second-endpoint labels at least four. The remaining four finite
coordinate domains have exactly 36 profiles. Independent symbolic
enumeration in reverse coordinate order leaves the same eight profiles
and rejects the same 28. Every cost and escape pairing matches exactly.

The other support blocks have inverse-entry sum 13. Ordinary positive
definite Cauchy–Schwarz bounds their escape contribution by less than
`1/1000`. Reverse Lorentz Cauchy–Schwarz bounds the frame contribution by
`-sqrt((9999/10000)u)`, where `u=-y^2`. For each of the eight finite profiles,
the review verifies the strictly positive squared margin at
`u=max(A10_cost,1/a_hi^2)`. Larger `u` only strengthens it. The two remaining
dotted endpoint cases use their separately checked affine-in-`u` margins.
This covers every support contact pattern under the stated symmetric arc
hypothesis, including unbounded finite labels and parallel endpoints.

Generation and fresh saved replay passed, as did source hashes and local
links. This strengthens a construction-family exclusion; it establishes
neither a general H20 nonexistence theorem nor a new dimension bound.

```sh
python3 scripts/fresh_h20_symmetric_arc_escape_gate_20261016_review.py
python3 scripts/fresh_h20_symmetric_arc_escape_gate_20261016_review.py --verify
```
