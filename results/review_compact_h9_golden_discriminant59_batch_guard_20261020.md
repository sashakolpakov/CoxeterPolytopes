# Independent paper and source review of the discriminant-59 route

20 October 2026. **PASS within the stated local construction scope.**
This review used source reading and independent hand algebra only. I did
not import or execute the worker, generate or replay a certificate, run a
search, compile Lean, or build TeX. No compact H9 example or full-host
reflectivity conclusion follows from this review.

Reviewed final sources:

* `results/compact_h9_golden_discriminant59_edge_route_20261020.md`,
  SHA256 `725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`;
* `scripts/check_compact_h9_golden_discriminant59_edges_20261020.py`,
  SHA256 `01a569b66b77c38056255a94dce765e4784d809dfe86735ac5bef869755edf3d`;
* the preparation note of the same date, the pinned Bugaenko fixture and
  the pure quadratic-field/matrix helper used by the worker.

At the parent's request I made two mathematical clarifications before
these final pins: the arbitrary-arrangement angle list includes 2-phi
and 3-phi in addition to the five nonzero fundamental-wall values; and
Section 4 explicitly obtains X=x-y in O and Z in O from E8 dual
coordinates. The worker's note pin and separate arrangement divisibility
list, and the preparation hashes, were updated accordingly. This review
is independent of execution and does not describe those edits as an
independent replay.

## Lattice, complete zero chamber, and actual base facet

I checked det(B)=-5-12phi, delta=phi*b, N(b)=59, the physical Lorentz
signature and positive conjugate binary block. The valuation-one
determinant precludes a proper integral overlattice. The primitive
reflector divisibility argument uses the actual even lattice, its
adjugate, coordinate primitivity and the stated PID/unit facts; it gives
the necessary norm representatives 2 and 2b modulo unit squares. It is
not a bounded root-shell inference.

The long-short and long-long exclusions work for arbitrary arrangement
mirrors after the corrected finite-angle list. For long-short pairs the
nonzero squared-cosine values have norms divisible only by 2,3,5. For
long-long pairs reduction in the one-dimensional radical modulo (b)
forces the squared-cosine value to be 4 modulo (b); none of the finite
values permits this. Parallel mirrors would give a nonzero K-rational
null vector, excluded by the positive conjugate form. Thus the entire
long-mirror separation argument survives the correction.

The normal r has Gr=(0,-b;0), so its reflection and its y=0 orthogonal
lattice are exactly as claimed. The ambient mirrors crossing it are
short and orthogonal, and the intrinsic norm-two reflectors extend
integrally. The identification with the full compact Bugaenko facet
uses the inherited complete arrangement/transcript, not the sixteen-wall
seed. At p the zero-root energy argument excludes every nonzero binary
part of a short root and every nonzero E8 part of a long root. The
remaining system is exactly ordinary E8 plus +/-r. This justifies the
actual full-chamber germ and the inherited 41 finite vertices.

## Both all-height outgoing-edge arguments

For the first transverse edge I rederived the unit determinant of (r,c)
and its two-by-two Gram. A first finite wall along the retained E8 edge
must be binary, since E8 has no connected finite rank-nine extension.
For an interior cut its coefficients lambda and mu have opposite signs.
The physical and positive-conjugate bounds on lambda, or on gamma after
mu=b*gamma in the long case, give the claimed nonzero integral norm
strictly below one. This excludes all intervening mirrors, independently
of height. The endpoint h is finite and c is an actual first wall.

For the next edge I independently reconstructed the E7+c finite
extension inventory: short isolated; short joined to c with label 3 or
5; the unique E7-to-E8 extension; or long orthogonal to every retained
wall. No finite rank-nine connected component contains E7 and c.
The root decomposition gives X=x-y integral. Its E8 dual coordinate is
Z/2 or (Z-3)/2, so the norm arguments really concern algebraic integers.
The five equations have R=4,3,3-phi,1,4b, and the collision parameter is
delta*X/Z. Contact by tau_star implies exactly the two displayed
embedding bounds.

For R=1 and 3-phi their product is less than one. For R=3,4 it makes X
a unit; the two strict embedding bounds and the golden unit group force
X=1. The R=4 square has rational norm 45, whereas R=3 gives the unique
positive Z=2phi^2 and precisely u. In the long case E8 reflection
integrality forces Z/(2b) integral and nonzero, with norm below one.
Consequently the second first-simple-wall assertion is complete and has
no unexamined Pell height range. Old-cap compatibility is checked only
after this selection, so it cannot improperly discard a first wall.

## Pentagon, reversal, and exact worker coverage

The E7 section has the claimed five cyclic sides and angles. Its four
already actual vertices and its positive-definite E7+I2(5) corner are
distinct. Projection onto the latter actual wall set is justified by
its inverse-positive acute Gram and strict feasibility against every
other actual wall. Convexity then proves equality with the compact
pentagon, rather than merely finding five feasible prefix vertices.
Exactly three of these vertices lie off the inherited r facet, giving
the lower bound 44.

I checked the reversal columns against the literal native node order,
the actual half-glue y=((phi-1)x-r)/2, and S^2=1. This is a chamber
involution, not a period. The r/c swap is blocked by the nonsquare b;
the displayed c/u swap is blocked by reflection integrality. The
nonorthogonal short mirrors of Section 7 are genuine norm-two integral
reflectors, so omitting them from an orthogonal c-lattice calculation
would be invalid.

The worker source reconstructs all sixteen normals and S without using
a saved result. It checks both embeddings, every reflection-matrix
entry, all 120 pairs, the exact S-fixed strict center, and 44 distinct
finite points by nine-wall positive Gram projection and all sixteen
inequalities. Its field arithmetic and sign tests are rational and
exact. The quadratic-profile interpolation evaluates six points that
determine every coefficient of a known degree-at-most-two polynomial;
this is an identity check, not sample evidence. Every five-profile
bound, the unique A2 root and the norm-45 nonsquare are tested directly.
The helper is imported only for definitions, and no helper producer or
model is called. Source hashes are streamed before and after; generation
refuses replacement and replay recomputes the canonical result.

The program deliberately does not certify the number-theoretic norm
pool, finite-diagram inventory, golden unit group, full Bugaenko
arrangement, or universal acute-projection geometry. Those are explicit
written or inherited inputs. Its 44-point verification is a lower
bound, not a complete vertex census or a finite-volume certificate.
The finite two-edge and pentagon assertions pass; full L59 reflectivity
and compactness remain open.
