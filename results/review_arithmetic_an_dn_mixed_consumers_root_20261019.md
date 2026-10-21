# Root review of the mixed A/D consumer lemma

19 October2026. **PASS for the stated whole-vertex, eight-field scope.**
This reviews `arithmetic_an_dn_mixed_end_refund_consumers_20261019.md`,
SHA256 `8714cb88ed454130b975d2441ffdf72b6518959869c65632625215c17e74bf0e`.
No scientific computation, enumeration, Lean or build was run. The
standalone A and D consumer proofs remain separate prerequisites; the
A proof has root and independent reviews, and the D proof has the new
independent full review. This note checks the missing composition
interface, not a global rank20 upper or a dimension theorem.

## Exhaustion of genuinely mixed faces

The white graph is A1+A_(n-3). Relative to its isolated node e, every
side participating in a selected corner is I (no white contact), L
(simple contact with e only), or F (simple contacts with e and one
long-path endpoint). Selected D corners are F/I with bond3, selected
A corners are F/L with product zero. Both endpoint placements of F
are retained. Two F roots cannot be adjacent in a finite corner:
their white contacts already produce a cycle. Consequently a mixed
face has one F, or two opposite F roots.

With one F the opposite side must have original-zero corners against
I and L, so it has no e contact and only the standard isolated or
single long-path continuation. If it were itself F it would be the
two-F case. With two F roots the remaining F/I corner is either zero
or simple; a larger label would be a nonfinite large branch. F/L
products must be zero because a nonzero product closes a cycle.
These are exactly the three cases in the source, including all
simultaneous selected corners. The connected-white D-whisker role
does not occur in a mixed face.

## Independent arithmetic checks

For opposite I and L walls let their product be -beta. The projected
squares1 and3/4 first give beta>sqrt(3)/2. The actual four-cycle
puts2beta in the integral cyclic ring. Direct elimination in the
old whole D_n gives

    q_d=n*beta^2/2+(n-2)*beta/2+n/8.

Conjugate positive definiteness yields the stated small negative
interval for every beta'. The independent-block determinant ratio
also gives K=Q(beta), as used in the proof. The degree-at-most-five
finite-angle list is complete; each surviving native finite label
has a conjugate outside that interval. Compactness excludes the
positive-null native beta=1 pair. Thus beta>1 and A=2beta+1 is an
integral primitive element with -2/5<A'<3/5 at all other embeddings.

I checked the short integral-interval lemma separately. Degree at
least four has norm below1. In degree three a negative small root
also makes the norm below1; with two positive small roots the norm
is1 and trace4, whereas r+t=4-A<1 contradicts r+t>=2/sqrt(A)>1.
In degree two the only allowed field is Q(sqrt3), with A=2+sqrt3.
The larger trace3 candidate list used in the source is harmless:
all those fields are excluded. In particular A<7/2 is impossible.

The source's preliminary x>1 exclusions use independent old-D-plus-c
blocks before any quadrilateral rank equation. The inherited near
and far finite-angle classifications therefore apply without the
old isolated-d condition. I checked that each use has the same
contacts, determinant-ratio field generation and cyclic integrality
interface; no desired mixed composition is assumed there.

## The three rank equations

Writing N=n-2 and s=(N+2)/(4N), I rederived the one-F Schur matrix
and its determinant identity

    (X^2/t-s)(beta^2-3/4)=(beta+1)/4.

Here t>0, t<=1 and X>1. At beta=5/4 the lower comparison is
25/1088>0. The expression increases thereafter, forcing beta<5/4,
hence the impossible A<7/2.

For two F roots and one zero F/I corner, elimination through the
old D block gives exactly the source's polynomial

    2N(beta^2-3/4)X^2-(N/2)(beta+1/2)X-1/16-sQ,
    Q=(N+2)beta^2/2+N*beta/2+(N-6)/8.

Its derivatives are positive for X>=1,beta>=5/4, and its value at
(1,5/4) is (47N^2-132N-52)/(128N)>0 for N17,18. Thus again beta<5/4.

For four selected corners, the symmetric/antisymmetric split gives

    (X-s)(beta^2-3/4)=(beta+1)/2.

The near placement forces beta<5/4. The far placement forces
beta<7/5 and hence A<19/5. The sole interval survivor has
beta=(1+sqrt3)/2 and x=(3+5sqrt3)/11. The separate actual a-b-c
triangle requires integral2x, whose trace is12/11. This last
candidate is therefore excluded by a distinct checked arithmetic
condition, not discarded as an unwanted real model.

## Composition boundary

The mixed faces cannot become unpaid. Together with the standalone
theorems this permits the A two-end and D three-unit deletions
simultaneously. B reductions have white graphs containing a4-bond,
so their decreasing roles do not mix with these faces. Their new
coupons are nonnegative.

At n20 this supplies simultaneous connected A20 and D20 costs35.
The balanced B20 and B19+A1 costs remain141/4, and outside-factor
and other-field coverage is still required. These local numbers do
not establish a global allowance35 or compact arithmetic ceiling19.
There is no script or Lean verification of this new composition yet.
