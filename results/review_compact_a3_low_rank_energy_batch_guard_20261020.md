# Independent review of A3 energy and low-rank A forests

20 October 2026. Full paper/source review; no execution, import,
enumeration, optimizer, build or formal verification.

**PASS within the stated scope.** I independently read and rederived
the completed source
`results/compact_a3_block_energy_and_low_rank_forest_exclusion_20261020.md`,
SHA256
`70ab9f7e9fad80976b772e09388746db79c7bd161fc856bdc4ae0d38618b5068`.
No mathematical flaw or correction was found.

The result is a pointwise restriction: at the SAME vertex farthest
from a fixed interior point, each A3 component has energy>1/4 when
the whole finite diagram uses only A1,A2,A3; consequently the whole
forest A1^a+A2^r+A3^s is excluded in dimension n>=14. The theorem is
geometric, with no arithmetic premise. It does not exclude higher
rank or other finite components and proves no universal dimension
ceiling by itself.

## 1. Same-wall propagation and the A3 component

I checked the Lorentz signs, finite-point strict norm ||x||<1, and
additivity of orthogonal component energies. The reciprocal-shadow
alternative retains the DISTINCT actual wall condition. Its consequence
r_i*r_j>=K_ij is legitimate only when the two actual opposites differ,
as stated. No integer surrogate for a real height or formal corner
is substituted here.

For adjacent A3 directions K_ij=k=1/sqrt3, an assumed component
energy<=1/4 gives tau_i+tau_j<=sqrt((1+k)/2). I independently
differentiated log H(t), H(t)=2t/(1+t^2), and confirmed its strict
concavity on(0,1). The resulting upper product is
32(1+k)/(9+k)^2. Its strict inequality below k reduces to74k>39,
which follows from74^2>3*39^2. Thus the SAME actual wall closes all
three directions. Individual constraints give r_i<=4/5 and the
height window A^2<=(16/25)v_i^2.

The adjacent multi-contact exclusion uses an ACTUAL retained white
triangle at the endpoint of the third edge; its Gram cannot be
positive definite once all three contacts are at least simple3.
This is not a general prohibition of arbitrary nonpositive diagrams.
For the two old endpoints I rederived

    q-(16/25)v_2^2
       =(59/50)(u^2+w^2)+(9/25)u*w.

The simple/simple case is affine. The next native values yield the
strict lower bound177/200+9sqrt2/100>1, so arbitrary spectator costs
cannot rescue a multi-contact wall.

## 2. Complete spectator costs and single-contact exclusions

All spectator costs in(7) follow from the displayed inverse diagonals
and the positive definite common-white extension. Adjacent double
contacts contain the affine simple triangle; the two A3 endpoints
contain the affine four-cycle. Thus no multi-contact spectator is
missing. I checked the sum gaps after1/4,1/2,5/8; three positive
summands already cost at least3/4, so the argument includes arbitrarily
many spectator factors.

For a single A3 endpoint contact the minimal possible tangent square
is greater than16/225. The reciprocal far-end ratio is exactly3,
contradicting r_far<=4/5. For a single middle contact the possible
tangent square is at least1/6; both endpoint tangents then have
square>=1/2. Their correlation1/3 requires energy at least
9/2-3sqrt2>1/4, with the exact comparison289>288. The isolated
label6 and middle label4 affine exclusions are correctly explicit.
These cases exhaust the contact vector after the multi-contact step.

## 3. All new A3 spectators in the A1/A2 accounting

I independently checked BOTH complete release-cost windows against
the new A3 costs, including sums with the old A1/A2 costs.

For an A1-short direction, the only surviving additional possibilities
are P, Q5 and R7..12 in(13). The other sums are boundary values or
contain the displayed positive affine path/star. The label13 cutoff
from sin(pi/13)<22/91<1/4 is valid. P has too small a height product
for any positive overlap. Q5's middle overlap with an A3 endpoint
costs1/4, larger than its possible height product. Two R walls at
opposite endpoints have overlap1/8 BELOW their height product, while
at equal endpoints it is ABOVE; neither can be orthogonal. Their
overlaps with all old A1-short profiles are also too large or absent.
Thus a new profile creates only one A1-short direction and coexists
with none of the other short walls. The old3/25 debit remains valid.

For an A2-short direction, the additional sums in the label3 window
are exactly3/4 and17/24, with the decompositions listed. I verified
the affine F4 path, three-arm E6, and double-fork witnesses and their
positive null coefficients. In the label4 window the new costs3/8
and1/2 have the stated affine witnesses; the latter five coefficients
(1,1/sqrt2,1,1/2,1/2) solve the unit Gram equations directly.
For label5 the sole new case is G. Its mate white norm is negative,
so its ACTUAL opposite differs, precisely as required to apply the
shadow inequality. That mate has tangent square>2/3 and contributes
energy>1/4. It therefore does not create a new deficient A2 block.

Finally all new A1-short profiles have A1/A3 support disjoint from
the isolated7 A2 wall. Their products lie in the forbidden gap(-1/2,0).
Thus the A1 debit and exceptional isolated7 debit still cannot occur
together. Every inherited consumer/payment is explicitly retained.

## 4. Composition and limitations

The rank-normalized costs satisfy e3/3=1/12>319/4000=e2/2 and
1/10>e2/2. With isolated7 the rank14 bound is2093/2000>1.
Without it and with at least one A1 factor the bound is4067/4000>1;
with no A1 factor it is at least2233/2000>1. These comparisons are
at the one farthest vertex, and strictness is supplied by present
components rather than fictional credits for absent factors.

The inherited geometric inputs remain compact simplicity, all-wall
acuteness, the native Coxeter gap, and compactness exclusion of
positive affine null combinations. The prior A1/A2 source is pinned
inside the reviewed note. I am a developer of the general farthest-
vertex/shadow mechanism, but not an author of this A3 continuation or
spectator audit; the independent derivations above check the newly
assembled proof. This review makes no claim of external refereeing
or an end-to-end geometric formalization.
