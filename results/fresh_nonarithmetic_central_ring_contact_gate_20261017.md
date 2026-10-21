# An all-label obstruction to the final central cap in twenty literal rings

None of the twenty exact ring recipes below admits one additional
norm-2 Coxeter wall inducing an actual S6 facet and terminating all four
consecutive ring edges at finite vertices. The obstruction covers every
finite contact label and arbitrary real coefficient fields. It already
holds before imposing the other ambient walls or a common center.

This concerns the twenty selected rings in the
[central-patch test](fresh_nonarithmetic_release47_central_patch_20261017.md).
It is not a theorem about every coherent fourth cap in the predecessor
bank, additional central walls, changed ring parameters, or a different
patch. The exact rings are reconstructed recipe objects. This package
does not certify their entire ambient native-contact lists or promote
them as new chambers. The frozen 64-wall control is unchanged.

## Exact ring reconstruction

The [source](../scripts/fresh_nonarithmetic_central_ring_contact_gate_20261017.py)
starts with the frozen exact residual Gram in the basis [8,30,44,52]
from the [control boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md).
It reconstructs the changed 12 and then L,P,Q,R from their three finite
native contacts and norm equations. White A2 components are retained in
every native-to-projected conversion.

For each three-root frame, the stored inverse gives the particular
solution. A recorded old or previously constructed normal supplies a
nonzero orthogonal kernel. Its square and the required positive square
root are certified with outward rational intervals. The saved sign
chooses one of the two exact norm branches. Numerical coordinates serve
only to identify which branch was meant by the earlier discovery recipe;
they do not accept a contact or an exclusion.

Products between constructed normals are propagated using the projection
formula, rather than repeatedly multiplying expanded residual-coordinate
vectors. For example, if a new root is p+sλk, its product with another
root t is computed as t·p+sλ(t·k), with both terms evaluated from the
literal frame Gram. This avoids interval dependency while enclosing the
same exact algebraic product. Every inverse pivot and all required kernel
signs are saved.

The four projected ring roots span the residual Lorentz four-space; the
checker verifies the inverse pivots of their Gram G. Thus four prescribed
products determine a unique residual vector.

## Complete contact reduction

The retained support is S6=H4+A2. As in the
[seven-profile proof](fresh_nonarithmetic_a2_escape_partition_20261017.md),
an actual new S6 facet must be right to H4. The inverse H4 diagonals
are 14+6√5, 21+9√5, 10+4√5 and 3+√5, all greater than two. Any
nonzero native Coxeter contact already costs more than the available
norm two.

The inverse A2 Gram is J=(1/3)[[2,1],[1,2]]. Two nonzero contacts cost
at least two; one label at least six also costs at least two. The only
possibilities are consequently the zero A2 profile or one label 3,4,5
at node 5 or 6. These are the seven profiles used here.

There is no missing high-label tail among contacts with the four ring
walls. In a finite norm-2 Gram matrix, let a≥√3 be a bond with label
at least six. If either endpoint has another nonright incident bond
b≥1, the corresponding three-by-three determinant is

\[
8-2(a^2+b^2+c^2)-2abc\le8-2(3+1)=0,
\]

where c≥0. This contradicts positive definiteness. Both endpoints of
that high bond must therefore be isolated within every required finite
corner.

For a central bond to ring vertex j to be high, j must have the zero A2
profile and be right to both adjacent ring vertices. The central cap
must also have zero A2 profile and be right to those two neighbors.
The literal ring labels identify all such eligible j. In each of the
twenty rings, there is no eligible opposite pair. Two adjacent high
bonds are excluded by the same three-by-three determinant. Thus at
most one central bond can be high; its remaining opposite contact has
label 2,3,4 or 5.

## Low labels: the inverse-Gram norm obstruction

Let c be one of the seven negative A2 contact vectors and c_i the A2
profile of the i-th ring root. For proposed native ring labels m_i,
the required residual products and residual norm are

\[
b_i=-2\cos(\pi/m_i)-c^TJc_i,
\qquad \rho=2-c^TJc.
\]

The unique residual vector with these four products has square
b^T G^{-1}b. Hence a necessary norm equation is

\[
b^TG^{-1}b-\rho=0.
\]

For each ring the checker tests all seven profiles and all 4^4 label
quadruples in {2,3,4,5}. Every one of the 35,840 resulting exact interval
expressions is strictly positive or strictly negative. The smallest
certified absolute norm error is greater than 0.00006064. The saved
data give sign counts, a positive rational lower bound and an ordered
interval digest for each ring/profile pair. No native-wall rejection
or numerical angle recognition is needed for this part.

## High labels: a finite norm calculation without a cutoff

For each eligible high contact, the other three contacts are two right
angles and one of labels 2,3,4,5. Those three products determine a line;
the norm equation gives both exact branches. All 64 real branches are
evaluated, without choosing or bounding the purported high label.

Their products with the high-contact ring root are strictly greater
than −√3. The smallest certified gap above −√3 is greater than 0.5225.
No branch can therefore have any finite label at least six. This
exhausts the remaining contact profiles and proves the stated
single-central-cap obstruction.

## Scope and verification

The theorem requires all four intended endpoints to be finite. It does
not silently treat a zero Schur residual as finite. The original S6
support contains a retained π/5 edge, so the previously proved cusp
obstruction makes compactness of this face necessary in a finite-volume
completion; that broader cusp theorem is not re-proved here. The result
still concerns this one central-cap topology and these twenty ring
recipes, not arbitrary subsequent completion.

The [certificate](fresh_nonarithmetic_central_ring_contact_gate_20261017.json)
binds the exact control, numerical recipe provenance and arithmetic
helpers. It shares the frozen outward-rational-interval engine and is
not described as an independent arithmetic implementation.
**FINAL/FROZEN:** generation 20552 and fresh saved-data replay 40176
both terminated successfully. Source bindings, syntax and local report
links pass. No new example or dimension bound is claimed.
