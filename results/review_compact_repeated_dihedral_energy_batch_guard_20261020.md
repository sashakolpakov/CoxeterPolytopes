# Independent review of repeated dihedral corner energy

20 October 2026. Paper/source review only; no execution, enumeration,
imports, optimizer, build or formal verification.

**PASS.** I read the complete source
`results/compact_repeated_dihedral_corner_energy_20261020.md` and found
no flaw in its stated pointwise theorem: a farthest actual vertex of a
compact Coxeter polytope cannot have whole diagram I2(m)^r, for any
m>=2 and n=2r>=14. Mixed labels and other finite factors are outside
its scope. There is no new universal dimension ceiling here.
The reviewed source SHA256 is
`8149d0f3ade29b3c8d7644847c972fccfd28c64cfe2153b5563684c6b2d0d9a8`.

I independently rederived the release formula

    delta=1-E-g^2,
    cosh L=sqrt(d)*(b+c*g)/sqrt(delta),
    d=csc^2(pi/m), c=cos(pi/m).

Both positive contacts give the stated cot(theta/2)/sqrt3 bound. A
one-contact finite wall has positive definite replacements at BOTH
directions. Acute projection proves both actual endpoints; the
white-Gram assertion is not treated merely as a formal continuation.
The mate formula and the fact that a one-contact dotted wall cannot
close its mate agree with the exact reciprocal-shadow alternative.

For m=5, the complete off-block inventory follows directly from
k=d/4=(5+sqrt5)/10: a larger single contact or any double contact in
another block already costs>1, and2k>1. Thus E is exactly0 or k.
For E=k I checked separately both-positive, g=0, and b=0 cases. In the
last case positivity forces g=1/2 and

    cosh^2 L=k*c^2/(3/4-k)>7.

There is no omitted zero-contact case because it would give negative
height square. For E=0 the only weak possibility is a single finite
contact, necessarily at least1/sqrt2; its ACTUAL mate has
cosh L>=sqrt(d-1)>4/3. Hence each block contributes energy>1/7.

For m>=6 every off-block contact costs at least d/4>=1, so the finite
endpoint white Gram excludes it. At m=6 the one-contact3 case has zero
height; all other cases have a sufficiently long direct edge or mate.
For m>=8 the minimum finite contact1/2 gives mate bound
cot(theta)/sqrt3>4/3.

At m=7, a single3 wall has squared height h=d/4-1 in(0,1/2).
Distinct exceptional BLOCKS require distinct walls with disjoint
support. Their product would be-h in(-1/2,0), forbidden for actual
outward Coxeter normals. Thus at most one block is exceptional, even
though its wall can close both directions of that block. Every other
block has an edge with cosh L>2, hence energy>1/3. This proves the
claimed stronger exclusion for r>=4 in this label.

The inherited m=2,3,4 sources give the displayed thresholds12,14,12.
Their valid shared linear tests are also valid local tests at the SAME
farthest vertex; no averaging substitution is used. Summing component
energies is legitimate because the old components are orthogonal.

The only external geometric inputs are compact simplicity, all-wall
acuteness, actual native Coxeter products, and the compactness exclusion
of positive affine null combinations. All profiles beyond those inputs
are derived in the source. No arithmetic field restriction, facet
reflection group, or facet-descent premise is used.

## Successor review: arbitrary mixed labels at least four

I subsequently read the ENTIRE consolidated successor of the same
source, SHA256
`baef1c51a99822cdb1e5930b4f23cfc5aaf9e6a386b1dc602a1c3e7a6fe2642c`.
The preceding paragraphs bind the earlier repeated-label edition;
they are retained as its review history. **The successor also PASSES.**
It proves the additional pointwise exclusion of EVERY whole forest
sum_j I2(m_j), with all m_j>=4 and rank2r>=14.

The full off-block inventory is exactly E=0,1/2,k5: a B2 spectator
permits only one simple3 contact of cost1/2; an I2(5) spectator permits
only one simple3 contact of cost k5; every higher label is sealed.
Two spectator contacts cannot coexist in the finite white set.
The zero-contact distinguished-block case has negative height square,
so it does not leave an unclassified outgoing edge.

I checked every new asymmetric case in Section4. For E=1/2 and b=0,
positivity forces g=1/2; at target m=4 the height is EXACTLY zero,
whereas m>=5 has direct cosh L=cot(theta)>11/8. For E=1/2 and g=0,
the target m=5 contact3 has an ACTUAL positive-definite mate, whose
cosh is cot(pi/5); the target m=4 contact3 again has zero height.
All larger contacts or labels have the displayed strong direct edge.

For E=k5 and target m=4, contact3 at the released node has a finite
mate with norm delta=3/4-k5 and cosh-square>5. A contact at least4
has the strong direct bound. With b=0 the same positive denominator
and d>=2 give direct cosh-square>5. The m>=5 direct cases are
minimized at m=5, giving cosh>=phi. These are actual-edge arguments;
the weak direct edge is never silently credited with its mate's length.

Every strong block contributes more than3/19. The ONLY weak case is
still the E=0 single3 wall in a7 block, so the old disjoint-support
argument permits at most one such block even among mixed labels.
Its squared mate half-tangent exceeds9/100. I checked the exact
2652 cross-product difference in its rational comparison. Thus seven
blocks have energy greater than18/19+9/100=1971/1900>1.

No correction was needed. The stronger statement is geometric, not
arithmetic, and still leaves A1/A2 mixtures with higher labels and
connected components of rank at least three outside its scope.
