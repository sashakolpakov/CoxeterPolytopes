# Outside-author review: golden D_n long-end-five invariant lattices

20 October 2026. **PASS**, for the precise whole-D_n, isolated long-end5,
compact arithmetic first-type profile in the source. Full paper audit;
no source edits, enumeration, imported worker, scientific execution,
formal verification or dimension promotion.

Reviewed in full:

* `results/golden_dn_long_end_five_invariant_lattices_20261020.md`,
  SHA256 `58e8553bd53f00c78aaa336b5798ea83abbe4148c14c8403ad14e1416f5f9f2b`.

I derived the algebra, complete residue alternatives and norm hulls before
using the conclusion. The earlier B-star interface is not substituted
for the D-star argument: individual spatial sign flips are unavailable
here, and the fourfold sandwich and odd-dimensional deep lifts matter.

## Algebra and every residue layer

With f=<-u>+I_n, the norm-two cap e0+phi*ei has reflection block
[[phi,-phi],[1,-phi]]. Even spatial sign changes give
(1-dij)(1-dik)=4Eii for three distinct spatial indices. Multiplication
by permutations gives every spatial4Eij; the displayed products with
Ri successively give4Ei0,4E0i,4E00 with unit coefficients. Thus the
algebra contains4End(L). Coordinate-content normalization in the PID
really gives4L<=Lambda<=L, including control at odd primes; it does
not justify replacing4 by2.

Over F4 the functional ell=phi*t+sum xi is invariant and v=(phi^2,1^n)
is fixed. A transposition difference from unequal spatial coordinates
generates the entire spatial augmentation space. Applying Ri supplies
the remaining direction of H=ker ell. If all spatial coordinates are
equal, Ri-invariance instead forces the fixed line<v>. Hence0,<v>,H,V
is the complete subspace list. This remains correct when n is odd and
<v> lies inside H.

The two residue layers satisfy W1<=W2 because twice a lattice lift is
in Lambda intersect2L. For W1=H, a lifted ei+ek and the difference
under dij put ei in W2, whose permutation/Ri orbit spans V. For
W1=<v>, the differences put all ei+ej and then H in W2. Thus even n
forces W2=V; odd n allows precisely H or V. In the deep case, choosing
one lift lambda=v+2z gives Lambda=O*lambda+2E, with no omitted lift
parameters. W1=V gives L by Nakayama on L/Lambda, annihilated by4.
This verifies all alternatives in(3), including odd n.

## The two norm-recovery equalities

The standard lattice L is unimodular because u is a unit. Modulo2,
f(x)=ell(x)^2, so E is exactly its even-norm sublattice. The bilinear
pairing B(v,x) reduces to ell(x), giving E*=L+O*v/2 and
f(v)=n-1-2phi. Its2-valuation is0 for even n,1 for odd n.

In x=y+a*v/2, y in L, the first norm term and the cross term aB(y,v)
are integral. If a is a dyadic unit, the last term has valuation-2
or-1 respectively; it cannot cancel with those integral terms.
If a is divisible by2, x belongs to L. This proves(5), not merely
one inclusion. There are no denominators at other primes.

For the odd deep lattice, lambda lies in E and
f(lambda)=f(v)+4B(v,z)+4f(z) still has valuation1. Writing an element
of Lambda/2 as a*lambda/2+e, e in E, repeats the same noncancellation
argument. Integral norm forces a divisible by2, and then the element
is in E; conversely E is present. Thus(6) is exact. Duality and these
intrinsic norm sets give O(f,Lambda)<=O(f,L) for every listed lattice.
No equality of all these orthogonal groups is required or inferred.

## Arithmetic and the final direction

The full-rank normalization is legitimate over the arithmetic field:
along the connected3-bonds, old root lengths have square ratios in K,
so K-rescalings make their norms a common totally positive value b.
The totally positive form scale2/b makes every old norm2. The label5
square gives phi in K. Its sole cap contact g then satisfies
g^2=q*phi^2/2, so a K-rescaling -phi/g makes that contact -phi and its
norm2. The displayed n+1 roots span, fixing exactly this model.
A proper totally real extension of Q(phi) would have a nonphysical
embedding fixing this full Lorentz Gram, contradicting admissibility.

After that K-rational realization, an arithmetic subgroup preserving
a lattice has finite index in Gamma. Summing its finitely many coset
translates supplies a full O-lattice stable under Gamma. Thus the
classification applies to the ENTIRE actual reflection group, not
only the old vertex stabilizer. Every facet reflection therefore
preserves L and belongs to its full integral reflection group W(L).
The unit time change t->u*t identifies this literal lattice with the
standard <-phi>+I_n lattice integrally.

I independently checked the published account in Kellerhals,
*Hyperbolic Orbifolds of Minimal Volume* (2014), section4.1.1,
printed p474, equation(4.3) and its following paragraph:
https://homeweb.unifr.ch/kellerha/pub/CMFT.pdf . It states the cutoff
n<=7 for this exact lattice and attributes it to Bugaenko1984. The
original Bugaenko proof was not re-proved or independently retrieved
here. This remains an explicit external input, rather than an
all-commensurability-class assertion.

The covolume implication has the correct direction. The full integral
orthogonal group is an arithmetic lattice; nonreflectivity for n>=8
means W(L) has infinite covolume. A discrete group containing the
finite-covolume Gamma cannot have infinite covolume. Equivalently a
full W(L)-chamber refining the compact Gamma chamber would itself
be compact. The inclusion Gamma<=W(L), not merely commensurability
of lattices, therefore gives the contradiction.

No gap found. Spectator components, different D_n contact positions,
nonarithmetic polytopes, universal ceilings and novelty remain outside
this verdict. No executable or formal-geometry coverage is claimed.
