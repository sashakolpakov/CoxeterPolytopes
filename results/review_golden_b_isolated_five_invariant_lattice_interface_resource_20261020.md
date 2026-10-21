# Outside audit of the golden B-star lattice interface

20 October 2026. **PASS.** Full outside-author paper audit with own
algebra, finite-field and lattice derivations; no scientific import,
enumeration, worker, Lean or TeX execution. Internal review only.

Bound source:
`results/golden_b_isolated_five_invariant_lattice_interface_20261020.md`,
SHA256 `e0e347a7642ca196ec416cc8889adf666bf4903a3ebfe865e9d65549257b4e10`.
The scope is isolated native5 at long nodes1 or5 of a WHOLE B_n
vertex, n>=8; it is not a whole-B_n or dimension exclusion.

## Algebra and every invariant lattice

Both displayed caps have norm2 and unit time/nonzero spatial
coordinates. Signed permutations give all spatial2E_ij. Multiplying
the cap reflection on either side isolates every mixed2E_i0/2E_0i
with unit coefficient; multiplying again isolates2E_00. Thus the
generated O-algebra really contains2End(L), without assuming a
particular lattice at odd primes. For the coordinate-content ideal I,
the resulting2IL<=Lambda<=IL and the PID assertion give the complete
global sandwich2L<=Lambda<=L after scalar normalization.

I checked the two cap reductions over O/2=F4 directly. For one or
five supported coordinates they are the same rank-one operator in the
source. For five, replacing one supported index by another produces
every spatial augmentation difference; n>5 is sufficient. The common
kernel is the stated line D and the image-orbit span is the stated
hyperplane E. If U is not contained in D, one rank-one image lies in U
and its permutation orbit spans E. Hence exactly0,D,E,V occur. Since
Lambda contains2L, it is its residue subspace's full inverse image;
there are no missing dyadic lifts or odd-prime alternatives.

For d=(phi,1,...,1), f(d)=n-1-2phi has2-adic valuation0 when n is
even and exactly1 when n is odd, because the residue of phi is not
in F2. Thus every nonintegral element of Lambda_D/2 has a genuinely
negative-valuation last term in the displayed norm expansion, while
the other terms are integral. The integral-norm hull is exactly L.
This proves O(Lambda_D)<=O(L). Conversely modulo2, B(d,y)=sqrt(f(y))
holds and determines d uniquely, so every O(L) isometry fixes its
residue and preserves Lambda_D. The equality of groups follows.
Unimodularity gives Lambda_E=2Lambda_D^*, and duality gives the same
group for E. All four possible lattices therefore have exactly O(L),
not just commensurable orthogonal groups.

## Full-rank arithmetic star and actual chamber application

Normalize the B_n short root to norm1 and long roots to norm2. A label5
contact allows cap normalization to norm2 over the actual field;
its spatial vector is phi times the first k coordinate vectors.
The resulting full Lorentz star has Gram entries in Q(phi). A proper
ground-field extension would provide another real embedding fixing
this entire indefinite Gram matrix, contrary to admissibility. Hence
the ground field is exactly Q(phi). At its other embedding the cap
condition is k*(phi')^2<2, giving k<=5. The time squares at k1,5
are phi^-1 and phi^5 and produce exactly the two source caps after
the stated time-axis normalization.

Arithmeticity supplies a full invariant O-lattice: sum the translates
of a lattice preserved by a finite-index arithmetic subgroup over a
finite coset transversal. The proved equality O(Lambda)=O(L) then
places EVERY wall reflection in O(L), and hence in its full reflection
subgroup. A full-arrangement chamber at a generic point of the assumed
compact polytope is contained in that polytope. This is actual chamber
containment, not a transfer through a similar form or a finite index
claim in the wrong group.

The precise external cutoff G is indispensable. I independently read
the PDF text of [Kellerhals2014](https://homeweb.unifr.ch/kellerha/pub/CMFT.pdf),
Section4.1.1, printedp474: equation(4.3) displays this literal golden
diagonal and the following sentence gives reflectivity only for n<=7,
with the resulting groups cocompact. The [Vinberg1985 survey](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng),
Introduction p.38 and Section5.6 p.64, independently identifies the
same family and Bugaenko attribution. The original Bugaenko1984 proof
was not retrieved or reconstructed here; G remains an explicit import.
With G, the contained full chamber cannot be compact in n>=8.

No flaw found. Nodes2/3/4, dotted or multicontact caps and other finite
types are outside this theorem. No arbitrary golden-form, uniform
dimension or commensurability-class exclusion follows; no new
manuscript, computation or formalization is implied.
