# Independent terminal and E8 audit for the affine-tree quadrilateral argument

The [independent checker](../scripts/fresh_simple_cusp_affine_tree_terminal_audit_20261017.py)
and [exact data](fresh_simple_cusp_affine_tree_terminal_audit_20261017.json)
verify the terminal normalization and the two exceptional E8 calculations
in the proposed affine B/C/D quadrilateral completion. The mathematical
case split passes this read after the corrections below. The all-finite
quadrilateral transfer remains a separate obligation; this note makes
no ambient dimension claim.

The review binds the finalized [completion source](../scripts/fresh_simple_cusp_affine_tree_completion_20261017.py)
and [data](fresh_simple_cusp_affine_tree_completion_20261017.json). Its
whole-source read covers the corrected position partition, profile lists,
Green-function derivation, D2 bounds and exceptional branches. It does
not claim a separate replay of the producer's 82,542 formula-pair probes.

Assume throughout an all-simple finite-volume Coxeter polytope of dimension
n>=18, with the finite/ideal weights in the
[frozen interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md).
The [proper-affine exclusion](fresh_simple_cusp_affine_exclusion_20261017.md)
and the [short-end branch](fresh_simple_cusp_affine_tree_short_ends_20261017.md)
are inherited. In particular an adjacent corner must be finite or an entire
connected affine diagram of order n. A proper affine component cannot be
left beside another component.

## Two omitted affine replacements, now accounted for

Deleting a doubled terminal a leaves a long finite B or D diagram. Its complete
one-node replacements have four profiles: isolated; bond3 at its long
endpoint; bond4 at that endpoint; and bond3 at the next node. The last is
a full affine B or D diagram, with the old long endpoint and new root as fork
leaves. It is not excluded merely for being affine.

Deleting a fork leaf leaves a long B or D diagram. There are again four
profiles: isolated; bond3 at its long endpoint (the other fork leaf);
bond4 at that endpoint; and bond3 at its next node (the old fork center).
The doubled endpoint replacement is full affine C or B, respectively.
The original draft omitted these two full-order affine profiles. The
checker reconstructs their full positive null vectors explicitly.

These lists are complete in every required rank. The retained connected
B/D block has rank n-1>=17, so finite exceptional extensions are impossible.
Deleting the added node from the classical finite/affine paths and forks
gives exactly the listed locations and labels. Cycles cannot contain this
retained B/D block. This is a classification argument, not an inference
from a finite label search.

For the doubled-terminal next-node profile, at x=-a.c=0 the five vertices
`a, long endpoint, next node, c, following path node` induce affine B4.
The center is the next node, its two short leaves are c and the following
node, and its length-two arm ends at the doubled terminal a. Thus x>=1.
For the fork-leaf doubled-mate profile, x=0 similarly gives affine B4 on
`c, mate, fork center, a, first path node`. This exclusion is available,
although its weaker nonnegative-x bound already suffices.

## All-rank normalization and separation

Let nu be the positive cusp kernel, with coefficients one on the main
path and fork centers, 1/sqrt(2) at a doubled terminal and 1/2 at fork
leaves. Then `L=diag(nu) G diag(nu)` is a graph Laplacian. Every main-path
or doubled-terminal edge has conductance one; each fork-leaf edge has
conductance one half. Thus a fork-leaf edge has resistance two.

Ground the two marked roots a,b. Dead branches without sources do not
change their effective resistance R. Schur elimination gives
`R^(-1)[[1,-1],[-1,1]]` on the scaled marked roots. On the surviving path,
the Dirichlet Green function between ordered positions xi<=eta is
`xi(R-eta)/R`. An off-path branch adds its branch resistance to self terms,
but not to cross terms with a separated support at the other end.

For a new contact row w use the Laplacian source `diag(nu)w`.
Its first path moment l gives residual cross contact l/R to the opposite
marked root and residual square `c0+l^2/R`. If the left and right supports
are separated, their projection cross term is exactly lm/R. This proves
the factorized determinant identity in every rank. It also explains why
the fork mate has an additional branch-resistance contribution; its
normalization cannot be copied from a main-path node.

The independently reconstructed terminal data are:

| Mark and replacement | c0 | l | height h |
|---|---:|---:|---:|
| doubled terminal, isolated |2|0|x/sqrt(2)|
| doubled terminal, endpoint3 |1|1|1+x/sqrt(2)|
| doubled terminal, endpoint4 |0|sqrt(2)|(x+2)/sqrt(2)|
| doubled terminal, next-node3 |0|2|1+x/sqrt(2)|
| fork leaf, isolated |2|0|x/2|
| fork leaf, mate3 |1|1|(x+1)/2|
| fork leaf, mate4 |0|sqrt(2)|(x+sqrt(2))/2|
| fork leaf, center3 |0|2|1+x/2|

For `g=2l/h+c0/h^2`, the newly added doubled-terminal profile has
`g<=8-4sqrt(2)<5/2`. Every doubled-terminal profile therefore has g<=5/2.
The newly added fork-mate4 profile has g<=4 even before excluding x=0;
the uniform fork bound remains four. The other bounds and proper-affine
low-contact exclusions in the corrected proof are valid.

An unpaid doubled terminal is at least six graph steps from a main-path
mark; an unpaid fork leaf is at least five. These distances follow from
the actual six-node end-segment weight, including both interchangeable
fork leaves. Their resistances are P and P+1, respectively. Opposite
terminal pairs have resistance n-1 in all three B/C/D shapes. The two
uniform terminal contributions are therefore at most5, 13/2 or8, below
R when n>=18. Local supports are separated in these cases.

## Main-path and fork-center completeness

The fork center is the D2 outer-end case: two orthogonal A1 roots.
Its inverse is I/2, old attachment (1,1), and cusp coefficients (1/2,1/2).
The arbitrary outer-only contact estimates and its single-contact affine
E8/B4 exclusions give g<=14/9. Its inner joins touch one leaf by3 or both
leaves by3; their bounds are 14/9 and1. A doubled single-leaf join creates
a proper affine component beside the other A1 and is not allowed.

For other main-path ends, the inherited short-end classification applies.
When one inner block is small, the opposite outer block has rank at least
n-9>=9. The only additional finite inner extensions are D6 to E7 and D7
to E8. B inner blocks have rank at least five, excluding F4, and D inner
blocks have rank at least six. D6 to E7 is paid by component size.
For the D2/P5 to E8 case, the marked old root is the branching vertex.
The special-vertex assertion is explicitly supported by Vinberg's
§6.3, Lemma3, printed p.81: the exceptional unpaid conclusion there is
of type D, not E8. See the
[original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).
This citation supplies only that local finite-weight input.

The ordinary bound g<=14/9 gives a strict determinant gap for main-path
distance at least four. The two remaining literal E8 configurations are
checked separately below; no unresolved small-rank tail is hidden in this
ordinary estimate.

## The main-path D3/P4 exception

Take the outer D3 with center0 and leaves L,L', marked a=main1 and
b=main5. The exceptional d meets L' by3. The relevant white matrix is
`D3 direct-sum A3`; the far end contributes its already proved unit
port cost. Independent elimination gives alpha=1/4,A=1/2,D=5/4.
At y=-b.d=0, the original old tree together with d contains the induced
affine E8 with arms1,2,5. Hence k=y+1/2>=3/2 and v=1/k<=2/3.

All four determinant expressions in the corrected proof are independently
expanded from the actual four-by-four Schur matrix. Their rectangle lower
bounds are 1/18,11/144,1/18,14/225 for isolated, inner-only, center-join and
same-leaf join. The other-leaf join increases the expression by uv/4.
Outer-only c is paid instead: c,d lie in one finite component of order
at most five in the opposite corner. A disconnected proper affine
component would violate the inherited simple-cusp exclusion.

## The fork-leaf P6 exception and the graph correction

Here the literal white block is the six-node path

`mate L -- center0 -- main1 -- main2 -- main3 -- main4`.

Its inverse gives alpha=1/7,A=5/14,D=8/7 for d attached to L by3.
The same proper affine E8 exclusion gives v<=2/3.

The opposite components must be identified correctly:

- If c meets L by3, the opposite component is **finite D8**, with c,d
  its two fork leaves. Their fork-invariant end segment is paid.
- If c meets L by4, the opposite corner contains the four-node affine
  star with bonds4,3,3, which is forbidden.
- If c meets center0 by3, the opposite component is **finite E8** with
  arms1,2,4. It is not an affine E8 obstruction.

For the last case, literal elimination gives B=2/7,C=4/7,E=5/7 and

`det Schur/(h^2 k^2)=1-4u/7-5v/7-v^2/28`.

At x=-a.c=0 there is an induced affine D4 centered0 with leaves a,c,L
and main1. Thus x>=1 and u=1/(1+x/2)<=2/3 already give the sufficient
lower bound8/63. The additional x=1 affine-A2 exclusion yields the stronger
bound in the corrected primary draft, but is unnecessary here.

For isolated c, the determinant is
`1-5v/7-v^2/28-2u^2/7`. The low-contact exclusions imply x>=phi>3/2,
so u=2/x<4/3. The expression decreases in both positive variables and
is zero at u=4/3,v=2/3. Its actual value is therefore strictly positive.

## Verification scope

The checker imports no producer arithmetic. It uses physical norm-two
SymPy matrices, exact inversion, literal positive kernels, and Sylvester
minors. It confirms 48 terminal normalizations at orders18,19,26; these
are checks of the all-rank network identities above, not a rank cutoff.
It independently expands all exceptional Schur determinants, verifies the
positive lower bounds, and distinguishes finite D8/E8 from the proper
affine stars used in the exclusions. Generation and fresh saved replay
pass (final bound-input run30394 completed both with exit0). The input
proper-affine exclusion and published finite weights are
inherited. No separate all-finite-corner transfer, whole mixed-rule proof,
or new dimension bound is asserted by this audit.
