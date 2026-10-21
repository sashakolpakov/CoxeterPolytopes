# Outside audit: complete D common-2/3 inventory

20 October 2026. **PASS**. Full independently reviewed source:
`results/compact_dn_low_index_common_pair_inventory_and_relaxation_20261020.md`,
SHA256 `2bd35ff5bd18599c5f4624128e055c3669373eea88d4b9edf8cb7f830c9a99f3`.
No scientific execution, import, optimizer, Lean or TeX job, accepted-source
edit, or status promotion. I previously independently audited the imported
three-profile source b03a8498; that dependency is disclosed here.

## Complete geometric inventory

Deleting2 retains D_(r-2), rank at least7. I checked its complete finite
one-node inventory directly: two contacts create a nonpositive cycle;
an internal3 has a second fork or degree-four obstruction; a bond at
least4 has inverse-Gram cost at least1. Thus only no contact, simple3
at its long leaf3, or simple3 at a fork leaf remains. Deleting3 retains
left A2, so its only possible single contacts have labels3,4,5; two
contacts form a nonpositive triangle and a6 already gives affine6,3.
Every spectator is retained in both actual finite replacements.

For the long-leaf3 case, own2 zero gives the exact two-fork affine
principal diagram, independently of old1 and spectators. Positive
own2 rules out any old1 contact by the left A2 triangle. Own3 gives
the affine triangle; own4 is the accepted three-profile branch. Own5
has white cost phi^2/3 and therefore no spectator, because every nonzero
spectator cost is at least1/4. I obtain delta=(2-phi)/3 and
cosh L3=(2+5phi)/sqrt2>7 exactly as stated.

For a fork3, r>=10 already gives the retained affine E8 arms1,2,5.
At r9 the right D7 plus c is finite E8. Any old1 or spectator contact
would extend this whole E8 in a finite replacement, so both are absent.
I checked every possible remaining own2 contact: zero gives full-star
affine E8;3 gives an all3 cycle;4 gives the affine3,3,4 star;5 gives
the nonpositive retained3,5,3 path with determinant
9/16-phi^2/4<0. No nonparabolic extension or hidden spectator case
is assumed away in this exceptional-rank argument.

With no right contact the six table entries exhaust the left A2
possibilities. The node1/4 terminal4 extension, node2/3 two-fork
diagram, and node2/4 star are exact affine principal obstructions.
The two5 entries are the accepted isolated branches. For node1/3,
the full height Q-1/2 and retained white cost1/3 imply
0<A^2<1/6. Some spectator must be touched; every such contact has
beta_j^2 C_jj>=1/4>A^2. The reciprocal shadow therefore makes the
same actual c close a spectator direction, contradicting the accepted
farthest orthogonal-component theorem. This works for arbitrary
multi-component spectators and does not use an edge-only substitute
for the farthest projection.

These checks prove completeness of the new common-2/3 classification.
Composing its surviving three branches with b03a8498 yields exactly
the ACTUAL projection conclusion ||x_D||^2>9/14 for r>=9. It does
not yield that bound for every abstract edge-covering vector.

## Exact relaxation and all-pairs implication

I independently checked the rational vector, all partial sums, and
its norm641213/1000000. The cleared norm gap is23018. Direct hand
expansion of the six adjacent-tree expressions gives exactly the six
positive numerators displayed in the source:

    1732252476, 188938803039, 227885305759,
    850862963679, 57148587900, 251045977500,

all over10^12. The direction1/2 numerator is exactly
-338772348164 over10^12.

There is no correlation-factorization gap. On the direction tree
2--3--4--5--6--7 with leaves8,9 at7, the long correlations are
sqrt(i/j), each long-to-fork correlation is sqrt(i)/3, and the
fork/fork correlation is7/9. The latter equals
(sqrt7/3)^2, so EVERY endpoint correlation is the product of the
tree-edge correlations on its path. Multiplying the established
pair inequalities gives

    F(h_a)F(h_b) product_(internal v) F(h_v)^2
      >= product_(path edges) K=K_ab.

Since0<F(h_v)<1, this implies the needed endpoint law. Thus the
table does establish all distinct-pair laws within directions2..9,
including the two forks. The closed-chamber ties and final zero
coordinate can be perturbed strictly: the norm and all required
inequalities have positive finite margins. No actual cap normals,
native configuration, or compact polytope is thereby constructed.

I found no correction to request. The theorem excludes a common
2/3 opposite under the stated actual projection upper; the optional
witness accurately exposes the remaining direction1 boundary.
No unconditional D rank/14 rate or ambient ceiling is proved.
