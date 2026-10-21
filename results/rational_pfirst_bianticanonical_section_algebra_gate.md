# The bi-anticanonical section algebra is a genuine Property-4 gate

**Status (2026-08-11): proved same-frame no-go and classification-free
carrier-relative criterion over \(\mathbf C\).**  The marked lattice,
canonical vector, proper-point proximity graph, complete labelled-exceptional
contraction fan, divisor degree/multiplicity arrays, and formal jet-matrix
descriptors do not determine whether the physical class \(-2K\) is a reduced
irreducible rational boundary or the nonreduced double of a genus-one
anticanonical curve.

The missing datum is already visible in the first two pieces of the
anticanonical section algebra: the actual ranks and kernels of the cubic and
double-point sextic evaluation maps, their multiplication, and the divisor
factorization type of the sextic kernel.  This gives a precise
classification-free gate for the adapted `P10` failure.  It does **not** prove
universal incidence nonemptiness or claim that the two surfaces below have
the same full nef cone.

The finite arithmetic replay is
[`certify_rational_pfirst_bianticanonical_section_algebra_gate.py`](../scripts/certify_rational_pfirst_bianticanonical_section_algebra_gate.py),
with its frozen
[`archive`](../data/enumerations/rational_pfirst_bianticanonical_section_algebra_gate.json)
and
[`tests`](../tests/test_rational_pfirst_bianticanonical_section_algebra_gate.py).

## 1. One marked blowup frame, two opposite boundary answers

Let \(p_1,\ldots,p_{10}\) be ten ordered distinct proper points of
\(\mathbf P^2\), and put

\[
                 S=\operatorname {Bl}_{p_1,\ldots,p_{10}}\mathbf P^2.
\]

In the total-transform basis \(H,E_1,\ldots,E_{10}\), every such blowup has

\[
 \operatorname {NS}(S)=I_{1,10},\qquad
 K=-3H+\sum_{i=1}^{10}E_i,\qquad
 F=-K=3H-\sum_{i=1}^{10}E_i,                         \tag{1.1}
\]

and

\[
 B=-2K=6H-2\sum_{i=1}^{10}E_i.                       \tag{1.2}
\]

Thus

\[
 F^2=-1,\quad p_a(F)=1,\qquad
 B^2=-4,\quad K.B=2,\quad p_a(B)=0.                  \tag{1.3}
\]

There are two configurations with this identical marked and proximity
data.

### Theorem 1.1 (same-frame bi-anticanonical dichotomy)

There exist two ordered ten-point configurations \(z_{\mathrm{ell}}\) and
\(z_{\mathrm{Cob}}\), all of whose points are proper and distinct, such that the
corresponding blowups satisfy

\[
\begin{array}{c|cc|c}
 &h^0(-K)&h^0(-2K)&\text{unique bi-anticanonical divisor}\\ \hline
 S_{\mathrm{ell}}&1&1&2F,\quad F\text{ smooth irreducible of genus }1,\\
 S_{\mathrm{Cob}}&0&1&B,\quad B\text{ reduced, irreducible, smooth rational}.
\end{array}                                               \tag{1.4}
\]

For both configurations:

1. the marked \((\operatorname {NS},K,F,B)\) tuple is (1.1)--(1.2);
2. every proximity-parent list is empty;
3. the ten labelled exceptional curves are disjoint `(-1)`-curves, so all
   \(2^{10}\) labelled contraction subsets and all \(10!\) complete
   contraction orders occur; and
4. the formal \(-K\) and \(-2K\) systems are respectively the same
   \(10\times10\) cubic matrix descriptor and \(30\times28\) sextic
   double-jet matrix descriptor.

Only the actual coordinate specialization, matrix rank/kernel, and kernel
factorization differ.

### Proof: the genus-one branch

Take a smooth plane cubic \(C\) and any ten distinct points on it.  Its strict
transform \(F\) has class \(-K\) and remains a smooth irreducible genus-one
curve.

Any other cubic \(Q\) through the ten points would satisfy

\[
                         C.Q=9<10,                         \tag{1.5}
\]

unless it contained \(C\).  Bézout therefore forces \(Q=C\), and
\(h^0(-K)=1\).

Now let \(G\) be a sextic having multiplicity at least two at every
\(p_i\).  If \(C\) were not a component, local intersection along the smooth
branch of \(C\) would give

\[
                         C.G=18<20.                        \tag{1.6}
\]

Hence \(G=C+Q\), where \(Q\) is a cubic through all ten points.  The preceding
paragraph gives \(Q=C\), so \(G=2C\).  After blowing up,

\[
                    |-K|=\{F\},\qquad |-2K|=\{2F\}.       \tag{1.7}
\]

The physical numerical class \(B=2F\) therefore exists only as a nonreduced
double in this system.

### Proof: the strong-Coble branch

Take an irreducible plane sextic \(\Gamma\) having exactly ten ordinary
nodes and no other singularity, and let the \(p_i\) be its nodes.  Such a
sextic exists over \(\mathbf C\).  One completely exact choice in this
project is the Winger sextic over \(\mathbf Q(\sqrt5)\): its frozen
certificate proves that its ten centers are distinct, its cubic matrix has
rank ten, its double-point sextic matrix has rank 27, the unique sextic is
absolutely irreducible, and its projective singular scheme consists of ten
reduced ordinary nodes.  See
[`h3_winger_arrangement_coble_natural_cone.md`](h3_winger_arrangement_coble_natural_cone.md).

There is also a construction-independent existence argument: a general
three-dimensional base-point-free subspace of
\(H^0(\mathbf P^1,\mathcal O(6))\) maps \(\mathbf P^1\) birationally to a
plane sextic.  The generic-projection theorem gives only transverse double
points; the genus formula then gives exactly

\[
             p_a(\Gamma)=\frac{(6-1)(6-2)}2=10             \tag{1.8}
\]

nodes.  The explicit Winger certificate removes any reliance on generic
language in the frozen positive control.

If a cubic \(Q\) passed through all ten nodes, the local intersection bound
at a node would give

\[
                    \Gamma.Q=18<10(2\cdot1)=20,            \tag{1.9}
\]

contradicting Bézout.  Thus \(h^0(-K)=0\).  If a second sextic \(G\) had
multiplicity at least two at all ten nodes, then

\[
                    \Gamma.G=36<10(2\cdot2)=40,            \tag{1.10}
\]

unless \(G\) contained \(\Gamma\).  Irreducibility and equality of degrees
then force \(G=\Gamma\).  Hence \(h^0(-2K)=1\).

Blowing up an ordinary node separates its two smooth branches.  The strict
transform \(B\) is therefore smooth and irreducible.  Its arithmetic genus
is zero by (1.3), so it is rational.  This proves (1.4).

## 2. What the numerical and proximity layers cannot see

For ten proper points, the proximity matrix is the identity and the
exceptional intersection matrix is \(-I_{10}\), independently of their
coordinates.  Contracting any subset of the labelled exceptional curves
does not change the self-intersection or mutual intersection of the rest.
Consequently the entire Boolean fan generated by these ten labelled
contractions is the same on \(S_{\mathrm{ell}}\) and \(S_{\mathrm{Cob}}\), not merely
its numerical shadow.

The Stage-B1-type formal divisor arrays are also identical:

\[
\begin{array}{c|c|c|c}
\text{class}&\text{plane degree}&\text{multiplicity vector}
 &\text{formal matrix shape}\\ \hline
-K&3&(1,\ldots,1)&10\times10,\\
-2K&6&(2,\ldots,2)&30\times28.
\end{array}                                               \tag{2.1}
\]

After coordinate specialization their ranks and kernels are

\[
\begin{array}{c|cc|c}
 &\operatorname {rk}M_3&\operatorname {rk}M_{6,2}
 &\ker M_{6,2}\\ \hline
z_{\mathrm{ell}}&9&27&\langle C^2\rangle,\\
z_{\mathrm{Cob}}&10&27&\langle\Gamma\rangle.
\end{array}                                               \tag{2.2}
\]

This proves the following no-go without a surface classification.

### Corollary 2.1 (enriched numerical/proximity no-go)

No predicate factoring only through

\[
 (\operatorname {NS},K),\quad
 \text{the proper-point proximity DAG},\quad
 \text{the labelled-exceptional contraction fan},\quad
 \text{degree/multiplicity arrays},\quad
 \text{formal jet descriptors}                            \tag{2.3}
\]

can decide whether \(-2K\) has a reduced irreducible representative.

This is deliberately not a same-nef-cone claim.  The complete stocks of
negative curves can differ between the two surfaces.  The conclusion is
exactly that enriching the raw lattice by numerical contraction and
proximity data still does not close the effectivity/reducedness gate exposed
by adapted `P10`.

## 3. The first exact discriminator: a fat-point Rees slice

Let \(Z=p_1+\cdots+p_{10}\) be the reduced ordered point scheme and write

\[
 \mathscr R_{\mathrm{fat}}(Z)
   =\bigoplus_{m,d\ge0}
      H^0\!\left(\mathbf P^2,\mathcal I_Z^{(m)}(d)\right)  \tag{3.1}
\]

for the relevant bigraded slice of its symbolic Rees algebra.  The two pieces

\[
 \mathscr R_{\mathrm{fat}}(Z)_{(1,3)}=H^0(-K),\qquad
 \mathscr R_{\mathrm{fat}}(Z)_{(2,6)}=H^0(-2K)            \tag{3.2}
\]

and the multiplication map

\[
 \operatorname {Sym}^2 H^0(-K)\longrightarrow H^0(-2K)   \tag{3.3}
\]

already separate the two branches.  On \(S_{\mathrm{ell}}\), (3.3) is an
isomorphism between one-dimensional spaces and its image is \(C^2\).  On
\(S_{\mathrm{Cob}}\), its source is zero and the degree-two generator is the
divisor-prime section \(\Gamma\).

For a general packet, absence from the image of (3.3) only excludes a square;
it does not exclude a product of two different divisor classes.  The correct
universal test is therefore divisor primality, equivalently geometric
integrality of the zero divisor, rather than merely “not a square.”

## 4. A classification-free Cox-prime criterion

The preceding discriminator has an intrinsic surface formulation.  No finite
generation of the Cox ring is required; only its indicated homogeneous
pieces and divisor multiplication are used.

### Proposition 4.1 (rank-eleven boundary criterion)

Let \(S/\mathbf C\) be a smooth projective rational surface with
\(K_S^2=-1\).  The following are equivalent.

1. There is a reduced irreducible divisor \(B\in|-2K_S|\).
2. One has

   \[
        H^0(-K_S)=0,\qquad H^0(-2K_S)=\mathbf C s_B,        \tag{4.1}
   \]

   and the zero divisor of \(s_B\) is geometrically integral.
3. In Cox-ring language, the degree \(-K_S\) piece is zero and the degree
   \(-2K_S\) piece is generated by a divisor-prime homogeneous element.

When these conditions hold, \(B\) is automatically a smooth rational
`(-4)`-curve and \(S\) is a one-boundary strong Coble surface.

### Proof

Assume (1).  If \(D\in|-2K_S|\), then

\[
                             D.B=B^2=-4.                    \tag{4.2}
\]

Two effective divisors without a common component have nonnegative
intersection.  Since \(B\) is irreducible, (4.2) forces \(B\) to be a
component of \(D\).  The residual effective divisor is linearly trivial and
is therefore zero.  Thus \(|-2K_S|=\{B\}\).  If
\(A\in|-K_S|\), then \(2A\in|-2K_S|\), so \(2A=B\), impossible because
\(B\) is reduced.  This proves (2).  The equivalence with (3) is the
divisor-section dictionary for the Picard-graded Cox ring.

Conversely, (2) supplies the required reduced irreducible divisor.  Its class
has square \(-4\), and adjunction gives

\[
                 p_a(B)=1+\frac{B^2+K_S.B}{2}=0.           \tag{4.3}
\]

A geometrically integral curve of arithmetic genus zero on a smooth surface
is smooth and rational.  The vanishing and uniqueness in (4.1) are exactly
the strong-Coble Hilbert conditions.  \(\square\)

### Corollary 4.2 (exact-nef handoff)

Let \(P\) be the ideal closure of a finite-sided finite-volume Coxeter
chamber in a chosen component of the positive cone of an integral Lorentzian
lattice \(L\).  Assume all of the following data from the maximal
strong-Coble return theorem.

1. The distinct oriented physical facet classes are
   \(B_1,\ldots,B_b,D_1,\ldots,D_f\), where \(b\ge1\), the \(B_i\) are
   pairwise orthogonal with \(B_i^2=-4\), and

   \[
                         \sum_{i=1}^bB_i=-2K.               \tag{4.4}
   \]

   Every \(D_j\) has the declared rational-curve adjunction row
   \((D_j^2,K.D_j)=(-1,-1)\) or \((-2,0)\), and these classes define all
   the displayed irredundant halfspaces of \(P\).
2. Every side reflection preserves the adapted host

   \[
      M=\left\langle L,B_1/2,\ldots,B_b/2\right\rangle,
      \qquad q_M=2q_L.                                    \tag{4.5}
   \]

3. The side group \(W(P)\) has no proper discrete reflection overgroup in
   \(O^+(M)\).

Let \(S\) be an oriented marked smooth projective rational surface carrying
\((\operatorname {NS}(S),K_S)\) to \((L,K)\).  Then \(S\) is an exact
physical full-nef carrier for this packet if and only if every displayed
physical facet class \(B_i,D_j\) has a divisor-prime homogeneous section on
this same \(S\).

Indeed, necessity is immediate from the meaning of a physical facet.  For
sufficiency, the numerical adjunction rows make the prime divisors smooth
rational curves.  Distinct effective integral curves have nonnegative
intersection, so the orthogonal boundary curves are disjoint; their reduced
sum is \(-2K\).  Intersecting any other member of \(|-2K|\) with each
\(B_i\) forces every \(B_i\) as a component, hence the displayed sum is the
unique member.  Its reducedness rules out a half in \(|-K|\), so the surface
is strong Coble.  The already proved
[`reflection-maximal strong-Coble return theorem`](reflection_maximal_coble_exact_nef_return.md)
then excludes every hidden negative wall and returns the literal complete nef
chamber.

This corollary is carrier-relative: it says exactly which geometric invariant
finishes Property 4 once a marked surface point has been constructed.  It
does not prove that the corresponding prime-section locus in the universal
cluster stack is nonempty for every raw packet.

## 5. Consequence for the POLYTOPE \(\longrightarrow\) SURFACE program

The adapted `P10` obstruction cannot be repaired by adding only a proximity
DAG or a numerical contraction fan to the raw chamber.  The smallest exact
plane-compiler target must retain the actual universal evaluation maps and
cut out a common locally closed locus on which

1. the \(-K\) matrix has the required rank;
2. the \(-2K\) kernel has the required dimension;
3. its generator has exact cluster multiplicities and is geometrically
   integral; and
4. every other proposed facet kernel is divisor-prime on the same point.

On a reflection-maximal standard-Coble branch, those four gates are enough:
the exact-nef return is then automatic.  What remains genuinely open is the
universal nonemptiness of this prime-section locus from raw Coxeter data.
This is a geometric rank/factorization problem, not a missing lattice
identity.

## 6. Exact scope

This theorem does not assert:

- that the two blowups have the same full nef cone or the same raw Coxeter
  chamber;
- that an arbitrary ten-point configuration lies on either branch;
- that proximity plus **actual coordinate-valued** evaluation matrices is
  insufficient—the actual matrices are precisely what distinguish the
  branches;
- finite generation of either Cox ring;
- a new complex Coble surface or a new abstract polytope; or
- global Property 4 or universal essential surjectivity.

It proves a sharper local no-go and a sufficient-and-necessary
carrier-relative invariant at the precise physical boundary where adapted
`P10` fails.

## Replay

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_rational_pfirst_bianticanonical_section_algebra_gate.py \
  --check-archive

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_rational_pfirst_bianticanonical_section_algebra_gate.py
```
