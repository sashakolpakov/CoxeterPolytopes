# The characteristic-five Type-VI Coble full-nef chamber

**Status (2026-08-09): proved exact positive `SURFACE -> POLYTOPE`
theorem.**  Its inverse-carrier score is **N/A/L0**: the construction starts
from Kondō's surface and proves no `POLYTOPE -> SURFACE` implication.

**Directional update.**  This report's proof remains surface-first, but a
separate raw-wall theorem now extracts the unique `P2` frame, reconstructs
the characteristic-five Desargues/cut design, builds the integral ten-node
sextic, and returns this exact 78-wall chamber.  That companion receives
structural P-first `L3` credit; see the
[Desargues inverse theorem](coble_type_vi_char5_desargues_inverse.md).

## The theorem

Let `S` be the unique strong Coble surface of Type VI over an algebraically
closed field of characteristic `5`.  The projectivized full
Néron--Severi nef cone of `S` is the finite-volume, noncompact arithmetic
Coxeter polytope

```text
COBLE-VI-p5-b1-H10-78  in H10.
```

Its exact census is

| facets | primitive wall norms `-1,-2` | finite vertices | ideal vertices | hyperideal rays | edges |
|---:|---:|---:|---:|---:|---:|
| 78 | `58,20` | 21,661 | 1,797 | 0 | 166,320 |

All 78 walls are irredundant.  Their primitive normals span

\[
                  \operatorname{NS}(S)\cong I_{1,10},
\]

all 78 reflections are integral involutive isometries, and every pair of
walls has a Coxeter angle or is parallel/ultraparallel.  This is the fourth
one-boundary Coble full-nef chamber in `H10` completed exactly in this
project, after Types I, II, and VII.

The reproducible artifacts are:

- [summary certificate](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json);
- [complete ray/incidence archive](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_rays.jsonl.gz);
- [complete edge/endpoint archive](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_edges.jsonl.gz);
- [certifier](../scripts/certify_coble_type_vi_char5_one_boundary_full_ns_h10.py), reusing the input-generic [finite-diagram enumerator](../scripts/certify_coble_type_vii_char5_one_boundary_diagrams.cpp);
- [regression and full-replay tests](../tests/test_coble_type_vi_char5_one_boundary_full_ns_h10.py).

## Why these are all nef walls

Kondō's Theorem 1.1/Table 2 gives a unique Type-VI one-boundary Coble
surface in characteristic `5`.  Example 4.12 constructs it as a quotient of
the supersingular `K3` surface of Artin invariant one and exhibits twenty
`(-2)`-curves with graph `Gamma_VI`; Proposition 4.13 gives

```text
Aut(S)=S5,       R(S)=(E6+A4,{0}).
```

Write

\[
                   B\in|-2K_S|,\qquad b=B/2=-K_S.
\]

The one boundary is a rational `(-4)`-curve, so `b^2=-1`.  Kondō's
Proposition 3.1 gives the integral splitting

\[
 \operatorname{CM}(S)=b^\perp\cong E_{10},\qquad
 \operatorname{NS}(S)=E_{10}\mathbin\perp\langle-1\rangle
                     \cong I_{1,10}.                         \tag{1}
\]

The transport from Kondō's labelled `Gamma_VI` to Shimada's exact `20E`
coordinates is not based on an unweighted graph resemblance.  Their full
twenty-by-twenty root Gram matrices agree entrywise.  Moreover the
twenty-root coordinate matrix in `E10` has rank ten and Smith diagonal

```text
1,1,1,1,1,1,1,1,1,1.
```

Its HNF basis has determinant one, and its Gram has determinant `-1` and
Smith diagonal `1^10`.  Thus the twenty roots already generate all of
`E10`; their Gram isometry extends integrally with no discriminant glue or
overlattice choice.  This is a stronger marking bridge than a graph match.

Kondō's finite-index root argument (Proposition 3.7) proves that these
twenty curves are the complete effective irreducible-root list.  Since
there is only one boundary, all twenty are actual `(-2)`-curves.

The `(-1)`-curves are recovered from every cusp, rather than sampled.  By
Kondō's Lemma 3.11, each primitive isotropic ray `f` of the complete
Coble--Mukai chamber gives a genus-one fibration.  Proposition 3.8 says that
the one-boundary case comes from an index-two Halphen surface by blowing up
the singularity of one reduced `I1` or `II` fibre.  Its exceptional curve
`E` obeys

\[
                         B+2E=2f,
                    \qquad E=f-b.                            \tag{2}
\]

Conversely, adjunction and the one-boundary Halphen description send every
`(-1)`-curve back to the primitive nef cusp `f=E+b`.  Hence (2) is a
bijection

```text
primitive cusps of the CM chamber  <-->  all (-1)-curves on S.
```

Exact double description of the known `20E` Coble--Mukai chamber gives 57
primitive cusps.  The complete full-nef wall list is consequently

\[
 \{b\}\ \cup\ \{R_1,\ldots,R_{20}\}\ \cup\
 \{f-b:f\text{ is one of the 57 primitive CM cusps}\},       \tag{3}
\]

with `1+20+57=78` members.  In geometric terms these are one boundary
`(-4)`-curve, twenty `(-2)`-curves, and fifty-seven `(-1)`-curves; the
boundary mirror uses the primitive norm-`-1` normal `b`.  Kondō's
negative-curve trichotomy says that boundary, `(-2)`, and `(-1)` curves are
the only negative curves on a Coble surface.  Formula (3) exhausts all
three classes, so the cone certified below is the full nef cone, not a
chosen reflection subcone.

## Coxeter and arithmetic checks

The complete primitive wall stock in the marking (1) has

```text
norm -1: 58       norm -2: 20       wall rank: 11
Smith span: 1^11  (all of I_(1,10)).
```

Its exact pair-intersection profile is

```text
(-1,-1): pairing 0/1/2/3/4 = 435/612/390/150/66
(-2,-1): pairing 0/1/2/3   = 600/330/210/20
(-2,-2): pairing 0/1/2     = 135/45/10.
```

Thus the normalized Coxeter census is

```text
m2=1,170; m3=45; m4=330; parallel=622; ultraparallel=836;
non-Coxeter=0.
```

Every reflection matrix is checked entrywise for integrality,
involutivity, and preservation of the lattice Gram.  Together with the
finite-volume result below, this makes the group an arithmetic reflection
subgroup of the classical `I_(1,10)` reflection group.

## Literal Vinberg finite-volume certificate

The finite-volume assertion is not inferred from noncompactness or from a
partial vertex list.  Both parts of Vinberg's criterion are checked
literally, including disconnected finite edge diagrams and disconnected
affine cusp diagrams.

First, `cddexec_gmp --repall` computes the complete rational
V-representation, incidence, and adjacency:

```text
23,458 extreme rays = 21,661 finite + 1,797 ideal + 0 hyperideal.
```

Every facet's incident rays span a rank-ten face.  Every finite vertex is
simple.  All cusps are nonsimple, with the complete incidence census

| incident walls | finite vertices | ideal vertices |
|---:|---:|---:|
| 10 | 21,661 | 0 |
| 11 | 0 | 480 |
| 12 | 0 | 772 |
| 13 | 0 | 460 |
| 14 | 0 | 75 |
| 15 | 0 | 10 |

The sum of the coherently oriented rays is a positive-square strict
interior point, with positive slack at every wall.  Its exact Lorentz
pairing with every causal extreme ray has the same strict sign; the minimum
is `362134>0`.  Thus all rays lie on one component of the positive cone.

Second, a compiled exhaustive classifier tests every pair-compatible
rank-nine and rank-ten mask componentwise against the full finite
crystallographic Coxeter list:

| order | compatible masks screened | elliptic diagrams |
|---:|---:|---:|
| 9 | 26,769,085 | 166,320 |
| 10 | 40,619,993 | 21,661 |

The 166,320 elliptic rank-nine masks equal the cdd geometric-edge supports
set-for-set.  The 21,661 elliptic rank-ten masks equal the finite endpoint
supports set-for-set.  Containment in all finite and ideal active sets gives
exactly two endpoints for every edge:

```text
FF=69,210; FI=78,190; II=18,920; bad edges=0.
```

Third, the parabolic pass attaches every compatible nonorthogonal node to
every connected finite component through order twelve.  Exact
determinant-zero and positive-principal-deletion tests find all 5,551
connected affine components:

```text
order 2/3/4/5/6/7
count 622/1,890/1,410/804/625/200.
```

Every one of the 1,797 cusps is disconnected.  The number of connected
affine components per cusp is

```text
2/3/4/5/6 components: 480/772/460/75/10 cusps.
```

For each connected affine component `J`, its complete orthogonal
zero-neighbourhood `N(J)` is one of these exact ideal active sets.
Conversely every ideal active set decomposes into affine components of
total Coxeter rank nine, and every component recovers the whole set as its
zero-neighbourhood.  There are zero bad parabolic neighbourhoods.  This is
precisely the check that distinguishes finite-covolume noncompactness from
an infinite-volume chamber with escaping ends.

Therefore the elliptic-edge endpoint criterion and the full parabolic
zero-neighbourhood criterion both pass, proving finite volume.  The 1,797
ideal vertices prove noncompactness.

## Euler characteristic, volume, and index

An independent connected-elliptic/Solomon computation gives

\[
                    \chi(W_{78})=-\frac{527}{6480}.
\]

The dimension-ten Gauss--Bonnet formula gives

\[
 \operatorname{vol}(P_{78})
   =-\frac{32\pi^5}{945}\,\chi(W_{78})
   =\frac{1054\pi^5}{382725}.                               \tag{4}
\]

Since `chi(W(P10))=-1/11147673600`, the exact reflection-subgroup index is

\[
                        [W(P10):W_{78}]=906608640.            \tag{5}
\]

The same number already occurs for Shimada's known Type-VI `20E` induced
chamber in `H9`: his published `volindex` is
`2^7*3^4*5=51840`, and
`46998591897600/51840=906608640`.  Equality of the seed and lifted indices
is structural evidence for a general one-boundary cusp-extension theorem;
it is not novelty evidence.

## Exact lineage boundary

The claims separate as follows.

| layer | status |
|---|---|
| Kondō Type-VI characteristic-five one-boundary surface | known |
| twenty-root `Gamma_VI`/Shimada `20E` CM chamber in `H9` | known |
| the 57-cusp count | exact invariant computed from that known seed; no standalone novelty claim |
| ambient `I_(1,10)` arithmetic class and `W(P10)` | classical |
| first finite-volume Coxeter polytope in dimension ten | no |
| complete 57-cusp-to-`(-1)` lift and full `H10/78` nef chamber | new exact positive result in this project ledger |
| literal abstract `H10/78` worldwide priority | open; not claimed |
| `POLYTOPE -> SURFACE` inverse theorem or dimension bound | none; `N/A/L0` |

At CoxIter commit `b800be48d8240fbabd2d1c9542e61054584ad7fe`, a bounded
header scan finds no `H10/78` entry among 1,214 ordinary records or 2,959
records in the doubling archive.  Exact-phrase searches likewise produced
the known Type-VI/`20E` material but no printed 78-wall chamber.  These are
only bounded discovery checks, not a worldwide nonexistence or priority
proof.

In particular, this theorem must not receive inverse-carrier credit.  It
starts from a classified surface, transports its complete `20E` CM cone,
and adds the algebro-geometrically forced cusp lifts.  A genuine
`POLYTOPE -> SURFACE` theorem would instead need a surface-free non-singleton
class of raw polytopes, an exhaustive carrier construction for that class,
and algebraic geometry pushed back to every accepted input.

## Sources

- Shigeyuki Kondō, [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591), especially Theorem 1.1/Table 2, Propositions 3.1, 3.4, 3.7, and 3.8, Lemma 3.11, Example 4.12, and Proposition 4.13.
- Simon Brandhorst and Ichiro Shimada, [*Borcherds' method for Enriques surfaces*](https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/L10L26/PaperL10L26/L10L26.pdf), especially the `20E` row of the induced-chamber table.
- E. B. Vinberg, [*Hyperbolic reflection groups*](https://www.mathnet.ru/eng/rm2140), Russian Mathematical Surveys 40:1 (1985), especially Proposition 4.2.
