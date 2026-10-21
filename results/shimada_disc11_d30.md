# The exact physical \(D_{30}\) K3 eigensection

**Status.** Complete local reduction and internal first-star traversal.  This
is a theorem about a divisor-defined section of the nef cone of Shimada's
discriminant-11 singular K3 surface.  It is not yet a global eigen-nef,
finite-volume, or compactness theorem.

## 1. A nontrivial two-path identity

There are two previously certified routes from the physical chamber
\(D_{68}\) to representative \(D_{30}\):

\[
\begin{array}{ccc}
&D_{68}&\\[-2mm]
{}^{(14,0)}\swarrow&&\searrow^{(11,1)}\\
D_{181}&&D_{86}\\
{}_{(10,0)}\searrow&&\swarrow_{(23,0)}\\
&D_{30}.&
\end{array}
\]

The two cumulative integral row-isometries are **entrywise equal**.  Thus
these are not merely isomorphic sections: both walks reach the same physical
copy of \(D_{30}\) in the marked Neron--Severi lattice.  Representative
\(D_{30}\) has chamber-stabilizer order one, so there is no hidden stabilizer
ambiguity.  The common transport has SHA-256 fingerprint

```text
89ee0f1b46508b99b776d83b38d9c70bc7293a142f05c84f7d5256ba21810a6c
```

This closes the first nontrivial square in the restricted induced-chamber
complex.  In algebro-geometric language, the local divisor cones obtained by
the two wall-crossing paths glue to the same marked cone.

## 2. Complete restricted half-space theorem

Fix the transported ten-dimensional golden eigenspace used in the
\(D_{68}\), \(D_{181}\), and \(D_{86}\) calculations.  Restricting every
primitive divisor wall of the common physical \(D_{30}\) gives

\[
\boxed{
34\text{ ambient divisor walls}
\longrightarrow 26\text{ positive rays}
\longrightarrow 12\text{ facets}.}
\]

The facets split as

\[
\boxed{8\text{ effective-root facets}+4\text{ internal Borcherds facets}.}
\]

The eight root classes are

\[
\begin{gathered}
\{(0,0),(2,0)\},\quad \{(1,0),(15,0)\},\quad
\{(3,0),(16,0)\},\quad \{(4,0),(9,0)\},\\
\{(5,0),(6,0)\},\quad \{(7,0),(17,0)\},\quad
\{(8,0),(14,0)\},\quad \{(11,0),(12,0)\}.
\end{gathered}
\]

The ratios of the second restriction to the first are, in the displayed
order,

\[
\phi-1,\ \phi,\ \phi,\ \phi-1,\ \phi-1,\ \phi,\ \phi-1,\ \phi-1.
\]

The four internal facets are precisely

\[
\boxed{(19,0),\quad(20,0),\quad(22,0),\quad(26,0).}
\]

An exact strict interior point has trace square

\[
-\frac{1988609}{7200}+\frac{107631}{400}\phi>0
\]

and minimum primitive-wall pairing

\[
\frac{481}{60}-\frac{99}{20}\phi>0.
\]

For each of the twelve asserted facets the certificate archives a point of
\(\mathbf Q(\phi)^{10}\) which is timelike, vanishes on exactly the target
positive-proportional class, and is strictly positive on every other
primitive wall.  Thus every displayed inequality is a genuine hyperbolic
facet; floating-point convex-hull output is not used in the proof.

## 3. Exact redundancy and the curve-wall distinction

The remaining fourteen restricted rays have exact positive
\(\mathbf Q(\phi)\)-linear expressions in the twelve facet representatives.
For example, writing \(w_i\) for wall \((i,0)\),

\[
w_{13}=(\phi-1)(w_{19}+w_{26}),                                      \tag{1}
\]

\[
w_{18}=\phi w_5+\phi w_{20}+(\phi-1)w_{22},                           \tag{2}
\]

and

\[
w_{31}=w_5+\phi w_{20}.                                               \tag{3}
\]

The executable certificate records and checks all fourteen identities
coordinate by coordinate.  All coefficients are strictly positive at the
distinguished real embedding.

There are 18 ambient `_minustwo` walls.  Exact recovery using the inverse
Neron--Severi Gram matrix proves that every one is an integral class of square
\(-2\).  Under the golden structure:

- 17 are integral field-reflection roots;
- their shells are eight roots of norm \(2\) and nine of norm
  \(4-2\phi\);
- they span field rank 9;
- \((10,0)\) is the unique non-field root; and
- \((13,0)\) is a genuine field root but is hidden behind the two internal
  walls in (1), so it is not a facet of this induced section.

The eight root facets account for 16 ambient curves in the expected paired
lift.  The last point is geometrically useful: being an effective
\((-2)\)-curve and being visible as a facet of one Borcherds induced chamber
section are different conditions.  Internal chamber walls can temporarily
shield a genuine root wall; it reappears only after the induced chambers are
pasted.

## 4. Crossing every internal facet

Exact opposite-side timelike points cross all four internal facets:

| \(D_{30}\) wall | adjacency | target representative | crossing scale |
|---:|:---:|---:|---:|
| 19 | `_isom` | \(D_{16}\) | \(1/10\) |
| 20 | `_isom` | \(D_{86}\) | \(1/3\) |
| 22 | `_isom` | \(D_{61}\) | \(1/3\) |
| 26 | `_newcham` | \(D_{181}\) | \(1/2\) |

At each shared face there is exactly one zero and no negative primitive-wall
pairing on either side.  Each opposite-side point is timelike and pairs
strictly positively with every primitive wall of its target chamber.

Moreover, the target at wall 20 is exactly the already certified physical
\(D_{86}\), and the target at wall 26 is exactly the already certified
physical \(D_{181}\).  Hence the square in Section 1 closes in both
directions.  The only new forward representatives are

\[
\boxed{D_{16}\quad\text{and}\quad D_{61}.}
\]

## 5. Algebro-geometric meaning and scope

This calculation stays entirely on the proposed algebro-geometric track.
The singular K3 supplies the Neron--Severi lattice and its positive cone;
effective \((-2)\)-divisors supply the genuine nef boundaries; Shimada's
induced chambers give a finite local atlas; and restriction to the fixed
golden eigenspace produces exact semialgebraic hyperbolic cells.  The
internal Borcherds facets are atlas seams, while the paired root facets are
the geometric curve walls that survive after gluing.

What remains open is global: one must continue past \(D_{16}\), \(D_{61}\),
and the other branch-frontier chambers until either the root-only pasted cell
closes with a finite-volume/compactness certificate or an exact obstruction
appears.  The present theorem proves one complete local cell and a commuting
gluing square, not termination of that process.

The external source is Shimada's `disc11.txt`, SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

Reproduce with

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d30.py \
  /path/to/disc11.txt
SHIMADA_DISC11=/path/to/disc11.txt PYTHONPATH=src:. \
  pytest -q tests/test_shimada_disc11_d30.py
```

The verifier pins the physical transport, all 34 transported wall rows, and
all exact restricted functionals by independent SHA-256 fingerprints.
