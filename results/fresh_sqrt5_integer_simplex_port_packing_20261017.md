# Integer capacities for disjoint maximal-simplex port blocks

The actual-provenance rank13 dual satisfies the old port capacities but
violates stronger integer packing inequalities at41 vertex forests.
Adding every three-port matching inequality and the full-thirteen-port
integer capacities leaves an exact repaired dual

\[
 \frac{2228350155976161360}{184610141382384827}
 =12.0705728260\ldots>12.
\]

This is a partial-model obstruction, not a dimension theorem. It uses
the frozen [actual provenance split](fresh_sqrt5_actual_provenance_dual_20261017.md).
It does not yet combine these new columns with later vertex-labelled
prism-base refinements.

## The global incidence inequality

At a finite vertex of a compact Coxeter polytope, the port sets C of
distinct maximal simplex faces are pairwise disjoint by the
[disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md).
For any set U of ports and integer t>=2,

\[
 \sum_C\left\lfloor\frac{|C\cap U|}{t}\right\rfloor
 \le \left\lfloor\frac{\sum_C|C\cap U|}{t}\right\rfloor
 \le \left\lfloor\frac{|U|}{t}\right\rfloor.          \tag{1}
\]

This is an integer packing condition, not merely the sum of the
one-port inequalities. Singleton maximal-edge blocks contribute zero.

For a three-port set U and t=2, a simplex block contributes one exactly
when it occupies at least two of U. Disjointness allows at most one such
block. Its indicator is

\[
 {\bf1}_{|C\cap U|\ge2}
 =\sum_{p\in\binom U2}{\bf1}_{p\subset C}
       -2{\bf1}_{U\subset C}.                    \tag{2}
\]

The values at intersection sizes0,1,2,3 are0,0,1,1. Thus all new
coefficients are nonnegative, despite the subtraction in this formula.
Quadrilateral/prism rows receive zero for this new charge family.
The existing mixed-three-port charges have a different indicator:
they require at least two occupied pairs; in particular they give
zero to a triangle corner. Equation(2) gives a triangle corner a
positive coefficient whenever its two ports lie in U.

The second family takes U to be all13 ports and t=2,3,4. Its respective
vertex costs are6,4,3. A d-dimensional simplex corner contributes
floor(d/t). In this compact domain d<=4, so higher t contributes zero.

These inequalities can be assigned separate nonnegative charge
coefficients and added to existing charges. Their separate vertex
upper bounds add; no disjointness between charge families is needed.

## Exact symmetry averaging and reconstruction

The [checker](../scripts/fresh_sqrt5_integer_simplex_port_packing_20261017.py)
reconstructs each complete13-node finite forest. It enumerates all78
literal port pairs and all286 literal triples and groups them by labelled
diagram automorphisms. This averaging makes no assumption that the
ambient chamber has those symmetries.

For a triple orbit O and pair orbit p, let h(O,p) count members of O
containing a fixed member of p. The count is constant on p. The program
reconstructs all incidences, checks exact divisibility by the pair-orbit
size and verifies

\[
 \sum_p |p|h(O,p)=3|O|.
\]

For each positive simplex row, its new orbit coefficient is the sum of
its old pair coefficients times h(O,p), minus twice its full-triple
coefficient. The latter is supplied by the frozen complete positive-row
feature export. Every resulting coefficient is checked nonnegative.
A corner with d ports contributes in total

\[
 \binom d2(13-d)+\binom d3
\]

over all literal triples; the code checks this identity for every
simplex row. This tests both the subtracted complete triples and the
new contributions from triangles.

The full-port family independently reconstructs the number of each
forest's vertices on a simplex from its one-port incidences: divide
the incidence total at that forest by d. The result must be integral,
and the total over all corners is d+1. All14,433 corner multiplicities
of the5,462 positive simplex rows are checked.

## The certificate

The [saved output](fresh_sqrt5_integer_simplex_port_packing_20261017.json)
exports62,353 new columns:59,962 triple-orbit matching columns and2,391
full-port columns, together with their coefficients on all5,462 positive
simplex rows. It detects66 violated columns at41 vertex forests.
Each vertex mass is increased to the maximum required by its new
columns; all face masses and all1,858 actual provenance masses stay
unchanged.

The checker then recounts all211,339 inherited positive-cost columns,
all59,962 existing mixed columns, every new column, all5,693 conservative
signed transfers, and all4,969 actual shared-base inequalities. The
unchanged actual provenance split remains valid. All checks pass after
the vertex-mass correction, giving the exact ratio above. No optimizer runs.

Outside the positive support every load is zero. The exported feature
matrix does not reconstruct zero-mass simplex rows for a full-family
primal solve. Restricting a new charge to the explicitly exported face
classes is itself a valid weaker charge, but must be stated explicitly
if used in optimization. Other subset sizes and the complete disjoint-set
packing polytope are not included.

**FINAL/FROZEN.** Generation66524 and fresh saved replay72709 both
terminated with exit0/PASS. Syntax, source hashes and local links were
checked. The frontier_literature agent independently reviewed the integer
inequality, orbit normalization, coefficient formula, code and stated
scope; that internal review passed. It is not a second replay of all
62,353 columns or external refereeing. The numerical scout that first tested
only the three full-port capacities found one violation and a weaker
ratio12.17024420; the present exact result includes every triple orbit.
No frozen predecessor was edited.

