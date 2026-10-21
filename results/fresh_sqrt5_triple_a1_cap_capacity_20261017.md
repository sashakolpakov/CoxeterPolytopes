# Three-A1 cap capacities over the fixed field Q(sqrt5)

The exact rank13 domain gives **93 cap constraints, including12 violations** of the current incidence witness. Fourteen columns have positive load. These are additional constraints on a partial incidence model; this package does not repair the witness or prove a dimension bound.

The [source](../scripts/fresh_sqrt5_triple_a1_cap_capacity_20261017.py) and [certificate](fresh_sqrt5_triple_a1_cap_capacity_20261017.json) reconstruct all93 selected actual dotted-prism base occurrences whose deleted-wall complement has three internal isolated A1 components. They bind the [previous final witness](fresh_sqrt5_multispectator_cap_joint_ports_20261017.md). The inherited full prism-to-base maps retain the literal quadrilateral provenance, coface root, complete containing-simplex star, all spectator copies, and the base vertex set.

For mutually orthogonal norm2 spectator roots q,p,r, write the deleted cap normal as

\[
b=v-\frac{\sqrt\tau}{2}q-\frac{p+r}{2},\qquad v^2=1-\frac\tau2.
\]

The source reconstructs each core projection and checks this norm identity exactly. There are eight Lorentz cores and five values of tau. Core automorphisms preserve the simplex black/white coloring; an orientation is the simultaneous orbit of the complete contact vector and the base vertex subset. The counts include73 one-orientation and20 two-orientation columns.

For a fixed orientation, two different placements with the same q have mutual product1+|P∩P'|/2>0 and are incompatible. For distinct q,q', with P={p,r}, their product is

\[
1-\tau/2+j/2+c\sqrt\tau/2,
\quad j=|P\cap P'|,\quad c=1_{q\in P'}+1_{q'\in P}.
\]

The seven possible (j,c) values are(0,0),(1,0),(2,0),(0,1),(1,1),(0,2),(1,2). All35 pair tests check acuteness, algebraic integrality of the squared product, and membership of that squared product in the **fixed arithmetic ground field K=Q(sqrt5)**. The last requirement follows from the Vinberg two-cycle criterion in norm2 normalization. An algebraic integer outside K is still forbidden. The [independent field-gate review](fresh_sqrt5_spectator_groundfield_gate_20261017.md) reproduces these35 tests and the44 inherited two-A1 tests using separate symbolic field embeddings and minimal polynomials.

For every actual spectator count k, the saved compatibility graph contains all k·binom(k−1,2) positions(q,{p,r}). The symmetric group acts transitively, so a maximum clique can be chosen to contain(0,{1,2}). A complete include/exclude recurrence on its neighbors proves the exact graph maximum. The certificate saves the full adjacency graph, transitivity maps, recurrence children and a matching clique witness; no heuristic truncation or optimizer is used. Multiplying this maximum by the number of core orientations is a valid upper bound even when different orientations have further incompatibilities. All actual rank13 values of k are covered; no assertion about arbitrary rank is needed.

Generation and a fresh saved reconstruction replay both pass. A separate narrow source review confirms the pair formula, exhaustive seven patterns, transitivity reduction, recurrence and orientation factor; it does not duplicate the raw arithmetic reconstruction. Source and data are frozen with SHA256 respectively `18d5aee2a06876ea5184539a1119636d22faf0003d701de638237575a3154b99` and `f08967d709a412e26cb7fa9b2616e6fa2f288513f13b53283faf67e111df7f59`.
