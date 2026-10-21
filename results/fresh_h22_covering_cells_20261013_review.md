# Review of the covering criterion and contact filter

Mathematical source review and a fresh `--verify` run of the [checker](../scripts/fresh_h22_covering_cells_20261013.py) pass. The [main report](fresh_h22_covering_cells_20261013.md) distinguishes its proved finite criterion from an unperformed H21 cell enumeration.

The cell inequalities have the correct sign: minimizing `H_i(v)=−(b_i,v)` imposes `(b_j−b_i,v)≤0`. Each cell is a pointed polyhedral subcone of the future causal cone. Its extreme-ray decomposition and the reverse Lorentz triangle inequality prove sufficiency of the ray inequalities. Necessity at null rays follows by continuity of the minimum over the finite family, including when a cell has no timelike interior. Lower-dimensional cells and null rays cannot be omitted from a complete implementation. The projective timelike-center slice is compact, so the proposed exact H/V verification protocol has the required bounded domain.

The contact-graph obstruction correctly excludes parallel pairs as graph edges. Distinct graph components can share only finitely many tangencies; removing those cannot disconnect the open convex old chamber in dimension at least two. Every retained cap here covers a nonempty open region, so an isolated component cannot be ignored. This proves necessity of finite-angle bridges, while supplying no sufficient coverage condition.

An independent root calculation of the frozen H22 graph reproduced 325 upper caps, 43,900 finite-contact edges and component sizes `305,3,1,...,1` with 17 singletons. The checker separately verifies the rational parent graph and the omitted-contact interval gates, while explicitly inheriting the frozen new-cap finite identities. This review does not mislabel that inheritance as another full arithmetic reconstruction.

The standalone H3 example is also replayed: all nine facets, two finite vertices, two bounded ideal vertices, eight spacelike cone vertices, seven distinct minimum-cell rays and four failed new rays. Its strict escape has norm `9/16`. It refutes coverage based only on old finite vertices and inherited cusp links within the acute Coxeter setting itself.

The large-T generalization's same-side dotted bound and central-gap bound follow from the displayed squared inequalities with positive right sides. Finite many old vertices permit one common T satisfying all bounds. The conclusion is failure of the proposed shortcut, not a new finite-volume example or a nonexistence theorem for extensions.
