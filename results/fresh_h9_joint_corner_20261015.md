# Connecting the joint H9 corner to an older finite-edge component

An exact new cap joins the joint H3+D6 corner to a previously separate two-vertex component. The resulting chamber has **33 actual facets, 92 finite vertices and 253 finite edges**. Its finite-edge graph has seven components, down from eight. The same strict escape remains, so the chamber still has infinite volume. This is a connection between existing finite corners, not a finite-volume example.

## Exact new-corner frontiers and bounded scouts

The frozen [32-wall chamber](fresh_h9_joint_pair_20261015_two_extension.json) has a joint corner requiring its last two walls. The [frontier checker](../scripts/fresh_h9_joint_corner_20261015_germs.py) reconstructs its nine open directions exactly and verifies that each passes the light cone before meeting another old wall. The [frontier certificate](fresh_h9_joint_corner_20261015_germs.json) records all direction vectors, norm polynomials, strict spacelike witnesses and complete spherical endpoint profiles.

The nine rank-eight supports have types A2+D6, A1+A1+D6, D6+I2(5), A1+A1+A3+H3, twice A5+H3, D5+H3, A1+D4+H3, and A2+A3+H3. They admit 122 finite contact profiles and five symbolic dihedral tails. No claim is made that the later numerical sampling exhausts those tails.

The [rescore](fresh_h9_joint_corner_20261015_rescore.json) finds no proposal among the preceding 4,566-record bank compatible with both joint caps. The [fresh endpoint scout](fresh_h9_joint_corner_20261015_scout.json) instead saves 2,267 numerical proposals on the new corner's nine directions. None cuts the current escape in that bounded test. Generation and saved replay pass for both calculations.

A [finite-edge graph screen](fresh_h9_joint_corner_20261015_graph_screen.json) identifies six records that would join the joint corner to an older component. The screen uses numerical positive-definite Gram tests and does not certify those six walls. Its selected label-3/label-5 candidate is reconstructed exactly below. Numerical records can duplicate a geometric proposal, and no arbitrary-cap exclusion follows from these counts.

## The exact connecting cap

Write the existing field as

    K = Q(√5, √D),    D = 6 + (16/5)√5.

The cap has square two, label 3 at wall 30, label 5 at wall 8, and right contacts at walls `[0,1,2,5,6,10,29]`. These nine linear contacts give a particular solution p orthogonal to their one-dimensional kernel v, with v's final coordinate normalized to one. The norm equation requires

    δ = 1988340828/1990805 + (4448872824/9954025)√5,
    normal = p ± √δ v.

The [exact checker](../scripts/fresh_h9_joint_corner_20261015_bridge.py) selects the unique compatible orientation and verifies all 32 new pair relations. It also proves δ is not a square in K: neither δ nor δ/D is a square in Q(√5), as their rational field norms are not rational squares. Hence the construction uses the genuine degree-eight field K(√δ). No general restriction on finite labels in this extension is inferred.

The [certificate](fresh_h9_joint_corner_20261015_bridge.json) checks strict-center projection witnesses for all 33 facets. All 90 old finite vertices remain strict against the new wall. Exactly nine old walls have finite angle with it, leaving nine possible eight-wall supports for a new finite vertex. Exact positive-definiteness and full chamber inequalities leave precisely

    [0,1,2,5,6,8,10,29,32],
    [0,1,2,5,6,10,29,30,32].

Together they add three finite edges, joining old vertices 19 and 20 to the joint corner, old vertex 21. The joined component has five vertices after the two additions. The full finite-vertex and finite-edge census is therefore 92 and 253, with 322 open finite-corner directions. This local enumeration is exhaustive for the new cap; no full census of exterior rays or ideal vertices is claimed.

The old strict spacelike escape has strictly negative pairing with the cap. Its inherited and new inequalities therefore prove infinite volume without needing a new exterior-ray census.

The [independent review](fresh_h9_joint_corner_20261015_bridge_review.md) is responsible for separate arithmetic and geometry verification; its status is recorded in that report and the phase checkpoint. The primary generation and saved replay already pass.

## Scope of the progress

The connection reduces the number of finite-edge components, which is relevant to the compact H9 target. It does not close the chamber or its open directions, and the principal 81-vertex component is still separate. The retained escape and explicit frontier geometry remain the next obstructions. Repeating bounded scouts or adding isolated facets does not establish compactness or finite volume.
