# Incidence-preserving bounded cluster recovery

The [source](../scripts/fresh_h20_near_closed_cluster_recovery_20261016.py) and [data](fresh_h20_near_closed_cluster_recovery_20261016.json) passed generation and a fresh saved replay. This is a numerical diagnostic; it proves no continuous exclusion.

All 100 previously saved acute clusters were screened using normals projected off both S16 and wall 151. Only nine have at least three potential finite neighbors for every new facet. The best acute cluster, step 62, fails: its family-0 cap at parameter −0.8244363976171009 has only old wall 95 and one new cap as potential finite neighbors.

Five eligible clusters (22, 26, 28, 63, 78) were tested. Each trial preserves every intrinsic finite adjacency, chooses from the two closest permitted Coxeter angles per finite pair, and tests all 121 accumulated sample points. The explicit limits are 32 angle assignments per cluster and 500 optimizer evaluations per assignment. None of the 160 attempts recovers a Coxeter cluster. These are local optimization failures, not infeasibility certificates.

The projected pair Gram being positive definite is only a necessary incidence condition. It does not certify a closed facet; a compact three-dimensional edge also needs two finite endpoint vertices. The next diagnostic therefore targets missing edge endpoints directly. The source/data/report of this prefix are frozen.
