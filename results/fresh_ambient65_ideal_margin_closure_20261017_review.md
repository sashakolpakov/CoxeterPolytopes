# Independent reconstruction of the ambient-65 cusp margins

All **1,650 terminal profiles** and every saved scalar field agree with an independent standard-library calculation. The unique minimum normalized margin is **61/1936**, at ranks (9,10,10,10,10,15). The unique minimum absolute margin is **445/23**, at the single rank64 factor.

The [independent checker](../scripts/fresh_ambient65_ideal_margin_closure_20261017_review.py) imports no producer or project arithmetic. Its [certificate](fresh_ambient65_ideal_margin_closure_20261017_review.json) records compact scalar rows, the two complete extremal rows, and matching SHA-256 digests of all complete primary and independently reconstructed rows. It compares the scalar data in the [primary certificate](fresh_ambient65_ideal_margin_closure_20261017.json).

The enumeration differs from the producer's sorted-partition recursion: it chooses an optional low rank1..9 and an optional high rank18..64 first, then enumerates all nonnegative multiplicity vectors at ranks10..17 with the remaining total. The resulting1,650 rank64 profiles are distinct, and their set equals the complete saved row set.

For each profile, integer sliding-window convolution computes all65 coefficients of H=∏(1+x+⋯+x^r). The checker verifies symmetry and V=H(1)=∏(r+1). It computes G using symmetric cumulative-distribution tails,

    G = Σ[k=1..32] (V − 2 Σ[j=0..k−1] H_j),

and separately checks the absolute-deviation identity

    2G = 64V − Σ[j=0..64] H_j |64−2j|.

The exact rational cost schedule is r/2 through10; then16/3,74/13,6,94/15,107/16,7 at ranks11..16, and163/23 at every rank≥17. The curvature is (5G−64V)/2. Every H coefficient, V, G, cost, curvature, normalized margin and absolute margin is compared with its primary row; all margins are strictly positive. Both stated minima and their unique attaining profiles are independently recovered.

Generation and a separate fresh `--verify` invocation both terminated exit0/PASS. This is an independent scalar arithmetic/interface check. It does not repeat geometric censuses, the all-rank concentration proof, the finite-coupon budget, or the changing additive global-constant argument. The whole primary source/report is therefore not bound as an independently reviewed global theorem; the complete scalar rows are bound by their canonical content digest. This distinction permits the separate global-constant proof to be finalized without implying that it was audited here.
