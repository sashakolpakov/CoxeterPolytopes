# A reciprocal-square sharpening of the opposition cap

Date: 2026-10-18. This separate packet extends the frozen opposition-moments
packet without changing it. It proves the universal necessary condition

    2B <= (516992/3125) f = 165.43744 f

for any finite-volume right-angled twelve-polytope. Here B counts canonical
three-faces P3 and f counts finite vertices. No dimension improvement or
Lean verification is claimed.

Use the notation and proved moment identities of
`right_angled_mixed_opposition_moments_research_20261018.md`. On each finite
P3-opposition component, the vector h=3q-sum n_i is constant and spacelike.
Choose a future null z perpendicular to h, put E(q)=-z.q>0, and for all
220 formal triple flips let r=E(q_S)/E(q)=mu+d with mu=5/4. We have

    r>=1/2,  average d=0,  average d^2=9/176,
    -average d^3 <=27/1760,  average d^4<=3/440.

The exact identity

    1/r^2 = 1/mu^2-2d/mu^3+3d^2/mu^4-4d^3/mu^5
            + d^4(4r+mu)/(mu^5 r^2)

follows by multiplying by mu^5 r^2. The checker verifies every polynomial
coefficient after this multiplication. Since 4/r+mu/r^2 decreases for
r>0, its maximum for r>=1/2 is 8+5=13. Therefore

    average r^(-2)
       <= (4/5)^2 + (3/mu^4)(9/176)
                   + (4/mu^5)(27/1760)+(13/mu^5)(3/440)
        = 129248/171875.

Take g(q)=E(q)^(-2). Actual opposition neighbors form a subset of the 220
formal triple flips, so their ratio sum sum g(q_S)/g(q) is at most 220
times this bound. Summing over all finite vertices, each of the B
undirected opposition edges contributes a reciprocal pair of positive
ratios, whose sum is at least two. This proves the displayed cap.

The exact scalar target from the separate triangle-minimum packet is
6620/59, approximately 112.20339. This reciprocal-square cap remains above
that target; it does not close the aggregate dimension-twelve route.

Replay with
`python3 scripts/check_right_angled_mixed_opposition_power_two_research_20261018.py`.
The checker uses standard-library rational arithmetic only. Its saved JSON
certifies the polynomial and rational coefficient calculations. The
universal Lorentz/opposition reasoning is a written proof.
