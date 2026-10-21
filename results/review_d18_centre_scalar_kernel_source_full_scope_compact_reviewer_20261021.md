# D18 centre scalar kernel: full source and interface review

Verdict: **SOURCE/MATHEMATICS PASS, uncompiled**. This review reads all89 Lean lines and all92 assumption-map lines. It does not certify elaboration, an object file, runtime or actual axiom output. No compiler, import execution, mathematical worker, test, Lean or TeX job was run by this reviewer.

The reviewed files are:

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `formal/CoxeterProofs/D18CentreScalarKernel.lean` | 3568 | `998ae907c17886fd0c72a7bce74ce36126d128b6e51453966982025410cc9eb6` |
| `docs/d18_centre_scalar_kernel_logic_20261021.md` | 4538 | `dc257720eb80b1fef4f6562698a09524dfc621deea13600920fe2b5e3f6cf846` |

Disclosure: I authored the written D18+actual-I geometric application, `results/compact19_d18_isolated_i_complete_centre_and_head_exclusion_candidate_20261021.md`, current SHA-256 `8b8e4b90f4059a85778fe0a9c033d80c10d02a932d877b0d993dc588973bed73`. I did not author this new Lean module or assumption map. This is an independent review of their source and scalar translation, not an outside-author re-review of my own geometric theorem.

The five declarations have the following complete mathematical scope.

1. **Exact identity, lines28–31.** Put A=96,874,800 and b=1,528,000. The multiplier is4A and the square coefficient is2A, so the quadratic and linear coefficients agree identically. The constant difference is
   `387499200*6100 - 1528000^2 = 2363745120000 - 2334784000000 = 28961120000`.
   Thus the displayed square completion is correct in a commutative ring, with no numerical assertion supplied as a hypothesis.
2. **Positivity, lines34–48.** The ordered-ring square theorem supplies nonnegativity of the actual square. The positive remainder makes the scaled quadratic strictly positive. If the quadratic were nonpositive, multiplication by the nonnegative multiplier would make that same product nonpositive, a contradiction. The code obtains square sign and multiplication sign from the ordered-ring laws; neither is a geometric input. Source inspection of the pinned v4.32.2 `Init/Grind/Ordered/Ring.lean` confirms the implicit argument of `OrderedRing.sq_nonneg` and the two argument signs of `mul_nonpos_of_nonneg_of_nonpos` match these uses. This is source inspection, not a tactic or API execution test.
3. **Continuous scalar bound, lines51–57.** Subtracting34,190,400 from the supplied analytic lower polynomial leaves exactly the positive quadratic, since34,196,500−34,190,400=6100. Hence `1000000*F>34190400`, equivalently `10000*F>341904`. No parameter interval or claimed uniform margin is a premise of `quadratic_positive`; the analytic inequality itself is explicitly primitive in this declaration.
4. **K10 margin, lines60–69.** The given coefficient bounds and the newly derived strict F bound imply
   `49*T > 49*.87866 + 34.1904 - 14*2.731 + 10 = 49.01074`.
   Multiplication by100,000 gives the exact intermediate bound `4900000*T>4901074`, which strictly exceeds the value at T=1. The height upper bound has the correct sign because its coefficient is−14.
5. **K11 margin, lines72–81.** The separate values give
   `49*T > 49*.90246 + 34.1904 - 14*2.8854 + 11 = 49.01534`,
   hence `4900000*T>4901534` and T>1. This branch uses its own base and height bounds; it does not silently reuse the K10 constants.

The module's type and order assumptions, lines21–25, describe an arbitrary linearly ordered commutative ring with compatible strict order and a positive unit. All arithmetic is in that type. There is no truncating natural subtraction, integer stand-in for real heights, division, square-root operation, extra axiom declaration, `sorry` or native evaluation. The polynomial identity is a ring-normalization obligation; subsequent steps are polynomial or linear ordered-ring consequences. I find no source-level mathematical or API blocker, while leaving actual elaboration and tactic performance to a separately admitted compilation.

The application map, especially lines51–92, correctly draws the conditional boundary. The written Section5 radical estimate supplies `analytic_lower` after x=1/b and z=x−7/6; the chord/tangent estimates, radical signs and rational square-root comparisons are not formalized here. The written paid reference (application lines98–107), multiplied by49, has exactly the Lean form `49*T>=49*base+F-14*f+K` with F=14G−R. For strong5, the once-counted joint-I gain gives base>.8244+.05426=.87866 or base>.8482+.05426=.90246. For weak5, the bare references>.8870 and>.9104 separately exceed those respective floors and do not require a second I charge. The application supplies the positive endpoint coefficient, actual all-wall feasible core, reference orthogonality and original-vector norm interpretation. None is derived by this module, and none has been replaced with a supplied T>1 premise.

There are exactly five named `#print axioms` requests, lines85–89, matching the five declarations. At this reviewed preparation stage there are **zero actual outputs** and no execution PASS. The map does not claim a D18 exclusion or dimension theorem has been formalized. Published first-type automaticity, Meyer or other field results remain separately recorded external mathematical inputs where used; this nonarithmetic scalar kernel neither needs nor proves them. No source change is requested.
