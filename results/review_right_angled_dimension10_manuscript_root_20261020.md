# Root review of the consolidated dimension-ten manuscript

20 October 2026. **PASS: complete mathematical translation, final
typesetting and all three figures.** This is one successor manuscript
in the mixed right-angled family. It is not an extra proposed paper
alongside the preserved ceiling-eleven predecessor.

Final source `manuscripts/right_angled_dimension10.tex`, SHA256
`2fd439454de989f34590681f11e0475c137effb9faa71a72011107bf5b0922c1`.
Final PDF `manuscripts/build/right_angled_dimension10.pdf`, SHA256
`9b1362cb77405084b83fc6798bb9a4b951d453e32c1eda1d1f632018195dc205`:
12 pages, three vector figures, no final-pass warnings.

I read every line of the893-line first mathematical draft, independently
checked the scalar/count derivations against their complete sources,
and reviewed every subsequent diff. The final source includes only
the documented wording/strictness clarifications and figure spacing
changes beyond that full read. Initial draft and first-layout source/PDF
are retained in byte-verified gzip archives under `source_archives`.

## Scope and complete arguments

The theorem excludes all finite-volume right-angled polytopes in
dimensions at least eleven. The new mathematical contradiction concerns
mixed dimension eleven. Its higher tail uses eleven-face heredity and
the published compact4/entirely-ideal6 bounds; it imports neither the
earlier internal ceiling eleven nor an unproved higher-dimensional case.
The manuscript correctly limits inherited right-angle statements to
faces of dimension at least three, with polygon conventions explicit.

The actual opposition proof checks every other ambient wall, proves
feasibility of all five bipyramid vertices, and closes the whole
three-face by convexity. It uses one invariant ambient h and one fixed
future null vector on each actual graph component. The frame equations
and reciprocal reversal are derived, with no sign assumption on p_i
or substitution of a maximum-height average.

The bad-triangle attachment, small-type counts, R7 degree pattern,
one-cusp ten-facet exclusion and direct H<=2E argument are all in the
paper. In particular the last argument avoids metric D3 recognition.
The signed defects are not assumed positive. All ordinary count
multiplicities refer to the same actual P3 count B as transport.
The strict DS-weighted row, balanced corner row, defect/cusp elimination
and lower3722B>159195f+2794E are derived explicitly. The external
eight-dimensional input is accurately stated as a lower bound82,
not an asserted attained minimum.

Every scalar subset size is covered. Sizes0--6 have an elementary
proof;7/8 have explicit rational interval upper envelopes and covering
grids;9 has the complete global analytic proof;10 has full joint moments,
corrected rational envelopes, cleared polynomials and Bernstein
conversion;11 has a direct rational moment comparison. The finite
positivity algorithms and table sizes are stated in the paper. The
certificate records sufficient exact polynomial data, gaps, digests and
interval bounds for full reconstruction; it is not a sampled search.

The k9 clipping comparison is weak at d=0, and its reciprocal coefficient
is at most r^-2+4r^-1 before the strict constant bound. Those harmless
strictness errors in the pinned research-source display are corrected
in this final TeX. The theorem and the exact identities are unchanged.
The cleared difference is displayed with its sign unambiguous.

Finally the manuscript enlarges only nonnegative formal terms and sums
only actual directed oppositions. Reverse-edge pairing gives an
ordinary upper4B<=168f+3E. Its exact comparison with the strict lower
would force5742f+5E<0. The separate conditional Lean module checks
this logic from its explicitly listed geometry/analytic interfaces;
the paper does not relabel them as formal geometric theorems.

## Build, figures and resource use

The final bounded build took1.578 seconds, with85.4MB peak process-tree
RSS, one shared lock, and automatic removal of its temporary directory.
The12-page PDF has zero warnings. The first layout's crowded captions
in figures2/3 were corrected by changing only node widths, fonts and
positions. I directly inspected the final figure pages5 and11;
figure1/page2 is byte-identical to its previously inspected rendering.
All labels, ideal/finite vertex distinctions, actual-versus-formal
domains and reciprocal arrows are legible without overlap.

The final visual ledger is
`results/mixed_dimension10_checked_figure_visual_review_20261020.json`.
After the final wording clarifications all three pages were rendered
again and proved byte-identical to their directly inspected versions.
All nine disposable PNGs and three owned preview directories were removed,
reclaiming2247379 bytes. The final source, PDF, exact certificates,
byte-verified archives and compact ledgers remain. Disk was about181GiB,
well above the user's cleanup trigger; no cleanup daemon was started.

Written mathematical review, exact replay, conditional kernel checking,
typesetting and visual inspection retain separate scopes. This PASS is
internal review, not external refereeing or end-to-end Lean geometry.
