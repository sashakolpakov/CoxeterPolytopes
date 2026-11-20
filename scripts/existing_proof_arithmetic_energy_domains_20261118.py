"""Pure recipe data for arithmetic15 recovery from EXISTING proved inputs.

No file access, solver, worker import, enumeration of geometric objects or
import-time catalogue construction.  Root may call catalogue() under its
bounded owner.  Every domain inclusion/local theorem/cover is an explicit
mathematical input, not something this data supplier verifies.

This is the pruned-cover initial recovery: 14 rank17 and 19 rank16 source
groups, with every displayed spectator alternative retained. Symbolic T2,
T3,T4 records denote the WHOLE stated uniform theorem family. They are not
sampled native labels. No actual forest is asserted realizable.
"""

SOURCE = {
    "path": "manuscripts/revisions/compact_arithmetic_dimension_bound_consolidated.tex",
    "sha256": "d0f040b253146f318bafd2391c8809c9d995f53fcd920477f78f42038d424295",
}

# Exact source labels, not certificates constructed by this module.
THEOREMS = {
    "GEO": ("original-vector geometry; disjoint original components",),
    "COVER16": ("prop:f16", "tab:f16", "app:arithranksixteenclosure"),
    "COVER17": ("prop:f17", "tab:f17", "app:arithrankseventeenclosure"),
    "ORD": ("subsec:shared-seals", "app:arithcoverinputs"),
    "DEF": ("actual I/J/N7/O7 definitions and compatibility",),
    "PROP": ("actual typed-defect propagation in prop:f16 / prop:f17",),
    "IJOINT": ("lem:Ijoint",),
    "A10d": ("app:arithA10consumers",),
    "A11": ("app:arithA11uniform",),
    "A11I": ("app:arithA11Ijoint",),
    "A12": ("app:arithA12universal",),
    "A13": ("app:arithA13universal",),
    "A13O": ("app:arithA13seven",),
    "A14": ("app:arithA14fourfifths",),
    "A14ORD": ("app:arithA14ordinarycore", "app:arithA14ordinarycomposition"),
    "A14DEF": ("app:arithA14typed", "app:arithA14typedcomposition"),
    "A15": ("app:arithsummedfeet", "app:arithgoldenend"),
    "A15DEF": ("app:arithA15typed",),
    "A3O": ("lem:A3point26",),
    "DIH": ("prop:ordinarydihedralthreenineteenths: every ordinary I2(m), m>=4",
            "subsec:shared-seals: separate A2=I2(3) floor 319/2000 > 3/19"),
    "BD14DEF": ("app:arithBD14typed", "app:arithBD14typedhead"),
    "BD15DEF": ("app:arithBD15typed",),
    "D13": ("app:arithD13rankthree",),
    "D14ORD": ("cor:D14ordinarydihedral",),
    "B17SURV": ("app:arithB17physical", "app:arithB17compatibility", "app:arithB17domain", "DATA3"),
    "D17SURV": ("app:arithD17geometry", "app:arithD17domain", "DATA3"),
    "B16SURV": ("app:arithB16edges", "app:arithB16domain", "DATA4"),
    "X-A15T2": ("app:arithA15ordinary",),
    "X-A16R1": ("app:arithA16rankone",),
    "X-A17": ("app:wholearithA17", "app:arithA17datainclusion", "DATA1", "DATA2"),
    "X-B16I": ("app:arithB16Ishort", "app:arithB16Ihead"),
    "X-D16R1": ("app:arithD16ordinary", "app:arithD16ordinaryfinish"),
    "X-D16I": ("app:arithD16Ifoundation", "app:arithD16Iquadraticfinish"),
    "X-A13I": ("app:arithA13I",),
    "X-A15R1": ("app:arithrankone", "app:arithA15nonisolated", "complete isolated A15 rank-one subcases"),
    "X-A16": ("app:wholearithA16", "complete actual finite-anchor/time-ratio proof"),
    "X-BD15R1": ("app:arithBD15ordinary", "separate B15/D15 nonisolated reductions"),
    "X-B15I": ("app:arithB15Iwindow", "app:arithB15Ipurepair"),
    "X-D15I": ("app:arithD15Icentre", "app:arithD15Iforknonisolated"),
    "X-D16DYADIC": ("sec:d16", "prop:d16reduce", "prop:noD16"),
    "DATA1": ("prop:pair-output(i): complete quadratic A17 empty domain",),
    "DATA2": ("prop:pair-output(ii): both complete cubic A17 finite-anchor empty domains",),
    "DATA3": ("prop:pair-output(iii): complete shared quadratic BD17 output, TWO survivors",),
    "DATA4": ("prop:pair-output(iv): complete whole-B16 output, ONE survivor",),
    "ARITH": ("actual common totally real field, original root lines, integral cycles, field squares, all nonphysical embeddings",),
    "LATTICE": ("literal invariant lattice and full reflection-group inclusion wherever the local source uses the golden-lattice theorem",),
}


def catalogue():
    """Return the complete finite recipe catalogue; no theorem is evaluated.

    Lower bounds are exact rational strings. Each slot is one ORIGINAL
    component. Rows on several slots are genuinely joint costs, not extra
    energy variables. Complete local X exits force ONLY that matched mass
    to zero. The caller must retain/validate all named mathematical inputs.
    """
    branches = []

    def add(key, dimension, group, components, rows=(), exit_key=None):
        slots = [{"key": "s%d" % j, "component": kind, "rank": rank,
                  "role": role} for j, (kind, rank, role) in enumerate(components)]
        domain = {"cover_group": "%d.%d" % (dimension, group),
                  "complete_original_forest": [dict(s) for s in slots],
                  "retained_hypotheses": ["GEO", "ARITH", "LATTICE",
                                           "COVER%d" % dimension],
                  "all_local_theorem_hypotheses_retained": True}
        branches.append({
            "key": key, "dimension": dimension, "domain": domain,
            "slots": slots,
            "rows": [{"key": key + ".cost%d" % j,
                      "support": ["s%d" % i for i in support],
                      "lower": lower, "theorem": theorem}
                     for j, (support, lower, theorem) in enumerate(rows)],
            "exit": None if exit_key is None else {
                "theorem": exit_key,
                "hypotheses": {"complete_domain": domain,
                               "no_exclusion_of_other_branches": True}},
        })

    def c(kind, rank, role="ordinary"):
        return (kind, rank, role)

    I = c("A1", 1, "actual_short_I: 0<h<1/2")
    O1 = c("A1", 1)
    T2 = c("I2(m)", 2, "ordinary, every integer m>=3")
    T3 = c("A3/B3/H3", 3, "ordinary complete family")
    T4 = c("A4/B4/D4/F4/H4", 4, "ordinary complete family")
    defects = (("N7", c("A2", 2, "actual_N7")),
               ("O7", c("I2(7)", 2, "actual_O7")))
    ordinary_five = (
        ("A5", (c("A5", 5),)), ("D5", (c("D5", 5),)),
        ("A4_A1", (c("A4", 4), O1)),
        ("H4_A1", (c("H4", 4), O1)),
        ("T3_T2", (T3, T2)), ("T2_T2_A1", (T2, T2, O1)),
    )

    # All fourteen rank-seventeen groups.
    for kind in ("A5", "D5"):
        add("a17.g01." + kind, 17, 1, (c("A11", 11), I, c(kind, 5)),
            (((0,), "13/20", "A11"), ((2,), "5/14", "ORD")))
    for tag, defect in defects:
        for kind in ("A4", "B4", "D4", "F4"):
            add("a17.g02." + tag + "." + kind, 17, 2,
                (c("A11", 11), defect, c(kind, 4)),
                (((0,), "13/20", "A11"), ((1,), "9/100", "DEF"),
                 ((2,), "2/7", "ORD")))
    for tag, tail in ordinary_five:
        add("a17.g03." + tag, 17, 3, (c("A12", 12),) + tail,
            (((0,), "5/7", "A12"),
             (tuple(range(1, len(tail) + 1)), "5/14", "ORD")))
    add("a17.g04.I_T4", 17, 4, (c("A12", 12), I, T4),
        (((0,), "5/7", "A12"), ((1, 2), "23/70", "IJOINT")))
    add("a17.g04.I_T3_A1", 17, 4, (c("A12", 12), I, T3, O1),
        (((0,), "5/7", "A12"), ((1, 2), "2/7", "IJOINT"), ((3,), "1/10", "ORD")))
    add("a17.g04.I_T2_T2", 17, 4, (c("A12", 12), I, T2, T2),
        (((0,), "5/7", "A12"), ((1, 2), "21/100", "IJOINT"), ((3,), "1/7", "ORD")))
    add("a17.g05.T3", 17, 5, (c("A14", 14), T3),
        (((0,), "4/5", "A14"), ((1,), "1/4", "ORD")))
    add("a17.g05.T2_A1", 17, 5, (c("A14", 14), T2, O1),
        (((0,), "4/5", "A14"), ((1,), "1/7", "ORD"), ((2,), "1/10", "ORD")))
    add("a17.g05.I_T2", 17, 5, (c("A14", 14), I, T2),
        (((0,), "4/5", "A14"), ((1, 2), "21/100", "IJOINT")))
    add("a17.g06.T2", 17, 6, (c("A15", 15), T2), exit_key="X-A15T2")
    add("a17.g06.A1_A1", 17, 6, (c("A15", 15), O1, O1),
        (((0,), "2087/2500", "A15"), ((1,), "1/10", "ORD"), ((2,), "1/10", "ORD")))
    add("a17.g06.I_A1", 17, 6, (c("A15", 15), I, O1),
        (((0,), "2087/2500", "A15"), ((1, 2), "5051/25000", "IJOINT")))
    for tag, defect in defects:
        add("a17.g06." + tag, 17, 6, (c("A15", 15), defect),
            (((0,), "574463/625000", "A15DEF"), ((1,), "9/100", "DEF")))
    for tag, rankone in (("ordinary", O1), ("I", I)):
        add("a17.g07." + tag, 17, 7, (c("A16", 16), rankone), exit_key="X-A16R1")
    add("a17.g08.whole", 17, 8, (c("A17", 17),), exit_key="X-A17")
    for tag, defect in defects:
        add("a17.g09." + tag, 17, 9, (c("B15", 15), defect),
            (((0,), "93859/100000", "BD15DEF"), ((1,), "9/100", "DEF")))
    add("a17.g10.I", 17, 10, (c("B16", 16), I), exit_key="X-B16I")
    add("a17.g11.whole", 17, 11, (c("B17", 17),),
        (((0,), "12523/12500", "B17SURV"),))
    for tag, defect in defects:
        add("a17.g12." + tag, 17, 12, (c("D15", 15), defect),
            (((0,), "91859/100000", "BD15DEF"), ((1,), "9/100", "DEF")))
    add("a17.g13.ordinary", 17, 13, (c("D16", 16), O1), exit_key="X-D16R1")
    add("a17.g13.I", 17, 13, (c("D16", 16), I), exit_key="X-D16I")
    add("a17.g14.whole", 17, 14, (c("D17", 17),),
        (((0,), "12523/12500", "D17SURV"),))

    # All nineteen rank-sixteen groups, independently of the rank17 cover.
    for kind in ("A5", "D5"):
        add("a16.g01." + kind, 16, 1, (c("A10", 10), I, c(kind, 5)),
            (((0,), "404/625", "A10d"), ((2,), "5/14", "ORD")))
    for tag, defect in defects:
        for kind in ("A4", "B4", "D4", "F4"):
            add("a16.g02." + tag + "." + kind, 16, 2,
                (c("A10", 10), defect, c(kind, 4)),
                (((0,), "404/625", "A10d"), ((1,), "9/100", "DEF"), ((2,), "2/7", "ORD")))
    for tag, tail in ordinary_five:
        actual_tail = (c("A3", 3), T2) if tag == "T3_T2" else tail
        add("a16.g03." + tag, 16, 3, (c("A11", 11),) + actual_tail,
            (((0,), "13/20", "A11"),
             (tuple(range(1, len(actual_tail) + 1)), "5/14", "ORD")))
    add("a16.g04.I_T4", 16, 4, (c("A11", 11), I, T4),
        (((0, 1), "18/25", "A11I"), ((2,), "2/7", "ORD")))
    add("a16.g05.O7_A3", 16, 5, (c("A11", 11), defects[1][1], c("A3", 3)),
        (((0,), "13/20", "A11"), ((1,), "9/100", "DEF"), ((2,), "13/50", "A3O")))
    for tag, tail in (("A4", (c("A4", 4),)), ("H4", (c("H4", 4),)),
                      ("T2_T2", (T2, T2))):
        add("a16.g06." + tag, 16, 6, (c("A12", 12),) + tail,
            (((0,), "5/7", "A12"), (tuple(range(1, len(tail) + 1)), "2/7", "ORD")))
    add("a16.g07.I_T3", 16, 7, (c("A12", 12), I, T3),
        (((0,), "5/7", "A12"), ((1, 2), "2/7", "IJOINT")))
    for tag, defect in defects:
        add("a16.g08." + tag, 16, 8, (c("A12", 12), defect, T2),
            (((0,), "5/7", "A12"), ((1,), "9/100", "DEF"), ((2,), "1/5", "PROP")))
    add("a16.g09.A3", 16, 9, (c("A13", 13), c("A3", 3)),
        (((0,), "53/70", "A13"), ((1,), "1/4", "ORD")))
    add("a16.g09.T2_A1", 16, 9, (c("A13", 13), T2, O1),
        (((0,), "53/70", "A13"), ((1,), "1/7", "ORD"), ((2,), "1/10", "ORD")))
    add("a16.g09.I_T2", 16, 9, (c("A13", 13), I, T2), exit_key="X-A13I")
    add("a16.g09.O7_A1", 16, 9, (c("A13", 13), defects[1][1], O1),
        (((0,), "233/300", "A13O"), ((1,), "9/100", "DEF"), ((2,), "2/15", "PROP")))
    add("a16.g10.T2", 16, 10, (c("A14", 14), T2),
        (((0,), "16/19", "A14ORD"), ((1,), "3/19", "DIH")))
    add("a16.g10.A1_A1", 16, 10, (c("A14", 14), O1, O1),
        (((0,), "4/5", "A14"), ((1,), "1/10", "ORD"), ((2,), "1/10", "ORD")))
    add("a16.g10.I_A1", 16, 10, (c("A14", 14), I, O1),
        (((0,), "4/5", "A14"), ((1, 2), "5051/25000", "IJOINT")))
    for tag, defect in defects:
        add("a16.g10." + tag, 16, 10, (c("A14", 14), defect),
            (((0,), "45643222599/50000000000", "A14DEF"), ((1,), "9/100", "DEF")))
    for tag, rankone in (("ordinary", O1), ("I", I)):
        add("a16.g11." + tag, 16, 11, (c("A15", 15), rankone), exit_key="X-A15R1")
    add("a16.g12.whole", 16, 12, (c("A16", 16),), exit_key="X-A16")
    for tag, defect in defects:
        add("a16.g13." + tag, 16, 13, (c("B14", 14), defect),
            (((0,), "91859/100000", "BD14DEF"), ((1,), "9/100", "DEF")))
    add("a16.g14.ordinary", 16, 14, (c("B15", 15), O1), exit_key="X-BD15R1")
    add("a16.g14.I", 16, 14, (c("B15", 15), I), exit_key="X-B15I")
    add("a16.g15.whole", 16, 15, (c("B16", 16),),
        (((0,), "4025629707/4000000000", "B16SURV"),))
    add("a16.g16.I_T2", 16, 16, (c("D13", 13), I, T2),
        (((0,), "20609/25000", "D13"), ((1, 2), "21/100", "IJOINT")))
    add("a16.g16.O7_A1", 16, 16, (c("D13", 13), defects[1][1], O1),
        (((0,), "20609/25000", "D13"), ((1, 2), "21/100", "D13")))
    add("a16.g17.ordinary", 16, 17,
        (c("D14", 14), c("I2(m)", 2, "ordinary, every integer m>=4")),
        (((0,), "5301/6250", "D14ORD"), ((1,), "3/19", "DIH")))
    for tag, defect in defects:
        add("a16.g17." + tag, 16, 17, (c("D14", 14), defect),
            (((0,), "91859/100000", "BD14DEF"), ((1,), "9/100", "DEF")))
    add("a16.g18.ordinary", 16, 18, (c("D15", 15), O1), exit_key="X-BD15R1")
    add("a16.g18.I", 16, 18, (c("D15", 15), I), exit_key="X-D15I")
    add("a16.g19.whole", 16, 19, (c("D16", 16),), exit_key="X-D16DYADIC")

    return {
        "family": "existing-proof-compact-arithmetic-energy",
        "dimensions": [16, 17],
        "sources": {SOURCE["path"]: SOURCE["sha256"]},
        "premises": {key: {"source": SOURCE["path"], "anchors": list(anchors),
                           "scope": "complete proved source input with every original hypothesis retained"}
                     for key, anchors in THEOREMS.items()},
        "coverage": {"kind": "proved_pruned_complete_terminal_covers",
                     "groups": {"16": list(range(1, 20)), "17": list(range(1, 15))},
                     "expected_recipe_counts_from_source": {"16": 49, "17": 39},
                     "source": dict(SOURCE), "theorems": ["COVER16", "COVER17"],
                     "all_pruning_baseline_proofs_retained": True,
                     "no_realization_of_formal_forests_claimed": True,
                     "append_unchanged_compact_18_through_29_catalogue": True,
                     "above29_requires_explicit_published_compact_ceiling": True,
                     "threshold_only_filters_fixed_actual_dimension_labels": True},
        "branches": branches,
        "hand_witness": {
            b["key"]: [{"row": r["key"], "weight": "1"} for r in b["rows"]]
            for b in branches},
    }
