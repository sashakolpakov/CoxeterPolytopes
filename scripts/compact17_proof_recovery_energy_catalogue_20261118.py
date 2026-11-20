"""Pure, unexecuted catalogue recipes for the existing compact-17 proof.

This module supplies theorem inputs, not new geometric proofs.  The caller
materializes iter_branches() under its own bounded scientific owner.  Nothing
is enumerated at import time, and there is no CLI, solver, file I/O or cache.
All costs refer to one original farthest vector.  A terminal contradiction
is reported as an exit, never disguised as an energy floor.
"""

from fractions import Fraction


SOURCE = "manuscripts/revisions/compact_farthest_vertex_bounds_consolidated.tex"
SOURCE_SHA256 = "6222f5c6cee8b30b2c1019d648bc31cde8bfcddf08f031a0355a3e7137c107ac"
DIMENSIONS = tuple(range(18, 30))

# Labels refer to the complete statements/proofs in SOURCE, not merely to a
# historical execution or to a final theorem excluding an entire dimension.
THEOREMS = {
    "original_geometry": {"anchors": ["sec:components", "sec:verification"], "lines": [362, 1360],
                          "scope": "One original farthest p,z,x; orthogonal incident-component energies sum to norm squared <1."},
    "classification": {"anchors": ["sec:verification"], "lines": [1360],
                       "scope": "Full spherical classification and actual defect compatibility, not realization of abstract forests."},
    "external_ceiling29": {"anchors": ["app:tail:prop:24to29"], "lines": [15124],
                          "scope": "The separately retained published Vinberg compact ceiling29, used only above29."},
    "nondominant_rate": {"anchors": ["sec:components", "sec:verification"], "lines": [362, 1360],
                         "scope": "Derived from all ordinary nondominant rank/14 costs and one maximum defect debit1/14: total rankS costs at least(S-1)/14. J owns two ranks; both typed debits37/700 are smaller."},
    "dominant_half": {"anchors": ["sec:components", "app:Achain", "app:prop:BDnine"],
                      "lines": [414, 2407, 2208],
                      "scope": "Each actual A_r,r>=7 or B_r/D_r,r>=9 has energy>1/2; arbitrary finite spectators."},
    "ordinary": {"anchors": ["sec:components", "app:ordinary"], "lines": [362, 1451],
                 "scope": "Unmarked ordinary irreducible; arbitrary finite spectators."},
    "defect": {"anchors": ["sec:components", "app:ordinary"], "lines": [380, 1451],
               "scope": "At most one actual I/J/N7/O7 profile; J owns two old A1 factors."},
    "I_joint": {"anchors": ["eq:joint", "app:Ijoint"], "lines": [397, 2034],
                "scope": "One actual I plus the stated ordinary A1/T2/T3 supports."},
    "typed_propagation": {"anchors": ["app:typedpropagation"], "lines": [2117],
                          "scope": "An actual N7/O7 forces the stated cost on each other ordinary support."},
    "A_chain": {"anchors": ["sec:components", "app:Achain"], "lines": [414, 2407],
                "scope": "Actual A component, rank at least7/8/13, and exact rank17, respectively; arbitrary spectators."},
    "A12": {"anchors": ["app:A12"], "lines": [4150], "scope": "Actual A12, arbitrary finite spectators."},
    "A15": {"anchors": ["app:A15"], "lines": [3939], "scope": "Actual A15, arbitrary finite spectators."},
    "A16": {"anchors": ["app:prop:A16"], "lines": [3466], "scope": "Actual A16, arbitrary finite spectators; energy >43/50."},
    "A18": {"anchors": ["app:Achain", "sec:components"], "lines": [2407, 430],
            "scope": "Actual A18, arbitrary finite spectators; energy >91/100."},
    "A19": {"anchors": ["app:tail:lem:A19ordinary"], "lines": [16801],
            "scope": "Actual A19, arbitrary finite spectators; energy >563/625."},
    "A20": {"anchors": ["app:tail:lem:A20A1"], "lines": [15995],
            "scope": "Actual A20, arbitrary finite spectators; energy >.92274."},
    "BD9": {"anchors": ["app:prop:BDnine"], "lines": [2208],
            "scope": "Actual B_r or D_r, r>=9; arbitrary finite spectators."},
    "prefix": {"anchors": ["prop:prefix", "app:exact18layers", "app:tail:prop:24to29"],
               "lines": [453, 4730, 15160],
               "scope": "6<=M<=16; B_r:r>=M+5, D_r:r>=M+6; all finite spectators. The global accounting paragraph is used only for its explicitly displayed rational prefix evaluations."},
    "prefix_gain": {"anchors": ["app:tail:lem:prefixgain"], "lines": [15774],
                    "scope": "Exactly the generic-prefix hypotheses. Replace its floor by floor+3969/800000; do not charge a separate support."},
    "BD16_ordinary_T2": {"anchors": ["app:exact18layers"], "lines": [4769, 4801, 4918],
                         "scope": "Entire B16 or D16 plus one ordinary T2 at actual dimension 18; fresh short-layer proof, not generic-prefix rank extension."},
    "B17_ordinary_A1": {"anchors": ["app:exact18layers"], "lines": [4923, 4968],
                        "scope": "Entire B17 plus one ordinary nonshort A1 at dimension 18."},
    "D19_ordinary_A1": {"anchors": ["app:tail:lem:D4proper", "app:tail:lem:prefixgain"],
                        "lines": [16026, 15774, 17930],
                        "scope": "Entire D19 plus ordinary nonshort A1; fresh D4 coverage through fifteen before the prefix-gain calculation, as explicitly consumed in the exact20 accounting."},
    "large_caps": {"anchors": ["app:tail:prop:caps"], "lines": [14306],
                   "scope": "Any A21..29, B23..29, or D24..29 component, arbitrary finite spectators."},
    "classical_I": {"anchors": ["app:tail:prop:classicalI"], "lines": [15189],
                    "scope": "Entire B_r+I,r>=20, or D_r+I,r>=21; NO further spectators."},
    "A16_typed": {"anchors": ["app:A16typed18"], "lines": [11583], "scope": "Entire A16+N7 or A16+O7, actual dimension18; both complete endpoint families."},
    "A17_rank1": {"anchors": ["app:A17rankone18"], "lines": [11941], "scope": "Entire A17 plus ordinary nonshort A1 or actual I, dimension18."},
    "BD16_typed": {"anchors": ["app:BD16typed18"], "lines": [12320], "scope": "Entire B16/D16 plus N7/O7, dimension18; both signed D forks and full B short end."},
    "B17_I": {"anchors": ["app:B17I18"], "lines": [12925], "scope": "Entire B17+I, actual dimension18."},
    "D17_A1": {"anchors": ["app:D17ordinary18"], "lines": [12735], "scope": "Entire D17 plus ordinary nonshort A1, dimension18."},
    "D17_I": {"anchors": ["app:D17I18"], "lines": [13814], "scope": "Entire D17+I, actual dimension18."},
    "whole_A18": {"anchors": ["app:prop:wholeA18", "app:wholeA18"], "lines": [5895, 4980], "scope": "Entire A18 at an original farthest vertex; all native/continuous endpoint roles."},
    "whole_B18": {"anchors": ["sec:b18", "app:wholeB18"], "lines": [638, 5918], "scope": "Entire B18; complete far/near/no-left and common/distinct terminal-cap families."},
    "whole_D18": {"anchors": ["sec:d18-setup", "app:D18foundationcomposition", "sec:composition"],
                  "lines": [662, 11533, 1255, 1340],
                  "scope": "Entire D18; complete initial/head/fork families, both signed forks, all own parameters. Stop at the local whole-D18 conclusion before the global theorem paragraph."},
    "A17_typed": {"anchors": ["app:tail:lem:A17typed"], "lines": [17951], "scope": "Entire A17 plus one actual typed rank-two remainder, dimension19."},
    "A18_I": {"anchors": ["app:tail:lem:A18I"], "lines": [18054], "scope": "A18 and actual I; arbitrary further finite spectators allowed. Generator uses only the entire two-component case."},
    "proper19": {"anchors": ["app:tail:lem:proper19standard"], "lines": [18691], "scope": "Entire D17+ordinary T2, B17+N7/O7, or B18/D18+ordinary A1, dimension19."},
    "D17_typed": {"anchors": ["app:tail:lem:D17typed"], "lines": [18265], "scope": "Entire D17+N7/O7, dimension19."},
    "B18_I": {"anchors": ["app:tail:lem:B18I"], "lines": [18188], "scope": "Entire B18+I, dimension19."},
    "D18_I": {"anchors": ["app:tail:lem:D18I"], "lines": [18452], "scope": "Entire D18+I, dimension19."},
    "whole_A19": {"anchors": ["app:exact19:prop:wholeA19"], "lines": [19733], "scope": "Entire A19, all endpoint patterns and native/continuous parameters."},
    "whole_B19": {"anchors": ["app:exact19:prop:wholeB19"], "lines": [21873], "scope": "Entire B19, full strong/helper and retained terminal families."},
    "whole_D19": {"anchors": ["app:exact19:prop:wholeD19"], "lines": [25061], "scope": "Entire D19, full native-head and terminal conjunction at dimension19."},
    "A19_I": {"anchors": ["app:tail:lem:A19I"], "lines": [16885], "scope": "A19 and actual I, arbitrary further finite spectators; generator uses entire two-component case."},
    "proper20": {"anchors": ["app:tail:lem:proper20classical"], "lines": [16687], "scope": "Entire B18/D18+N7/O7 or B19/D19+I, dimension20."},
    "whole_A20": {"anchors": ["app:tail:lem:wholeA20"], "lines": [17042], "scope": "Entire A20 at actual dimension20."},
    "whole_B20": {"anchors": ["app:tail:lem:wholeB20"], "lines": [17207], "scope": "Entire B20 at actual dimension20."},
    "whole_D20": {"anchors": ["app:tail:lem:wholeD20"], "lines": [17527], "scope": "Entire D20 at actual dimension20."},
    "A19_typed": {"anchors": ["app:tail:lem:A19typed"], "lines": [15970], "scope": "A19 plus N7/O7; generator uses entire dimension21 forest."},
    "A20_I": {"anchors": ["app:tail:lem:A20A1"], "lines": [15995], "scope": "Actual A20 cannot coexist with I; generator uses entire dimension21 forest."},
    "proper21": {"anchors": ["app:tail:lem:proper21"], "lines": [16082], "scope": "Entire D19+N7/O7 or D20+I, dimension21."},
    "whole_B21": {"anchors": ["app:tail:prop:wholeB21"], "lines": [16194], "scope": "Entire B21 at actual dimension21."},
    "whole_D21": {"anchors": ["app:tail:prop:wholeD21"], "lines": [16270], "scope": "Entire D21 at actual dimension21."},
    "whole_B22_D22_D23": {"anchors": ["app:tail:prop:whole58"], "lines": [15556], "scope": "Entire B22, D22, or D23 at its actual rank, not a whole-rank exclusion."},
}

ORDINARY = {
    ("A", 1): "1/10", ("A", 2): "319/2000", ("I", 2): "1/7",
    ("A", 3): "1/4", ("B", 3): "2/7", ("H", 3): "2/7",
    ("A", 4): "2/7", ("B", 4): "1/3", ("D", 4): "1/3",
    ("F", 4): "1/3", ("H", 4): "11/36", ("A", 5): "3/8",
    ("A", 6): "3/7", ("B", 5): "3/7", ("B", 6): "3/7",
    ("D", 5): "5/14", ("D", 6): "3/7", ("B", 7): "1/2",
    ("D", 7): "1/2", ("B", 8): "4/7", ("D", 8): "4/7",
    ("E", 6): "9/20", ("E", 7): "5/9", ("E", 8): "2/3",
}
PREFIX_DECIMALS = {
    7: ".7349975", 8: ".7683275", 9: ".797425", 10: ".823735",
    11: ".8472225", 12: ".8689525", 13: ".88864", 14: ".90715",
    15: ".9240775", 16: ".9401975",
}
WHOLE_EXITS = {
    ("A", 18): "whole_A18", ("B", 18): "whole_B18", ("D", 18): "whole_D18",
    ("A", 19): "whole_A19", ("B", 19): "whole_B19", ("D", 19): "whole_D19",
    ("A", 20): "whole_A20", ("B", 20): "whole_B20", ("D", 20): "whole_D20",
    ("B", 21): "whole_B21", ("D", 21): "whole_D21",
    ("B", 22): "whole_B22_D22_D23", ("D", 22): "whole_B22_D22_D23",
    ("D", 23): "whole_B22_D22_D23",
}


def _rational(value):
    value = Fraction(value)
    return f"{value.numerator}/{value.denominator}"


def _types():
    """All nondominant types that can fit spectator rank at most7."""
    out = [("A", r) for r in range(1, 7)]
    out += [("I", 2)]  # I2(m), all integers m>=4; A2 is m=3.
    out += [("B", r) for r in range(3, 8)]
    out += [("D", r) for r in range(4, 8)]
    out += [("E", 6), ("E", 7), ("F", 4), ("H", 3), ("H", 4)]
    return sorted(out, key=lambda item: (item[1], item[0]))


def _forests(types, remaining, start=0, prefix=()):
    """Canonical multiset recursion; no energy or dominant-count pruning."""
    if remaining == 0:
        yield prefix
        return
    for i in range(start, len(types)):
        item = types[i]
        if item[1] > remaining:
            break
        yield from _forests(types, remaining - item[1], i, prefix + (item,))


def _defects(forest):
    """Exchangeable equal-type slots are relabelled, not removed geometrically."""
    yield "none", ()
    a1 = [i for i, item in enumerate(forest) if item == ("A", 1)]
    if a1:
        yield "I", (a1[0],)
    if len(a1) >= 2:
        yield "J", tuple(a1[:2])
    for kind, old in (("N7", ("A", 2)), ("O7", ("I", 2))):
        for i, item in enumerate(forest):
            if item == old:
                yield kind, (i,)
                break


def _slots(forest, defect, owned):
    out = []
    for i, (family, rank) in enumerate(forest):
        role = defect if i in owned else "ordinary"
        parameter = {"m": 3} if (family, rank) == ("A", 2) else None
        if family == "I":
            parameter = {"m": 7} if role == "O7" else {"integer_m_min": 4, "integer_m_max": None}
        out.append({"key": f"c{i:03d}", "family": family, "rank": rank,
                    "role": role, "native_parameter": parameter})
    return out


def _exit(forest, defect, owned):
    for i, (family, rank) in enumerate(forest):
        if (family == "A" and rank >= 21) or (family == "B" and rank >= 23) or (family == "D" and rank >= 24):
            return "large_caps", {"component": i, "arbitrary_finite_spectators": True}
    if len(forest) == 1 and defect == "none":
        theorem = WHOLE_EXITS.get(forest[0])
        if theorem:
            return theorem, {"entire_forest": True, "all_native_and_continuous_roles": True}
    if len(forest) != 2:
        return None
    for i, (family, rank) in enumerate(forest):
        j = 1 - i
        other = forest[j]
        role = defect if j in owned else "ordinary"
        typed = role in ("N7", "O7")
        ordinary_a1 = role == "ordinary" and other == ("A", 1)
        ordinary_t2 = role == "ordinary" and other[1] == 2
        theorem = None
        if role == "I" and ((family == "B" and rank >= 20) or (family == "D" and rank >= 21)):
            theorem = "classical_I"
        elif family == "A":
            if rank == 16 and typed:
                theorem = "A16_typed"
            elif rank == 17 and (ordinary_a1 or role == "I"):
                theorem = "A17_rank1"
            elif rank == 17 and typed:
                theorem = "A17_typed"
            elif rank == 18 and role == "I":
                theorem = "A18_I"
            elif rank == 19 and role == "I":
                theorem = "A19_I"
            elif rank == 19 and typed:
                theorem = "A19_typed"
            elif rank == 20 and role == "I":
                theorem = "A20_I"
        elif family in ("B", "D"):
            if rank == 16 and typed:
                theorem = "BD16_typed"
            elif rank == 17 and role == "I":
                theorem = "B17_I" if family == "B" else "D17_I"
            elif family == "D" and rank == 17 and ordinary_a1:
                theorem = "D17_A1"
            elif rank == 17 and typed:
                theorem = "proper19" if family == "B" else "D17_typed"
            elif family == "D" and rank == 17 and ordinary_t2:
                theorem = "proper19"
            elif rank == 18 and ordinary_a1:
                theorem = "proper19"
            elif rank == 18 and role == "I":
                theorem = "B18_I" if family == "B" else "D18_I"
            elif (rank == 18 and typed) or (rank == 19 and role == "I"):
                theorem = "proper20"
            elif family == "D" and ((rank == 19 and typed) or (rank == 20 and role == "I")):
                theorem = "proper21"
        if theorem:
            return theorem, {"entire_forest": True, "component": i, "spectator": j,
                             "spectator_role": role, "all_native_and_continuous_roles": True}
    return None


def _rows(forest, defect, owned):
    rows = []

    def add(support, lower, theorem, **parameters):
        support = tuple(sorted(support))
        assert support and len(support) == len(set(support)) and theorem in THEOREMS
        rows.append({"key": f"r{len(rows):03d}", "support": [f"c{i:03d}" for i in support],
                     "lower": _rational(lower), "theorem": "compact." + theorem,
                     "parameters": parameters, "source_inequality": "strict"})

    if defect == "J":
        add(owned, "17/50", "defect", profile="J")
    elif defect in ("N7", "O7"):
        add(owned, "9/100", "defect", profile=defect)
    for i, (family, rank) in enumerate(forest):
        if i in owned:
            continue  # An actual defect never also gets an ordinary floor.
        if (family, rank) in ORDINARY:
            add((i,), ORDINARY[family, rank], "ordinary", family=family, rank=rank)
        if family == "A":
            for threshold, lower in ((7, "1/2"), (8, "4/7"), (13, "5/7")):
                if rank >= threshold:
                    add((i,), lower, "A_chain", rank=rank, threshold=threshold)
            if rank == 17:
                add((i,), "6/7", "A_chain", rank=17)
            special = {12: ("13/20", "A12"), 15: ("413/500", "A15"),
                       16: ("43/50", "A16"), 18: ("91/100", "A18"),
                       19: ("563/625", "A19"), 20: (".92274", "A20")}.get(rank)
            if special:
                add((i,), special[0], special[1], rank=rank)
        if family in ("B", "D") and rank >= 9:
            add((i,), "9/14", "BD9", family=family, rank=rank)
            for m, decimal in PREFIX_DECIMALS.items():
                if rank >= m + (5 if family == "B" else 6):
                    add((i,), decimal, "prefix", M=m, family=family, rank=rank)
                    add((i,), Fraction(decimal) + Fraction(3969, 800000), "prefix_gain",
                        M=m, family=family, rank=rank, replaces="prefix")
    ordinary = [i for i in range(len(forest)) if i not in owned]
    a1 = [i for i in ordinary if forest[i] == ("A", 1)]
    if defect == "I":
        for i in ordinary:
            rank = forest[i][1]
            lower = "1/5" if rank == 1 else "21/100" if rank == 2 else "2/7" if rank == 3 else None
            if lower:
                add(owned + (i,), lower, "I_joint", ordinary_rank=rank)
        for p, i in enumerate(a1):
            for j in a1[p + 1:]:
                add(owned + (i, j), "3/10", "I_joint", ordinary_A1_count=2)
    if defect in ("N7", "O7"):
        for i in ordinary:
            if forest[i][1] == 2:
                add((i,), "1/5", "typed_propagation", profile=defect, other_rank=2)
            elif forest[i] == ("A", 1):
                add((i,), "1/5" if defect == "N7" else "2/15", "typed_propagation",
                    profile=defect, other_rank=1)
    # The three exact-rank continuations below are genuinely distinct from
    # the generic prefix theorem.  No ordinary charge is assigned to I.
    if len(forest) == 2 and defect == "none":
        for i, (family, rank) in enumerate(forest):
            other = forest[1 - i]
            if family in ("B", "D") and rank == 16 and other[1] == 2:
                add((i,), ".8689525", "BD16_ordinary_T2", family=family, rank=rank)
            elif family == "B" and rank == 17 and other == ("A", 1):
                add((i,), ".90715", "B17_ordinary_A1")
            elif family == "D" and rank == 19 and other == ("A", 1):
                add((i,), ".91211125", "D19_ordinary_A1")
    return rows


def catalogue_header():
    return {
        "family": "existing_compact17_proof_energy_recovery_v1",
        "dimensions": list(DIMENSIONS),
        "sources": {SOURCE: SOURCE_SHA256},
        "premises": {"compact." + key: value for key, value in THEOREMS.items()},
        "coverage": {"ranks": "Every actual rank18..29 independently; no rank deletion or face inheritance.",
                     "hierarchy": "Number of dominant components0, >=2, or1. For exactly1 retain component-cap exits and any complete aggregate context whose uniform cost proves>=1. Otherwise expand every nondominant spectator partition. Only the source-proved S<=4 residual contexts reach a nontrivial recursion.",
                     "dominant": "A_r:r>=7; B_r,D_r:r>=9.",
                     "forests": "Every noncoarsened context keeps every nondominant spherical partition and every allowed defect profile. Coarse branches record their complete uniform aggregate rows or local component-cap exit, never silently omit a context.",
                     "native_bucket": "A2=I2(3); I2(m) for all integer m>=4 is a uniform theorem domain, not a finite label cutoff. O7 fixes its owned slot to m=7.",
                     "defects": "none, I, J, N7, O7; at most one, complete actual profiles; exchangeable equal-type slots may be relabelled.",
                     "branch_parameters": "Every actual own/white/native/parallel/dotted/foreign-contact parameter remains quantified in its complete local provider.",
                     "non_energy_exits": "Only explicit local component/entire-forest contradictions; no final dimension-exclusion theorem is a row."},
        "status": "Theorem-input catalogue; this object itself is not an execution or LP certificate.",
    }


def iter_branches():
    """Root-owned bounded execution only; every branch is independently usable."""
    types = _types()
    for dimension in DIMENSIONS:
        yield {"key": f"compact.n{dimension}:no_dominant", "dimension": dimension,
               "domain": {"dominant_count": 0, "complete_all_types_and_defects": True,
                          "component_partition_declared": False},
               "slots": [{"key": "all", "rank": dimension, "kind": "aggregate_orthogonal_energy"}],
               "rows": [{"key": "baseline", "support": ["all"],
                         "lower": _rational(Fraction(dimension - 1, 14)), "theorem": "compact.nondominant_rate"}],
               "exit": None, "hand_witness": {"baseline": "1/1"}}
        yield {"key": f"compact.n{dimension}:at_least_two_dominants", "dimension": dimension,
               "domain": {"dominant_count_min": 2, "component_partition_declared": False,
                          "rank_accounting": "D1 and D2 are any two distinct actual dominant components; rest is the orthogonal sum of every remaining component and may be empty."},
               "slots": [{"key": "D1", "rank": None, "kind": "actual_dominant_energy"},
                         {"key": "D2", "rank": None, "kind": "actual_dominant_energy"},
                         {"key": "rest", "rank": None, "kind": "aggregate_orthogonal_energy"}],
               "rows": [{"key": "half1", "support": ["D1"], "lower": "1/2", "theorem": "compact.dominant_half"},
                        {"key": "half2", "support": ["D2"], "lower": "1/2", "theorem": "compact.dominant_half"}],
               "exit": None, "hand_witness": {"half1": "1/1", "half2": "1/1"}}
        for family, minimum in (("A", 7), ("B", 9), ("D", 9)):
            for rank in range(minimum, dimension + 1):
                spectator_rank = dimension - rank
                cap = ((family == "A" and rank >= 21) or
                       (family == "B" and rank >= 23) or
                       (family == "D" and rank >= 24))
                if cap:
                    slots = [{"key": "D", "rank": rank, "kind": "actual_dominant_energy"}]
                    if spectator_rank:
                        slots.append({"key": "S", "rank": spectator_rank, "kind": "aggregate_orthogonal_energy"})
                    yield {"key": f"compact.n{dimension}:one_{family}{rank}:component_cap",
                           "dimension": dimension,
                           "domain": {"dominant_count": 1, "dominant_family": family,
                                      "spectator_rank": spectator_rank,
                                      "complete_all_nondominant_spectators_and_defects": True,
                                      "component_partition_declared": False},
                           "slots": slots, "rows": [],
                           "exit": {"theorem": "compact.large_caps",
                                    "hypotheses": {"component": "D", "family": family, "rank": rank,
                                                   "arbitrary_finite_spectators": True}}}
                    continue
                # Coarsening retains a complete branch and explicit valid LP
                # rows. It never deletes a context because of an unrecorded sum.
                lower = Fraction(1, 2)
                provider = {"theorem": "compact.dominant_half", "parameters": {}}
                if 0 < spectator_rank < 8:
                    provider = max(_rows(((family, rank),), "none", ()),
                                   key=lambda row: Fraction(row["lower"]))
                    lower = Fraction(provider["lower"])
                if spectator_rank > 0 and lower + Fraction(spectator_rank - 1, 14) >= 1:
                    yield {"key": f"compact.n{dimension}:one_{family}{rank}:S{spectator_rank}_coarse",
                           "dimension": dimension,
                           "domain": {"dominant_count": 1, "dominant_family": family,
                                      "spectator_rank": spectator_rank, "complete_all_nondominant_spectators_and_defects": True,
                                      "component_partition_declared": False},
                           "slots": [{"key": "D", "rank": rank, "kind": "actual_dominant_energy"},
                                     {"key": "S", "rank": spectator_rank, "kind": "aggregate_orthogonal_energy"}],
                           "rows": [{"key": "dominant", "support": ["D"], "lower": _rational(lower),
                                     "theorem": provider["theorem"], "parameters": provider["parameters"]},
                                    {"key": "spectator", "support": ["S"], "lower": _rational(Fraction(spectator_rank - 1, 14)), "theorem": "compact.nondominant_rate"}],
                           "exit": None, "hand_witness": {"dominant": "1/1", "spectator": "1/1"}}
                    continue
                for spectators in _forests(types, spectator_rank):
                    # Keep the unique dominant first; the remaining components
                    # retain the complete canonical spectator partition.
                    forest = ((family, rank),) + spectators
                    forest_key = "+".join("I2ge4" if f == "I" else f"{f}{r}" for f, r in forest)
                    for defect, owned in _defects(forest):
                        terminal = _exit(forest, defect, owned)
                        yield {"key": f"compact.n{dimension}:{forest_key}:{defect}", "dimension": dimension,
                               "domain": {"entire_forest": True, "dominant_count": 1,
                                          "spectator_rank": spectator_rank, "defect": defect,
                                          "owned_slots": [f"c{i:03d}" for i in owned],
                                          "one_original_farthest_vector": True,
                                          "component_partition_declared": True,
                                          "parameter_scope": "complete-provider-domains"},
                               "slots": _slots(forest, defect, owned),
                               "rows": _rows(forest, defect, owned),
                               "exit": None if terminal is None else {"theorem": "compact." + terminal[0], "hypotheses": terminal[1]}}


def catalogue():
    """Materialize only inside the coordinating agent's bounded owner."""
    result = catalogue_header()
    result["branches"] = list(iter_branches())
    return result
