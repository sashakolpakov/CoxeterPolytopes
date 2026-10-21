import json

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from scripts.screen_h9_14_timeout_blockers import (
    FACET_ARCHIVE,
    RESOLVED_TIMEOUTS,
    SUPPORT_ARCHIVE,
    inclusion_exclusion_profile,
)


EXPECTED_UNDECIDED = (
    0,
    1,
    3,
    4,
    5,
    6,
    7,
    9,
    12,
    14,
    16,
    17,
    18,
    31,
    32,
    34,
    35,
    37,
    43,
    44,
    46,
    47,
    74,
)


def test_timeout_blocker_screen_has_exact_23_record_coverage() -> None:
    completion = json.loads(FACET_ARCHIVE.read_text(encoding="utf-8"))
    records = tuple(
        record
        for record in completion["records"]
        if record["solver"]["status"] == "unknown"
        and record["index"] not in RESOLVED_TIMEOUTS
    )
    assert tuple(record["index"] for record in records) == EXPECTED_UNDECIDED
    for record in records:
        covers = minimal_graph_vertex_covers(record["graph6"])
        assert 2 <= len(covers) <= 6
        profile = inclusion_exclusion_profile(covers)
        assert profile
        assert max(degree for _, degree in profile) <= 4


def test_high_edge_cluster_has_two_strict_link_schedules_and_small_terms() -> None:
    support = json.loads(SUPPORT_ARCHIVE.read_text(encoding="utf-8"))
    source = {
        record["graph6"]: record
        for record in support["records"]
        if record["survives_second_layer"]
    }
    completion = json.loads(FACET_ARCHIVE.read_text(encoding="utf-8"))
    by_index = {record["index"]: record for record in completion["records"]}
    for index in (43, 44, 46, 47, 74):
        record = by_index[index]
        schedules = tuple(
            schedule
            for schedule in source[record["graph6"]]["dehn_sommerville_schedule"]
            if summed_vertex_link_g4(schedule["f_vector"]) >= 14
        )
        assert tuple(schedule["missing_four_count"] for schedule in schedules) == (0, 1)
        assert tuple(schedule["missing_five_count_offset"] for schedule in schedules) == (4, 10)
        assert tuple(schedule["f_vector"][8] for schedule in schedules) == (90, 88)
        covers = minimal_graph_vertex_covers(record["graph6"])
        assert len(covers) == 4
        assert {cover.bit_count() for cover in covers} == {3}
        assert max(degree for _, degree in inclusion_exclusion_profile(covers)) <= 2
