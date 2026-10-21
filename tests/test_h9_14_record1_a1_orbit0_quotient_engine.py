from __future__ import annotations

from hashlib import sha256
import json

from scripts.assemble_h9_14_record1_a1_orbit0_quotient_logs import (
    independent_target_audit,
)
from scripts.independent_h9_14_record1_a1_orbit0_quotient_audit import (
    ARCHIVE,
    CPP,
    TARGET0,
    build_audit,
    independent_candidates,
    mask,
)


def test_record1_quotient_engine_has_an_unforgeable_terminal_policy() -> None:
    source = CPP.read_text(encoding="utf-8")
    assert "S2({1,2}) x S3({3,4,5}) x" in source
    assert "distinguished seed is not unique in family" in source
    assert "packed key failed dense adjacency replay" in source
    assert "whole_census_complete" in source
    assert "complete_selected_classes_only" in source
    assert "cutoff_max_nodes" in source
    assert "complete_exhaustive" in source


def test_independent_bounded_quotient_replay_matches_archive() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert build_audit() == archived
    assert archived["engine_sha256"] == sha256(CPP.read_bytes()).hexdigest()
    assert archived["source_group"]["explicit_permutation_count"] == 60_480
    assert archived["source_group"]["candidate_orbit_count"] == 21
    assert archived["seed_stabilizer_depth_two_replay"] == {
        "direct_quotient_count": 32,
        "equivalence_relations_identical": True,
        "literal_family_count": 1277,
        "production_quotient_count": 32,
        "seed_class_index": 0,
        "seed_mask": 31,
        "seed_orbit_size": 3,
        "seed_stabilizer_order": 20160,
    }
    assert archived["cutoff_protocol"]["whole_census_complete"] is False


def test_known_target_passes_independent_assembly_gate() -> None:
    candidates = frozenset(independent_candidates())
    additions = sorted(mask(row) for row in TARGET0)
    # The assembler checks the key encoding separately from its nauty value;
    # this syntactically valid placeholder lets this test isolate the literal
    # local/blocker/ridge replay.
    target = {
        "source_A4_mask": 57,
        "additional_A5_masks": additions,
        "scheduled_facets": 158,
        "active_ridges": 711,
        "canonical_full_key": "17" + "00" * 41,
    }
    seed = next(value for value in additions if value == 410)
    assert independent_target_audit(target, candidates, seed) == {
        "active_ridge_count": 711,
        "facet_count": 158,
    }
