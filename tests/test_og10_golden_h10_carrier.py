import json
from pathlib import Path

from scripts.certify_og10_golden_h10_carrier import certify


def test_og10_golden_h10_carrier_certificate():
    result = certify()

    archived = json.loads(
        Path("data/candidates/og10_golden_h10_carrier.json").read_text()
    )
    assert archived == result

    assert result["verified"]
    assert result["field_lattice_rank_over_Z_phi"] == 11
    assert result["field_signatures"] == {
        "identity": [10, 1, 0],
        "conjugate": [11, 0, 0],
    }
    assert result["trace_lattice"]["rank"] == 22
    assert result["trace_lattice"]["signature"] == [1, 21]
    assert result["trace_lattice"]["phi_eigenspace_signature"] == [1, 10]
    assert result["trace_lattice"]["determinant"] == -44

    embedding = result["primitive_embedding"]
    assert embedding["complement_gram"] == [[2, 0], [0, 66]]
    assert embedding["gluing_graph_order"] == 44
    assert embedding["glued_rank6_determinant"] == 3
    assert embedding["glued_rank6_smith_invariants"] == [3]
    assert embedding["glued_rank6_isometry_type"] == "2U + A2(-1)"
    assert embedding["ambient_isometry_type"] == "3U + 2E8(-1) + A2(-1)"

    assert result["geometric_realization"]["conclusion"] == (
        "there_exists_a_projective_OG10_type_IHS_with_this_"
        "Neron_Severi_lattice_and_golden_H10_carrier"
    )
    assert result["geometric_realization"]["coxeter_wall_status"] == "not_yet_certified"

    wall_lift = result["root_wall_lift"]
    assert wall_lift["field_root_norm"] == 2
    assert wall_lift["divisor_squares"] == [-2, -4]
    assert wall_lift["divisor_pairing"] == -2
    assert wall_lift["restricted_hyperplanes_coincide_by_self_adjointness"]
