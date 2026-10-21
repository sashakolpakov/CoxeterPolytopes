from pathlib import Path

from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, identity_matrix


def _write_tiny_atlas(path: Path) -> None:
    path.write_text(
        """GramMatS:=[[0, 1], [1, 0]];
prMatSdual:=[[1, 0], [0, 1]];
maxchamnumb:=0;
chamberdata[0][weyl]:=[1, 1];
chamberdata[0][Deltaw]:=[[1, 0], [0, 1]];
chamberdata[0][orderaut]:=1;
chamberdata[0][generatorsaut]:=[[[1, 0], [0, 1]]];
chamberdata[0][numborbitswalls]:=2;
chamberdata[0][wallorbit][0]:=[[1, 0]];
chamberdata[0][adjacent][0]:=[_minustwo];
chamberdata[0][wallorbit][1]:=[[0, 1]];
chamberdata[0][adjacent][1]:=[_isom, 0, _via, [[1, 0], [0, 1]]];
chamberdata[0][adjacentweyl][1]:=[1, 1];
chamberdata[0][adjacentDeltaw][1]:=[[1, 0], [0, 1]];
""",
        encoding="utf-8",
    )


def test_lazy_shimada_atlas_parser_and_dual_convention(tmp_path: Path) -> None:
    source = tmp_path / "disc11-tiny.txt"
    _write_tiny_atlas(source)
    atlas = ShimadaDisc11Atlas(source)
    chamber = atlas.chamber(0)

    assert atlas.chamber_count == 1
    assert chamber.wall_count == 2
    assert chamber.adjacencies[0].kind == "_minustwo"
    assert chamber.adjacencies[1].kind == "_isom"
    assert chamber.adjacencies[1].target == 0
    assert atlas.projected_delta_w(0) == ((1, 0), (0, 1))
    assert atlas.projected_delta_w(0, adjacent_orbit=1) == (
        (1, 0),
        (0, 1),
    )
    assert atlas.dual_action(identity_matrix(2)) == identity_matrix(2)
    assert atlas.stabilizer_group(0) == (identity_matrix(2),)
    assert atlas.wall_member_transports(0, 0) == (identity_matrix(2),)
    state = atlas.representative_state(0)
    assert atlas.neighboring_state(state, 0) is None
    assert atlas.neighboring_state(state, 1) == state
    assert atlas.state_wall_orbits(state) == chamber.wall_orbits


def test_shimada_adjacency_iteration_is_orbit_level(tmp_path: Path) -> None:
    source = tmp_path / "disc11-tiny.txt"
    _write_tiny_atlas(source)
    atlas = ShimadaDisc11Atlas(source)

    assert [
        (index, orbit, adjacency.kind, adjacency.target)
        for index, orbit, adjacency in atlas.iter_adjacencies()
    ] == [
        (0, 0, "_minustwo", None),
        (0, 1, "_isom", 0),
    ]
