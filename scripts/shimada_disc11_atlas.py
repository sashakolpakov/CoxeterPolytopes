"""Lazy reader for Shimada's discriminant-11 induced-chamber atlas.

The 47 MB ``disc11.txt`` supplement is deliberately not vendored.  This
module indexes the assignments needed for chamber traversal and parses a
record only when it is requested.  In particular it omits the much larger
global automorphism data unless a caller asks for it through another parser.

Shimada uses row vectors.  A row isometry ``M`` acts on ``S`` from the right,
and on dual-basis wall rows through ``G^-1 M G``.  See
``AlgoAutcompdata.pdf``, pp. 2--3, for the data format and conventions.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from functools import cached_property, lru_cache
from pathlib import Path
from typing import Iterator

from sympy import Matrix


_CHAMBER_NAME = re.compile(
    r"^chamberdata\[(?P<chamber>\d+)\]\[(?P<field>[^]]+)\]"
    r"(?:\[(?P<orbit>\d+)\])?$"
)
_MAPLE_SYMBOL = re.compile(r"\b(_minustwo|_backto|_newcham|_isom|_via)\b")


def _tuplify(value):
    if isinstance(value, list):
        return tuple(_tuplify(item) for item in value)
    return value


def _parse_maple_literal(text: str):
    """Parse the integer/list subset of Maple syntax used by the atlas."""

    value = text.replace("\\", "").strip()
    if value.endswith(";"):
        value = value[:-1]
    value = _MAPLE_SYMBOL.sub(lambda match: repr(match.group(1)), value)
    return _tuplify(ast.literal_eval(value))


def _matrix_product(left, right):
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("incompatible matrix dimensions")
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def _row_product(row, matrix):
    if len(row) != len(matrix):
        raise ValueError("incompatible row and matrix dimensions")
    return tuple(
        sum(row[i] * matrix[i][j] for i in range(len(row)))
        for j in range(len(matrix[0]))
    )


def identity_matrix(rank: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(i == j) for j in range(rank)) for i in range(rank))


@dataclass(frozen=True)
class ShimadaAdjacency:
    """The neighbor across the representative member of a wall orbit."""

    kind: str
    target: int | None = None
    isometry: tuple[tuple[int, ...], ...] | None = None


@dataclass(frozen=True)
class ShimadaChamber:
    index: int
    weyl: tuple[int, ...]
    delta_w: tuple[tuple[int, ...], ...]
    order_aut: int
    generators_aut: tuple[tuple[tuple[int, ...], ...], ...]
    wall_orbits: tuple[tuple[tuple[int, ...], ...], ...]
    adjacencies: tuple[ShimadaAdjacency, ...]

    @property
    def wall_count(self) -> int:
        return sum(len(orbit) for orbit in self.wall_orbits)


@dataclass(frozen=True)
class ShimadaChamberState:
    """A physical chamber ``D_index * row_isometry`` in the nef cone."""

    index: int
    row_isometry: tuple[tuple[int, ...], ...]


class ShimadaDisc11Atlas:
    """Seekable, lazy interface to an external ``disc11.txt`` file."""

    _GLOBAL_FIELDS = {"GramMatS", "prMatSdual", "maxchamnumb"}
    _CHAMBER_FIELDS = {
        "weyl",
        "Deltaw",
        "orderaut",
        "generatorsaut",
        "numborbitswalls",
        "wallorbit",
        "adjacent",
        "adjacentweyl",
        "adjacentDeltaw",
    }

    def __init__(self, path: Path | str):
        self.path = Path(path)
        self._offsets: dict[str, int] = {}
        self._build_index()

    def _wanted_name(self, name: str) -> bool:
        if name in self._GLOBAL_FIELDS:
            return True
        match = _CHAMBER_NAME.match(name)
        return bool(match and match.group("field") in self._CHAMBER_FIELDS)

    def _build_index(self) -> None:
        with self.path.open(encoding="utf-8") as handle:
            while True:
                offset = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if ":=" not in line:
                    continue
                name = line.split(":=", 1)[0].strip()
                if self._wanted_name(name):
                    self._offsets[name] = offset

        missing = self._GLOBAL_FIELDS - self._offsets.keys()
        if missing:
            raise ValueError(f"missing Shimada assignments: {sorted(missing)}")

    def _assignment_text(self, name: str) -> str:
        try:
            offset = self._offsets[name]
        except KeyError as error:
            raise KeyError(f"assignment {name!r} is absent from Shimada data") from error
        pieces: list[str] = []
        with self.path.open(encoding="utf-8") as handle:
            handle.seek(offset)
            first = handle.readline()
            found_name, value = first.split(":=", 1)
            if found_name.strip() != name:
                raise ValueError("the Shimada assignment index is inconsistent")
            pieces.append(value.strip())
            while not pieces[-1].rstrip().endswith(";"):
                continuation = handle.readline()
                if not continuation:
                    raise ValueError(f"unterminated Shimada assignment {name!r}")
                pieces.append(continuation.strip())
        return "".join(pieces)

    def assignment(self, name: str):
        return _parse_maple_literal(self._assignment_text(name))

    @cached_property
    def gram(self) -> tuple[tuple[int, ...], ...]:
        return self.assignment("GramMatS")

    @cached_property
    def projection_to_s_dual(self) -> tuple[tuple[int, ...], ...]:
        return self.assignment("prMatSdual")

    @cached_property
    def max_chamber_index(self) -> int:
        return self.assignment("maxchamnumb")

    @property
    def chamber_count(self) -> int:
        return self.max_chamber_index + 1

    def _chamber_name(self, index: int, field: str, orbit: int | None = None) -> str:
        if not 0 <= index <= self.max_chamber_index:
            raise IndexError(f"chamber index {index} is outside the atlas")
        suffix = "" if orbit is None else f"[{orbit}]"
        return f"chamberdata[{index}][{field}]{suffix}"

    @lru_cache(maxsize=None)
    def adjacency(self, index: int, orbit: int) -> ShimadaAdjacency:
        value = self.assignment(self._chamber_name(index, "adjacent", orbit))
        kind = value[0]
        if kind == "_minustwo":
            return ShimadaAdjacency(kind=kind)
        if kind in {"_backto", "_newcham"}:
            return ShimadaAdjacency(kind=kind, target=int(value[1]))
        if kind == "_isom":
            if len(value) != 4 or value[2] != "_via":
                raise ValueError("malformed _isom adjacency")
            return ShimadaAdjacency(
                kind=kind,
                target=int(value[1]),
                isometry=value[3],
            )
        raise ValueError(f"unknown Shimada adjacency tag {kind!r}")

    @lru_cache(maxsize=None)
    def chamber(self, index: int) -> ShimadaChamber:
        orbit_count = int(
            self.assignment(self._chamber_name(index, "numborbitswalls"))
        )
        wall_orbits = tuple(
            self.assignment(self._chamber_name(index, "wallorbit", orbit))
            for orbit in range(orbit_count)
        )
        return ShimadaChamber(
            index=index,
            weyl=self.assignment(self._chamber_name(index, "weyl")),
            delta_w=self.assignment(self._chamber_name(index, "Deltaw")),
            order_aut=int(self.assignment(self._chamber_name(index, "orderaut"))),
            generators_aut=self.assignment(
                self._chamber_name(index, "generatorsaut")
            ),
            wall_orbits=wall_orbits,
            adjacencies=tuple(
                self.adjacency(index, orbit) for orbit in range(orbit_count)
            ),
        )

    def projected_delta_w(
        self,
        index: int,
        *,
        adjacent_orbit: int | None = None,
    ) -> tuple[tuple[int, ...], ...]:
        """Return ``pr_S(Delta_w)`` in the dual basis of ``S``.

        When ``adjacent_orbit`` is given, the explicitly archived physical
        neighbor of an ``_isom`` adjacency is used.
        """

        field = "Deltaw" if adjacent_orbit is None else "adjacentDeltaw"
        delta = self.assignment(
            self._chamber_name(index, field, adjacent_orbit)
            if adjacent_orbit is not None
            else self._chamber_name(index, field)
        )
        return tuple(
            _row_product(vector, self.projection_to_s_dual) for vector in delta
        )

    @cached_property
    def gram_inverse(self) -> Matrix:
        return Matrix(self.gram).inv()

    @lru_cache(maxsize=None)
    def dual_action(
        self, row_isometry: tuple[tuple[int, ...], ...]
    ) -> tuple[tuple[int, ...], ...]:
        """Return the dual-basis row action ``G^-1 M G``."""

        action = self.gram_inverse * Matrix(row_isometry) * Matrix(self.gram)
        if any(not entry.is_Integer for entry in action):
            raise ValueError("a Shimada isometry induced a nonintegral dual action")
        return tuple(tuple(int(entry) for entry in row) for row in action.tolist())

    def transform_dual_rows(
        self,
        rows: tuple[tuple[int, ...], ...],
        row_isometry: tuple[tuple[int, ...], ...],
    ) -> tuple[tuple[int, ...], ...]:
        dual = self.dual_action(row_isometry)
        return tuple(_row_product(row, dual) for row in rows)

    @lru_cache(maxsize=None)
    def stabilizer_group(
        self, index: int
    ) -> tuple[tuple[tuple[int, ...], ...], ...]:
        """Enumerate the finite chamber stabilizer from its generators."""

        chamber = self.chamber(index)
        identity = identity_matrix(len(self.gram))
        seen = {identity}
        frontier = [identity]
        while frontier:
            current = frontier.pop()
            for generator in chamber.generators_aut:
                product = _matrix_product(current, generator)
                if product not in seen:
                    seen.add(product)
                    frontier.append(product)
                    if len(seen) > chamber.order_aut:
                        raise ValueError("stabilizer generators exceed the archived order")
        if len(seen) != chamber.order_aut:
            raise ValueError("stabilizer generators do not attain the archived order")
        return tuple(sorted(seen))

    @lru_cache(maxsize=None)
    def wall_member_transports(
        self, index: int, orbit: int
    ) -> tuple[tuple[tuple[int, ...], ...], ...]:
        """Find stabilizer elements taking the first wall to every member."""

        chamber = self.chamber(index)
        walls = chamber.wall_orbits[orbit]
        first = walls[0]
        actions = self.stabilizer_group(index)
        result = []
        for member_index, wall in enumerate(walls):
            if member_index == 0:
                result.append(identity_matrix(len(self.gram)))
                continue
            matches = [
                action
                for action in actions
                if _row_product(first, self.dual_action(action)) == wall
            ]
            if not matches:
                raise ValueError("wall orbit is not generated by the chamber stabilizer")
            result.append(matches[0])
        return tuple(result)

    def representative_state(self, index: int) -> ShimadaChamberState:
        return ShimadaChamberState(index, identity_matrix(len(self.gram)))

    @lru_cache(maxsize=None)
    def state_wall_orbits(
        self, state: ShimadaChamberState
    ) -> tuple[tuple[tuple[int, ...], ...], ...]:
        """Return primitive dual walls of a physical chamber state."""

        dual = self.dual_action(state.row_isometry)
        return tuple(
            tuple(_row_product(wall, dual) for wall in orbit)
            for orbit in self.chamber(state.index).wall_orbits
        )

    def neighboring_state(
        self,
        state: ShimadaChamberState,
        orbit: int,
        member: int = 0,
    ) -> ShimadaChamberState | None:
        """Cross a specified physical wall, or return ``None`` at a K3 wall.

        If ``a`` takes the representative member of the orbit to ``member``,
        the representative neighbor is ``D_j h a`` for an ``_isom`` record
        and ``D_j a`` for ``_newcham``/``_backto``.  Applying the current
        physical transport ``g`` gives the cumulative right action ``h a g``.
        """

        chamber = self.chamber(state.index)
        if not 0 <= orbit < len(chamber.wall_orbits):
            raise IndexError("wall-orbit index is outside the chamber")
        if not 0 <= member < len(chamber.wall_orbits[orbit]):
            raise IndexError("wall-member index is outside the orbit")
        adjacency = chamber.adjacencies[orbit]
        if adjacency.kind == "_minustwo":
            return None
        if adjacency.target is None:
            raise ValueError("a non-boundary adjacency has no target")
        member_transport = self.wall_member_transports(state.index, orbit)[member]
        local_transport = (
            adjacency.isometry
            if adjacency.kind == "_isom"
            else identity_matrix(len(self.gram))
        )
        cumulative = _matrix_product(
            _matrix_product(local_transport, member_transport),
            state.row_isometry,
        )
        return ShimadaChamberState(adjacency.target, cumulative)

    def iter_adjacencies(self) -> Iterator[tuple[int, int, ShimadaAdjacency]]:
        """Iterate representative wall-orbit adjacencies without loading chambers."""

        for index in range(self.chamber_count):
            count = int(
                self.assignment(self._chamber_name(index, "numborbitswalls"))
            )
            for orbit in range(count):
                yield index, orbit, self.adjacency(index, orbit)


__all__ = [
    "ShimadaAdjacency",
    "ShimadaChamber",
    "ShimadaChamberState",
    "ShimadaDisc11Atlas",
    "identity_matrix",
]
