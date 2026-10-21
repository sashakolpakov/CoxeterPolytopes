"""Standard exact lattices used by calibration cases."""

from __future__ import annotations


def e8_cartan_gram() -> list[list[int]]:
    """The positive-definite E8 root lattice Gram matrix."""

    rank = 8
    gram = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    # Arms of lengths 2, 4, and 1 at vertex 2.
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]:
        gram[i][j] = gram[j][i] = -1
    return gram


def block_diagonal(*blocks: list[list[int]]) -> list[list[int]]:
    size = sum(len(block) for block in blocks)
    result = [[0 for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        if any(len(row) != len(block) for row in block):
            raise ValueError("blocks must be square")
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[offset + i][offset + j] = entry
        offset += len(block)
    return result


def hyperbolic_plane() -> list[list[int]]:
    return [[0, 1], [1, 0]]


def k3_lattice_gram() -> list[list[int]]:
    """The K3 lattice U^3 + E8(-1)^2, with signature (3, 19)."""

    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    return block_diagonal(
        hyperbolic_plane(),
        hyperbolic_plane(),
        hyperbolic_plane(),
        e8_negative,
        e8_negative,
    )
