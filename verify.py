"""Replay the finite cases and check residual placements by a second method."""

import argparse
from collections import deque
from functools import cache
from itertools import combinations_with_replacement
from pathlib import Path
from subprocess import PIPE, Popen


# Flocco, Pulaj, and Yerger (arXiv:2312.12618), Figure 1, in their numbering.
SECOND_FIGURE = ("12", "13", "24", "35", "36", "37", "45", "46", "47", "48", "58", "68", "78")
FIRST_TO_SECOND = (2, 1, 3, 5, 8, 6, 4, 7)
ROOTS = (9, 12, 36)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def square_data():
    inverse = {name: i for i, name in enumerate(FIRST_TO_SECOND)}
    factor = [set() for _ in range(8)]
    for edge in SECOND_FIGURE:
        a, b = (inverse[int(name)] for name in edge)
        factor[a].add(b)
        factor[b].add(a)
    neighbors = tuple(
        tuple(sorted([8 * candidate + b for candidate in factor[a]] +
                     [8 * a + candidate for candidate in factor[b]]))
        for a in range(8) for b in range(8)
    )
    distances = []
    for root in range(64):
        distance = {root: 0}
        pending = deque([root])
        while pending:
            vertex = pending.popleft()
            for child in neighbors[vertex]:
                if child not in distance:
                    distance[child] = 1 + distance[vertex]
                    pending.append(child)
        require(len(distance) == 64, "Disconnected square")
        distances.append(tuple(distance[vertex] for vertex in range(64)))
    return factor, neighbors, distances


def check_geometry(factor, distances):
    factor_distance = [[distances[a * 8][b * 8] for b in range(8)] for a in range(8)]
    require([(a, b) for a in range(8) for b in range(a) if factor_distance[a][b] == 3] == [(4, 1)],
            "Unexpected distance-three pairs")
    for root in range(8):
        signatures = [factor[root] & factor[v] for v in range(8) if factor_distance[root][v] == 2]
        require(any(all({a, b} & choices for choices in signatures)
                    for a in factor[root] for b in factor[root]), "No two-neighbor cover")
    for root, opposite in ((1, 4), (4, 1)):
        require(all(factor_distance[opposite][v] == 2 for v in factor[root]), "Nongeodesic neighbor")
        require(factor[opposite] == {v for v in range(8) if factor_distance[root][v] == 2},
                "Antipode neighbors differ from the distance-two layer")
    require({root for root, row in enumerate(distances) if max(row) == 6} == {9, 12, 33, 36},
            "Unexpected eccentricity-six roots")


def can_descend(root, vertices, capacity, neighbors, distances):
    """Work from the root recursively with ordinary integer tuples."""
    unit = {vertex: tuple(int(i == j) for j in range(len(vertices)))
            for i, vertex in enumerate(vertices)}
    distance = distances[root]

    def sums(choices):
        for a, b in combinations_with_replacement(choices, 2):
            total = tuple(x + y for x, y in zip(a, b))
            if all(value <= limit for value, limit in zip(total, capacity)):
                yield total

    @cache
    def ancestry(vertex):
        # A set with entries sorted by mass suffices to remove dominated vectors.
        candidates = {unit[vertex]} if vertex in unit else set()
        for parent in neighbors[vertex]:
            if distance[parent] == distance[vertex] + 1:
                candidates.update(sums(ancestry(parent)))
        minimal = []
        for candidate in sorted(candidates, key=sum):
            if not any(all(x <= y for x, y in zip(previous, candidate)) for previous in minimal):
                minimal.append(candidate)
        return tuple(minimal)

    return any(any(sums(ancestry(parent))) for parent in neighbors[root])


def replay(executable, support, stride, neighbors, distances):
    hard = {root: 0 for root in ROOTS}
    reviewed = {root: 0 for root in ROOTS}
    seen = set()
    command = [str(executable), "1", "64", "-1", str(support), "--emit-hard"]
    with Popen(command, stdout=PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.split()
            require(len(fields) >= 2, f"Malformed output {line!r}")
            root = int(fields[1])
            require(root in ROOTS, f"Unexpected root {root}")
            if fields[0] == "hard":
                vertices = tuple(map(int, fields[2::2]))
                amounts = tuple(map(int, fields[3::2]))
                require(len(set(vertices)) == support and len(amounts) == support
                        and all(0 <= v < 64 and v != root for v in vertices)
                        and sum(amounts) == 64 and min(amounts) > 0, f"Invalid placement {line!r}")
                if hard[root] % stride == 0:
                    require(can_descend(root, vertices, amounts, neighbors, distances),
                            f"Independent check failed for {line.strip()}")
                    reviewed[root] += 1
                hard[root] += 1
            else:
                require(fields[0] == "root" and root not in seen, f"Unexpected output {line!r}")
                counts = {key: int(value) for key, value in zip(fields[2::2], fields[3::2])}
                require(counts['support'] == support and counts['unresolved'] == 0,
                        f"Incomplete support check {line.strip()}")
                require(counts['leaves'] == counts['funnel'] + hard[root]
                        and hard[root] == counts['sweep'] + counts['ancestry'], f"Miscount {line.strip()}")
                require(not (support == 5 and counts['partitions'] != 6747), "Wrong partition count")
                seen.add(root)
                print(f"{line.strip()} independently_checked {reviewed[root]}", flush=True)
        require(process.wait() == 0 and seen == set(ROOTS), "Missing or failed run")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--audit-stride", type=int, default=1)
    arguments = parser.parse_args()
    require(arguments.audit_stride >= 1, "The audit stride must be positive")
    factor, neighbors, distances = square_data()
    check_geometry(factor, distances)
    require(sum(sum(p) == 64 for p in combinations_with_replacement(range(1, 61), 5)) == 6747,
            "Incorrect five-part control count")
    for support in range(1, 6):
        replay(arguments.executable.resolve(), support, arguments.audit_stride, neighbors, distances)


if __name__ == "__main__":
    main()
