# Five-vertex support in the original Lemke square

The support-five pebbling number of the original Lemke square is exactly **64**. In the same graph, every 64-pebble configuration with at least 60 occupied vertices solves every root. Any counterexample to Graham’s conjecture on this square must occupy between 6 and 59 vertices.

[Paper](paper.pdf) · [Source](paper.tex) · [Replay record](REPLAY.txt)

The proof treats roots of eccentricity four and five by shortest-path pooling. At eccentricity six it enumerates every possible assignment of up to five piles, up to symmetry of the root and ordering of equal piles. It constructs a legal solution for every assignment through distance tests, move sweeps or a consumption-vector recurrence. A second implementation checks every placement that survives the distance tests. The lower bound comes from placing 63 pebbles at distance six.

## Replay

The tools need a C++20 compiler and Python 3.10 or newer. The replay uses integer arithmetic and no third-party Python packages.

```sh
g++ -O3 -Wall -Wextra -Wpedantic -std=c++20 search.cpp -o search
python verify.py ./search
```

On Windows, name the executable `search.exe` and pass its path to the same verifier. The full check verifies every emitted placement. Use `--audit-stride 113` for a shorter replay that runs the entire enumeration and independently checks every 113th emitted placement.

## Prior work

Pulaj, Wood and Yerger proved the support-four bound for all products of eight-vertex Lemke graphs. The factor edges here follow their Figure 1, and the independent checker reconstructs the graph from Flocco, Pulaj and Yerger’s figure. The five-pile theorem here treats the square of the original factor. [Bilevel computation](https://arxiv.org/html/2411.19314v1#S4) · [Second graph source](https://arxiv.org/html/2312.12618v1#S2.F1)
