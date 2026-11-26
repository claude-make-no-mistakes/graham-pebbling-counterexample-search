#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int order = 64;
constexpr int max_piles = 5;
using Mask = std::uint64_t;
using Row = std::array<int, order>;
using Distances = std::array<Row, order>;
using Placement = std::array<int, max_piles>;

// Figure 1 of Pulaj, Wood, and Yerger (arXiv:2411.19314), zero-based here.
constexpr std::array<std::pair<int, int>, 13> edges{{
    {0, 1}, {0, 6}, {1, 2}, {2, 3}, {2, 5}, {2, 7}, {3, 4},
    {3, 6}, {4, 5}, {4, 6}, {4, 7}, {5, 6}, {6, 7}
}};

struct Square {
    std::array<std::vector<int>, order> neighbors;
    Distances shift{};

    Square() {
        for (auto [a, b] : edges) {
            for (int i = 0; i < 8; ++i) {
                add(a * 8 + i, b * 8 + i);
                add(i * 8 + a, i * 8 + b);
            }
        }
        for (int source = 0; source < order; ++source) {
            auto& distances = shift[source];
            distances.fill(-1);
            distances[source] = 0;
            std::queue<int> pending;
            pending.push(source);
            while (!pending.empty()) {
                const int vertex = pending.front();
                pending.pop();
                for (int next : neighbors[vertex]) {
                    if (distances[next] == -1) {
                        distances[next] = distances[vertex] + 1;
                        pending.push(next);
                    }
                }
            }
        }
    }

    void add(int a, int b) {
        neighbors[a].push_back(b);
        neighbors[b].push_back(a);
    }
};

struct Search {
    const Square& graph;
    int root;
    int piles;
    bool emit_hard;
    Placement amount{}, vertex{};
    std::array<Mask, order + 1> single{};
    // [amount at u][amount at v][u] indexes the mask of admissible v.
    std::vector<Mask> pair = std::vector<Mask>((order + 1) * (order + 1) * order);
    std::uint64_t partitions = 0, leaves = 0, unresolved = 0;
    std::uint64_t funnel_solved = 0, sweep_solved = 0, ancestry_solved = 0;

    // Each of the five lanes has six value bits and one borrow guard.
    static constexpr Mask guards = (Mask{1} << 6) | (Mask{1} << 13) |
                                   (Mask{1} << 20) | (Mask{1} << 27) | (Mask{1} << 34);

    static bool covers(Mask upper, Mask lower) {
        return (((upper | guards) - lower) & guards) == guards;
    }

    static void insert(std::vector<Mask>& antichain, Mask candidate) {
        for (Mask existing : antichain) {
            if (covers(candidate, existing)) return;
        }
        std::erase_if(antichain, [candidate](Mask existing) { return covers(existing, candidate); });
        antichain.push_back(candidate);
    }

    Search(const Square& square, int target, int support, bool list_hard)
        : graph(square), root(target), piles(support), emit_hard(list_hard) {
        const auto& droot = graph.shift[root];
        for (int u = 0; u < order; ++u) {
            for (int a = 1; a < (1 << droot[u]); ++a) {
                single[a] |= Mask{1} << u;
                for (int v = 0; v < order; ++v) {
                    if (u == v) continue;
                    int threshold = 1 << droot[v];
                    for (int w = 0; w < order; ++w) {
                        const int remaining = (1 << droot[w]) - (a >> graph.shift[u][w]);
                        threshold = std::min(threshold, std::max(0, remaining) << graph.shift[v][w]);
                    }
                    for (int b = 1; b < threshold; ++b) {
                        pair[index(a, b, u)] |= Mask{1} << v;
                    }
                }
            }
        }
    }

    static int index(int a, int b, int u) {
        return (a * (order + 1) + b) * order + u;
    }

    bool funnel() const {
        for (int w = 0; w < order; ++w) {
            int delivered = 0;
            for (int i = 0; i < piles; ++i) {
                delivered += amount[i] >> graph.shift[vertex[i]][w];
            }
            if (delivered >= (1 << graph.shift[w][root])) return true;
        }
        return false;
    }

    bool sweep(bool reverse, bool prefer_odd) const {
        Row counts{};
        for (int i = 0; i < piles; ++i) counts[vertex[i]] = amount[i];
        for (int layer = 6; layer > 0; --layer) {
            for (int i = 0; i < order; ++i) {
                const int u = reverse ? order - 1 - i : i;
                if (graph.shift[root][u] != layer) continue;
                while (counts[u] >= 2) {
                    int best = -1;
                    int best_score = -1;
                    for (int next : graph.neighbors[u]) {
                        if (graph.shift[root][next] != layer - 1) continue;
                        const int score = counts[next] + (prefer_odd && counts[next] % 2 ? order : 0);
                        if (score > best_score) {
                            best = next;
                            best_score = score;
                        }
                    }
                    counts[u] -= 2;
                    ++counts[best];
                }
            }
        }
        return counts[root] > 0;
    }

    bool descending_ancestry() const {
        Mask capacity = 0;
        std::array<std::vector<Mask>, order> possible;
        for (int i = 0; i < piles; ++i) {
            capacity |= Mask(amount[i]) << (7 * i);
            possible[vertex[i]].push_back(Mask{1} << (7 * i));
        }
        for (int layer = 6; layer > 0; --layer) {
            for (int u = 0; u < order; ++u) {
                if (graph.shift[root][u] != layer) continue;
                const auto& choices = possible[u];
                std::vector<Mask> merged;
                for (std::size_t i = 0; i < choices.size(); ++i) {
                    for (std::size_t j = i; j < choices.size(); ++j) {
                        const Mask consumption = choices[i] + choices[j];
                        if ((consumption & guards) == 0 && covers(capacity, consumption)) {
                            if (layer == 1) return true;
                            insert(merged, consumption);
                        }
                    }
                }
                for (int next : graph.neighbors[u]) {
                    if (graph.shift[root][next] != layer - 1) continue;
                    for (Mask consumption : merged) insert(possible[next], consumption);
                }
            }
        }
        return false;
    }

    void check_leaf() {
        ++leaves;
        if (funnel()) {
            ++funnel_solved;
            return;
        }
        if (emit_hard) print_placement("hard");
        for (bool reverse : {false, true}) {
            for (bool prefer_odd : {false, true}) {
                if (sweep(reverse, prefer_odd)) {
                    ++sweep_solved;
                    return;
                }
            }
        }
        if (descending_ancestry()) {
            ++ancestry_solved;
            return;
        }
        if (++unresolved <= 100) {
            print_placement("pending");
        }
    }

    void print_placement(const char* label) const {
        std::cout << label << ' ' << root;
        for (int i = 0; i < piles; ++i) std::cout << ' ' << vertex[i] << ' ' << amount[i];
        std::cout << '\n';
    }

    void place(int depth, std::array<Mask, max_piles> allowed) {
        if (depth >= piles || depth >= max_piles) {
            check_leaf();
            return;
        }
        while (allowed[depth]) {
            const int u = std::countr_zero(allowed[depth]);
            allowed[depth] &= allowed[depth] - 1;
            vertex[depth] = u;
            auto next = allowed;
            bool viable = true;
            for (int j = depth + 1; j < piles; ++j) {
                auto& choices = next.at(j);
                const int weight = amount.at(j);
                choices &= pair[index(amount[depth], weight, u)];
                if (amount[depth] == weight) choices &= ~(Mask{1} << u) & (~Mask{0} << u);
                if (!choices) {
                    viable = false;
                    break;
                }
            }
            if (viable) place(depth + 1, next);
        }
    }

    void partition(int depth, int remaining, int upper) {
        if (depth >= piles || depth >= max_piles) {
            if (remaining != 0) return;
            ++partitions;
            std::array<Mask, max_piles> allowed{};
            for (int i = 0; i < piles; ++i) allowed[i] = single[amount[i]];
            place(0, allowed);
            return;
        }
        const int slots = piles - depth;
        const int lower = (remaining + slots - 1) / slots;
        for (int a = lower; a <= std::min(upper, remaining - slots + 1); ++a) {
            amount[depth] = a;
            partition(depth + 1, remaining - a, a);
        }
    }

    void run(int minimum, int maximum) {
        for (int a = std::max(minimum, (order + piles - 1) / piles);
             a <= std::min(maximum, order - piles + 1); ++a) {
            amount[0] = a;
            if (single[a]) partition(1, order - a, a);
        }
        std::cout << "root " << root << " support " << piles
                  << " partitions " << partitions << " leaves " << leaves
                  << " funnel " << funnel_solved << " sweep " << sweep_solved
                  << " ancestry " << ancestry_solved
                  << " unresolved " << unresolved << '\n';
    }
};

int main(int argc, char** argv) {
    if (argc < 3 || argc > 6) {
        std::cerr << "usage: search MIN_MAX_PILE MAX_MAX_PILE [ROOT_0_TO_63|-1 [SUPPORT [--emit-hard]]]\n";
        return 2;
    }
    const int lower = std::stoi(argv[1]);
    const int upper = std::stoi(argv[2]);
    const int root = argc >= 4 ? std::stoi(argv[3]) : -1;
    const int support = argc >= 5 ? std::stoi(argv[4]) : max_piles;
    const bool emit_hard = argc == 6;
    if (lower < 1 || upper > order || upper < lower || root < -1 || root >= order ||
        support < 1 || support > max_piles || (emit_hard && std::string(argv[5]) != "--emit-hard")) {
        throw std::invalid_argument("invalid pile range, root, support, or output flag");
    }
    const Square square;
    if (root >= 0) {
        Search(square, root, support, emit_hard).run(lower, upper);
        return 0;
    }
    // Factor orbits under the permutation group on Figure 1 labels 4, 6, 8.
    const std::array representatives{0, 1, 2, 3, 4, 6};
    for (int i : representatives) {
        for (int j : representatives) {
            if (j >= i && (i == 1 || i == 4) && (j == 1 || j == 4)) {
                Search(square, i * 8 + j, support, emit_hard).run(lower, upper);
            }
        }
    }
}
