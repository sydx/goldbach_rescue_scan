// Exhaustive two-shift Goldbach rescue search.
// Build: g++ -O3 -march=native -std=c++17 -pthread goldbach_rescue_scan.cpp -o goldbach_rescue_scan
// Run:   ./goldbach_rescue_scan                    (default n <= 1,000,000,000)
//        ./goldbach_rescue_scan 1000000 --threads 4 --progress 100000
//        ./goldbach_rescue_scan --mode residue 1000000000
//
// For n == 1 (mod 3), the strong conjecture says that, if 2n-3 is
// composite, 2n-2 = t+c for primes t,t+2,c,c+4.  In mode "both"
// (the default), the program also tests the general conditional claim:
// whenever old Goldbach pairs exist and none has a direct +2 continuation,
// such a quartet exists.  A small-prime search is followed by an exhaustive
// fallback, so a reported counterexample is never inferred from a cutoff.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using u64 = std::uint64_t;
using Clock = std::chrono::steady_clock;

struct PrimeBits {
    u64 maximum;
    std::vector<u64> composite;

    explicit PrimeBits(u64 max_value)
        : maximum(max_value), composite(((max_value / 2 + 1) + 63) / 64, 0) {
        for (u64 p = 3; p <= maximum / p; p += 2) {
            if (!prime(p)) continue;
            for (u64 j = p * p / 2; j <= maximum / 2; j += p)
                composite[j >> 6] |= u64{1} << (j & 63);
        }
    }

    bool prime(u64 x) const {
        if (x == 2) return true;
        if (x < 3 || x > maximum || (x & 1) == 0) return false;
        const u64 j = x / 2;
        return !(composite[j >> 6] & (u64{1} << (j & 63)));
    }
};

struct Options {
    u64 limit = 1'000'000'000;
    u64 progress = 10'000'000;
    u64 small_bound = 1'000'000;
    u64 block_size = 10'000;
    unsigned threads = 1;
    unsigned samples = 8;
    std::string mode = "both";
};

static u64 number(const std::string& text) {
    std::size_t used = 0;
    const u64 value = std::stoull(text, &used);
    if (used != text.size()) throw std::invalid_argument("invalid integer: " + text);
    return value;
}

static Options parse(int argc, char** argv) {
    Options o;
    bool positional_seen = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [MAX_N] [options]\n"
                      << "  --mode both|general|residue   default: both\n"
                      << "  --threads N                   default: 1\n"
                      << "  --progress N                  print after approximately N more pivots; 0 disables\n"
                      << "  --small-bound N               small-pattern search radius; default: 1000000\n"
                      << "  --samples N                   number of example witnesses; default: 8\n"
                      << "  --block-size N                pivots per worker block; default: 10000\n"
                      << "MAX_N defaults to 1000000000. Search and primality tests are exhaustive.\n";
            std::exit(0);
        }
        auto value = [&]() -> std::string {
            if (i + 1 == argc) throw std::invalid_argument("missing value for " + arg);
            return argv[++i];
        };
        if (arg == "--mode") o.mode = value();
        else if (arg == "--threads") o.threads = static_cast<unsigned>(number(value()));
        else if (arg == "--progress") o.progress = number(value());
        else if (arg == "--small-bound") o.small_bound = number(value());
        else if (arg == "--samples") o.samples = static_cast<unsigned>(number(value()));
        else if (arg == "--block-size") o.block_size = number(value());
        else if (arg[0] != '-' && !positional_seen) {
            o.limit = number(arg);
            positional_seen = true;
        } else throw std::invalid_argument("unknown or repeated argument: " + arg);
    }
    if (o.limit < 4 || o.limit > (std::numeric_limits<std::uint32_t>::max() - 2ULL) / 2)
        throw std::invalid_argument("MAX_N must be between 4 and 2147483646");
    if (o.threads == 0 || o.block_size == 0 || o.small_bound < 7)
        throw std::invalid_argument("threads, block-size must be positive; small-bound >= 7");
    if (o.mode != "both" && o.mode != "general" && o.mode != "residue")
        throw std::invalid_argument("mode must be both, general, or residue");
    return o;
}

struct Witness {
    u64 t = 0, c = 0, probes = 0;
    bool fallback = false;
    explicit operator bool() const { return t != 0; }
};

struct Patterns {
    std::vector<std::uint32_t> twins, cousins;
    u64 bound;
    Patterns(const PrimeBits& prime, u64 requested)
        : bound(std::min(requested, prime.maximum - 4)) {
        for (u64 x = 3; x <= bound; x += 2) {
            if (!prime.prime(x)) continue;
            if (prime.prime(x + 2)) twins.push_back(static_cast<std::uint32_t>(x));
            if (prime.prime(x + 4)) cousins.push_back(static_cast<std::uint32_t>(x));
        }
    }
};

// t+c=target. A small twin or a small cousin is tried first. If neither
// works, the middle is scanned completely (no heuristic stopping rule).
static Witness quartet(u64 target, const PrimeBits& p, const Patterns& pat) {
    Witness w;
    if (target < 6 || (target & 1)) return w;
    const u64 upper = target - 3;
    const std::size_t count = std::max(pat.twins.size(), pat.cousins.size());
    for (std::size_t i = 0; i < count; ++i) {
        if (i < pat.twins.size() && pat.twins[i] <= upper) {
            const u64 t = pat.twins[i], c = target - t;
            ++w.probes;
            if (p.prime(c) && p.prime(c + 4)) { w.t = t; w.c = c; return w; }
        }
        if (i < pat.cousins.size() && pat.cousins[i] <= upper) {
            const u64 c = pat.cousins[i], t = target - c;
            ++w.probes;
            if (p.prime(t) && p.prime(t + 2)) { w.t = t; w.c = c; return w; }
        }
    }
    // The small searches covered t <= bound and c <= bound.
    w.fallback = true;
    const u64 first = (pat.bound & 1) ? pat.bound + 2 : pat.bound + 1;
    if (target > pat.bound + 1) {
        const u64 last = std::min(upper, target - pat.bound - 1);
        for (u64 t = first; t <= last; t += 2) {
            if (p.prime(t) && p.prime(t + 2)) {
                const u64 c = target - t;
                ++w.probes;
                if (p.prime(c) && p.prime(c + 4)) {
                    w.t = t; w.c = c; return w;
                }
            }
        }
    }
    return w;
}

// Existence of any old Goldbach pair with a forward twin member.
static bool direct_survival(u64 n, const PrimeBits& p, const Patterns& pat) {
    const u64 sum = 2 * n;
    const u64 upper = sum - 3;
    for (std::uint32_t t : pat.twins) {
        if (t > upper) break;
        if (p.prime(sum - t)) return true;
    }
    const u64 first = (pat.bound & 1) ? pat.bound + 2 : pat.bound + 1;
    for (u64 t = first; t <= upper; t += 2)
        if (p.prime(t) && p.prime(t + 2) && p.prime(sum - t)) return true;
    return false;
}

static bool has_old_pair(u64 n, const PrimeBits& p) {
    for (u64 a = 3; a <= n; a += 2)
        if (p.prime(a) && p.prime(2 * n - a)) return true;
    return false;
}

static std::pair<u64,u64> first_new_pair(u64 n, const PrimeBits& p) {
    const u64 sum = 2 * n + 2;
    for (u64 a = 3; a <= n + 1; a += 2)
        if (p.prime(a) && p.prime(sum - a)) return {a, sum - a};
    return {0,0};
}

struct Stats {
    u64 visited = 0, residue_targets = 0, other_targets = 0;
    u64 direct_survivors = 0, empty_old = 0, quartets = 0;
    u64 left_twin = 0, right_twin = 0, fallback = 0;
    u64 max_probes = 0, failures = 0;
    void add(const Stats& s) {
        visited += s.visited; residue_targets += s.residue_targets;
        other_targets += s.other_targets; direct_survivors += s.direct_survivors;
        empty_old += s.empty_old; quartets += s.quartets;
        left_twin += s.left_twin; right_twin += s.right_twin;
        fallback += s.fallback; failures += s.failures;
        max_probes = std::max(max_probes, s.max_probes);
    }
};

struct Example { u64 n, t, c, probes; bool fallback; };

int main(int argc, char** argv) {
    try {
        std::cout << std::unitbuf;
        const Options o = parse(argc, argv);
        const auto started = Clock::now();
        const u64 maximum = 2 * o.limit + 2;
        const u64 bytes = ((((maximum / 2 + 1) + 63) / 64) * 8);
        std::cout << "MAX_N=" << o.limit << " mode=" << o.mode
                  << " threads=" << o.threads << " sieve_MB="
                  << std::fixed << std::setprecision(1) << bytes / 1048576.0 << '\n';
        const PrimeBits prime(maximum);
        const Patterns pat(prime, o.small_bound);
        std::cout << "sieve ready in "
                  << std::chrono::duration<double>(Clock::now() - started).count()
                  << "s; small twins=" << pat.twins.size()
                  << ", small cousins=" << pat.cousins.size() << '\n';

        std::atomic<u64> next{4}, completed{0}, next_report{o.progress};
        std::atomic<bool> stop{false};
        std::mutex output_mutex;
        Stats totals;
        std::vector<Example> examples;
        auto worker = [&]() {
            while (!stop.load(std::memory_order_relaxed)) {
                const u64 first = next.fetch_add(o.block_size);
                if (first > o.limit) break;
                const u64 last = std::min(o.limit, first + o.block_size - 1);
                Stats local;
                for (u64 n = first; n <= last && !stop.load(std::memory_order_relaxed); ++n) {
                    ++local.visited;
                    const bool residue = (n % 3 == 1);
                    if (residue) {
                        // For n == 1 mod 3, every old pair without 3 has
                        // both endpoints 1 mod 3 and loses its direct +2 shift.
                        if (prime.prime(2 * n - 3)) { ++local.direct_survivors; continue; }
                        ++local.residue_targets;
                    } else {
                        if (o.mode == "residue") continue;
                        if (direct_survival(n, prime, pat)) {
                            ++local.direct_survivors; continue;
                        }
                        if (!has_old_pair(n, prime)) { ++local.empty_old; continue; }
                        ++local.other_targets;
                    }
                    const Witness w = quartet(2 * n - 2, prime, pat);
                    local.fallback += w.fallback;
                    local.max_probes = std::max(local.max_probes, w.probes);
                    if (w) {
                        ++local.quartets;
                        if (w.t + 2 <= w.c) ++local.left_twin;
                        else ++local.right_twin;
                        if (o.samples) {
                            std::lock_guard<std::mutex> lock(output_mutex);
                            if (examples.size() < o.samples)
                                examples.push_back({n,w.t,w.c,w.probes,w.fallback});
                        }
                    } else {
                        // The general statement is conditional on an old
                        // Goldbach pair. The residue statement is not.
                        if (residue && o.mode == "general" && !has_old_pair(n, prime)) {
                            --local.residue_targets;
                            ++local.empty_old;
                            continue;
                        }
                        ++local.failures;
                        const auto new_pair = first_new_pair(n, prime);
                        {
                            std::lock_guard<std::mutex> lock(output_mutex);
                            std::cout << "COUNTEREXAMPLE n=" << n
                                      << " residue=" << n%3
                                      << " 2n-3=" << 2*n-3
                                      << " quartet_probes=" << w.probes
                                      << " old_pair=" << has_old_pair(n, prime)
                                      << " next_pair=" << new_pair.first << '+' << new_pair.second
                                      << std::endl;
                        }
                        stop.store(true, std::memory_order_relaxed);
                        break;
                    }
                }
                {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    totals.add(local);
                }
                const u64 done = completed.fetch_add(local.visited) + local.visited;
                if (o.progress && !stop.load(std::memory_order_relaxed)) {
                    u64 milestone = next_report.load();
                    while (done >= milestone &&
                           next_report.compare_exchange_weak(milestone, milestone + o.progress)) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
                        std::cout << "progress visited=" << done << '/' << (o.limit - 3)
                                  << " strong_targets=" << totals.residue_targets
                                  << " other_extinctions=" << totals.other_targets
                                  << " rate=" << static_cast<u64>(done / std::max(seconds,0.001))
                                  << " pivots/s elapsed=" << std::setprecision(1) << seconds << "s\n";
                    }
                }
            }
        };

        std::vector<std::thread> threads;
        for (unsigned i = 1; i < o.threads; ++i) threads.emplace_back(worker);
        worker();
        for (auto& t : threads) t.join();
        std::sort(examples.begin(), examples.end(),
                  [](const Example& a, const Example& b) { return a.n < b.n; });
        for (const Example& e : examples)
            std::cout << "example n=" << e.n << " twin=" << e.t << ',' << e.t+2
                      << " cousin=" << e.c << ',' << e.c+4
                      << " shift=" << (e.t+2<=e.c ? "-1" : "2")
                      << " probes=" << e.probes << '\n';
        const double seconds = std::chrono::duration<double>(Clock::now() - started).count();
        std::cout << "SUMMARY " << (totals.failures ? "COUNTEREXAMPLE FOUND" : "NO COUNTEREXAMPLE")
                  << " visited=" << totals.visited
                  << " residue_targets=" << totals.residue_targets
                  << " other_extinctions=" << totals.other_targets
                  << " direct_survivors=" << totals.direct_survivors
                  << " empty_old=" << totals.empty_old
                  << " quartets=" << totals.quartets
                  << " first_witness_left=" << totals.left_twin
                  << " first_witness_right=" << totals.right_twin
                  << " exhaustive_fallbacks=" << totals.fallback
                  << " max_candidate_probes=" << totals.max_probes
                  << " failures=" << totals.failures
                  << " elapsed_s=" << std::setprecision(2) << seconds << '\n';
        return totals.failures ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
