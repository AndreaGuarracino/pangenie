# PanGenie Stage 13 Novel Optimizations

Date: 2026-07-13

Branch: `aou-perf-stage13-novel-optimizations`

Base: the uncommitted Stage 9 optimized implementation described in
`STAGE9_OPTIMIZATION_AUDIT_AND_RESULTS.md`.

## Confirmed default: automatic Jellyfish capacity

Stage 9 requested three billion Jellyfish slots for every input. Jellyfish rounds
that request to a power of two, producing a 4,294,967,296-slot table and the
measured 24 GB peak even for a reduced chr22 graph with only 43.8 million distinct
31-mers.

Stage 13 makes `-e 0` the default. In graph-only mode, the graph FASTA is the only
input allowed to create keys, so its byte size is used as a conservative initial
capacity request. Jellyfish rounds it to a power of two and retains its built-in
cooperative size doubling. A sizing underestimate can therefore cost time and
temporary memory but cannot silently discard keys or change counts. A nonzero
`-e` remains an explicit override. Stdin, file errors, and inputs at least 3 GB
fall back to the previous three-billion request.

For the v4.2.1 reduced chr22 graph:

| Metric | Stage 9 optimized | Stage 13 auto | Change |
|---|---:|---:|---:|
| Requested slots | 3,000,000,000 | 182,719,605 | -93.9% |
| Allocated slots | 4,294,967,296 | 268,435,456 | -93.8% |
| Count reads | 38.7383 s | 24.5200 s | -36.7% |
| End-to-end wall | 58.31 s | 41.10 s | -29.5% |
| Peak RSS | 24,338,424 KiB | 1,752,364 KiB | -92.8% |

Relative to the original Stage 9 baseline, Stage 13 is 45.3% faster
(`75.08 -> 41.10 s`) and uses 92.8% less peak RSS
(`24,368,540 -> 1,752,364 KiB`).

Correctness gates match:

- VCF body plus column header MD5: `39faa95b84f4d2d9b1d1a0827c2c8ced`
- VCF records: `331933`
- Histogram MD5: `568d5d9868d2e51b508d228faf0f694b`
- Histogram peaks: `31 (2220980), 63 (44563)`
- Production-linked unit tests: passed

One tighter trial requested 90,000,000 slots and allocated 134,217,728. It used
943,372 KiB but took 41.94 seconds, so it is not the runtime-default choice.

## Full-graph constraint

The old full chr22 graph is qualitatively different from the reduced v4.2.1 graph:

| Graph | FASTA bytes | Distinct canonical 31-mers |
|---|---:|---:|
| Reduced v4.2.1 chr22 | 182,719,605 | at most 43,782,833 (non-canonical KMC audit upper bound) |
| Full Stage 9 chr22 | 3,349,381,242 | 2,516,559,577 |

The full graph genuinely needs the existing 4.29-billion-slot Jellyfish capacity
at a reasonable load factor. Auto-sizing deliberately falls back to the old
request for this case. Lowering `-e` alone cannot solve full-HPRCv2 memory.

## KMC replacement experiment

KMC 3.2.4 was measured as an external exact counter:

| Workload | Wall | Peak RSS | Temporary data |
|---|---:|---:|---:|
| chr22 reads, 3.3 GB FASTQ | 9.13 s | 7,623,776 KiB | 1,245 MB |
| full chr22 graph, strict 8 GB mode, no database output | 19.07 s | 7,786,840 KiB | 3,053 MB |

KMC is substantially faster than the embedded Jellyfish read count on chr22, but
it is not yet the default because a complete equivalent backend must also:

1. Intersect read counts with the graph database while retaining read-side counts.
2. Return zero for every graph key absent from the reads.
3. Reproduce the graph-only histogram and count-width behavior.
4. Avoid materializing the 2.5-billion-key graph in an `unordered_set`.
5. Use KMC strict-memory mode and enforce temporary-disk limits on full WGS.

The previous KMC memory failure is consistent with duplicating the graph keys in a
node-based set and/or using non-strict KMC memory. It is not evidence that KMC's
disk-backed counter itself requires the old 24 GB peak.

## Next architectural backend

The best full-HPRCv2 design remains an exact static graph dictionary with stable
dense IDs:

- SSHash is exact and returns a dense ID or not-found for a queried k-mer.
- A dense count array preserves explicit zero-count graph keys and exact histogram
  semantics.
- KMC can produce reusable sample counts, with a streaming dictionary intersection.
- A bare minimal perfect hash or Bloom filter is insufficient without exact key
  verification.

For 2.516 billion graph keys, a four-byte dense count array alone is about 9.4 GiB.
This can still fit materially below the old 24 GB table if the exact dictionary is
compact, but build peak, lookup throughput, counter width, and full-WGS wall time
must be measured before adoption.

## CPU-target experiment

`PANGENIE_CPU_ARCH` is now a CMake cache option. The default remains `haswell`;
deployment-specific builds can request `native`. On the Ryzen 7 PRO 7840U test
host, one native run took 47.24 seconds versus the 41.10-second default run. The
runs were not temperature-controlled, so native remains opt-in and is not claimed
as an improvement.
