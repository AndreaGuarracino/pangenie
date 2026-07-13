# PanGenie Stage 9 Optimization Audit and Measured Results

Date: 2026-07-13

Repository: `pangenie-stage9-optimized`

Base commit: `e2eb0fba52ae3fbedce5e6b8c86f551717d72557`

Base branch: `aou-perf-stage9-parallel-index-read`

This report supersedes the conclusions in
`/home/guarracino/Desktop/pangenie-perf/PANGENIE_STAGE9_OPTIMIZATION_REPORT.md`.
That file remains untouched as provenance. The implementation described here is in
this repository and has not been committed.

The optimization contract is strict:

1. Reduce end-to-end runtime, not just an isolated microbenchmark.
2. Preserve genotype and serialized-index behavior.
3. Keep peak memory at or below the Stage 9 baseline.
4. Do not trade exact k-mer membership or counts for an approximation.

## 1. Result

One cold-cache chr22 run of the complete patch was compared with the recorded cold
Stage 9 run using the same v4.2.1 index, FASTQ, sample name, 16 counting threads,
16 PanGenie threads, and Jellyfish `-e 3000000000` setting.

| Phase | Stage 9 baseline | Optimized | Change |
|---|---:|---:|---:|
| Read `UniqueKmersMap` | 0.3967 s | 0.3762 s | -5.2% |
| Count read k-mers | 45.5654 s | 38.7383 s | -15.0% |
| Update unique k-mers | 7.1748 s | 6.9111 s | -3.7% |
| Select paths | 1.2666 s | 1.2587 s | -0.6% |
| HMM | 16.9352 s | 8.6875 s | -48.7% |
| Write VCF | 3.5612 s | 2.2790 s | -36.0% |
| **External wall clock** | **75.08 s** | **58.31 s** | **-22.3%, 1.288x faster** |
| **Peak RSS** | **24,368,540 KiB** | **24,338,424 KiB** | **-30,116 KiB, -0.12%** |

The peak remains dominated by the fixed Jellyfish hash. Memory before Jellyfish is
materially lower: PanGenie's reported RSS immediately after loading the index fell
from 0.137016 GB to 0.101008 GB, a 26.3% reduction.

These runtime values are one cold run per binary, not a statistically controlled
benchmark series. The 22.3% result is large enough to justify the patch, but the
phase-level percentages should be treated as directional until repeated in
alternating order and on full WGS data.

## 2. Correctness and compatibility gates

The optimized binary consumed the unchanged v4.2.1 cereal index and read all
219,553 indexed variants.

| Gate | Result |
|---|---|
| VCF records | 331,933 in both runs |
| VCF body plus `#CHROM` MD5, excluding `##` metadata | `39faa95b84f4d2d9b1d1a0827c2c8ced` in both runs |
| Histogram MD5 | `568d5d9868d2e51b508d228faf0f694b` in both runs |
| Histogram peaks | 31 (2,220,980), 63 (44,563) in both runs |
| Existing cereal index load | Passed |
| Unit tests | 2,947 assertions in 218 test cases passed |
| CTest registration | One production-linked test target, passed |

The body hash deliberately excludes `##fileDate` and other metadata lines. It
includes the column header and every data byte, so it is stronger than GT-only
intersection comparisons.

The allele-container and likelihood-container changes use custom cereal save/load
functions that retain the legacy map-shaped wire representation. Old indexes load,
and new indexes remain readable by code expecting that representation.

## 3. Corrections to the previous report

Several claims in the previous report were either incorrect or insufficiently
supported by its artifacts.

| Previous claim or recommendation | Corrected finding |
|---|---|
| The approximately 24-25 GB peak came from `CopyNumber` and allele maps. | The peak is the fixed Jellyfish hash selected by `-e 3e9`, rounded to its internal power-of-two capacity. The maps matter before and after counting, but they are not the 24 GB peak. |
| b2s had identical histogram semantics. | b2s built a histogram over selected query keys, whereas the Jellyfish path includes the complete graph-k-mer key universe. Matching genotypes on chr22/chrX does not establish semantic equivalence. b2s also saturates counts at 65,535. |
| PanGenie uses OpenMP scheduling. | This fork uses its own `ThreadPool`; OpenMP scheduling clauses are not a drop-in optimization here. |
| zlib-ng/libdeflate would speed the measured VCF writer. | The current genotype writer emits plain VCF through `std::ofstream`. Deflate libraries do not affect this path unless compressed output is added. |
| Current floating output is universally `%g` precision 6. | The writer uses field-specific stream precision: GL uses precision 4 and AF uses precision 6. A replacement formatter must reproduce each field exactly. |
| The current HMM recurrence was already fully optimal. | Forward/backward was reduced previously, but still stored dead matrix halves, allocated per column, repeated lookups, and recomputed weighted cells. Viterbi remained O(P^4). There was substantial exact work left. |
| `unordered_map` is automatically safe for unique-k-mer selection. | Selection iterates sorted keys and stops after a cap. Changing iteration order can change which k-mers are retained and therefore the model. A keyed lookup can use hashing; capped selection cannot do so without restoring a deterministic ordering rule. |
| `vector<bool>` to `uint8_t` reduces memory. | A byte per flag is approximately 8x the packed storage of `vector<bool>`. It can be faster, but it violates a memory-reduction claim unless the allocation is negligible and measured. This patch retains packed flags where appropriate. |
| Dense vectors are always appropriate for multiallelic IDs. | Allele IDs can be sparse. The safe representation is a sorted flat vector of present `(allele_id, info)` pairs, not a vector sized by the maximum ID. |
| The reported WGS split of 5,160/1,938/22/18 seconds was established by raw logs. | No immutable raw artifact for those numbers was available in the audited tree. They may be operational observations, but they are not used as measured evidence here. |

The original b2s experiment remains useful as a negative result: its cold runs were
slower and its WGS memory projection was unfavorable. It should not be promoted as
an exact replacement on the strength of two chromosome-level genotype comparisons.

## 4. Implemented optimizations

### 4.1 HMM core

Files: `src/hmm.cpp`, `src/hmm.hpp`, `src/columnindexer.cpp`,
`src/columnindexer.hpp`, `src/emissionprobabilitycomputer.*`, and
`src/transitionprobabilitycomputer.*`.

- Store only the upper triangle of symmetric forward/backward state matrices.
  Dead lower-triangle writes and normalization passes are removed.
- Pool and reuse HMM columns instead of allocating their vectors for every graph
  column.
- Remove backward row sums that are computed but never consumed.
- Compute each weighted backward cell once and reuse it.
- Replace repeated path scans in `ColumnIndexer` with direct aligned indexing and a
  checked fallback.
- Expose allele spans and inline the hottest allele/emission/transition accessors.
- Build transition probabilities on the stack, use `std::array<double, 3>`, and
  evaluate the shared exponential once.
- Pass the emission computer by reference and retain a non-owning pointer in the
  hot object rather than incrementing a shared-pointer reference count.
- Rewrite Viterbi from O(P^4) to exact O(P^2) per column using row/column maxima and
  top-two exclusion. The original largest-flat-index tie rule is preserved.
- Pass likelihood results by reference instead of deep-copying small containers.

The measured HMM time fell from 16.94 to 8.69 seconds. `objdump` inspection found
no remaining PLT calls for the hottest emission, transition, or allele accessor
symbols.

### 4.2 Resident data structures

Files: `src/copynumber.*`, `src/biallelicuniquekmers.*`,
`src/multiallelicuniquekmers.*`, and `src/genotypingresult.*`.

- Replace `CopyNumber`'s heap-backed three-element vector with
  `std::array<double, 3>` and fix its incorrect equality comparison.
- Replace the biallelic `std::map<bool, AlleleInfo16>` with a two-slot array plus
  presence bits.
- Replace the multiallelic red-black tree with a sorted flat vector, preserving
  sparse allele IDs and deterministic iteration.
- Replace the tiny genotype-likelihood tree with a sorted flat vector and reuse an
  output buffer.
- Preserve the legacy cereal map representation at the serialization boundary.

These changes remove millions of tiny node allocations and pointer chases without
changing the on-disk index contract.

### 4.3 Jellyfish counting and histogram

Files: `src/jellyfishcounter.*`, `src/jellyfishreader.*`, and
`src/kmercounter.hpp`.

- Compute the histogram with disjoint Jellyfish region slices in parallel.
- Use one private integer histogram per worker and reduce bins exactly in a fixed
  integer sum. There is no floating reduction-order change.
- Query the batch through one Jellyfish array handle instead of repeating the
  virtual/interface path per k-mer.
- Parse query strings through `string_view` and reuse a thread-local `mer_dna`.
- Increase parser blocks from 4 KiB to 64 KiB. At 16 workers this adds only a few
  MiB and did not increase peak RSS.
- Avoid `std::endl` in histogram output.

This is the main source of the measured 6.83-second counting reduction. It does not
replace Jellyfish or change its exact counts.

### 4.4 K-mer sidecar parsing and lookup

Files: `src/kmerparser.*` and the genotype path in `src/commands.cpp`.

- Replace per-line `istringstream` and token-vector construction with a manual
  `string_view` scanner and `from_chars` for integers.
- Give zlib a 1 MiB input buffer and read 64 KiB chunks.
- Reserve and reuse parse and query buffers.
- Batch the variant and flank queries while retaining backing storage until every
  view is consumed.
- Avoid repeated variant-map lookups and unnecessary `CopyNumber` temporaries.

The text gzip sidecar is still an avoidable boundary. Section 7 describes the
versioned binary replacement rather than hiding this cost behind more parser
micro-optimizations.

### 4.5 VCF writer

Files: primarily `src/graph.cpp` and `src/variant.*`.

- Give the output stream a 1 MiB buffer.
- Replace record-level `std::endl` with newline characters.
- Write INFO and GL fields directly instead of constructing a per-record
  `ostringstream` and temporary ALT string.
- Avoid deep copies of `Variant`, `GenotypingResult`, and sampled-panel objects on
  the common no-filter path.
- Add a non-mutating allele-without-flanks accessor. Copy and filter only when
  missing alleles require it.
- Preserve the existing stream formatting behavior for AF and GL fields.

The measured writer time fell from 3.56 to 2.28 seconds. This is a plain-text
buffering and allocation win, not a compression-library result.

### 4.6 General parsing and allocation cleanup

Files: `src/dnasequence.*`, `src/fastareader.*`, `src/graphbuilder.cpp`,
`src/variantreader.cpp`, `src/stepwiseuniquekmercomputer.cpp`, and
`src/threadpool.cpp`.

- Pass strings and sequence objects by const reference where ownership is not
  transferred.
- Decode directly into resized sequence storage and compare sequences without
  temporary strings.
- Use one map lookup in FASTA accessors and allocate shared objects directly with
  `make_shared`.
- Replace per-variant regular expressions with direct base validation.
- Parse positions with `from_chars`.
- Buffer graph/sidecar writes, render each sidecar line once, and avoid flushing on
  every record.
- Use explicit lambda captures in the thread pool.

### 4.7 Build and validation hygiene

Files: `CMakeLists.txt`, `tests/CMakeLists.txt`, and `.gitignore`.

- Enable CTest at the project level and register the test binary.
- Link tests against the actual production `PanGenieLib` instead of recompiling a
  second private copy of implementation files.
- Copy fixtures into the build tree so tests do not modify source fixtures.
- Compile tests as C++20, matching the production target.
- Add `-fno-semantic-interposition` for GNU/Clang so internal shared-library hot
  calls can be bound and inlined safely.
- Ignore the out-of-source `build/` directory.

Before this repair, `ctest` discovered zero tests, so a green CTest invocation was
not evidence that PanGenie had been exercised.

## 5. Memory model

The `-e 3000000000` Jellyfish sizing parameter produces a fixed-capacity hash that
dominates this test's 24.3 GB peak. Optimizing PanGenie's resident C++ objects cannot
substantially lower that peak while the same hash sizing remains in use.

That does not make the resident reductions irrelevant:

- They reduce index-load RSS by approximately 36 MB on chr22.
- They reduce allocation traffic and improve HMM time.
- They create headroom for a future exact static graph-k-mer dictionary.
- They matter for phases and deployments that do not retain the large Jellyfish
  hash.

An honest next memory target must therefore change the counting architecture while
retaining exact graph-key membership, integer counts, zero-count graph keys, and the
same histogram/selection rules. Swapping allocators or serializers alone cannot
remove the fixed hash.

## 6. Research findings

### 6.1 Current PanGenie direction

The 2026 PanGenie scaling preprint describes aggressive reference-path sampling and
reports large improvements over older chunking strategies. This fork already
contains the corresponding reduced-path workflow, so those paper-level gains are
context, not an additional optimization that can be claimed for this patch. The
current chr22 run sampled one subset of 89 paths.

Primary source:
[PanGenie scaling preprint](https://www.biorxiv.org/content/10.64898/2026.06.29.735275v1.full.pdf).

Upstream users have explicitly requested reuse of precomputed KMC/KFF counts because
k-mer counting dominates large workflows. That is the right architectural
direction, but the earlier attempt's large `unordered_set` of graph k-mers simply
moved the memory failure into membership filtering.

Primary sources:
[PanGenie issue 62](https://github.com/eblerjana/pangenie/issues/62),
[PanGenie repository](https://github.com/eblerjana/pangenie), and
[KFF specification](https://github.com/Kmer-File-Format/kff-reference).

### 6.2 Exact static graph-k-mer dictionary

The strongest next design is an immutable dictionary of every canonical graph
k-mer with a stable dense ID:

1. Build the dictionary once with the PanGenie index.
2. Initialize a dense count array to zero so absent graph k-mers retain exact zero
   counts and histogram semantics.
3. Stream a KFF/KMC count file and update only verified dictionary members.
4. Store selected k-mers in the index as dense IDs rather than 31-base text.

SSHash is a credible research candidate because it is designed for compact static
k-mer dictionaries. Its published HPRC experiment reports 3.718 billion canonical
31-mers at 11.93 bits/k-mer, 5.54 GB, with a 4 minute 45 second build. Those are
external results, not PanGenie measurements, and the exact non-member behavior and
build memory must be validated before adoption.

Primary sources:
[SSHash paper](https://jermp.github.io/assets/pdf/papers/2026.01.21.700884v1.full.pdf) and
[SSHash implementation](https://github.com/jermp/sshash).

PtrHash or another minimal perfect hash can provide dense IDs at still lower index
cost, but an MPHF alone is not an exact membership structure: a non-key also maps
to an ID. It therefore requires full-key verification or a collision-free exact
fingerprint scheme. A Bloom filter alone is also unacceptable because false
positives can alter counts and selection.

Primary source:
[PtrHash paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2025.21).

KMC remains a strong producer for reusable count databases, but it should feed the
exact dictionary rather than materialize all graph k-mers in a node-based hash set.

Primary sources:
[KMC repository](https://github.com/refresh-bio/KMC),
[KMC 3 paper](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399),
and [KFF paper](https://academic.oup.com/bioinformatics/article/38/18/4423/6651834).

### 6.3 Compression and input delivery

libdeflate is a high-performance DEFLATE implementation, but it intentionally does
not provide a streaming gzip API. It is useful if PanGenie adopts block-compressed
output or whole-block sidecars; it is not a drop-in answer for the current plain
VCF writer.

Primary source:
[libdeflate implementation and API constraints](https://github.com/ebiggers/libdeflate).

For large compressed FASTQ, decompression should occur outside the fixed Jellyfish
worker pool and be streamed through stdin, using a bounded producer such as `pigz`
when the existing CLI/pipeline supports it. This avoids materializing a decompressed
FASTQ while allowing parallel inflate. It must be benchmarked on the production
storage path because local chr22 counting is CPU-heavy, whereas remote/localized WGS
can be I/O-heavy.

### 6.4 Serialization and output templates

A faster serializer helps only the approximately 0.4-second chr22 load measured
here. The larger win is to stop serializing pointer-rich objects:

- Store allele data in flat, versioned arrays.
- Store k-mers as dense IDs.
- Memory-map immutable arrays when this reduces, rather than duplicates, resident
  pages.
- Keep a compatibility reader for the legacy cereal format.

Similarly, much of each VCF record is index-invariant. A versioned index can store
pre-rendered CHROM/POS/REF/ALT/INFO/FORMAT templates and append only sample-specific
GT/GQ/GL/KC fields. This attacks both CPU formatting and repeated string storage,
but must be evaluated against index size and page-cache pressure.

## 7. Prioritized next work

### P0: establish a statistically defensible WGS gate

- Freeze binaries, inputs, index hashes, compiler, CPU model, and command line.
- Run at least five cold repetitions per binary in alternating order.
- Pin CPU affinity and record frequency/power state if the environment permits.
- Report median, range, CPU time, page faults, filesystem bytes, and peak RSS.
- Compare the complete decompressed VCF body and histogram, not only intersecting
  GT fields.
- Run the same gate on a real WGS sample before merging.

### P1: remove avoidable input materialization

- Benchmark bounded `pigz -dc | PanGenie ... -i -` delivery where supported.
- Add a first-class KFF/KMC count-input option with an explicit count-width and
  canonicalization contract.
- Keep the current Jellyfish path as the reference implementation and exactness
  oracle.

### P2: build the exact static dictionary

- Prototype SSHash and PtrHash plus full-key verification on the actual graph-k-mer
  universe.
- Measure build peak, lookup throughput, serialized size, and false-member behavior.
- Preserve all graph keys, including zero-count keys, when computing the histogram.
- Use a dense integer count array with enough width to match Jellyfish values; do
  not silently saturate at 65,535.

This is the only researched route that can plausibly lower the approximately 24 GB
peak substantially while also eliminating repeated FASTQ counting across samples.

### P3: version the k-mer sidecar

- Replace `*_kmers.tsv.gz` text with a packed binary sidecar containing variant
  boundaries, dense IDs, flank metadata, and checksums.
- Keep a legacy text reader for existing indexes.
- Batch lookup by ID, eliminating gzip, decimal parsing, DNA text parsing, and
  repeated canonicalization at genotype time.

### P4: flatten the remaining index and VCF invariants

- Move `UniqueKmersMap` toward CSR/offset arrays suitable for read-only mapping.
- Pre-render index-invariant VCF record fragments.
- Measure index-size and page-cache effects; reject changes that only move heap RSS
  into duplicated mapped pages.

### P5: controlled code-generation experiments

- Train PGO on representative chr22, chrX, and WGS inputs, then validate on a held-
  out sample.
- Compare portable AVX2, host-native, and production-VM tuning separately.
- Test jemalloc/mimalloc only after the node-heavy allocations removed here; retain
  an allocator only if both wall time and peak RSS meet the gate.

## 8. Approaches not to ship

- b2s as an exact Jellyfish replacement without restoring the complete graph-key
  histogram and full count range.
- A Bloom filter as the final membership decision.
- A bare MPHF lookup without exact key verification.
- An unordered container in any capped, iteration-order-sensitive selection step.
- `vector<bool>` to one-byte flags under the claim that it reduces memory.
- zlib-ng or libdeflate as an explanation for plain-VCF writer speed.
- A serializer-only rewrite before flattening the data layout.
- Float HMM state or reassociated floating reductions without a full WGS output
  identity gate.
- `-ffast-math`; it permits transformations that violate exact likelihood behavior.

## 9. Reproduction

Build and test:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPKG_CONFIG_EXECUTABLE=/usr/bin/pkg-config
cmake --build build -j16
ctest --test-dir build --output-on-failure
```

Optimized chr22 command:

```bash
/usr/bin/time -v build/src/PanGenie \
  -f /home/guarracino/Desktop/pangenie-perf/builds/v4.2.1/build/index_chr22/idx \
  -i /home/guarracino/Desktop/pangenie-perf/data/chr22/reads.fq \
  -s hg002_chr22 \
  -j 16 -t 16 -e 3000000000 \
  -o /tmp/pangenie-opt-final
```

Body and record gates:

```bash
grep -v '^##' /tmp/pangenie-opt-final_genotyping.vcf | md5sum
grep -vc '^#' /tmp/pangenie-opt-final_genotyping.vcf
md5sum /tmp/pangenie-opt-final_histogram.histo
```

Expected outputs:

```text
39faa95b84f4d2d9b1d1a0827c2c8ced  -
331933
568d5d9868d2e51b508d228faf0f694b  /tmp/pangenie-opt-final_histogram.histo
```

The recorded Stage 9 baseline is in
`/home/guarracino/Desktop/pangenie-perf/tool-eval/cold_s9.log`; the original cold
comparison summary and driver are `tool-eval/cold_chr22.txt` and
`tool-eval/cold_chr22.sh` under that external performance workspace.
