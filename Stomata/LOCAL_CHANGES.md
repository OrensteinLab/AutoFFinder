# Local Stomata Changes for the AutoOffTarget Stage-1 Benchmark

This directory is based on Stomata v0.11.0, upstream commit
`f63d2cc536fe166728e4108cb4d57537c06dac56`. The changes below make Stomata
suitable for generating the raw end positions consumed by
`PostAutoFFinder` as an alternative to the FPGA-based first stage.

## Batch `--no-deduplicate` support

### `src/main.cpp`

The CLI's `disable_deduplication` value is now copied into
`BatchSearchConfig`.

**Why:** The benchmark uses `--spacer-file` and `--no-deduplicate` to request
every raw Myers end-position. Previously, the CLI accepted the option in batch
mode but did not pass it into the batch search, so halo deduplication remained
enabled silently.

### `src/search_pipeline.cpp`

The batch configuration's `disable_deduplication` value is now copied into the
per-spacer `SearchConfig` in both sequential and multithreaded CPU paths.

The multithreaded CPU path also now copies:

- `distance_mode`
- `search_start`
- `search_end`

**Why:** These values must survive a GPU-to-CPU fallback and explicit
`--cpu-only` execution. Without the propagation, CPU batch searches could use
default behavior instead of the CLI options selected for the benchmark.

## Explicit GPU hit-buffer exhaustion failure

### `src/gpu_engine.cu`

The two Myers GPU batch implementations now set a separate device overflow
flag and throw an error when their sparse hit buffer is exhausted instead of
clamping the hit count and returning the buffer's prefix:

- `myers_gpu_batch_threshold`
- `myers_gpu_batch_sparse`

**Why:** Silent truncation would produce incomplete candidate-position files
that look successful. The benchmark requires exhaustive results, so an
explicit failure is safer and makes the affected run identifiable.

The separate flag remains reliable even if the 32-bit hit counter itself wraps.
This does not increase the existing per-pattern capacity. It changes an
overflow from silent data loss into a visible error.

## Integration test

### `tests/integration/test_cli.cpp`

An integration test compares a normal batch search with the same search using
`--no-deduplicate` and verifies that the raw search emits more positions. The
test:

- Uses batch mode through `--spacer-file`.
- Forces the CPU path so it can run without a GPU.
- Disables mismatch annotation and scoring, matching the benchmark mode.
- Uses the process ID in its temporary filename to avoid collisions.

**Why:** This protects against regression of the original batch
`--no-deduplicate` bug.

## Direct AutoOffTarget candidate output

### `src/main.cpp`

The new `--autoofftarget-output-dir` batch option writes Stomata's in-memory,
per-spacer sorted hits directly as per-chromosome `endPosition:guideID` files.
`--autoofftarget-trim-trailing` removes synthetic trailing separator bases
from consideration.

Direct mode requires `--strand plus`; the benchmark invokes it once with a
forward index and once with a precomputed reverse-complement index. It creates
empty chromosome files when there are no hits, removes exact duplicate
endpoints without collapsing the nearby Myers halo, and reports direct file
writing time separately.

**Why:** The previous benchmark path formatted large in-memory BED text, wrote
it to disk, read it back with `awk`, and then wrote the final candidate files.
At high edit distances this conversion dominated total stage-1 time. Direct
output removes the intermediate BED file and text reparsing while preserving
the existing file boundary between stages 1 and 2.

### Benchmark runner

`run_stomata_stage1.sh` uses direct mode by default and records both the total
Stomata process time and `candidate_write_elapsed_seconds`. The derived
`search_without_candidate_write_seconds` shows the process time excluding the
direct candidate writer. Set `STOMATA_OUTPUT_MODE=bed` to retain the legacy
BED-plus-adapter path for comparison.

## Related runner design

The benchmark runner uses separate forward and reverse-complement genome
indexes. This is required because GPU sparse hits store the position relative
to a search slice in a 32-bit field. A combined forward-plus-reverse human
genome exceeds `2^32` bases and would wrap positions. Each separate hg38 and
CHM13 strand index is approximately 3.1 billion bases and remains below that
limit.

The separate-index workaround is implemented in `run_stomata_stage1.sh`.

## Validation

The modified GPU workflow was used to generate the hg38 and T2T-CHM13
AutoFFinder benchmark candidates. The following development checks were also
completed:

- Shell syntax validation for the benchmark runner.
- Synthetic forward/reverse FASTA preparation.
- Synthetic BED-to-candidate conversion.
- Synthetic two-pass benchmark orchestration.
- CSV column-count validation.
- Verification that each strand genome is below the 32-bit position limit.
- `git diff --check`.
