# AutoFFinder

AutoFFinder is a runtime-reconfigurable hardware-software pipeline for
multi-guide CRISPR/Cas9 off-target search:

1. **ReLev** scans a reference genome on an FPGA with parallel Levenshtein
   automata and reports candidate end positions.
   **Stomata** can alternatively generate compatible candidate files on an
   NVIDIA GPU with Myers' bit-parallel algorithm.
2. **PostAutoFFinder** validates those candidates on the CPU, reconstructs exact
   alignments, and enforces separate mismatch, bulge, PAM, and total-edit
   constraints.

The two stages are stored in this repository but are built and run separately.

## Repository layout

- `ReLev/`: FPGA overlay and host code for candidate generation. This directory
  is a Git submodule.
- `Stomata/`: Stomata v0.11.0 source with the modifications used to generate
  AutoFFinder-compatible stage-1 candidates, benchmark guide sets, and the
  stage-1 benchmark runner.
- `PostAutoFFinder/`: Java CPU post-processing implementation.
- `sgRNAs.txt`: example guide file.

## Clone

Clone the repository with its ReLev submodule:

```bash
git clone --recurse-submodules https://github.com/OrensteinLab/AutoFFinder.git
cd AutoFFinder
```

If the repository was cloned without submodules:

```bash
git submodule update --init --recursive
```

## Requirements

- Java 17 or later for PostAutoFFinder
- The ReLev prerequisites documented in
  [`ReLev/README.md`](ReLev/README.md) for FPGA candidate generation
- CUDA 11.8+, CMake 3.20+, and a C++17 compiler when using the optional Stomata
  stage 1; see [`Stomata/README.md`](Stomata/README.md)

No external Java libraries are required.

## How the pipeline works

ReLev reports candidate locations that satisfy a global Levenshtein-distance
threshold. A candidate identifies a guide and the end position of a possible
off-target alignment; it does not contain the reconstructed alignment or
separate mismatch and bulge counts.

PostAutoFFinder extracts a short genomic window ending at each candidate
position and performs semiglobal alignment reconstruction. It applies:

- a maximum total-edit threshold;
- a mismatch threshold for alignments without bulges;
- a mismatch threshold for alignments with bulges;
- a maximum bulge count;
- PAM matching, optionally allowing edits within the PAM; and
- optional best-alignment selection within a genomic window.

The general reconstruction algorithm builds a dynamic-programming matrix and
uses a memoized constrained traceback. Traceback states include the target and
text coordinates and the mismatch and bulge counts, avoiding repeated
evaluation of the same subproblem. Thread-local primitive workspaces reduce
allocation and garbage-collection overhead.

For configurations allowing at most one bulge, a fast path first evaluates
mismatch-only and single-bulge cases in linear time. It falls
back to the memoized algorithm whenever the shortcut cannot preserve the
general algorithm's result. This validated fast path is always used when
applicable.

## Build PostAutoFFinder

From the repository root:

```bash
mkdir -p bin
javac -d bin PostAutoFFinder/*.java
```

## PostAutoFFinder command

```bash
java -cp bin PostAutoFFinder.AutoOffTargetSearchAlign \
  <genome-fasta> \
  <guide-file-or-sequence> \
  <output-prefix> \
  <maxE> <maxM> <maxMB> <maxB> \
  <threads> \
  <best-in-window> <best-window-size> \
  <PAM> <allow-PAM-edits> \
  <candidate-source-path>
```

### Positional arguments

1. `genome-fasta`: reference genome in FASTA format.
2. `guide-file-or-sequence`: guide file, or one guide sequence supplied
   directly.
3. `output-prefix`: output path without the `.csv` suffix.
4. `maxE`: maximum total edits.
5. `maxM`: maximum mismatches for an alignment without bulges.
6. `maxMB`: maximum mismatches for an alignment containing bulges.
7. `maxB`: maximum number of bulges.
8. `threads`: number of CPU post-processing threads.
9. `best-in-window`: `true` to retain the best alignment in each locus window,
   otherwise `false`.
10. `best-window-size`: locus-window size used when `best-in-window` is true.
11. `PAM`: PAM suffix, for example `NGG`.
12. `allow-PAM-edits`: `true` to allow edits in the PAM, otherwise `false`.
13. `candidate-source-path`: path interpreted according to the selected
    candidate source.

Boolean arguments are case-sensitive and should be written as `true` or
`false`.

## Inputs

### Reference FASTA

Multi-FASTA input is supported. On the first run, PostAutoFFinder creates:

- `<genome-name>_split/`: one forward sequence file per FASTA record;
- `<genome-name>_split_rc/`: the corresponding reverse-complement files.

Existing split directories are reused on later runs. FASTA headers determine
the chromosome filenames, so candidate filenames must use the same names.

### Guide file

Provide one guide per line, including its PAM suffix. For example:

```text
GAGTCCGAGCAGAAGAAGAANGG
```

Guide IDs in candidate files are zero-based line indexes. DNA bases and `N` are
supported.

## Candidate input modes

Select the input mode with the Java system property
`autoffinder.candidateSource`. The default is `file`.

### Text candidate files (default)

Use candidate files produced by ReLev or by another compatible first stage:

```bash
java -cp bin PostAutoFFinder.AutoOffTargetSearchAlign \
  genome.fa sgRNAs.txt results/run \
  6 6 4 2 32 false 50 NGG false candidates
```

For a split chromosome file named `chr1.txt`, the candidate directory must
contain:

```text
candidates/chr1_fw.txt
candidates/chr1_rc.txt
```

Each line is:

```text
<end-position>:<guide-id>
```

For example:

```text
1048576:3
```

Missing chromosome/strand files are skipped. Candidate positions are relative
to the corresponding forward or reverse-complement split sequence.

### Stomata GPU candidate generation

The vendored Stomata source includes the changes used for the GPU stage-1
experiments:

- batch-mode propagation of `--no-deduplicate`;
- explicit failure instead of silent truncation when a Myers/Levenshtein GPU
  hit buffer fills;
- direct per-chromosome `endPosition:guideID` output for PostAutoFFinder; and
- removal of synthetic trailing separator bases from candidate output.

Build it from the repository root:

```bash
cmake -S Stomata -B Stomata/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=<gpu-architecture>
cmake --build Stomata/build -j"$(nproc)"
```

For one batch search, create separate Stomata indexes for the forward and
reverse-complement chromosome FASTAs, then run the binary once per index:

```bash
Stomata/build/src/stomata \
  --genome forward.fa.st \
  --spacer-file guides.tsv \
  --threshold 6 \
  --distance-mode levenshtein \
  --strand plus \
  --no-deduplicate \
  --no-compute-mismatches \
  --no-scores \
  --no-treat-u-as-t \
  --autoofftarget-output-dir candidates \
  --autoofftarget-trim-trailing 7
```

The full benchmark runner is
[`Stomata/run_stomata_stage1.sh`](Stomata/run_stomata_stage1.sh). It prepares
the forward and reverse-complement indexes, runs the edit-distance and
guide-count experiments, records timing and memory measurements, and writes
candidate directories consumable by PostAutoFFinder. The exact guide subsets
are in [`Stomata/benchmark_guides/`](Stomata/benchmark_guides/).

Configure its input chromosome directories through environment variables:

```bash
HG38_FORWARD_DIR=/data/hg38_split \
HG38_REVERSE_DIR=/data/hg38_split_rc \
STOMATA_OUTPUT_PARENT=/data/results \
Stomata/run_stomata_stage1.sh hg38
```

For CHM13, use `CHM13_FORWARD_DIR` and `CHM13_REVERSE_DIR`. Pass `all` to run
both assemblies. See [`Stomata/LOCAL_CHANGES.md`](Stomata/LOCAL_CHANGES.md)
for the implementation rationale and upstream provenance.

### Raw ReLev binary captures

Binary mode decodes the raw 8-byte records emitted by ReLev without converting
them to intermediate text files:

```bash
java \
  -Dautoffinder.candidateSource=binary \
  -Drelev.binary.forwardDir=/data/relev/forward \
  -Drelev.binary.reverseDir=/data/relev/reverse \
  -Drelev.binary.editDistance=6 \
  -cp bin PostAutoFFinder.AutoOffTargetSearchAlign \
  genome.fa sgRNAs.txt results/run \
  6 6 4 2 32 false 50 NGG false /data/relev
```

The forward and reverse directories must contain files named:

```text
<split-sequence-filename>_ed<distance>.out
```

For example, `chr1.txt_ed6.out` corresponds to `chr1.txt`.

If the directory properties are omitted, binary mode interprets argument 13 as
a root directory and defaults to:

```text
<root>/hg38_only_chrs_split/
<root>/hg38_only_chrs_split_rc/
```

The explicit directory properties are recommended for other genomes or naming
schemes.

Each raw record is little-endian and contains:

| Bytes | Field |
|---:|---|
| 4 | signed candidate end position |
| 2 | unsigned 16-bit match mask |
| 1 | unsigned automata group ID |
| 1 | valid flag |

Each group contains 16 guide lanes. A record with `valid = 0` terminates the
stream. The decoder rejects malformed records, missing terminators, guide IDs
outside the supplied guide count, and decreasing positions for the same guide.
Exact duplicate positions for a guide are removed.

Binary captures support 1--128 guides. Their active guide lanes must correspond
to the first lines of the supplied guide file.

When binary mode is used, PostAutoFFinder prints separate timings for:

- mapping and loading the raw buffers;
- Java binary decoding;
- candidate alignment post-processing; and
- decoding plus post-processing.

### Live FPGA adapter (advanced)

The `fpga` candidate source invokes ReLev through the `relev_jni` native
interface:

```bash
export JAVA_HOME=/path/to/jdk
make -C ReLev/fpga jni
```

This creates `ReLev/fpga/librelev_jni.so`. Then run PostAutoFFinder with the
JNI library and the ReLev FPGA image:

```bash
java \
  -Dautoffinder.candidateSource=fpga \
  -Drelev.xclbin=/absolute/path/to/ReLev/fpga/automata.hw.xclbin \
  -Drelev.nativeLibrary=/absolute/path/to/ReLev/fpga/librelev_jni.so \
  -cp bin PostAutoFFinder.AutoOffTargetSearchAlign \
  genome.fa sgRNAs.txt results/run \
  6 6 4 2 32 false 50 NGG false unused
```

This mode requires exactly 128 guides and an edit-distance threshold from 0 to
6. ReLev targets the AMD Alveo U280 by default and requires the XRT and Vitis
environment described in [`ReLev/README.md`](ReLev/README.md). Other FPGA
platforms require an appropriate platform path and memory mapping. Omit
`relev.nativeLibrary` only when `librelev_jni` is already available through
`java.library.path`.

## Output

PostAutoFFinder writes `<output-prefix>.csv` with these columns:

- `Chromosome`
- `Strand`
- `EndPosition`
- `Target`
- `SiteSeqPlusMaxEditsBefore`
- `#Edit`
- `AlignedTarget`
- `AlignedText`
- `#Mismatches`
- `#Bulges`

## Decoder checks

Compile the source and run the self-contained raw-decoder checks:

```bash
java -cp bin PostAutoFFinder.RelevBinaryDecoderTest
```

The checks cover interleaved automata groups, duplicate removal, equivalence to
the text parser, decreasing-position rejection, and missing-terminator
rejection.

To inspect one real binary capture:

```bash
java -cp bin PostAutoFFinder.RawBinaryCandidateSourceTest \
  <raw-root> <chromosome-file-name> <+|-> <edit-distance>
```

`RawBinaryCandidateSourceTest` expects the default hg38 subdirectory names
described above and decodes 128 guide lanes.

To summarize all captures for one edit distance:

```bash
java -cp bin PostAutoFFinder.RawBinaryBatchBenchmark <raw-root> <edit-distance>
```

## Notes

- PostAutoFFinder currently writes one combined CSV and processes chromosome
  and strand files sequentially while parallelizing candidates within each
  file.
- Best-in-window mode processes one chromosome/strand candidate stream in one
  worker so that a locus window cannot be split across independent workers.
- The optional live-FPGA adapter allocates an in-memory direct output buffer;
  its native JNI implementation must return ReLev's raw record format.
- The text and raw-binary modes are fully implemented in Java.
