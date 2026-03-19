# AutoFFinder

AutoFFinder is a runtime-reconfigurable hardware-software co-design for multi-gRNA CRISPR/Cas9 off-target search:

- **ReLev** (FPGA stage): high-throughput Levenshtein automata candidate generation.
- **PostAutoFFinder** (CPU stage): exact alignment reconstruction and biological filtering.

This repository stores the project code and is organized as a two-stage pipeline. It is **not** a single one-command tool.

## Repository structure

- `ReLev/` – FPGA code and host code for candidate generation (submodule).
- `PostAutoFFinder/` – Java post-processing code that consumes ReLev outputs and produces final off-target CSV results.

## Stage 1: Run ReLev separately

Use the ReLev documentation for build and execution details:

- [ReLev README](ReLev/README.md)

If cloning from scratch, make sure the submodule is available:

```bash
git submodule update --init --recursive
```

## Stage 2: Run PostAutoFFinder

### Requirements

- Java 17 (tested with OpenJDK 17)

### Compile

```bash
javac -d bin PostAutoFFinder/*.java
```

### Run

```bash
java -cp bin PostAutoFFinder.AutoOffTargetSearchAlign \
	<Genome reference FASTA> \
	<gRNA file path> \
	<Output prefix> \
	<maxE> <maxM> <maxMB> <maxB> \
	<Threads> \
	<Best-in-window> <Best-window-size> \
	<PAM> <Allow PAM edits> \
	<ReLev output folder>
```

> `AutoOffTargetSearchAlign` expects **13 positional arguments** in exactly this order.

## Input expectations and formats

### 1) Genome reference FASTA

- Multi-FASTA is supported.
- Each record is split into per-chromosome text files internally by `PostAutoFFinder`.

### 2) gRNA file (`<gRNA file path>`)

- Plain text file, one guide per line.
- Each line should include the guide **with PAM suffix** (for example, `NNNNNNNNNNNNNNNNNNNNNGG`).
- Allowed characters are DNA bases and `N` (see [sgRNAs.txt](sgRNAs.txt) for an example).

### 3) ReLev output folder (`<ReLev output folder>`)

This folder must contain ReLev candidate matches per chromosome and strand in text files named:

- `<chromosome>_fw.txt` for forward strand
- `<chromosome>_rc.txt` for reverse-complement strand

where `<chromosome>` must match the chromosome filenames generated from the FASTA headers.

Each file must contain one candidate match per line in this format:

```text
<end_position>:<target_id>
```

- `<end_position>`: integer genomic end position reported by ReLev.
- `<target_id>`: zero-based index of the gRNA in the input gRNA file.

## Arguments

1. **Genome reference FASTA**: Path to the FASTA genome.
2. **gRNA file path**: Path to text file containing guides (one per line, including PAM).
3. **Output prefix**: Output prefix. Final file is written as `<Output prefix>.csv`.
4. **maxE**: Maximum total edits (integer).
5. **maxM**: Maximum mismatches when no bulges are used (integer).
6. **maxMB**: Maximum mismatches when bulges are used (integer).
7. **maxB**: Maximum bulges (integer).
8. **Threads**: Number of CPU threads for post-processing (integer).
9. **Best-in-window**: `true` or `false`.
10. **Best-window-size**: Window size used when argument 9 is `true` (integer).
11. **PAM**: PAM sequence, e.g. `NGG`.
12. **Allow PAM edits**: `true` or `false`.
13. **ReLev output folder**: Folder containing ReLev candidate files (`*_fw.txt`, `*_rc.txt`).

## Output format

PostAutoFFinder writes one CSV file (`<Output prefix>.csv`) with columns:

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