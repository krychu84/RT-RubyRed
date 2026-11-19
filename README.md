# RT-RubyRed

**Real-time rRNA amplicon classification for Oxford Nanopore sequencing runs.**

RT-RubyRed runs the [RubyRed](https://github.com/reillyeo/RubyRed) pipeline *while the sequencer is still running*. It watches the folder MinKNOW writes `*.fastq.gz` files into, analyses each new batch of files as it arrives, and keeps a merged feature table and an interactive HTML report up to date for the whole run (typically 70–110 hours).

```
██████╗ ████████╗    ██████╗ ██╗   ██╗██████╗ ██╗   ██╗██████╗ ███████╗██████╗
██╔══██╗╚══██╔══╝    ██╔══██╗██║   ██║██╔══██╗╚██╗ ██╔╝██╔══██╗██╔════╝██╔══██╗
██████╔╝   ██║ █████╗██████╔╝██║   ██║██████╔╝ ╚████╔╝ ██████╔╝█████╗  ██║  ██║
██╔══██╗   ██║ ╚════╝██╔══██╗██║   ██║██╔══██╗  ╚██╔╝  ██╔══██╗██╔══╝  ██║  ██║
██║  ██║   ██║       ██║  ██║╚██████╔╝██████╔╝   ██║   ██║  ██║███████╗██████╔╝
╚═╝  ╚═╝   ╚═╝       ╚═╝  ╚═╝ ╚═════╝ ╚═════╝    ╚═╝   ╚═╝  ╚═╝╚══════╝╚═════╝
```

![RT-RubyRed report](docs/report.png)

## Contents

- [Features](#features)
- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [The HTML report](#the-html-report)
- [Species-level classifier (database)](#species-level-classifier-database)
- [Resource limits](#resource-limits)
- [Testing without a sequencer](#testing-without-a-sequencer)
- [Citation, license and author](#citation-license-and-author)

## Features

- **Watches a MinKNOW output folder** and picks up new `*.fastq.gz` files as they appear. The folder is only read, never written: MinKNOW folders are owned by `minknow:minknow`.
- **Batch processing:** new files → [torchlex](https://nanoporetech.com/resource-centre/torchlex-method-real-time-demultiplexing-barcoded-oxford-nanopore-reads) demultiplexing on the inner barcodes (`BRK01`–`BRK192`) → RubyRed classification → one feature table per batch.
- **Merged feature table** of all batches, updated after every batch.
- **Interactive HTML report** in a single offline file:
  - raw vs classified reads for the whole run in the page header, with the reads removed at each RubyRed step
  - two 96-well barcode plates (BRK01–BRK192) coloured by the read threshold, and raw reads per barcode
  - species-level relative abundance
  - rarefaction curves
  - live sliders for the abundance and read-count filters
  - copy-ready Materials and methods (wet-lab protocol and a data-analysis text with this run's tool versions and settings)
  - downloads of all feature tables (ZIP) and printing on A4
- **Per-taxon read files** (`rep_seqs/<taxon>.fasta.gz`) with read IDs that are unique across batches.
- **Uses RubyRed's own reference database** (the curated species-level MIMt 16S classifier distributed with RubyRed), so RT-RubyRed and plain RubyRed give the same classification. See [below](#species-level-classifier-database).
- **Safe to run next to MinKNOW:** it uses at most 50% of CPU cores (hard pinning) and 50% of RAM (cgroup limit), at low CPU and I/O priority.
- **Stop at any time with Ctrl-C:** the per-batch `BRKxx.fastq` files are joined into one file per barcode and the batch folders are removed. `--resume` continues an interrupted run.

## How it works

```
MinKNOW  ──>  watch folder (*.fastq.gz, read-only)
                    │  every 60 s: files not yet in input_file_index.tsv
                    ▼
            torchlex demultiplexing ──> fastq_demultiplexed/batchNN/BRKxx.fastq
                    ▼
            RubyRed (QIIME 2): length/quality filter, chimera removal,
            orientation, sklearn classification (species-level classifier)
                    ▼
            feature_tables/feature_table_batchNN.tsv
            rep_seqs/<taxon>.fasta.gz  (reads appended, IDs BRKxx_bNN_n)
                    ▼   (from batch 2 on)
            <run>_<dir>_<date>_feature_table.tsv   (all batches merged)
            <run>_<dir>_<date>_report.html         (rebuilt every batch)
```

When a batch finishes, the tool checks the folder again and processes only the files that arrived in the meantime. A file is added to the index only after its batch has completed. If a batch is interrupted, its files are therefore picked up again later, and no file is analysed twice.

Example terminal output:

```
  Run          ONT_250925_16S
  Watching     /data/ONT_250925_16S/.../fastq_pass/barcode94  (owned by minknow:minknow; read-only use)
  Results      /home/user/Desktop/RT-RubyRed_ONT_250925_16S_barcode94_20260925_0737
  Resources    16/32 cores (pinned 16-31), RAM cap 62.5G of 124.9G (cgroup enforced), nice 10 / ionice idle-ish
  Classifier   MIMt_16s_311025_species_classifier.qza
[07:37:17] ▶ batch01: 10 new file(s) (84.0M) → torchlex demultiplexing
[07:37:22]    torchlex: 86,575 reads assigned to 37 BRK barcode(s)
[07:37:22]    RubyRed classification (16 threads)…
[07:43:36] ✔ batch01 done in 6.3 min: 30,821 reads added to the feature table (run total 30,821)
             input : /data/ONT_250925_16S/.../fastq_pass/barcode94
             output: /home/user/Desktop/RT-RubyRed_ONT_250925_16S_barcode94_20260925_0737
             table : feature_tables/feature_table_batch01.tsv
```

## Requirements

- Linux (tested on Ubuntu 22.04).
- [RubyRed](https://github.com/reillyeo/RubyRed) with its `qiime2` conda environment and its reference database (`install.sh` sets both up). That environment must contain `qiime` (tested with QIIME 2 2024.5 and 2024.10, scikit-learn 1.4.2), `torchlex`, `cutadapt`, `chopper`, `vsearch`, `seqkit`, GNU `parallel` and `biom`, plus Python packages `pandas` and `psutil`.
- `systemd-run` (standard on Ubuntu) to enforce the RAM limit. Without it the tool still runs, but the RAM limit is not enforced.

## Installation

RT-RubyRed runs on top of the original **RubyRed** pipeline ([reillyeo/RubyRed](https://github.com/reillyeo/RubyRed)) and uses RubyRed's own reference database. The installer clones RubyRed from its repository, fetches the database (stored with [Git LFS](https://git-lfs.com)), and puts RT-RubyRed on your `PATH`:

```bash
mkdir -p ~/my_scripts && cd ~/my_scripts
git clone https://github.com/krychu84/RT-RubyRed.git
cd RT-RubyRed && ./install.sh
source ~/.bashrc
RT-RubyRed -v
```

What `install.sh` does:

1. **Conda environment `qiime2`**: used as is if it exists. Otherwise it is created exactly as in the RubyRed README (QIIME 2 2024.5 amplicon distribution, plus `chopper`, `seqkit` and `parallel`).
2. **RubyRed**: cloned from `https://github.com/reillyeo/RubyRed` into `~/my_scripts/RubyRed` (or `$RUBYRED_HOME`). An existing RubyRed is kept unchanged; `./install.sh --update-rubyred` pulls the latest version.
3. **RubyRed database**: `git-lfs` is installed into the environment if needed, then `git lfs pull` downloads the real `.qza` files. The installer stops if any file is still a Git LFS pointer (a ~130-byte text stub).
4. **torchlex** (from RubyRed's `scripts/torchlex`), `pandas` and `psutil`, and a check that every required command is present.
5. **PATH**: adds RT-RubyRed to `~/.bashrc` (and `RUBYRED_HOME` if RubyRed is not in the default location).

The launcher activates the `qiime2` environment itself. It finds conda in `~/miniforge3`, `~/miniconda3` or `~/anaconda3`; for any other location, set `CONDA_BASE=/path/to/conda`.

> **Database check at start-up.** RT-RubyRed reads the database paths from the `RubyRed` script and refuses to start if a file is missing or is a Git LFS pointer. Fix it with `cd ~/my_scripts/RubyRed && git lfs pull` (or run `./install.sh` again).

## Usage

Run it inside `tmux` or `screen`, so it keeps running if your SSH connection drops:

```bash
tmux new -s rt
RT-RubyRed -w /data/<experiment>/<sample>/<flowcell_run>/fastq_pass/barcodeNN
```

- Detach: `Ctrl-b` then `d`. Re-attach: `tmux attach -t rt`.
- **Stop and finalise:** `Ctrl-C`.

Examples:

```bash
# set the run name yourself (default: MinKNOW experiment name from the path)
RT-RubyRed -w .../fastq_pass/barcode94 -n GUT_STUDY_01

# finalise automatically 3 h after the last new file appeared
RT-RubyRed -w .../fastq_pass/barcode94 --idle-stop 3

# continue a run that was stopped or crashed (same results folder)
RT-RubyRed -w .../fastq_pass/barcode94 --resume ~/Desktop/RT-RubyRed_<...>

# your own RubyRed filters, or another classifier folder
RT-RubyRed -w .../fastq_pass/barcode94 --rubyred-args "-q 12 -l 900 -x 1500"
RT-RubyRed -w .../fastq_pass/barcode94 --classifier /path/to/MIMt_16s_311025_species
```

| Option | Default | Description |
|---|---|---|
| `-w, --watch DIR` | *(required)* | Folder with the incoming `*.fastq.gz` files |
| `-n, --name NAME` | from path | Run name |
| `-o, --outdir DIR` | `~/Desktop` | Where the results folder is created |
| `--resume DIR` | | Continue an earlier RT-RubyRed results folder |
| `--poll SEC` | 60 | Seconds between checks for new files |
| `--settle SEC` | 60 | Skip files modified in the last SEC seconds (still being written) |
| `--min-files N` | 1 | Minimum number of new files to start a batch |
| `--idle-stop HOURS` | 0 (never) | Finalise after HOURS without new files |
| `-R, --recursive` | off | Also look in sub-folders of the watch folder |
| `--cpu-frac F` | 0.5 | Share of CPU cores to use |
| `--mem-frac F` | 0.5 | Share of RAM to allow |
| `--classifier DIR` | RubyRed's database | Classifier folder `<dir>/<name>_{classifier,seqs,tax}.qza`; by default RubyRed's own database is used |
| `--subsample N` | no subsampling | RubyRed `-t` per batch |
| `--rubyred PATH` | `$RUBYRED_HOME/RubyRed` | RubyRed script |
| `--rubyred-args "..."` | | Extra RubyRed options (`-q`, `-l`, `-x`, `-c`, `-z` …) |
| `--barcodes FASTA` | RubyRed's | torchlex barcode file |

A batch takes about 2–6 minutes on a 32-core machine. Most of that time is spent in RubyRed's sklearn classification step.

## Output

```
~/Desktop/RT-RubyRed_<run>_<dir>_<YYYYMMDD_HHMM>/
├── <run>_<dir>_<date>_feature_table.tsv   merged feature table (all batches)
├── <run>_<dir>_<date>_report.html         interactive report (single offline file)
├── input_file_index.tsv                   every input file analysed, with its batch
├── batch_stats.json                       raw reads per BRK barcode and per batch
├── RT-RubyRed.log                         what was shown in the terminal
├── run_meta.json                          run settings
├── fastq_demultiplexed/
│   ├── batch01/ batch02/ …                during the run
│   └── BRK01.fastq BRK02.fastq …          after Ctrl-C (batches joined)
├── feature_tables/feature_table_batchNN.tsv
├── rep_seqs/<taxon>.fasta.gz              classified reads per taxon, all batches
├── rubyred_batches/batchNN/               RubyRed output per batch (log, rep_seqs, qzv)
└── logs/batchNN.log, batchNN_files.txt    full torchlex/RubyRed output, input list
```

- `<run>` is the MinKNOW experiment name (or `-n`); `<dir>` is the watched folder (e.g. `barcode94`).
- RubyRed numbers reads from 1 in every batch. In `rep_seqs/`, read IDs therefore carry the batch number (`>BRK185_b02_89`), so they are unique across batches. The part before the first `_` is the sample.
- Low-frequency taxa are not removed per batch (RubyRed `-f 1`). Use the abundance slider in the report instead.

## The HTML report

A single self-contained HTML file (plotly.js is embedded, so it works offline). It is rebuilt after every batch; re-open it to see the latest data.

- **Summary tiles:** batches, input files, samples, raw and classified reads, taxa.
- **Minimum relative abundance** (slider plus number box): taxa below the threshold in every shown sample are grouped as *Other* in the barchart and left out of the rarefaction curves.
- **Minimum classified reads per sample** (slider plus number box): hides samples with too few reads.
- **Raw reads per BRK barcode:**
  - *Barcode plates:* two 96-well plates, BRK01–BRK96 on plate 1 and BRK97–BRK192 on plate 2, filled column-wise (A1 = BRK01, B1 = BRK02 … H12 = BRK96; A1 = BRK97 … H12 = BRK192). Each well shows the BRK number and its classified reads. Wells turn **green** at or above the *minimum classified reads* value and **red** below it; grey wells have no reads. Hover a well for raw/classified reads and its top taxon.
  - *Reads per barcode:* stacked by batch, with classified reads as markers.
- **Relative abundance at species level:** taxa not resolved to species are shown at their lowest classified rank, e.g. *Staphylococcus sp.* or *Enterobacteriaceae (family)*.
- **Rarefaction curves per sample:** expected richness computed analytically (Hurlbert), without random subsampling.

All plots update when either slider moves, and every plot can be downloaded as SVG.

- **Downloads** (built in the browser, no server needed): *Download all feature tables (ZIP)* holds the merged feature table (counts), relative abundance (%), a table filtered with the current thresholds, a sample/well summary (plate, well, raw and classified reads, pass/below, top taxon) and every per-batch table (`feature_tables/feature_table_batchNN.tsv`). Each table can also be downloaded on its own.
- **Print panel / Print full report (A4):** each panel has its own *Print panel* button that prints only that panel (with the header and the current thresholds); *Print full report* sits at the bottom with the downloads. Both print A4 portrait, with both plates side by side and the species legend printed under the barchart.

## Species-level classifier (database)

By default RT-RubyRed uses **the database that ships with RubyRed** (`RubyRed/resources/classifiers/`), exactly as configured in the `RubyRed` script. With RubyRed of October 2026 this is `MIMt_16s_311025_species_amp341F1390R`:

- **Reference:** MIMt 16S rRNA database (Cabezas et al., 2024; release M2c_25_10).
- **Curated to species level**, then **retrained on the 341F–1390R amplicon** that RubyRed extracts with cutadapt (QIIME 2 2024.5, scikit-learn 1.4.2). The training command is in the classifier folder's `README.txt` in RubyRed.

RT-RubyRed passes `-w <RubyRed>/resources` to RubyRed, and the HTML report names the classifier that was used. To use another classifier folder (`<dir>/<name>_{classifier,seqs,tax}.qza`), add `--classifier DIR`.

### Species-level curation

Most species labels in the MIMt 16S reference database carry a strain name (e.g. `S__Escherichia_coli_strain_JCM_1649`), and some taxa appear under two names or two lineages. QIIME's `classify-sklearn` works down the ranks and stops at the first rank where no single name reaches the confidence threshold (0.9). Both problems split the probability and make reads stop too early. For example, *Akkermansia* reads stopped at **phylum** because the database had the class both as *Verrucomicrobiae* and as *Verrucomicrobiia*.

`database/make_species_classifier.py` corrected the taxonomy that RubyRed's classifier is now built on:

1. **Lineage clean-up at all ranks:**
   - synonyms and spelling variants are unified (e.g. *Euryarchaeota* → *Methanobacteriota*, *Verrucomicrobiae* → *Verrucomicrobiia*)
   - species names that use a misspelled or older genus name are corrected
   - every taxon is given a single parent lineage
2. **Species instead of strains:**
   - strain names, culture-collection numbers, *subsp.* and isolate codes are removed
   - where a species has several references, the strain moves to an 8th rank (`…;S__Escherichia_coli;T__JCM_1649`), so the classifier adds up the probability of all its strains
   - RubyRed collapses the table at level 7, so strains never appear in the results
3. **Species 16S cannot separate:**
   - *Akkermansia massiliensis* and *A. biwaensis* are reported as *A. muciniphila*
   - *Escherichia* and *Shigella* become one genus, *Escherichia-Shigella* (as in SILVA); *Shigella* species are reported as *E. coli*

Effect on 49,775 reads from a 16S run (original → cleaned full-length classifier, same reads):

| Classified to | original | cleaned |
|---|---|---|
| species | 8.3 % | **41.5 %** |
| genus | 58.8 % | 44.7 % |
| above genus | 32.9 % | **13.8 %** |
| reads stopping at phylum | 9,544 | **1** |

No read became less specific. Reads that remain at a high rank (e.g. class *Bacilli*) were mostly PCR chimeras: one half of the read from one genus, the other half from another.

Details, the full list of changes and how to add your own species merges: [database/README.md](database/README.md). The full-length (untrimmed) cleaned classifier from RT-RubyRed v1.0.0 is still attached to the [v1.0.0 release](https://github.com/krychu84/RT-RubyRed/releases/tag/v1.0.0) for use with `--classifier`.

## Resource limits

RT-RubyRed is meant to run on the sequencing computer during the run, so it limits itself:

- **CPU:** the process and all its children are pinned to half of the cores (`--cpu-frac`), which is a hard limit.
- **RAM:** each torchlex/RubyRed step runs in a `systemd-run --user --scope` with `MemoryMax` = 50% of RAM (`--mem-frac`) and no swap.
- **Priority:** `nice 10` and `ionice -c2 -n7`.
- **Monitoring:** CPU and RAM use are shown in the terminal every 30 s while a batch runs.

## Testing without a sequencer

`test/simulate_run.sh` copies fastq.gz files into a folder in batches, the way MinKNOW writes them. Run it in one terminal and RT-RubyRed in another:

```bash
# terminal 1: 10 files every 2 minutes
test/simulate_run.sh /path/to/existing/fastq_pass/barcodeNN /tmp/rt_watch 10 120
# terminal 2
RT-RubyRed -w /tmp/rt_watch -n TEST -o /tmp/rt_out --poll 20 --settle 5
```

## Citation, license and author

If you use RT-RubyRed, please cite **RubyRed**, on which it is built, and the protocol paper:

- **RubyRed:** Reilly, E. T., & Krych, Ł. RubyRed: a pipeline for read-by-read analysis of ONT-generated amplicon sequencing data. https://github.com/reillyeo/RubyRed
- **Protocol and RubyRed paper:** McPhillips, C. H., Reilly, E. T., Stolberg-Mathieu, G., Nielsen, K., Gottlieb, A. D., Madjarov, G., Roager, H. M., Nielsen, D. S., & Krych, L. (2026). Ultra-high multiplexing enables near-full-length 16S rRNA gene amplicon sequencing of over 1,200 gut microbiome samples on a single nanopore flow cell. *bioRxiv* (preprint). https://doi.org/10.64898/2026.08.29.747698 *(will be updated when the paper is published)*
- **RT-RubyRed:** Krych, L. (2026). RT-RubyRed: real-time rRNA amplicon classification for Oxford Nanopore sequencing runs. https://github.com/krychu84/RT-RubyRed
- Please also cite the tools RubyRed uses: QIIME 2 (Bolyen et al., 2019), q2-feature-classifier (Bokulich et al., 2018), RESCRIPt (Robeson et al., 2021), the MIMt database (Cabezas et al., 2024), torchlex (Schack et al., 2024), cutadapt (Martin, 2011), chopper (De Coster & Rademakers, 2023), VSEARCH (Rognes et al., 2016) and SeqKit (Shen et al., 2024). The report's *Materials and methods* panel lists them with full references and the versions used in your run.

`CITATION.cff` contains the same references in machine-readable form (GitHub shows it as "Cite this repository").

- **Author:** Lukasz Krych, University of Copenhagen, Department of Food Science.
- **License:** RT-RubyRed is released under the [MIT License](LICENSE).
- **Built on:**
  - [RubyRed](https://github.com/reillyeo/RubyRed) (Eoghan T. Reilly and Łukasz Krych)
  - [QIIME 2](https://qiime2.org)
  - [torchlex](https://nanoporetech.com/resource-centre/torchlex-method-real-time-demultiplexing-barcoded-oxford-nanopore-reads)
  - [scikit-learn](https://scikit-learn.org)
  - [plotly.js](https://plotly.com/javascript/) (MIT, bundled in `assets/`)
- The University of Copenhagen seal in the report belongs to the University of Copenhagen.
