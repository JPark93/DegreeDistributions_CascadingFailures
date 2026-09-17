# Network Cascade Simulation (Self-Contained)

A single-file R script that simulates how **graph family** (Erdős–Rényi, Barabási–Albert, Watts–Strogatz), **network density** (20 vs. 40 edges on 20 vertices), **centrality-based seed selection** (degree, strength, harmonic closeness, betweenness × High/Average/Low targets), and **activation thresholds** (four regimes) affect final cascade size and conditional duration.

The script is fully self-contained: graph generation, edge weighting, target selection, simulation, analysis, and plotting all live in one file. No project data or other R scripts are required. The full design is 2 densities × 3 families × 4 metrics × 3 seed classes × 4 threshold regimes = **288 cells × 1,000 replications = 288,000 cascades**.

> The complete user guide and data dictionary are embedded as comments at the top of the script. Run `Rscript Cascade_Simulation_Self_Contained.R --help` to print them without running the simulation.

---

## Requirements

- **R** (exact numerical replay was validated with R 4.4.2; other versions run but may change results slightly and will trigger a warning)
- **igraph** (validated with 2.1.1) and **ggplot2** — installed with:
  ```r
  install.packages(c("igraph", "ggplot2"), repos = "https://cloud.r-project.org")
  ```
- `parallel` ships with R — nothing extra to install.
- Internet access is only needed to install packages; the simulation itself runs locally.

---

## Quick start

1. Save the script (e.g., `Cascade_Simulation_Self_Contained.R`) in a writable folder.
2. Install the two packages above (once).
3. Run the **whole file** from a terminal:
  ```sh
  Rscript Cascade_Simulation_Self_Contained.R
  ```
   In RStudio, use the **Source** button (do not run selected lines — the script locates itself when sourced). You can also use `source("path/to/Cascade_Simulation_Self_Contained.R")`.
4. Wait for the final success message, then open the reported output folder. Start with:
   - `cascade_summary_all_conditions.csv` — the broadest results overview (288 rows for a full run)
   - `table2_change_in_failure_rate.csv` and `table3_steps_to_complete_failure.csv` — the two manuscript analyses
   - `figure04`–`figure09` PNG files — the figures

A full run typically takes several minutes, depending on hardware and worker count.

---

## Changing the settings

There are two equivalent ways to configure a run. **Command-line options take precedence** over the settings in the script.

### Option A: Edit the USER CONFIGURATION block (top of the script)

```r
STUDY_DESIGN <- 1L                    # 1: WS beta=.10; 2: WS beta=.05
THRESHOLD_CONDITION <- 0L             # 0: all; 1..4: one threshold regime
REPLICATIONS_OVERRIDE <- NA_integer_  # NA: 1,000; e.g., 10L for a quick check
N_WORKERS <- 3L                       # Automatically capped by cores/replications
OUTPUT_DIRECTORY_OVERRIDE <- NULL     # NULL: fresh folder beside this script
MAKE_FIGURES <- TRUE                  # Figures 5-9
MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE <- TRUE  # Figure 4
```

### Option B: Command-line options

| Option | Values | What it does |
|---|---|---|
| `--help` | — | Print the embedded user guide/data dictionary; does not run the simulation |
| `--smoke` | — | Run 10 replications per selected threshold regime (execution check only) |
| `--study=` | `1` or `2` | `1` = primary study (WS rewiring 0.10); `2` = sensitivity study (WS rewiring 0.05) |
| `--threshold=` | `0`–`4` | `0` = all four regimes; `1`–`4` = one regime (see table below) |
| `--replications=` | positive integer | Number of replications per regime (default 1,000) |
| `--workers=` | positive integer | Requested parallel worker processes |
| `--output=` | path | Output destination; relative paths resolve against the script's folder |
| `--no-figures` | — | Disable both figure-generation options (data only) |

### What each setting controls

| Setting | Effect |
|---|---|
| `STUDY_DESIGN` / `--study` | Chooses the Watts–Strogatz rewiring probability (0.10 primary vs. 0.05 sensitivity). Everything else (factors, base seed, 1,000 replications) is the same. |
| `THRESHOLD_CONDITION` / `--threshold` | Which activation-threshold regimes to run. Each vertex gets a threshold drawn Uniform on the regime's interval. |
| `REPLICATIONS_OVERRIDE` / `--replications` / `--smoke` | Replications per regime. `NA` = 1,000 (manuscript setting). Small values (e.g., 10) are for execution checks only, not manuscript data. |
| `N_WORKERS` / `--workers` | Parallel processes. Automatically capped at the replication count and one fewer than detected logical CPUs (minimum 1). Worker count does not change results. |
| `OUTPUT_DIRECTORY_OVERRIDE` / `--output` | `NULL` creates a fresh timestamped folder beside the script (e.g., `Cascade_Results_Study1_20250101_120000`). An absolute path is used as-is; a relative path resolves against the script's directory. **The destination must be new or empty.** |
| `MAKE_FIGURES` | Generates Figures 5–9 (final-size panels per metric + High−Low contrast figure). |
| `MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE` | Generates Figure 4 (six representative graphs) plus `representative_topology_selection.csv`. |

### Threshold regime reference

| Value | Regime | Vertex threshold draw |
|---|---|---|
| `0` | All four regimes | — |
| `1` | Sensitive | Uniform(0.10, 0.20) |
| `2` | Moderate | Uniform(0.20, 0.40) |
| `3` | Elevated | Uniform(0.40, 0.60) |
| `4` | Resilient | Uniform(0.60, 0.80) |

---

## Common workflows

```sh
# Fast execution check (~10 replications, one worker, no figures)
Rscript Cascade_Simulation_Self_Contained.R --smoke --workers=1 --no-figures

# Full primary study (all regimes, 1,000 replications, fresh output folder)
Rscript Cascade_Simulation_Self_Contained.R

# Sensitivity study (WS rewiring 0.05)
Rscript Cascade_Simulation_Self_Contained.R --study=2

# One threshold regime only, custom destination
Rscript Cascade_Simulation_Self_Contained.R --threshold=1 --output=Sensitive
```

Notes:
- Quote the entire `--output=PATH` argument if the path contains spaces.
- Each run starts from replication 1; there is **no resume/checkpoint** facility.
- Choose a fresh output folder for every run — the script refuses to write into a non-empty directory to protect previous results.

---

## What you get

All results are written to the output folder. CSVs have headers, no row-number column, and use `NA` for undefined values. Percentages are 0–100; final-active proportions are 0–1. Column-level details for every file are in the script header (or `--help`).

| File | Contents |
|---|---|
| `cascade_summary_all_conditions.csv` | **Start here.** One row per density/family/metric/target/threshold (288 rows full run): means, SDs, Monte Carlo SEs, 95% intervals, completion counts. |
| `table2_change_in_failure_rate.csv` | Manuscript Table 2: paired High−Low change in failure rate per replication, summarized (96 rows). |
| `table3_steps_to_complete_failure.csv` | Manuscript Table 3: conditional steps-to-complete-failure for Low and High targets (192 rows). |
| `cascade_replication_level.csv` | Raw per-cascade results (288,000 rows full run) — the auditable source for every summary. |
| `high_low_difference_replication_level.csv` | Raw paired High−Low differences (96,000 rows). |
| `density_change_in_high_low_difference.csv`, `target_position_failure_rate_gaps.csv`, `conditional_duration_density_comparison.csv` | Derived density and target-position comparisons. |
| `topology_diagnostics_*.csv`, `target_selection_*.csv`, `target_agreement_*.csv`, `edge_weight_diagnostics*.csv` | Replication-level and summary diagnostics for graphs, target choices, metric agreement, and edge weights. |
| `manuscript_simulation_results.rds` | Bundle of all summary tables; read with `readRDS()`. |
| `figure04_representative_topologies.png`, `figure05`–`figure08` (one per metric), `figure09_change_in_failure_rate.png` | Manuscript figures (if enabled). |
| `reproducibility_record.txt` | Configuration, seeds, RNG settings, timestamp, actual R/package session, source checksum. |
| `output_file_md5.csv` | MD5 checksums of every generated file (integrity check). |
| *(copy of the script)* | The script copies itself into the output folder so each result set is self-documenting. |

---

## Reproducibility

- Replication *r* uses `set.seed(3887891 + r)`; with 1,000 replications the seeds are 3887892–3888891 regardless of worker count or scheduling. Each threshold regime restarts the same seeds, so graphs, weights, and targets are matched across regimes.
- The script pins the RNG (`RNGversion("4.4.2")`, Mersenne-Twister / Inversion / Rejection).
- Exact numerical replay was validated with **R 4.4.2** and **igraph 2.1.1**. The script warns when versions differ and records the actual session in `reproducibility_record.txt`. Different package versions can change graph generation even with identical seeds.
- Figures were produced with ggplot2 4.0.0 in the archived environment; figure appearance may vary with your graphics-package versions.

## Troubleshooting

| Problem | Fix |
|---|---|
| `Missing R package(s)` error | Run `install.packages(c("igraph", "ggplot2"))`, then rerun. |
| `Rscript` not found (Windows) | Use the full path to `Rscript.exe`, or use RStudio's Source button. |
| "Run this file with Rscript or source(...)" error | You ran selected lines in RStudio; use **Source** so the script can locate itself. |
| `Output directory is not empty` | Pick a new or empty destination with `--output=PATH` (or leave `OUTPUT_DIRECTORY_OVERRIDE` as `NULL` for an automatic fresh folder). |
| Warning about R/igraph version | Expected if you are not on the validated versions; results may differ numerically. The actual session is recorded in `reproducibility_record.txt`. |
| Figure device fails | Rerun with `--no-figures` (data only), then rerun with a fresh output destination. |
| `NA` values for SD/SE | Normal for tiny replication counts (e.g., single-replication checks have no Monte Carlo SD). |

## Scope and limitations

The study is fixed at **20 vertices** with the configured densities, families, metrics, seed classes, and threshold intervals. This is **not** a general arbitrary-size simulator: changing graph size, factors, metric definitions, or generator order requires revising dependent checks, summaries, and figure labels elsewhere in the script (graph construction lives in Section 2). Interpretation notes (e.g., partial correlations are a modeling convention for cascade coupling, not identified causal transmission coefficients; Table 2/3 flags are descriptive, not significance tests) are documented in the script header.
