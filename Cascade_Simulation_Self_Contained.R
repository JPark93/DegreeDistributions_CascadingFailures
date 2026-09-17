# BEGIN USER GUIDE
# =============================================================================
# SELF-CONTAINED MANUSCRIPT CASCADE SIMULATION: USER GUIDE AND DATA DICTIONARY
# =============================================================================
# Purpose: simulate how graph family, density, centrality-based seed selection,
# and activation thresholds affect final cascade size and conditional duration.
# All graph generation, weighting, target selection, simulation, analysis, and
# plotting code is in this file. No project data or other R scripts are needed.
# R and the external packages listed below must be installed separately.
#
# QUICK START
# 1. Save this file in a writable folder. Install R if it is not installed.
# 2. In an R console, install missing dependencies:
#      install.packages(c("igraph", "ggplot2"),
#                       repos = "https://cloud.r-project.org")
#    parallel is included with R. Installation needs Internet access; the
#    simulation itself runs locally and does not need a network connection.
# 3. Review the USER CONFIGURATION block below, then run the whole saved file:
#      Rscript Cascade_Simulation_Self_Contained.R
#    On Windows, use the full path to Rscript.exe if it is not on PATH.
#    You may also use RStudio's Source button or R's source() on the saved file.
#    The file can be renamed or moved; its path is discovered when it runs.
# 4. Wait for the final success message and open the reported output directory.
#    cascade_summary_all_conditions.csv is the broadest results overview.
#    table2_change_in_failure_rate.csv and table3_steps_to_complete_failure.csv
#    contain the two manuscript analyses. PNG files contain the figures.
#
# REPRODUCIBLE ENVIRONMENT
# Exact numerical replay of the source simulation was validated using R 4.4.2
# and igraph 2.1.1; the archived graphics environment used ggplot2 4.0.0.
# Installing current package releases does not pin those historical versions.
# The script warns when R or igraph differs from the validated versions and
# records the actual session in reproducibility_record.txt. An alternative
# package version can change graph generation even when the seed is identical.
# The explicitly configured RNG is Mersenne-Twister / Inversion / Rejection,
# using RNGversion("4.4.2"). Figures may vary with graphics-package versions.
#
# CONFIGURATION
# STUDY_DESIGN: 1 = primary study with WS rewiring probability 0.10;
#               2 = sensitivity study with WS rewiring probability 0.05.
# Both designs otherwise use the same configured factors and base seed.
# THRESHOLD_CONDITION: 0 = all four regimes; 1 = Sensitive; 2 = Moderate;
#                     3 = Elevated; 4 = Resilient.
# REPLICATIONS_OVERRIDE: NA_integer_ = 1,000 replications per regime;
#                        10L is a small execution check, not manuscript data.
# N_WORKERS: requested parallel processes, capped at the replication count and
#            one fewer than detected logical CPUs, with a minimum of one.
# OUTPUT_DIRECTORY_OVERRIDE: NULL selects a fresh timestamped output folder;
#            an absolute path selects another destination, while a relative
#            path is resolved relative to this saved script's directory.
# MAKE_FIGURES: generate the four final-size figures and contrast figure.
# MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE: generate a separate six-graph illustration.
# The default folder is created beside this script. Explicit destinations must
# be new or empty, protecting previous runs. Choose a fresh folder for each run.
#
# COMMAND-LINE OPTIONS (override the configuration block)
#   --help             print usage without running the simulation
#   --smoke            run 10 replications per selected threshold regime
#   --study=1|2        choose primary or sensitivity study
#   --threshold=0..4   all regimes (0) or one specified regime (1 through 4)
#   --replications=N   choose a positive integer number of replications
#   --workers=N        request a positive integer number of worker processes
#   --output=PATH      select a new or empty output destination
#   --no-figures       disable both figure-generation options
# Example execution check:
#   Rscript Cascade_Simulation_Self_Contained.R --smoke --workers=1 --no-figures
# Example full primary study, with automatic fresh output destination:
#   Rscript Cascade_Simulation_Self_Contained.R
# Example sensitive-threshold subset with a custom destination:
#   Rscript Cascade_Simulation_Self_Contained.R --threshold=1 --output=Sensitive
# In a shell, quote the entire --output=PATH argument if its path has spaces.
#
# The supplied study is fixed at 20 vertices. Changing graph size, factors,
# metric definitions, or generator order requires revising dependent checks,
# summaries, and figure labels; this is not a general arbitrary-size simulator.
#
# STUDY DESIGN AND GRAPH CONSTRUCTION
# Factors: 2 densities x 3 graph families x 4 centrality metrics x 3 seed classes
#          x 4 threshold regimes = 288 cells, each with 1,000 replications.
# Full design: 288,000 target-level cascades, including the forced active seed.
# Sparse: 20 edges among 20 vertices, density 20/190 = 0.1052632.
# Denser: 40 edges among 20 vertices, density 40/190 = 0.2105263.
# All graphs are simple, undirected, and retain all 20 vertices, even isolates.
# ER = Erdos-Renyi G(n,m), with unrestricted connectivity and exact edge count.
# BA = Barabasi-Albert linear preferential attachment (power = 1), m = 1 or 2.
#      At n = 20 these initially give 19 or 37 edges, then receive 1 or 3 random
#      missing edges to match 20 or 40. No pruning is needed in these conditions.
# WS = one-dimensional Watts-Strogatz, nei = 1 or 2, with rewiring probability
#      0.10 (primary) or 0.05 (sensitivity). Its edge count is checked directly;
#      a mismatch stops the run. No connectivity restriction is imposed.
# Graph generation order is BA, ER, WS within each density; analysis order is
# S_ER, S_BA, S_WS, D_ER, D_BA, D_WS. These orders preserve RNG consumption.
# The comparison varies connectivity and other topology features as well as
# degree distribution. A 20-vertex BA graph is not evidence of a population
# power law; the WS generator label does not guarantee a small-world signature.
#
# EDGE WEIGHTS
# Every present edge receives a target positive partial correlation from
# Uniform(0.10, 0.20). Absent edges remain zero. Form a unit-diagonal precision
# matrix Omega with negative target values at present off-diagonal entries.
# If necessary, multiply all off-diagonal entries by one common positive factor
# so lambda_min(Omega) >= 0.05. Recover partial correlations as
#   rho_ij = -Omega_ij / sqrt(Omega_ii * Omega_jj).
# Cascade coupling is abs(rho_ij); all couplings here are positive. Uniform
# scaling preserves sparsity and relative weights, but can lower realized
# magnitudes below the requested interval. Scaling is recorded in diagnostics.
# Graph families receive independent edge draws from the same requested range;
# their realized edge-weight multisets are not constrained to be identical.
# Multiplying every weight in one graph by a positive constant leaves its
# cascade activation fractions unchanged. No observation sample is generated
# and no empirical network estimation is performed by this cascade script.
#
# CENTRALITY AND SEED SELECTION
# degree: binary number of incident edges, ignoring the simulated weights.
# strength: sum of absolute partial-correlation weights incident to a vertex.
# closeness: normalized harmonic centrality, sum of inverse finite shortest-path
#            distances divided by n - 1; unreachable vertices contribute zero.
# betweenness: normalized undirected weighted shortest-path betweenness;
#              normalization divides raw betweenness by choose(n - 1, 2).
# Both path metrics use edge lengths 1 / max(abs(rho_ij), machine epsilon).
# Nonfinite computed centralities are replaced with zero.
# Seed candidates exclude isolates but those vertices remain in the cascade.
# High = maximum eligible centrality; Low = minimum eligible centrality;
# Average = closest to the mean centrality among eligible vertices.
# Exact ties are sampled uniformly, in High, Low, Average order. A class can
# select the same vertex as another class; distinct seeds are not enforced.
# All metrics start tie selection from the same saved RNG state. Centrality
# selects the initial seed only; every metric uses the identical weighted
# cascade update. Seed eligibility is not restricted to the largest component.
#
# THRESHOLDS, MATCHING, AND SEEDS
# Sensitive Uniform(0.10, 0.20); Moderate Uniform(0.20, 0.40);
# Elevated Uniform(0.40, 0.60); Resilient Uniform(0.60, 0.80).
# Each graph receives one vector of 20 independently drawn vertex thresholds.
# All four metrics and all three seed classes share this vector and weighted
# graph within a replication. Different graphs receive separate threshold draws.
# Replication r starts with set.seed(3887891 + r). With 1,000 replications the
# seeds are 3887892 through 3888891, regardless of worker count or scheduling.
# Each threshold regime restarts those same replication seeds. Consequently,
# graphs, weights, seed selections, and underlying uniform threshold draws are
# matched across regimes; only the interval used to rescale thresholds changes.
# Graph/weight/target diagnostics are therefore saved once per replication,
# using the first requested threshold regime, and do not have a threshold column.
# Sparse and denser graphs are separate realizations, not nested edge sets.
# Study 2 changes the WS generator's random-number consumption; the shared base
# seed alone does not establish full cross-study pairing of weights and targets.
#
# CASCADE AND OUTCOMES
# At time zero only the selected seed is active. At each synchronous update,
# q_v = sum_j(weight_vj * active_j) / sum_j(weight_vj). An inactive vertex
# activates when q_v >= threshold_v, including equality. All vertices use the
# previous state for the update, and activation is absorbing. For isolates q_v
# is defined as zero, so a nonseed isolate cannot activate with these thresholds.
# The process stops when no vertex changes state. At most 19 state-changing
# updates are possible after one forced seed in a 20-vertex graph.
# final_active_count: number of active vertices at the stopping point, 1 to 20.
# final_active: final_active_count / 20, including the seed; range 0.05 to 1.00.
# steps_to_stopping_point: synchronous updates adding at least one active
#   vertex; excludes time-zero seeding and the final unchanged checking update.
#   Zero means the cascade never spread beyond its seed.
# complete_failure: integer 1 when final_active_count == 20, otherwise 0.
# Conditional steps to complete failure uses steps_to_stopping_point only for
# complete_failure == 1. It is undefined (NA) if no replication completes, and
# is an update-layer count rather than clock time. A disconnected graph cannot
# reach all 20 vertices from one seed under these positive thresholds.
#
# SUMMARY STATISTICS AND INTERPRETATION
# SD describes variation across replications; Monte Carlo SE = SD / sqrt(n).
# 95% Monte Carlo intervals use mean +/- 1.96 * SE. Final-proportion intervals
# are truncated to [0, 1]; conditional-duration lower bounds are truncated at 0.
# SD and SE are NA when fewer than two contributing observations exist. A
# conditional mean with one completing replication exists but has no SD/SE/CI.
# Completion percentages use all replications as denominator. Conditional
# duration summaries use only completers and may compare different selected
# subsets across targets, densities, and graph families. Always read duration
# alongside n_complete_failures and percent_complete_failures.
# Paired change in failure rate = High final_active - Low final_active within
# the same replication, graph, metric, and threshold regime. It can be negative.
# Table 2 descriptive flags identify maximum unrounded mean change across graph
# families for each density/metric/threshold; exact ties are all flagged.
# Table 3 descriptive flags identify minimum defined unrounded conditional mean
# across graph families for each density/metric/target/threshold; ties are all
# flagged, undefined cells never are. Neither flag is a significance test.
# Partial correlations serve as a modeling convention for cascade coupling;
# they are not identified causal transmission coefficients. These simulations
# do not establish a clinical intervention effect or an optimal treatment target.
#
# OUTPUT DATA DICTIONARY
# CSV files have headers, no row-number column, and use NA for undefined values.
# Percentages are in 0-100 units; final-active proportions and their differences
# are in proportion units. Vertex identifiers are local 1-based indices within
# a simulated graph, not persistent vertices shared across graph realizations.
# Common keys: replication (1...R), density_condition (Sparse/Denser), network
# (Erdos-Renyi/Barabasi-Albert/Watts-Strogatz), metric
# (degree/strength/closeness/betweenness), target (Low/Average/High), threshold
# (0.10-0.20/0.20-0.40/0.40-0.60/0.60-0.80). graph_key encodes density and family,
# e.g., S_ER or D_BA. Table 3 calls its target column target_centrality.
#
# cascade_replication_level.csv
#   One row per replication x density x family x metric x target x threshold.
#   Columns: the six common keys, target_vertex, final_active_count,
#   final_active, steps_to_stopping_point, complete_failure. This is the raw
#   auditable source for every cascade summary; 288,000 rows for a full run.
# cascade_summary_all_conditions.csv
#   One row per density/family/metric/target/threshold (288 for a full run).
#   mean_/sd_/se_final_active and ci95_lo_/ci95_hi_final_active describe size.
#   mean_/sd_/se_steps_unconditional describe stopping time over ALL cascades.
#   n_complete_failures, percent_complete_failures, and n_replications provide
#   completion counts and denominator. mean_/sd_/se_steps_to_complete_failure
#   and ci95_lo_/ci95_hi_steps_to_complete_failure describe completers only.
# high_low_difference_replication_level.csv
#   Paired raw comparisons (96,000 rows): keys omit target and include
#   final_active_high, final_active_low, high_low_difference = high minus low.
# table2_change_in_failure_rate.csv
#   Paired-comparison summaries (96 rows): mean/median/SD/SE/95% CI of
#   high_low_difference; percent_high_greater, percent_equal, percent_high_lower;
#   n_replications; is_table2_max. Percentages classify each paired difference.
# table3_steps_to_complete_failure.csv
#   Low and High only (192 rows): conditional mean/SD, n_complete_failures,
#   percent_complete_failures, n_replications, and is_table3_fastest.
#   Average targets' conditional summaries are in the all-conditions file.
# density_change_in_high_low_difference.csv
#   One row per family/metric/threshold: mean_high_low_difference_sparse,
#   mean_high_low_difference_denser, and denser_minus_sparse_difference.
# target_position_failure_rate_gaps.csv
#   Mean final proportions for most/average/least targets; most_minus_average
#   and average_minus_least are differences of the corresponding cell means.
# conditional_duration_density_comparison.csv
#   Low/High conditional means, completion counts, and completion percentages
#   for each density; sparse_minus_denser_steps is the difference of conditional
#   means and is not a comparison limited to jointly completing replications.
#
# topology_diagnostics_replication_level.csv
#   One row per graph/replication (6,000 full-run rows): n_nodes, n_edges,
#   density, edges_before_matching, edges_added, edges_removed; mean_degree,
#   degree_variance, maximum_degree, global_clustering, mean_path_length,
#   n_components, n_isolates, largest_component, connected (0/1).
#   mean_path_length is UNWEIGHTED and averages reachable pairs only; it is
#   different from the weighted distances used to select centrality targets.
# topology_diagnostics_summary.csv
#   Six family/density summaries of those topology measures, with replication
#   count and percent_connected; sd_degree_variance is across-graph SD.
# target_selection_replication_level.csv
#   One row per graph/metric/replication (24,000 rows): high/average/low_target,
#   high/average/low_ties (numbers of candidates), high_average_same,
#   high_low_same, average_low_same (0/1), eligible_n, centrality_sd (SD across
#   eligible vertices in that realization).
# target_selection_summary.csv
#   Mean tie counts, percentages with more than one candidate, percentages of
#   same-vertex selections, and mean_centrality_sd by graph family/density/metric.
# target_agreement_replication_level.csv and target_agreement_summary.csv
#   Compare every pair of metrics (metric_1, metric_2). Raw high/average/low_agreement
#   indicators identify identical selected vertex IDs; summaries give percentages.
#   Raw full-run count: 36,000. Agreement is not a correlation of metric scores.
# edge_weight_diagnostics_replication_level.csv
#   One row per graph/replication: n_weights, weight_sum, weight_sum_squares,
#   minimum_weight, maximum_weight, within_graph_sd, scale_factor.
# edge_weight_diagnostics.csv
#   Six family/density summaries. mean_weight and sd_weight pool individual
#   edge weights, whereas mean_within_graph_sd averages graph-level SDs.
#   percent_scaled, mean_scale_factor, minimum_scale_factor describe precision
#   adjustments. Each family contributes 20,000 sparse and 40,000 denser weights
#   in the full 1,000-replication design, retained once across threshold regimes.
#
# manuscript_simulation_results.rds
#   Read with results <- readRDS("path/to/manuscript_simulation_results.rds").
#   Named elements: configuration, cascade_summary, table2, table3,
#   density_change_summary, target_position_gaps,
#   conditional_duration_density_comparison, topology_summary, target_summary,
#   target_agreement_summary, edge_weight_summary. This is a summary bundle;
#   raw replications are in the CSVs. Graph adjacency matrices, per-vertex
#   centrality vectors, and threshold vectors are not separately serialized;
#   the saved code, configuration, environment, and seeds specify regeneration.
# figure05_failure_rate_degree.png through figure08_failure_rate_betweenness.png
#   Final-size means and normal-approximation Monte Carlo intervals, by seed
#   class, graph family, density, and threshold regime.
# figure09_change_in_failure_rate.png
#   Paired High-minus-Low means and Monte Carlo intervals.
# figure04_representative_topologies.png and representative_topology_selection.csv
#   Six illustrative graphs, selected independently within each family/density.
#   Selection restricts to the modal connectivity status and minimizes squared
#   standardized deviation from medians of seven topology features; MAD scales
#   use SD and then 1 as fallbacks. Ties select the lower replication index.
#   The CSV adds representativeness_score and modal_connected to diagnostics.
#   Illustration layouts have their own seed and do not alter simulation data.
# reproducibility_record.txt
#   Configuration, seeds, RNG settings, timestamp, actual R/package session,
#   source checksum, and identifiers of the archived numerical provenance.
# output_file_md5.csv
#   file and md5 columns for generated data, figures, record, and copied script.
#   The script copies itself into the output folder. Checksums detect file
#   changes; timestamp and session-dependent files can differ across reruns.
#
# PRACTICAL CHECKS
# Confirm the success message, expected row counts, configuration, and actual
# versions. Use a small replication override for execution checks, and restore
# NA_integer_ for the manuscript run. Each run starts from replication 1;
# there is no resume/checkpoint facility or automatic comparison with archives.
# The code checks edge counts, cascade bounds, row counts, representative-graph
# reconstruction, and presence of expected outputs. Single-replication checks
# naturally have undefined Monte Carlo SD/SE values. If a figure device fails,
# use --no-figures (or disable the corresponding configuration option), then
# rerun with a fresh output destination to generate data only.
# Full numerical reruns can take several minutes; duration depends on hardware.
# =============================================================================
# END USER GUIDE

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================
# Edit these settings OR use the command-line options documented above.
# Command-line options take precedence. The simulation starts when sourced.

STUDY_DESIGN <- 1L                    # 1: WS beta=.10; 2: WS beta=.05
THRESHOLD_CONDITION <- 0L             # 0: all; 1..4: one threshold regime
REPLICATIONS_OVERRIDE <- NA_integer_  # NA: 1,000; e.g., 10L for a quick check
N_WORKERS <- 3L                       # Automatically capped by cores/replications
OUTPUT_DIRECTORY_OVERRIDE <- NULL    # NULL: fresh folder beside this script
MAKE_FIGURES <- TRUE                  # Figures 5-9
MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE <- TRUE  # Figure 4

# =============================================================================
# SETUP: FILE LOCATION, COMMAND-LINE OPTIONS, AND DEPENDENCIES
# =============================================================================

# Source frames take precedence over --file, which may refer to a caller script.
# This also supports RStudio's Source command and moving or renaming this file.
source_files <- vapply(sys.frames(), function(frame) {
  candidate <- frame$ofile
  if (!is.character(candidate) || length(candidate) != 1L || is.na(candidate)) {
    return("")
  }
  candidate <- path.expand(candidate)
  # source(chdir=TRUE) keeps the caller's directory in owd and its original,
  # possibly relative filename in ofile.
  if (!grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", candidate) &&
      is.character(frame$owd) && length(frame$owd) == 1L) {
    candidate <- file.path(frame$owd, candidate)
  }
  candidate
}, character(1L))
source_files <- source_files[nzchar(source_files)]
command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
script_path <- if (length(source_files)) {
  tail(source_files, 1L)
} else if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  stop(
    "Run this file with Rscript or source(\"path/to/this_file.R\"). ",
    "In RStudio, use Source rather than running selected lines.",
    call. = FALSE
  )
}
script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
project_directory <- dirname(script_path)
cli_arguments <- if (length(source_files)) character() else {
  commandArgs(trailingOnly = TRUE)
}

if ("--help" %in% cli_arguments) {
  guide <- readLines(script_path, warn = FALSE)
  guide_start <- match("# BEGIN USER GUIDE", guide)
  guide_end <- match("# END USER GUIDE", guide)
  cat(paste(sub("^# ?", "", guide[seq.int(guide_start + 1L, guide_end - 1L)]),
            collapse = "\n"), "\n")
  quit(save = "no", status = 0L)
}

integer_option <- function(value, name, minimum, maximum = .Machine$integer.max) {
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || is.na(number) || !is.finite(number) ||
      number != floor(number) || number < minimum || number > maximum) {
    stop(name, " must be a whole number between ", minimum, " and ", maximum,
         ".", call. = FALSE)
  }
  as.integer(number)
}

if ("--smoke" %in% cli_arguments) REPLICATIONS_OVERRIDE <- 10L
for (argument in cli_arguments) {
  if (argument == "--smoke") next
  if (argument == "--no-figures") {
    MAKE_FIGURES <- FALSE
    MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE <- FALSE
    next
  }
  if (!grepl("^--(study|threshold|replications|workers|output)=.+$", argument)) {
    stop("Unknown or incomplete option: ", argument, ". Use --help.",
         call. = FALSE)
  }
  option <- sub("=.*$", "", argument)
  value <- sub("^[^=]*=", "", argument)
  switch(option,
    "--study" = { STUDY_DESIGN <- integer_option(value, option, 1L, 2L) },
    "--threshold" = {
      THRESHOLD_CONDITION <- integer_option(value, option, 0L, 4L)
    },
    "--replications" = {
      REPLICATIONS_OVERRIDE <- integer_option(value, option, 1L)
    },
    "--workers" = { N_WORKERS <- integer_option(value, option, 1L) },
    "--output" = { OUTPUT_DIRECTORY_OVERRIDE <- value }
  )
}
STUDY_DESIGN <- integer_option(STUDY_DESIGN, "STUDY_DESIGN", 1L, 2L)
THRESHOLD_CONDITION <- integer_option(
  THRESHOLD_CONDITION, "THRESHOLD_CONDITION", 0L, 4L
)
N_WORKERS <- integer_option(N_WORKERS, "N_WORKERS", 1L)
if (length(REPLICATIONS_OVERRIDE) != 1L) {
  stop("REPLICATIONS_OVERRIDE must be NA or one positive integer.", call. = FALSE)
}
if (!is.na(REPLICATIONS_OVERRIDE)) {
  REPLICATIONS_OVERRIDE <- integer_option(
    REPLICATIONS_OVERRIDE, "REPLICATIONS_OVERRIDE", 1L,
    .Machine$integer.max - 3887891L
  )
}
for (setting in c("MAKE_FIGURES", "MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE")) {
  value <- get(setting)
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(setting, " must be TRUE or FALSE.", call. = FALSE)
  }
}
if (!is.null(OUTPUT_DIRECTORY_OVERRIDE) &&
    (!is.character(OUTPUT_DIRECTORY_OVERRIDE) ||
     length(OUTPUT_DIRECTORY_OVERRIDE) != 1L ||
     is.na(OUTPUT_DIRECTORY_OVERRIDE) || !nzchar(OUTPUT_DIRECTORY_OVERRIDE))) {
  stop("OUTPUT_DIRECTORY_OVERRIDE must be NULL or one nonempty path.",
       call. = FALSE)
}

required_packages <- c("igraph", "ggplot2")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nRun install.packages(c(\"igraph\", \"ggplot2\")) in R, then rerun this file.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(igraph))
library(parallel)
suppressPackageStartupMessages(library(ggplot2))

study_conditions <- list(
  `1` = list(
    label = "Primary manuscript study: beta = .10",
    n_nodes = 20L,
    edge_counts = c(Sparse = 20L, Denser = 40L),
    ws_rewiring = 0.10,
    n_replications = 1000L,
    base_seed = 3887891L,
    output_directory =
      "Cascade_Results_Study1"
  ),
  `2` = list(
    label = "Watts-Strogatz sensitivity study: beta = .05",
    n_nodes = 20L,
    edge_counts = c(Sparse = 20L, Denser = 40L),
    ws_rewiring = 0.05,
    n_replications = 1000L,
    base_seed = 3887891L,
    output_directory =
      "Cascade_Results_Study2"
  )
)

study_key <- as.character(as.integer(STUDY_DESIGN))
if (!study_key %in% names(study_conditions)) {
  stop("STUDY_DESIGN must be 1 or 2.")
}
design <- study_conditions[[study_key]]

if (!is.na(REPLICATIONS_OVERRIDE)) {
  design$n_replications <- as.integer(REPLICATIONS_OVERRIDE)
}
if (design$n_replications < 1L) stop("At least one replication is required.")

output_directory_name <- if (THRESHOLD_CONDITION == 0L) {
  design$output_directory
} else {
  paste0(design$output_directory, "_Threshold", THRESHOLD_CONDITION)
}
if (!is.na(REPLICATIONS_OVERRIDE)) {
  output_directory_name <- paste0(
    output_directory_name,
    "_Smoke_n",
    design$n_replications
  )
}

output_directory <- if (is.null(OUTPUT_DIRECTORY_OVERRIDE)) {
  run_name <- paste0(output_directory_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  destination <- file.path(project_directory, run_name)
  duplicate <- 0L
  while (file.exists(destination)) {
    duplicate <- duplicate + 1L
    destination <- file.path(project_directory, sprintf("%s_%03d", run_name, duplicate))
  }
  destination
} else {
  requested_output <- path.expand(OUTPUT_DIRECTORY_OVERRIDE)
  if (grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", requested_output)) {
    requested_output
  } else {
    file.path(project_directory, requested_output)
  }
}
if (file.exists(output_directory) && !dir.exists(output_directory)) {
  stop("The output path is a file. Choose a directory with --output=PATH.",
       call. = FALSE)
}
if (dir.exists(output_directory) && length(list.files(
    output_directory, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty: ", output_directory,
       "\nChoose a new or empty directory with --output=PATH.", call. = FALSE)
}

metrics <- c("degree", "strength", "closeness", "betweenness")
target_order <- c("Low", "Average", "High")
graph_keys <- c("S_ER", "S_BA", "S_WS", "D_ER", "D_BA", "D_WS")
threshold_ranges <- list(
  Sensitive = c(0.10, 0.20),
  Moderate = c(0.20, 0.40),
  Elevated = c(0.40, 0.60),
  Resilient = c(0.60, 0.80)
)
if (!THRESHOLD_CONDITION %in% 0:4) {
  stop("THRESHOLD_CONDITION must be 0, 1, 2, 3, or 4.")
}
active_threshold_ranges <- if (THRESHOLD_CONDITION == 0L) {
  threshold_ranges
} else {
  threshold_ranges[as.integer(THRESHOLD_CONDITION)]
}
threshold_codes <- vapply(
  active_threshold_ranges,
  function(x) sprintf("%.2f-%.2f", x[1L], x[2L]),
  character(1L)
)
threshold_display <- setNames(
  sprintf("%s (%s)", names(active_threshold_ranges), threshold_codes),
  threshold_codes
)

network_names <- c(
  ER = "Erdos-Renyi",
  BA = "Barabasi-Albert",
  WS = "Watts-Strogatz"
)
density_names <- c(S = "Sparse", D = "Denser")

# The simulated partial correlations are positive and uniformly distributed.
weight_lower <- 0.10
weight_upper <- 0.20
minimum_precision_eigenvalue <- 0.05
exclude_isolates_as_targets <- TRUE

# Explicit version and RNG settings reproduce the archived R 4.4.2 stream.
RNGversion("4.4.2")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

expected_versions <- c(R = "4.4.2", igraph = "2.1.1")
if (as.character(getRversion()) != expected_versions[["R"]]) {
  warning("Exact replay was validated with R 4.4.2; current R is ", getRversion(), ".")
}
if (as.character(packageVersion("igraph")) != expected_versions[["igraph"]]) {
  warning(
    "Exact replay was validated with igraph 2.1.1; current igraph is ",
    packageVersion("igraph"),
    "."
  )
}

stopifnot(
  design$n_nodes == 20L,
  identical(names(design$edge_counts), c("Sparse", "Denser")),
  all(design$edge_counts >= design$n_nodes - 1L),
  all(design$edge_counts <= choose(design$n_nodes, 2L)),
  N_WORKERS >= 1L
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_directory) || file.access(output_directory, 2L) != 0L) {
  stop("Cannot write to output directory: ", output_directory, call. = FALSE)
}

# =============================================================================
# 2. GRAPH GENERATION AND EDGE WEIGHTS
# =============================================================================

density_from_key <- function(key) {
  density_names[[strsplit(key, "_", fixed = TRUE)[[1L]][1L]]]
}

network_from_key <- function(key) {
  network_names[[strsplit(key, "_", fixed = TRUE)[[1L]][2L]]]
}

adjust_edge_count <- function(graph, target_edges) {
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
  target_edges <- as.integer(target_edges)

  while (gsize(graph) > target_edges) {
    removable <- setdiff(seq_len(gsize(graph)), as.integer(bridges(graph)))
    if (!length(removable)) {
      stop("No non-bridge edge was available for edge-count matching.")
    }
    graph <- delete_edges(
      graph,
      removable[sample.int(length(removable), 1L)]
    )
  }

  while (gsize(graph) < target_edges) {
    missing_edges <- as_edgelist(
      complementer(graph, loops = FALSE),
      names = FALSE
    )
    selected <- missing_edges[sample.int(nrow(missing_edges), 1L), ]
    graph <- add_edges(graph, as.integer(selected))
  }

  graph
}

generate_topologies <- function(n_nodes, target_edges, ws_rewiring) {
  target_edges <- as.integer(target_edges)

  # Barabasi-Albert graphs are generated first to preserve the study's random
  # number sequence, then matched to the requested edge count.
  edges_per_new_vertex <- max(1L, as.integer(round(target_edges / n_nodes)))
  graph_ba <- sample_pa(
    n = n_nodes,
    power = 1,
    m = edges_per_new_vertex,
    directed = FALSE
  )
  graph_ba <- simplify(graph_ba, remove.multiple = TRUE, remove.loops = TRUE)
  ba_edges_before_matching <- gsize(graph_ba)
  graph_ba <- adjust_edge_count(graph_ba, target_edges)
  graph_ba <- set_graph_attr(
    graph_ba,
    "edges_before_matching",
    value = ba_edges_before_matching
  )

  graph_er <- sample_gnm(
    n = n_nodes,
    m = target_edges,
    directed = FALSE,
    loops = FALSE
  )
  graph_er <- set_graph_attr(
    graph_er,
    "edges_before_matching",
    value = gsize(graph_er)
  )

  lattice_neighbors <- max(
    1L,
    as.integer(round(target_edges / n_nodes))
  )
  graph_ws <- sample_smallworld(
    dim = 1,
    size = n_nodes,
    nei = lattice_neighbors,
    p = ws_rewiring,
    loops = FALSE,
    multiple = FALSE
  )
  ws_edges_before_matching <- gsize(graph_ws)
  if (ws_edges_before_matching != target_edges) {
    stop("The configured Watts-Strogatz graph did not have the target edge count.")
  }
  graph_ws <- set_graph_attr(
    graph_ws,
    "edges_before_matching",
    value = ws_edges_before_matching
  )

  stopifnot(
    gsize(graph_er) == target_edges,
    gsize(graph_ba) == target_edges,
    gsize(graph_ws) == target_edges
  )

  # The return order matches the original analysis after all three graph
  # generators have consumed their random numbers.
  list(ER = graph_er, BA = graph_ba, WS = graph_ws)
}

add_partial_correlation_weights <- function(
    graph,
    lower = 0.10,
    upper = 0.20,
    minimum_eigenvalue = 0.05
) {
  edge_list <- as_edgelist(graph, names = FALSE)
  n_edges <- nrow(edge_list)
  magnitudes <- runif(n_edges, lower, upper)

  off_diagonal <- matrix(0, vcount(graph), vcount(graph))
  if (n_edges) {
    off_diagonal[edge_list] <- -magnitudes
    off_diagonal[cbind(edge_list[, 2L], edge_list[, 1L])] <- -magnitudes
  }

  minimum_off_diagonal_eigenvalue <- min(
    eigen(off_diagonal, symmetric = TRUE, only.values = TRUE)$values
  )
  scale_factor <- 1
  if (1 + minimum_off_diagonal_eigenvalue < minimum_eigenvalue) {
    scale_factor <-
      (1 - minimum_eigenvalue) / abs(minimum_off_diagonal_eigenvalue)
  }

  precision <- diag(vcount(graph)) + scale_factor * off_diagonal
  diagonal_scale <- 1 / sqrt(diag(precision))
  partial_correlation <-
    -outer(diagonal_scale, diagonal_scale) * precision
  diag(partial_correlation) <- 1

  if (n_edges) {
    E(graph)$pcor <- partial_correlation[edge_list]
    E(graph)$weight <- abs(E(graph)$pcor)
  }
  graph <- set_graph_attr(graph, "pcor_scale_factor", value = scale_factor)
  graph
}

prepare_graphs <- function(configuration) {
  density_codes <- c(Sparse = "S", Denser = "D")
  graphs <- unlist(lapply(names(configuration$edge_counts), function(density) {
    generated <- generate_topologies(
      configuration$n_nodes,
      configuration$edge_counts[[density]],
      configuration$ws_rewiring
    )
    names(generated) <- paste(density_codes[[density]], names(generated), sep = "_")
    generated
  }), recursive = FALSE)

  graphs <- lapply(graphs, add_partial_correlation_weights,
    lower = configuration$weight_lower,
    upper = configuration$weight_upper,
    minimum_eigenvalue = configuration$minimum_precision_eigenvalue
  )
  graphs[configuration$graph_keys]
}

# =============================================================================
# 3. CENTRALITY, TARGET SELECTION, AND CASCADE PROCESS
# =============================================================================

compute_centrality <- function(graph, metric) {
  coupling <- abs(E(graph)$pcor)
  distance <- 1 / pmax(coupling, .Machine$double.eps)

  value <- switch(
    metric,
    degree = degree(graph, mode = "all", loops = FALSE),
    strength = strength(
      graph,
      mode = "all",
      loops = FALSE,
      weights = coupling
    ),
    closeness = harmonic_centrality(
      graph,
      mode = "all",
      weights = distance,
      normalized = TRUE
    ),
    betweenness = betweenness(
      graph,
      directed = FALSE,
      weights = distance,
      normalized = TRUE
    ),
    stop("Unknown centrality metric.")
  )

  value[!is.finite(value)] <- 0
  as.numeric(value)
}

select_targets <- function(graph, centrality, exclude_isolates = TRUE) {
  eligible <- seq_len(vcount(graph))
  if (exclude_isolates) {
    nonisolates <- which(degree(graph, loops = FALSE) > 0)
    if (length(nonisolates)) eligible <- nonisolates
  }

  values <- centrality[eligible]
  if (!length(values) || any(!is.finite(values))) {
    stop("No valid target candidates were available.")
  }

  select_one <- function(candidates) {
    candidates[sample.int(length(candidates), 1L)]
  }

  high_candidates <- eligible[values == max(values)]
  low_candidates <- eligible[values == min(values)]
  distance_from_mean <- abs(values - mean(values))
  average_candidates <- eligible[distance_from_mean == min(distance_from_mean)]

  # This draw order is retained for exact reproduction of the original study.
  high_target <- select_one(high_candidates)
  low_target <- select_one(low_candidates)
  average_target <- select_one(average_candidates)

  list(
    targets = c(
      High = high_target,
      Average = average_target,
      Low = low_target
    ),
    tie_counts = c(
      High = length(high_candidates),
      Average = length(average_candidates),
      Low = length(low_candidates)
    ),
    centrality_sd = stats::sd(values),
    eligible_n = length(eligible)
  )
}

run_cascade <- function(graph, target, thresholds) {
  adjacency <- as.matrix(
    as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)
  )
  incident_weight <- rowSums(adjacency)
  state <- integer(vcount(graph))
  state[target] <- 1L
  state_changing_steps <- 0L

  for (step in seq_len(vcount(graph))) {
    active_incident_weight <- as.numeric(adjacency %*% state)
    active_fraction <- numeric(vcount(graph))
    defined <- incident_weight > 0
    active_fraction[defined] <-
      active_incident_weight[defined] / incident_weight[defined]

    updated_state <- pmax(
      state,
      as.integer(active_fraction >= thresholds)
    )
    if (identical(updated_state, state)) break
    state <- updated_state
    state_changing_steps <- state_changing_steps + 1L
  }

  c(
    final_active_count = sum(state),
    final_active = mean(state),
    steps_to_stopping_point = state_changing_steps
  )
}

topology_diagnostics <- function(graph, replication, graph_key) {
  graph_degree <- degree(graph, loops = FALSE)
  graph_components <- components(graph, mode = "weak")
  edges_before <- as.integer(graph_attr(graph, "edges_before_matching"))

  data.frame(
    replication = replication,
    graph_key = graph_key,
    density_condition = density_from_key(graph_key),
    network = network_from_key(graph_key),
    n_nodes = vcount(graph),
    n_edges = gsize(graph),
    density = edge_density(graph, loops = FALSE),
    edges_before_matching = edges_before,
    edges_added = max(0L, gsize(graph) - edges_before),
    edges_removed = max(0L, edges_before - gsize(graph)),
    mean_degree = mean(graph_degree),
    degree_variance = stats::var(graph_degree),
    maximum_degree = max(graph_degree),
    global_clustering = transitivity(
      graph,
      type = "global",
      isolates = "zero"
    ),
    mean_path_length = mean_distance(
      graph,
      directed = FALSE,
      unconnected = TRUE,
      weights = NA
    ),
    n_components = graph_components$no,
    n_isolates = sum(graph_degree == 0),
    largest_component = max(graph_components$csize),
    connected = as.integer(graph_components$no == 1L),
    stringsAsFactors = FALSE
  )
}

weight_diagnostics <- function(graph, replication, graph_key) {
  weights <- E(graph)$weight
  data.frame(
    replication = replication,
    graph_key = graph_key,
    density_condition = density_from_key(graph_key),
    network = network_from_key(graph_key),
    n_weights = length(weights),
    weight_sum = sum(weights),
    weight_sum_squares = sum(weights^2),
    minimum_weight = min(weights),
    maximum_weight = max(weights),
    within_graph_sd = stats::sd(weights),
    scale_factor = as.numeric(graph_attr(graph, "pcor_scale_factor")),
    stringsAsFactors = FALSE
  )
}

simulate_replication <- function(replication, threshold_limits, threshold_code, cfg) {
  set.seed(cfg$base_seed + as.integer(replication))
  graphs <- prepare_graphs(cfg)
  metrics_for_run <- cfg$metrics

  cascade_rows <- list()
  topology_rows <- list()
  target_rows <- list()
  agreement_rows <- list()
  weight_rows <- list()

  for (graph_key in names(graphs)) {
    graph <- graphs[[graph_key]]
    thresholds <- runif(
      vcount(graph),
      min = threshold_limits[1L],
      max = threshold_limits[2L]
    )

    topology_rows[[graph_key]] <-
      topology_diagnostics(graph, replication, graph_key)
    weight_rows[[graph_key]] <-
      weight_diagnostics(graph, replication, graph_key)

    centralities <- setNames(
      lapply(metrics_for_run, function(metric) compute_centrality(graph, metric)),
      metrics_for_run
    )

    # Every metric begins target tie-breaking from the same random-number state.
    target_rng_state <- .Random.seed
    target_information <- setNames(lapply(metrics_for_run, function(metric) {
      assign(".Random.seed", target_rng_state, envir = .GlobalEnv)
      select_targets(
        graph,
        centralities[[metric]],
        exclude_isolates = cfg$exclude_isolates_as_targets
      )
    }), metrics_for_run)

    for (metric in metrics_for_run) {
      target_info <- target_information[[metric]]
      targets <- target_info$targets

      target_rows[[paste(graph_key, metric, sep = "__")]] <- data.frame(
        replication = replication,
        graph_key = graph_key,
        density_condition = density_from_key(graph_key),
        network = network_from_key(graph_key),
        metric = metric,
        high_target = unname(targets["High"]),
        average_target = unname(targets["Average"]),
        low_target = unname(targets["Low"]),
        high_ties = unname(target_info$tie_counts["High"]),
        average_ties = unname(target_info$tie_counts["Average"]),
        low_ties = unname(target_info$tie_counts["Low"]),
        high_average_same = as.integer(targets["High"] == targets["Average"]),
        high_low_same = as.integer(targets["High"] == targets["Low"]),
        average_low_same = as.integer(targets["Average"] == targets["Low"]),
        centrality_sd = target_info$centrality_sd,
        eligible_n = target_info$eligible_n,
        stringsAsFactors = FALSE
      )

      # Retain the original cascade evaluation order.
      for (target_class in c("High", "Average", "Low")) {
        result <- run_cascade(
          graph,
          unname(targets[target_class]),
          thresholds
        )
        cascade_rows[[paste(graph_key, metric, target_class, sep = "__")]] <-
          data.frame(
            replication = replication,
            density_condition = density_from_key(graph_key),
            network = network_from_key(graph_key),
            metric = metric,
            target = target_class,
            target_vertex = as.integer(unname(targets[target_class])),
            threshold = threshold_code,
            final_active_count = as.integer(result["final_active_count"]),
            final_active = unname(result["final_active"]),
            steps_to_stopping_point = as.integer(
              result["steps_to_stopping_point"]
            ),
            complete_failure = as.integer(
              result["final_active_count"] == cfg$n_nodes
            ),
            stringsAsFactors = FALSE
          )
      }
    }

    metric_pairs <- combn(metrics_for_run, 2L, simplify = FALSE)
    for (pair in metric_pairs) {
      first <- target_information[[pair[1L]]]$targets
      second <- target_information[[pair[2L]]]$targets
      agreement_rows[[paste(graph_key, pair, collapse = "__")]] <- data.frame(
        replication = replication,
        graph_key = graph_key,
        density_condition = density_from_key(graph_key),
        network = network_from_key(graph_key),
        metric_1 = pair[1L],
        metric_2 = pair[2L],
        high_agreement = as.integer(first["High"] == second["High"]),
        average_agreement = as.integer(
          first["Average"] == second["Average"]
        ),
        low_agreement = as.integer(first["Low"] == second["Low"]),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    cascades = do.call(rbind, cascade_rows),
    topology = do.call(rbind, topology_rows),
    targets = do.call(rbind, target_rows),
    target_agreement = do.call(rbind, agreement_rows),
    weights = do.call(rbind, weight_rows)
  )
}

# =============================================================================
# 4. RUN THE MONTE CARLO STUDY
# =============================================================================

detected_workers <- suppressWarnings(detectCores(logical = TRUE))
if (!is.finite(detected_workers)) detected_workers <- 1L
workers <- max(
  1L,
  min(
    as.integer(N_WORKERS),
    design$n_replications,
    max(1L, detected_workers - 1L)
  )
)

configuration_for_workers <- c(
  design,
  list(
    metrics = metrics,
    graph_keys = graph_keys,
    weight_lower = weight_lower,
    weight_upper = weight_upper,
    minimum_precision_eigenvalue = minimum_precision_eigenvalue,
    exclude_isolates_as_targets = exclude_isolates_as_targets
  )
)

cat(sprintf("Study design: %d\n", STUDY_DESIGN))
cat(sprintf("Threshold condition: %d\n", THRESHOLD_CONDITION))
cat(sprintf("Design: %s\n", design$label))
cat(sprintf("Replications per threshold range: %d\n", design$n_replications))
cat(sprintf("Workers: %d\n", workers))
cat(sprintf("Output directory: %s\n", output_directory))

cluster <- makeCluster(workers)
invisible(clusterEvalQ(cluster, {
  RNGversion("4.4.2")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  library(igraph)
}))
clusterExport(
  cluster,
  c(
    "metrics", "graph_keys", "network_names", "density_names",
    "weight_lower", "weight_upper", "minimum_precision_eigenvalue",
    "exclude_isolates_as_targets", "density_from_key", "network_from_key",
    "adjust_edge_count", "generate_topologies",
    "add_partial_correlation_weights", "prepare_graphs",
    "compute_centrality", "select_targets", "run_cascade",
    "topology_diagnostics", "weight_diagnostics", "simulate_replication"
  ),
  envir = environment()
)

cascade_by_threshold <- vector("list", length(active_threshold_ranges))
topology_raw <- NULL
target_raw <- NULL
target_agreement_raw <- NULL
weight_raw <- NULL

tryCatch(
  {
    for (threshold_index in seq_along(active_threshold_ranges)) {
      threshold_name <- names(active_threshold_ranges)[threshold_index]
      threshold_limits <- active_threshold_ranges[[threshold_index]]
      threshold_code <- threshold_codes[threshold_index]

      cat(sprintf(
        "\nRunning %s thresholds (%s) ...\n",
        threshold_name,
        threshold_code
      ))

      replications <- parLapply(
        cluster,
        seq_len(design$n_replications),
        simulate_replication,
        threshold_limits = threshold_limits,
        threshold_code = threshold_code,
        cfg = configuration_for_workers
      )

      cascade_by_threshold[[threshold_index]] <- do.call(
        rbind,
        lapply(replications, `[[`, "cascades")
      )

      # The same seeded graph and target realizations are reused across threshold
      # ranges, so diagnostics are retained once.
      if (threshold_index == 1L) {
        topology_raw <- do.call(rbind, lapply(replications, `[[`, "topology"))
        target_raw <- do.call(rbind, lapply(replications, `[[`, "targets"))
        target_agreement_raw <- do.call(
          rbind,
          lapply(replications, `[[`, "target_agreement")
        )
        weight_raw <- do.call(rbind, lapply(replications, `[[`, "weights"))
      }

      cat(sprintf(
        "Completed %s: %s cascade results.\n",
        threshold_name,
        format(nrow(cascade_by_threshold[[threshold_index]]), big.mark = ",")
      ))
    }
  },
  finally = {
    try(stopCluster(cluster), silent = TRUE)
  }
)
cluster <- NULL

cascade_raw <- do.call(rbind, cascade_by_threshold)
rownames(cascade_raw) <- NULL
rownames(topology_raw) <- NULL
rownames(target_raw) <- NULL
rownames(target_agreement_raw) <- NULL
rownames(weight_raw) <- NULL

expected_cascade_rows <-
  design$n_replications *
  length(active_threshold_ranges) *
  length(graph_keys) *
  length(metrics) *
  length(target_order)

stopifnot(
  nrow(cascade_raw) == expected_cascade_rows,
  nrow(topology_raw) == design$n_replications * length(graph_keys),
  all(cascade_raw$final_active >= 1 / design$n_nodes),
  all(cascade_raw$final_active <= 1),
  all(cascade_raw$final_active_count == design$n_nodes * cascade_raw$final_active),
  all(cascade_raw$steps_to_stopping_point >= 0L),
  all(cascade_raw$steps_to_stopping_point <= design$n_nodes - 1L),
  all(
    cascade_raw$complete_failure ==
      as.integer(cascade_raw$final_active_count == design$n_nodes)
  )
)

# =============================================================================
# 5. CELL SUMMARIES FOR TABLES AND FIGURES
# =============================================================================

safe_sd <- function(x) {
  if (length(x) >= 2L) stats::sd(x) else NA_real_
}

safe_se <- function(x) {
  if (length(x) >= 2L) stats::sd(x) / sqrt(length(x)) else NA_real_
}

split_by_columns <- function(data, columns) {
  split(
    data,
    interaction(data[columns], drop = TRUE, lex.order = TRUE)
  )
}

cascade_summary <- do.call(rbind, lapply(
  split_by_columns(
    cascade_raw,
    c("density_condition", "network", "metric", "target", "threshold")
  ),
  function(cell) {
    final_active <- cell$final_active
    all_steps <- cell$steps_to_stopping_point
    complete <- cell$final_active_count == design$n_nodes
    complete_steps <- all_steps[complete]
    n_total <- length(final_active)
    n_complete <- length(complete_steps)
    final_se <- safe_se(final_active)
    all_duration_se <- safe_se(all_steps)
    complete_duration_se <- safe_se(complete_steps)

    data.frame(
      density_condition = cell$density_condition[1L],
      network = cell$network[1L],
      metric = cell$metric[1L],
      target = cell$target[1L],
      threshold = cell$threshold[1L],
      mean_final_active = mean(final_active),
      sd_final_active = safe_sd(final_active),
      se_final_active = final_se,
      ci95_lo_final_active = max(0, mean(final_active) - 1.96 * final_se),
      ci95_hi_final_active = min(1, mean(final_active) + 1.96 * final_se),
      mean_steps_unconditional = mean(all_steps),
      sd_steps_unconditional = safe_sd(all_steps),
      se_steps_unconditional = all_duration_se,
      n_complete_failures = n_complete,
      percent_complete_failures = 100 * n_complete / n_total,
      mean_steps_to_complete_failure = if (n_complete) {
        mean(complete_steps)
      } else {
        NA_real_
      },
      sd_steps_to_complete_failure = safe_sd(complete_steps),
      se_steps_to_complete_failure = complete_duration_se,
      ci95_lo_steps_to_complete_failure = if (n_complete >= 2L) {
        max(0, mean(complete_steps) - 1.96 * complete_duration_se)
      } else {
        NA_real_
      },
      ci95_hi_steps_to_complete_failure = if (n_complete >= 2L) {
        mean(complete_steps) + 1.96 * complete_duration_se
      } else {
        NA_real_
      },
      n_replications = n_total,
      stringsAsFactors = FALSE
    )
  }
))
rownames(cascade_summary) <- NULL

high_results <- cascade_raw[cascade_raw$target == "High", ]
low_results <- cascade_raw[cascade_raw$target == "Low", ]
paired_results <- merge(
  high_results[
    , c(
      "replication", "density_condition", "network", "metric", "threshold",
      "final_active"
    )
  ],
  low_results[
    , c(
      "replication", "density_condition", "network", "metric", "threshold",
      "final_active"
    )
  ],
  by = c(
    "replication", "density_condition", "network", "metric", "threshold"
  ),
  suffixes = c("_high", "_low"),
  sort = FALSE
)
paired_results$high_low_difference <-
  paired_results$final_active_high - paired_results$final_active_low

contrast_summary <- do.call(rbind, lapply(
  split_by_columns(
    paired_results,
    c("density_condition", "network", "metric", "threshold")
  ),
  function(cell) {
    difference <- cell$high_low_difference
    difference_se <- safe_se(difference)
    data.frame(
      density_condition = cell$density_condition[1L],
      network = cell$network[1L],
      metric = cell$metric[1L],
      threshold = cell$threshold[1L],
      mean_high_low_difference = mean(difference),
      median_high_low_difference = median(difference),
      sd_high_low_difference = safe_sd(difference),
      se_high_low_difference = difference_se,
      ci95_lo_high_low_difference = mean(difference) - 1.96 * difference_se,
      ci95_hi_high_low_difference = mean(difference) + 1.96 * difference_se,
      percent_high_greater = 100 * mean(difference > 0),
      percent_equal = 100 * mean(difference == 0),
      percent_high_lower = 100 * mean(difference < 0),
      n_replications = length(difference),
      stringsAsFactors = FALSE
    )
  }
))
rownames(contrast_summary) <- NULL

# Table 2 shading: maximum unrounded mean change across graph families within
# each metric-by-density-by-threshold comparison. Exact ties are all flagged.
contrast_summary$is_table2_max <- FALSE
table2_groups <- split_by_columns(
  contrast_summary,
  c("density_condition", "metric", "threshold")
)
for (group in table2_groups) {
  group_key <- with(
    group,
    paste(density_condition, metric, threshold, sep = "__")
  )[1L]
  full_key <- with(
    contrast_summary,
    paste(density_condition, metric, threshold, sep = "__")
  )
  candidates <- which(full_key == group_key)
  contrast_summary$is_table2_max[candidates] <-
    contrast_summary$mean_high_low_difference[candidates] ==
      max(contrast_summary$mean_high_low_difference[candidates])
}

# Derived changes reported in the Results: denser minus sparse Delta_F.
sparse_contrast <- contrast_summary[
  contrast_summary$density_condition == "Sparse",
  c("network", "metric", "threshold", "mean_high_low_difference")
]
names(sparse_contrast)[4L] <- "mean_high_low_difference_sparse"
denser_contrast <- contrast_summary[
  contrast_summary$density_condition == "Denser",
  c("network", "metric", "threshold", "mean_high_low_difference")
]
names(denser_contrast)[4L] <- "mean_high_low_difference_denser"
density_change_summary <- merge(
  sparse_contrast,
  denser_contrast,
  by = c("network", "metric", "threshold"),
  sort = FALSE
)
density_change_summary$denser_minus_sparse_difference <-
  density_change_summary$mean_high_low_difference_denser -
    density_change_summary$mean_high_low_difference_sparse

# Derived most-minus-average and average-minus-least final-cascade gaps.
target_mean_columns <- c(
  "density_condition", "network", "metric", "threshold",
  "target", "mean_final_active"
)
target_means <- cascade_summary[, target_mean_columns]
high_means <- target_means[target_means$target == "High", -5L]
average_means <- target_means[target_means$target == "Average", -5L]
low_means <- target_means[target_means$target == "Low", -5L]
names(high_means)[5L] <- "mean_final_active_most"
names(average_means)[5L] <- "mean_final_active_average"
names(low_means)[5L] <- "mean_final_active_least"
target_position_gaps <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c("density_condition", "network", "metric", "threshold"),
    sort = FALSE
  ),
  list(high_means, average_means, low_means)
)
target_position_gaps$most_minus_average <-
  target_position_gaps$mean_final_active_most -
    target_position_gaps$mean_final_active_average
target_position_gaps$average_minus_least <-
  target_position_gaps$mean_final_active_average -
    target_position_gaps$mean_final_active_least

topology_summary <- do.call(rbind, lapply(
  split_by_columns(topology_raw, c("density_condition", "network")),
  function(cell) data.frame(
    density_condition = cell$density_condition[1L],
    network = cell$network[1L],
    n_replications = nrow(cell),
    mean_edges = mean(cell$n_edges),
    mean_density = mean(cell$density),
    mean_degree = mean(cell$mean_degree),
    mean_degree_variance = mean(cell$degree_variance),
    sd_degree_variance = safe_sd(cell$degree_variance),
    mean_maximum_degree = mean(cell$maximum_degree),
    mean_global_clustering = mean(cell$global_clustering),
    mean_path_length = mean(cell$mean_path_length),
    percent_connected = 100 * mean(cell$connected),
    mean_components = mean(cell$n_components),
    mean_isolates = mean(cell$n_isolates),
    mean_largest_component = mean(cell$largest_component),
    stringsAsFactors = FALSE
  )
))
rownames(topology_summary) <- NULL

target_summary <- do.call(rbind, lapply(
  split_by_columns(target_raw, c("density_condition", "network", "metric")),
  function(cell) data.frame(
    density_condition = cell$density_condition[1L],
    network = cell$network[1L],
    metric = cell$metric[1L],
    n_replications = nrow(cell),
    mean_high_ties = mean(cell$high_ties),
    mean_average_ties = mean(cell$average_ties),
    mean_low_ties = mean(cell$low_ties),
    percent_high_tied = 100 * mean(cell$high_ties > 1L),
    percent_average_tied = 100 * mean(cell$average_ties > 1L),
    percent_low_tied = 100 * mean(cell$low_ties > 1L),
    percent_high_average_same = 100 * mean(cell$high_average_same),
    percent_high_low_same = 100 * mean(cell$high_low_same),
    percent_average_low_same = 100 * mean(cell$average_low_same),
    mean_centrality_sd = mean(cell$centrality_sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
))
rownames(target_summary) <- NULL

target_agreement_summary <- do.call(rbind, lapply(
  split_by_columns(
    target_agreement_raw,
    c("density_condition", "network", "metric_1", "metric_2")
  ),
  function(cell) data.frame(
    density_condition = cell$density_condition[1L],
    network = cell$network[1L],
    metric_1 = cell$metric_1[1L],
    metric_2 = cell$metric_2[1L],
    n_replications = nrow(cell),
    percent_high_agreement = 100 * mean(cell$high_agreement),
    percent_average_agreement = 100 * mean(cell$average_agreement),
    percent_low_agreement = 100 * mean(cell$low_agreement),
    stringsAsFactors = FALSE
  )
))
rownames(target_agreement_summary) <- NULL

weight_summary <- do.call(rbind, lapply(
  split_by_columns(weight_raw, c("density_condition", "network")),
  function(cell) {
    total_n <- sum(cell$n_weights)
    total_sum <- sum(cell$weight_sum)
    total_sum_squares <- sum(cell$weight_sum_squares)
    pooled_mean <- total_sum / total_n
    pooled_variance <-
      (total_sum_squares - total_sum^2 / total_n) / (total_n - 1L)

    data.frame(
      density_condition = cell$density_condition[1L],
      network = cell$network[1L],
      n_weights = total_n,
      mean_weight = pooled_mean,
      sd_weight = sqrt(max(0, pooled_variance)),
      minimum_weight = min(cell$minimum_weight),
      maximum_weight = max(cell$maximum_weight),
      mean_within_graph_sd = mean(cell$within_graph_sd),
      percent_scaled = 100 * mean(cell$scale_factor < 1),
      mean_scale_factor = mean(cell$scale_factor),
      minimum_scale_factor = min(cell$scale_factor),
      stringsAsFactors = FALSE
    )
  }
))
rownames(weight_summary) <- NULL

# Table 3 intentionally contains only the least and most central targets.
table3_complete_failure_duration <- cascade_summary[
  cascade_summary$target %in% c("Low", "High"),
  c(
    "density_condition", "network", "metric", "target", "threshold",
    "mean_steps_to_complete_failure", "sd_steps_to_complete_failure",
    "n_complete_failures", "percent_complete_failures", "n_replications"
  )
]

names(table3_complete_failure_duration)[
  names(table3_complete_failure_duration) == "target"
] <- "target_centrality"

# Table 3 shading: minimum unrounded conditional mean across graph families.
# Cells without a complete cascade are undefined and cannot be flagged.
table3_complete_failure_duration$is_table3_fastest <- FALSE
table3_groups <- split_by_columns(
  table3_complete_failure_duration,
  c("density_condition", "metric", "target_centrality", "threshold")
)
for (group in table3_groups) {
  group_key <- with(
    group,
    paste(density_condition, metric, target_centrality, threshold, sep = "__")
  )[1L]
  full_key <- with(
    table3_complete_failure_duration,
    paste(density_condition, metric, target_centrality, threshold, sep = "__")
  )
  candidates <- which(full_key == group_key)
  estimable <- candidates[
    !is.na(
      table3_complete_failure_duration$mean_steps_to_complete_failure[candidates]
    )
  ]
  if (length(estimable)) {
    fastest <- min(
      table3_complete_failure_duration$mean_steps_to_complete_failure[estimable]
    )
    table3_complete_failure_duration$is_table3_fastest[estimable] <-
      table3_complete_failure_duration$mean_steps_to_complete_failure[estimable] ==
        fastest
  }
}

# Density comparison for the conditional Table 3 outcome. These means can be
# based on different completed-cascade subsets, so counts and percentages are
# retained for both densities.
table3_density_fields <- c(
  "network", "metric", "target_centrality", "threshold",
  "mean_steps_to_complete_failure", "n_complete_failures",
  "percent_complete_failures"
)
table3_sparse <- table3_complete_failure_duration[
  table3_complete_failure_duration$density_condition == "Sparse",
  table3_density_fields
]
table3_denser <- table3_complete_failure_duration[
  table3_complete_failure_duration$density_condition == "Denser",
  table3_density_fields
]
names(table3_sparse)[5:7] <- paste0(names(table3_sparse)[5:7], "_sparse")
names(table3_denser)[5:7] <- paste0(names(table3_denser)[5:7], "_denser")
conditional_duration_density_comparison <- merge(
  table3_sparse,
  table3_denser,
  by = c("network", "metric", "target_centrality", "threshold"),
  sort = FALSE
)
conditional_duration_density_comparison$sparse_minus_denser_steps <-
  conditional_duration_density_comparison$mean_steps_to_complete_failure_sparse -
    conditional_duration_density_comparison$mean_steps_to_complete_failure_denser

# =============================================================================
# 6. WRITE AUDITABLE DATA FILES
# =============================================================================

write.csv(
  cascade_raw,
  file.path(output_directory, "cascade_replication_level.csv"),
  row.names = FALSE
)

write.csv(
  cascade_summary,
  file.path(output_directory, "cascade_summary_all_conditions.csv"),
  row.names = FALSE
)
write.csv(
  paired_results,
  file.path(output_directory, "high_low_difference_replication_level.csv"),
  row.names = FALSE
)
write.csv(
  contrast_summary,
  file.path(output_directory, "table2_change_in_failure_rate.csv"),
  row.names = FALSE
)
write.csv(
  table3_complete_failure_duration,
  file.path(output_directory, "table3_steps_to_complete_failure.csv"),
  row.names = FALSE
)
write.csv(
  density_change_summary,
  file.path(output_directory, "density_change_in_high_low_difference.csv"),
  row.names = FALSE
)
write.csv(
  target_position_gaps,
  file.path(output_directory, "target_position_failure_rate_gaps.csv"),
  row.names = FALSE
)
write.csv(
  conditional_duration_density_comparison,
  file.path(output_directory, "conditional_duration_density_comparison.csv"),
  row.names = FALSE
)
write.csv(
  topology_raw,
  file.path(output_directory, "topology_diagnostics_replication_level.csv"),
  row.names = FALSE
)
write.csv(
  topology_summary,
  file.path(output_directory, "topology_diagnostics_summary.csv"),
  row.names = FALSE
)
write.csv(
  target_raw,
  file.path(output_directory, "target_selection_replication_level.csv"),
  row.names = FALSE
)
write.csv(
  target_summary,
  file.path(output_directory, "target_selection_summary.csv"),
  row.names = FALSE
)
write.csv(
  target_agreement_raw,
  file.path(output_directory, "target_agreement_replication_level.csv"),
  row.names = FALSE
)
write.csv(
  target_agreement_summary,
  file.path(output_directory, "target_agreement_summary.csv"),
  row.names = FALSE
)
write.csv(
  weight_summary,
  file.path(output_directory, "edge_weight_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  weight_raw,
  file.path(output_directory, "edge_weight_diagnostics_replication_level.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    configuration = configuration_for_workers,
    cascade_summary = cascade_summary,
    table2 = contrast_summary,
    table3 = table3_complete_failure_duration,
    density_change_summary = density_change_summary,
    target_position_gaps = target_position_gaps,
    conditional_duration_density_comparison =
      conditional_duration_density_comparison,
    topology_summary = topology_summary,
    target_summary = target_summary,
    target_agreement_summary = target_agreement_summary,
    edge_weight_summary = weight_summary
  ),
  file.path(output_directory, "manuscript_simulation_results.rds")
)

# =============================================================================
# 7. MANUSCRIPT FIGURES 5-9
# =============================================================================

if (MAKE_FIGURES) {
  plot_summary <- cascade_summary
  plot_summary$target <- factor(
    plot_summary$target,
    levels = target_order
  )
  plot_summary$network <- factor(
    plot_summary$network,
    levels = c("Erdos-Renyi", "Barabasi-Albert", "Watts-Strogatz")
  )
  plot_summary$density_label <- factor(
    plot_summary$density_condition,
    levels = c("Sparse", "Denser"),
    labels = c(
      "Sparse (20 edges; ~10% density)",
      "Denser (40 edges; ~20% density)"
    )
  )
  plot_summary$threshold_label <- factor(
    plot_summary$threshold,
    levels = threshold_codes,
    labels = unname(threshold_display[threshold_codes])
  )

  metric_titles <- c(
    degree = "degree",
    strength = "strength",
    closeness = "harmonic closeness",
    betweenness = "betweenness"
  )
  metric_figure_numbers <- c(
    degree = "05",
    strength = "06",
    closeness = "07",
    betweenness = "08"
  )

  for (metric in metrics) {
    metric_data <- droplevels(plot_summary[plot_summary$metric == metric, ])
    plot_metric <- ggplot(
      metric_data,
      aes(
        target,
        mean_final_active,
        shape = network,
        color = network,
        group = network
      )
    ) +
      geom_errorbar(
        aes(
          ymin = ci95_lo_final_active,
          ymax = ci95_hi_final_active
        ),
        width = 0.15,
        linewidth = 0.4,
        position = position_dodge(width = 0.4)
      ) +
      geom_line(
        linewidth = 0.6,
        position = position_dodge(width = 0.4)
      ) +
      geom_point(
        size = 3,
        position = position_dodge(width = 0.4)
      ) +
      scale_y_continuous(
        "Mean final active proportion",
        limits = c(0, 1)
      ) +
      scale_x_discrete(NULL) +
      scale_color_brewer("Graph family", palette = "Set1") +
      scale_shape_manual(
        "Graph family",
        values = c(
          "Erdos-Renyi" = 16,
          "Barabasi-Albert" = 17,
          "Watts-Strogatz" = 15
        )
      ) +
      facet_grid(density_label ~ threshold_label) +
      labs(title = paste("Cascade size by", metric_titles[[metric]], "target selection")) +
      theme_minimal(base_size = 13) +
      theme(
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )

    metric_file <- sprintf(
      "figure%s_failure_rate_%s.png",
      metric_figure_numbers[[metric]],
      metric
    )
    ggsave(
      file.path(output_directory, metric_file),
      plot_metric,
      width = 14,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }

  contrast_plot_data <- contrast_summary
  contrast_plot_data$threshold_label <- factor(
    contrast_plot_data$threshold,
    levels = threshold_codes,
    labels = unname(threshold_display[threshold_codes])
  )
  contrast_plot_data$density_label <- factor(
    contrast_plot_data$density_condition,
    levels = c("Sparse", "Denser"),
    labels = c(
      "Sparse (20 edges; ~10% density)",
      "Denser (40 edges; ~20% density)"
    )
  )
  contrast_plot_data$metric_label <- factor(
    contrast_plot_data$metric,
    levels = metrics,
    labels = c("Degree", "Strength", "Harmonic closeness", "Betweenness")
  )
  contrast_plot_data$network <- factor(
    contrast_plot_data$network,
    levels = c("Erdos-Renyi", "Barabasi-Albert", "Watts-Strogatz")
  )

  plot_difference <- ggplot(
    contrast_plot_data,
    aes(
      threshold_label,
      mean_high_low_difference,
      color = network,
      group = network
    )
  ) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_errorbar(
      aes(
        ymin = ci95_lo_high_low_difference,
        ymax = ci95_hi_high_low_difference
      ),
      width = 0.10,
      linewidth = 0.4,
      position = position_dodge(width = 0.25)
    ) +
    geom_line(
      linewidth = 0.7,
      position = position_dodge(width = 0.25)
    ) +
    geom_point(
      size = 2.4,
      position = position_dodge(width = 0.25)
    ) +
    scale_color_brewer("Graph family", palette = "Set1") +
    scale_y_continuous(
      expression("Mean change in failure rate (" * Delta[F] * ")")
    ) +
    scale_x_discrete(NULL) +
    facet_grid(density_label ~ metric_label) +
    labs(
      title = expression(
        "Change in failure rate (" * Delta[F] * ") by density"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )

  ggsave(
    file.path(output_directory, "figure09_change_in_failure_rate.png"),
    plot_difference,
    width = 15,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# =============================================================================
# 8. REPRESENTATIVE TOPOLOGY FIGURE
# =============================================================================

if (MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE) {
  representative_features <- c(
    "degree_variance",
    "maximum_degree",
    "global_clustering",
    "mean_path_length",
    "n_components",
    "n_isolates",
    "largest_component"
  )

  choose_representative <- function(cell) {
    modal_connected <- as.integer(mean(cell$connected) >= 0.5)
    candidates <- cell[cell$connected == modal_connected, , drop = FALSE]
    values <- candidates[, representative_features, drop = FALSE]
    centers <- vapply(values, median, numeric(1L), na.rm = TRUE)
    scales <- vapply(values, stats::mad, numeric(1L), na.rm = TRUE)
    fallback <- vapply(values, stats::sd, numeric(1L), na.rm = TRUE)
    invalid_scale <- !is.finite(scales) | scales <= 0
    scales[invalid_scale] <- fallback[invalid_scale]
    scales[!is.finite(scales) | scales <= 0] <- 1
    standardized <- sweep(as.matrix(values), 2L, centers, "-")
    standardized <- sweep(standardized, 2L, scales, "/")
    candidates$representativeness_score <- rowSums(standardized^2)
    candidates$modal_connected <- modal_connected
    candidates <- candidates[order(
      candidates$representativeness_score,
      candidates$replication
    ), , drop = FALSE]
    candidates[1L, , drop = FALSE]
  }

  selected_representatives <- do.call(rbind, lapply(
    split_by_columns(topology_raw, c("density_condition", "network")),
    choose_representative
  ))
  rownames(selected_representatives) <- NULL

  graph_diagnostics_for_check <- function(graph) {
    graph_degree <- degree(graph)
    graph_components <- components(graph)
    c(
      n_edges = gsize(graph),
      degree_variance = stats::var(graph_degree),
      maximum_degree = max(graph_degree),
      global_clustering = transitivity(
        graph,
        type = "global",
        isolates = "zero"
      ),
      mean_path_length = mean_distance(
        graph,
        directed = FALSE,
        unconnected = TRUE,
        weights = NA
      ),
      n_components = graph_components$no,
      n_isolates = sum(graph_degree == 0),
      largest_component = max(graph_components$csize),
      connected = as.integer(graph_components$no == 1L)
    )
  }

  regenerate_topologies <- function(replication) {
    set.seed(design$base_seed + as.integer(replication))
    list(
      Sparse = generate_topologies(
        design$n_nodes,
        design$edge_counts[["Sparse"]],
        design$ws_rewiring
      ),
      Denser = generate_topologies(
        design$n_nodes,
        design$edge_counts[["Denser"]],
        design$ws_rewiring
      )
    )
  }

  network_codes <- c(
    "Erdos-Renyi" = "ER",
    "Barabasi-Albert" = "BA",
    "Watts-Strogatz" = "WS"
  )
  representative_graphs <- list()

  for (row_index in seq_len(nrow(selected_representatives))) {
    selected_row <- selected_representatives[row_index, , drop = FALSE]
    regenerated <- regenerate_topologies(selected_row$replication)
    network_code <- network_codes[[selected_row$network]]
    graph <- regenerated[[selected_row$density_condition]][[network_code]]
    observed <- graph_diagnostics_for_check(graph)
    expected <- unlist(selected_row[names(observed)], use.names = TRUE)
    if (!isTRUE(all.equal(
      as.numeric(observed),
      as.numeric(expected),
      tolerance = 1e-12
    ))) {
      stop("A representative graph did not match its stored diagnostics.")
    }
    representative_graphs[[paste(
      selected_row$density_condition,
      network_code,
      sep = "_"
    )]] <- graph
  }

  write.csv(
    selected_representatives,
    file.path(output_directory, "representative_topology_selection.csv"),
    row.names = FALSE
  )

  topology_colors <- c(
    "Barabasi-Albert" = "#D95F59",
    "Watts-Strogatz" = "#4C9F70",
    "Erdos-Renyi" = "#4C78A8"
  )

  graph_layout <- function(graph, network, replication) {
    if (network == "Watts-Strogatz") {
      return(layout_in_circle(graph, order = seq_len(vcount(graph))))
    }
    set.seed(900000L + as.integer(replication))
    layout_with_fr(graph, niter = 1500L, grid = "nogrid")
  }

  plot_representative_graph <- function(graph, selected_row) {
    network <- selected_row$network
    density <- selected_row$density_condition
    coordinates <- graph_layout(
      graph,
      network,
      selected_row$replication
    )
    coordinates <- norm_coords(
      coordinates,
      xmin = -1,
      xmax = 1,
      ymin = -1,
      ymax = 1
    )
    graph_degree <- degree(graph)

    if (network == "Watts-Strogatz") {
      edge_list <- as_edgelist(graph, names = FALSE)
      ordinary_distance <- abs(edge_list[, 1L] - edge_list[, 2L])
      circular_distance <- pmin(
        ordinary_distance,
        vcount(graph) - ordinary_distance
      )
      lattice_radius <- if (density == "Sparse") 1L else 2L
      long_range <- circular_distance > lattice_radius
      edge_color <- ifelse(long_range, "#E6862A", "#A7AFB8A6")
      edge_width <- ifelse(long_range, 2.25, 1.15)
    } else {
      edge_color <- rep("#9DA7B3A6", gsize(graph))
      edge_width <- rep(1.15, gsize(graph))
    }

    plot(
      graph,
      layout = coordinates,
      vertex.size = 9 + 2 * sqrt(graph_degree),
      vertex.color = adjustcolor(topology_colors[[network]], alpha.f = 0.92),
      vertex.frame.color = "white",
      vertex.label = NA,
      edge.color = edge_color,
      edge.width = edge_width,
      edge.curved = 0,
      margin = 0.17,
      rescale = FALSE,
      asp = 1
    )
    title(main = network, cex.main = 1.15, font.main = 2, col.main = "#17324D")
    mtext(
      sprintf(
        "rep %d | m=%d | C=%.3f | L=%.3f",
        selected_row$replication,
        selected_row$n_edges,
        selected_row$global_clustering,
        selected_row$mean_path_length
      ),
      side = 1,
      line = 0.35,
      cex = 0.70,
      col = "#5F6B76"
    )
  }

  draw_topology_panel <- function() {
    network_order <- c("Barabasi-Albert", "Watts-Strogatz", "Erdos-Renyi")
    density_order <- c("Sparse", "Denser")
    layout(matrix(seq_len(6L), nrow = 2L, byrow = TRUE))
    par(
      mar = c(2.5, 2.2, 3.0, 1.2),
      oma = c(4.2, 4.4, 3.3, 1.0),
      xpd = NA,
      bg = "white",
      family = "sans"
    )

    for (density in density_order) {
      for (network in network_order) {
        selected_row <- selected_representatives[
          selected_representatives$density_condition == density &
            selected_representatives$network == network,
          , drop = FALSE
        ]
        graph_key <- paste(density, network_codes[[network]], sep = "_")
        plot_representative_graph(
          representative_graphs[[graph_key]],
          selected_row
        )
        if (network == network_order[1L]) {
          density_label <- if (density == "Sparse") {
            "20 edges\ndensity 0.105"
          } else {
            "40 edges\ndensity 0.211"
          }
          mtext(
            density_label,
            side = 2,
            line = 1.2,
            cex = 0.87,
            font = 2,
            col = "#17324D"
          )
        }
      }
    }

    mtext(
      "Representative simulated graph topologies",
      side = 3,
      outer = TRUE,
      line = 1.45,
      cex = 1.45,
      font = 2,
      col = "#17324D"
    )
    mtext(
      paste(
        "Node size is proportional to degree. Each graph is nearest its cell",
        "median across recorded topology diagnostics."
      ),
      side = 1,
      outer = TRUE,
      line = 2.0,
      cex = 0.82,
      col = "#5F6B76"
    )
  }

  png(
    file.path(output_directory, "figure04_representative_topologies.png"),
    width = 3900,
    height = 2600,
    res = 320,
    bg = "white",
    type = if (capabilities("cairo")) "cairo-png" else getOption("bitmapType")
  )
  draw_topology_panel()
  dev.off()
}

# =============================================================================
# 9. REPRODUCIBILITY RECORD
# =============================================================================

configuration_lines <- c(
  sprintf("Study design: %d", STUDY_DESIGN),
  sprintf("Threshold condition: %d", THRESHOLD_CONDITION),
  sprintf("Study label: %s", design$label),
  sprintf("Completed (UTC): %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("Base seed: %d", design$base_seed),
  sprintf("Replications per threshold range: %d", design$n_replications),
  sprintf("Workers used: %d", workers),
  sprintf("Nodes: %d", design$n_nodes),
  sprintf(
    "Edges: Sparse=%d; Denser=%d",
    design$edge_counts[["Sparse"]],
    design$edge_counts[["Denser"]]
  ),
  sprintf("Watts-Strogatz rewiring probability: %.2f", design$ws_rewiring),
  sprintf("Metrics: %s", paste(metrics, collapse = ", ")),
  sprintf("Threshold ranges: %s", paste(threshold_codes, collapse = ", ")),
  sprintf("Validated R version: %s", expected_versions[["R"]]),
  sprintf("Validated igraph version: %s", expected_versions[["igraph"]]),
  sprintf("RNG kind: %s", paste(RNGkind(), collapse = ", ")),
  sprintf("Source path: %s", script_path),
  sprintf("Output directory: %s", normalizePath(output_directory, winslash = "/")),
  sprintf("Command-line options: %s", paste(cli_arguments, collapse = " ")),
  sprintf("MAKE_FIGURES: %s", MAKE_FIGURES),
  sprintf("MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE: %s",
          MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE),
  sprintf("Run script MD5: %s", unname(tools::md5sum(script_path))),
  "Simulation source: Manuscript_Simulation_Reproduction_CompleteFailure_20260817.R",
  paste0(
    "Simulation source SHA256: ",
    "3763C45645370F691758B15A00729B265A85DAF3FDAA6123758BB88993B06CC0"
  ),
  paste0(
    "Primary-study archived backbone SHA256: ",
    "E9CC833FB8D35F399BBA578C9B33D812EDE6AC108D29DF9EA544C927EDF59D74"
  ),
  paste0(
    "Primary-study archived beta=.10 runner SHA256: ",
    "9E5EA48D317CE8F801BF602A7EA11CF4C4C8E72F32EE6E7E1D755FA9E7BED52D"
  ),
  paste0(
    "Primary-study archived result RDS SHA256: ",
    "22E42EB105DDB4709793D854F4C085D6A50E0770570B5604E0A98CCD536B4EF6"
  ),
  "",
  capture.output(sessionInfo())
)
writeLines(
  configuration_lines,
  file.path(output_directory, "reproducibility_record.txt")
)

file.copy(
  script_path,
  file.path(output_directory, basename(script_path)),
  overwrite = TRUE
)

generated_output_names <- c(
  "cascade_replication_level.csv",
  "cascade_summary_all_conditions.csv",
  "high_low_difference_replication_level.csv",
  "table2_change_in_failure_rate.csv",
  "table3_steps_to_complete_failure.csv",
  "density_change_in_high_low_difference.csv",
  "target_position_failure_rate_gaps.csv",
  "conditional_duration_density_comparison.csv",
  "topology_diagnostics_replication_level.csv",
  "topology_diagnostics_summary.csv",
  "target_selection_replication_level.csv",
  "target_selection_summary.csv",
  "target_agreement_replication_level.csv",
  "target_agreement_summary.csv",
  "edge_weight_diagnostics.csv",
  "edge_weight_diagnostics_replication_level.csv",
  "manuscript_simulation_results.rds",
  "reproducibility_record.txt",
  basename(script_path)
)
if (MAKE_FIGURES) {
  generated_output_names <- c(
    generated_output_names,
    sprintf("figure%02d_failure_rate_%s.png", 5:8, metrics),
    "figure09_change_in_failure_rate.png"
  )
}
if (MAKE_REPRESENTATIVE_TOPOLOGY_FIGURE) {
  generated_output_names <- c(
    generated_output_names,
    "representative_topology_selection.csv",
    "figure04_representative_topologies.png"
  )
}
output_files <- file.path(output_directory, generated_output_names)
if (any(!file.exists(output_files))) {
  stop(
    "Expected output files were not generated: ",
    paste(basename(output_files[!file.exists(output_files)]), collapse = ", ")
  )
}
output_hashes <- data.frame(
  file = basename(output_files),
  md5 = unname(tools::md5sum(output_files)),
  stringsAsFactors = FALSE
)
write.csv(
  output_hashes,
  file.path(output_directory, "output_file_md5.csv"),
  row.names = FALSE
)

cat("\nSimulation and all requested summaries completed successfully.\n")
cat(sprintf("Primary output: %s\n", normalizePath(output_directory)))
