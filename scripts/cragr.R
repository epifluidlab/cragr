

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# example
# script_args <- list(
#   input = "~/Downloads/Pilot2-1.hg19.frag.bed.gz",
#   output_dir = here::here("sandbox"),
#   sample_id = "S1",
#   gc = NULL,
#   mappability = here("inst/extdata/mappability.subset.bedGraph.gz"),
#   mappability_threshold = 0.9,
#   chrom = NULL,
#   exclude_chrom = NULL,
#   min_mapq = 30L,
#   min_fraglen = 50L,
#   max_fraglen = 1000L,
#   chrom_sizes = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr"),
#   exclude_region = system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed", package = "cragr"),
#   exclude_soft_clipping = FALSE,
#   window_size = 200L,
#   step_size = 20L
# )

if (interactive()) {
  if (is.null(get0("script_args")) || is.null(get0("subcommand")))
    stop_quietly()

  # Build a string to record all parameters
  param_str <-
    paste(names(script_args),
          script_args,
          sep = "=",
          collapse = "; ")
} else {
  # The first argument is subcommand. Should be one of the following:
  # * main: run the full pipeline, i.e. calculate IFS scores and then call hotspots
  # * ifs: only calculate IFS scores
  # * hotspot: only call hotspots from a existing IFS score track
  args <- commandArgs(trailingOnly = TRUE)
  subcommand <- args[1]
  args <- args[-1]

  if (!subcommand %in% c("main", "ifs", "hotspot"))
    stop("Subcommand should be one of the following: main, ifs, and hotspot")

  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("--input"), help = "Path to the inpuut file"),
      # Can be either "frag" or "ifs"
      optparse::make_option(
        c("--input-type"),
        type = "character",
        default = "frag",
        help = "Type of the input. 'frag' for fragment data file, and 'ifs' for pre-calculated IFS scores. Default is 'frag'"
      ),
      # optparse::make_option(c("--output-dir"), type = "character", help = "Directory for output files"),
      optparse::make_option(c("--prefix"), type = "character", help = "Prefix of the output file names"),
      optparse::make_option(
        c("--gc"),
        type = "character",
        default = NULL,
        help = "Path to the GC% file. Default is NULL, i.e. do NOT perform GC correction"
      ),
      optparse::make_option(
        c("--mappability"),
        type = "character",
        default = NULL,
        help = "Path to the mappability file. Default is NULL, i.e. do NOT exclude fragments from low-mappability regions"
      ),
      optparse::make_option(
        c("--mappability-threshold"),
        default = 0.9,
        help = "Threshold for low-mappability filtering"
      ),
      optparse::make_option(
        c("--chrom"),
        type = "character",
        default = NULL,
        help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used"
      ),
      optparse::make_option(
        c("--exclude-chrom"),
        type = "character",
        default = NULL,
        help = "Exclude chromosomes from the analysis. Separated by colons, such as 12:16:X"
      ),
      optparse::make_option(c("--min-mapq"),
                            default = 30L,
                            help = "Minimal MAPQ for fragments included in the analysis"),
      optparse::make_option(c("--min-fraglen"),
                            default = 50L,
                            help = "Minimal length for fragments included in the analysis"),
      optparse::make_option(c("--max-fraglen"),
                            default = 1000L,
                            help = "Maximal length for fragments included in the analysis"),
      optparse::make_option(
        c("--exclude-region"),
        type = "character",
        default = system.file(
          "extdata",
          "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed",
          package = "cragr"
        ),
        help = "BED files defining regions to be excluded from the analysis, separated by colon. Default is the ENCODE Blacklist: https://www.nature.com/articles/s41598-019-45839-z, which is included in this R package"
      ),
      optparse::make_option(
        c("--chrom-sizes"),
        type = "character",
        default = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr"),
        help = "Path to the chrom_sizes file. Default is GRCh37 hg19, which is included in this R package"
      ),
      optparse::make_option(
        c("--exclude-soft-clipping"),
        action = "store_true",
        default = FALSE,
        help = "Exclude fragments with leading soft-clipping from the analysis"
      ),
      optparse::make_option(c("--window-size"), default = 200L, help = "Size of the sliding window. Default is 200"),
      optparse::make_option(c("--step-size"), default = 20L, help = "Step size of the sliding window. Default is 20"),
      optparse::make_option(
        c("--cpois"),
        action = "store_true",
        default = FALSE,
        help = "Use continuous Poisson model to call hotspots"
      ),
      optparse::make_option(c("--fdr"), default = 0.01, help = "FDR cut-off value used in hotspot calling. Default is 0.01"),
      optparse::make_option(
        c("--merge-distance"),
        action = "store",
        type = "integer",
        default = NULL,
        help = "During hotspot calling, two hotspots with distance smaller than this will be merged together. If not specified, the sliding window size will by used, i.e. 200bp by default"
      )
    )
  )
  script_args <-
    optparse::parse_args(
      parser,
      args = args,
      convert_hyphens_to_underscores = TRUE
    )

  library(tidyverse)
  library(magrittr)
  library(cragr)

  # Process arguments
  if (!("chrom" %in% names(script_args)))
    script_args$chrom <- NULL
  # if (is_null(script_args$chrom_sizes))
  #   script_args$chrom_sizes <- system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr")
  # if (is_null(script_args$exclude_region))
  #   script_args$exclude_region <- system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed", package = "cragr")

  c("chrom", "exclude_chrom", "exclude_region") %>%
    walk(function(arg) {
      if (!is.null(script_args[[arg]])) {
        script_args[[arg]] <<-
          str_split(script_args[[arg]], pattern = ":")[[1]]
      }
    })

  c("min_mapq", "min_fraglen", "max_fraglen", "window_size", "step_size") %>%
    walk(function(arg) {
      if (!is_null(script_args))
        script_args[[arg]] <<- as.numeric(script_args[[arg]])
    })

  if (!is_null(script_args$mappability_threshold))
    script_args$mappability_threshold <- as.numeric(script_args$mappability_threshold)

  if (!script_args$input_type %in% c("frag", "ifs")) {
    stop("--input-type should be either frag or ifs.")
  }

  if (is_null(script_args$merge_distance))
    script_args$merge_distance <- script_args$window_size

  if (script_args$window_size %% script_args$step_size != 0)
    stop("window_size must be multiples of step_size")

  # Build a string to record all parameters
  param_str <- paste0("Parameters: ", paste(commandArgs(trailingOnly = TRUE), collapse = " "))
}

logging::loginfo(str_interp("Argument summary:"))
names(script_args) %>%
  walk(function(arg) {
    full_arg <- str_replace_all(arg, "_", "-")
    logging::loginfo(str_interp("--${full_arg}: ${script_args[[arg]]}"))
  })


if (script_args$input_type == "frag") {
  logging::loginfo("Input: fragment data")

  logging::loginfo(str_interp("Loading fragments: ${script_args$input}"))
  if (!is_null(script_args$chrom)) {
    frag <-
      load_fragments(script_args$input, region = script_args[["chrom"]])
  } else {
    frag <- load_fragments(script_args$input)
  }

  if (!is_null(script_args$exclude_chrom)) {
    frag <- frag[!(chrom %in% script_args$exclude_chrom)]
  }



  logging::loginfo("Fragments summary:")
  print(rbind(frag[, .(count = length(start)), by = chrom], list(
    chrom = "Total", count = nrow(frag)
  )))

  logging::loginfo("Calculating IFS scores ...")
  ifs <- ifs_score(
    frag,
    window_size = script_args$window_size,
    step_size = script_args$step_size,
    gc = script_args$gc,
    blacklist_region = script_args$exclude_region,
    mappability_region = script_args$mappability,
    mappability_threshold = script_args$mappability_threshold,
    chrom_sizes = script_args$chrom_sizes,
    min_mapq = script_args$min_mapq,
    min_fraglen = script_args$min_fraglen,
    max_fraglen = script_args$max_fraglen,
    exclude_soft_clipping = script_args$exclude_soft_clipping
  )
  rm(frag)
} else if (script_args$input_type == "ifs") {
  logging::loginfo("Input: IFS scores")

  ifs <- load_bed(script_args$input)
} else {
  stop(str_interp("Invalid input type: ${script_args$input_type}"))
}

logging::loginfo("Calculating global p-values ...")
ifs <- calc_pois_pval(ifs, cpois = script_args$cpois)

logging::loginfo("Calculating local p-values ...")
ifs <-
  calc_pois_pval_local(
    ifs,
    cpois = script_args$cpois,
    window_size = script_args$window_size,
    step_size = script_args$step_size,
    local_layout = list(`5k` = 5e3L, `10k` = 10e3L)
  )

logging::loginfo("Calling hotspots ...")

pval_local_cutoff <- 1e-5

if (script_args$cpois) {
  hotspot_cpois <-
    call_hotspot(
      ifs,
      use_cpois = TRUE,
      fdr_cutoff = script_args$fdr,
      local_pval_cutoff = pval_local_cutoff,
      merge_distance = script_args$merge_distance
    )
}

hotspot_standard <-
  call_hotspot(
    ifs,
    use_cpois = FALSE,
    fdr_cutoff = script_args$fdr,
    local_pval_cutoff = pval_local_cutoff,
    merge_distance = script_args$merge_distance
  )

logging::loginfo("Writing results to disk ...")
output_dir <- dirname(script_args$prefix)
if (!dir.exists(output_dir))
  dir.create(output_dir)


# This is ifs. Each row corresponds to a 200bp region, and they are overlapping.
# We're more interested in a BEDGRAPH-style non-overlapping track
#    chrom    start      end    score cov mappability    gc   score0 cov_corrected      pval pval_adjust pval_cpois pval_cpois_adjust
# 1:    21  9412520  9412720 27.58913   9   0.9050833 0.295 17.90432      13.78962 0.9599079           1  0.9611725                 1
# 2:    21  9412540  9412740 29.36389  10   0.9750000 0.290 19.67908      14.78962 0.9841276           1  0.9828842                 1
# 3:    21  9412560  9412760 27.48691   9   0.9417500 0.270 17.80209      13.78962 0.9599079           1  0.9594159                 1
# 4:    21 10872940 10873140 71.61589  31   0.9162698 0.470 61.93108      35.78962 1.0000000           1  1.0000000                 1
# 5:    21 10873160 10873360 73.56209  33   0.9200000 0.510 63.87727      37.78962 1.0000000           1  1.0000000                 1
#
# After this transformation, the IFS track looks like this:
#    chrom    start      end    score cov mappability    gc   score0 cov_corrected      pval pval_adjust pval_cpois pval_cpois_adjust
# 1:    21  9412610  9412630 27.58913   9   0.9050833 0.295 17.90432      13.78962 0.9599079           1  0.9611725                 1
# 2:    21  9412630  9412650 29.36389  10   0.9750000 0.290 19.67908      14.78962 0.9841276           1  0.9828842                 1
# 3:    21  9412650  9412670 27.48691   9   0.9417500 0.270 17.80209      13.78962 0.9599079           1  0.9594159                 1
# 4:    21 10873030 10873050 71.61589  31   0.9162698 0.470 61.93108      35.78962 1.0000000           1  1.0000000                 1
# 5:    21 10873250 10873270 73.56209  33   0.9200000 0.510 63.87727      37.78962 1.0000000           1  1.0000000                 1
ws_ratio <- script_args$window_size %/% script_args$step_size
ifs[, start := as.integer(round(start + (ws_ratio - 1) / 2 * script_args$step_size))][, end := as.integer(start + script_args$step_size)]
write_bed(
  ifs,
  file_path = str_interp(
    "${script_args$prefix}.ifs.bedGraph.gz"
  ),
  create_index = TRUE
)
write_bed(
  hotspot_standard,
  file_path = str_interp(
    "${script_args$prefix}.hotspot.bed.gz"
  ),
  create_index = TRUE
)

if (script_args$cpois) {
  write_bed(
    hotspot_cpois,
    file_path = str_interp("${script_args$prefix}.hotspot.cpois.bed.gz"),
    create_index = TRUE
  )
}
