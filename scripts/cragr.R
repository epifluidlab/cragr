

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
  if (is.null(get0("script_args")))
    stop_quietly()

  # Build a string to record all parameters
  param_str <-
    paste(names(script_args),
          script_args,
          sep = "=",
          collapse = "; ")
} else {
  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("--input")),
      optparse::make_option(c("--output-dir")),
      optparse::make_option(c("--sample-id")),
      optparse::make_option(c("--gc")),
      optparse::make_option(c("--mappability")),
      optparse::make_option(c("--mappability-threshold"), default = "0.9"),
      optparse::make_option(c("--chrom"),
                            help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used"),
      optparse::make_option(c("--exclude-chrom"),
                            help = "Exclude chromosomes from the analysis. Separated by colons, such as 12:16:X"),
      optparse::make_option(
        c("--min-mapq"),
        type = "integer",
        default = 30L,
        help = "Minimal MAPQ for fragments to take part in the analysis"
      ),
      optparse::make_option(
        c("--min-fraglen"),
        type = "integer",
        default = 50L,
        help = "Minimal length for fragments to take part in the analysis"
      ),
      optparse::make_option(
        c("--max-fraglen"),
        type = "integer",
        default = 1000L,
        help = "Maximal length for fragments to take part in the analysis"
      ),
      optparse::make_option(
        c("--exclude-region"),
        type = "character",
        default = NULL,
        help = "BED files defining regions to be excluded from the analysis, separated by colon"
      ),
      optparse::make_option(
        c("--chrom-sizes"),
        type = "character",
        default = NULL,
        help = "Path to the chrom_sizes file"
      ),
      optparse::make_option(c("--exclude-soft-clipping"), default = FALSE),
      optparse::make_option(c("--window-size"), default = 200L),
      optparse::make_option(c("--step-size"), default = 20L),
      optparse::make_option(c("--cpois"), default = FALSE)
    )
  )
  script_args <-
    optparse::parse_args(
      parser,
      args = commandArgs(trailingOnly = TRUE),
      convert_hyphens_to_underscores = TRUE
    )

  library(tidyverse)
  library(magrittr)
  library(cragr)

  # Process arguments
  if (!("chrom" %in% names(script_args)))
    script_args$chrom <- NULL
  if (is_null(script_args$chrom_sizes))
    script_args$chrom_sizes <- system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr")
  if (is_null(script_args$exclude_region))
    script_args$exclude_region <- system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed", package = "cragr")

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

  # Build a string to record all parameters
  param_str <- paste0("Parameters: ", paste(commandArgs(trailingOnly = TRUE), collapse = " "))
}

logging::loginfo(str_interp("Argument summary:"))
names(script_args) %>%
  walk(function(arg) {
    full_arg <- str_replace_all(arg, "_", "-")
    logging::loginfo(str_interp("--${full_arg}: ${script_args[[arg]]}"))
  })

logging::loginfo(str_interp("Loading fragments: ${script_args$input}"))
if (!is_null(script_args$chrom)) {
  frag <- load_fragments(script_args$input, region = script_args[["chrom"]])
} else {
  frag <- load_fragments(script_args$input)
}

if (!is_null(script_args$exclude_chrom)) {
  frag <- frag[!(chrom %in% script_args$exclude_chrom)]
}

logging::loginfo("Fragments summary:")
print(rbind(frag[, .(count = length(start)), by = chrom], list(chrom = "Total", count = nrow(frag))))

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

logging::loginfo("Calling hotspots ...")
call_peak(ifs, cpois = script_args$cpois)

logging::loginfo("Writing results to disk ...")
if (!dir.exists(script_args$output_dir))
  dir.create(script_args$output_dir)

data.table::setnames(ifs, old = "chrom", new = "#chrom")
data.table::fwrite(
  ifs,
  file = str_interp(
    "${script_args$output_dir}/${script_args$sample_id}.bedGraph.gz"
  ),
  quote = FALSE,
  sep = "\t",
  na = ".",
  scipen = 999,
  col.names = TRUE
)


