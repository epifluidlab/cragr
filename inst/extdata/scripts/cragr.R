stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


parse_script_args <- function() {
  if (interactive()) {
    script_args <- get0("script_args")
    subcommand <- get0("subcommand")

    if (is.null(script_args) || is.null(subcommand))
      stop()
    else {
      return(list(subcommand, script_args))
    }
  } else {
    # The first argument is subcommand. Should be one of the following:
    # * main: run the full pipeline, i.e. calculate IFS scores and then call hotspots
    # * ifs: only calculate IFS scores
    # * hotspot: only call hotspots from a existing IFS score track
    args <- commandArgs(trailingOnly = TRUE)
    subcommand <- args[1]
    args <- args[-1]

    if (!subcommand %in% c("ifs", "peak", "hotspot"))
      stop("Subcommand should be one of the following: stage1, stage2")

    # Run in CLI script mode
    parser <- optparse::OptionParser(
      option_list = list(
        optparse::make_option(c("-i", "--input"), help = "Path to the inpuut file. If there are multiple input files, they should be separated by colons"),
        optparse::make_option(c("--output-ifs"), type = "character", help = "Directory for IFS output"),
        optparse::make_option(c("--output-hotspot"), type = "character", help = "Directory for hotspot output"),
        optparse::make_option(c("--genome"), default = "hs37-1kg", help = "Genome of the input"),
        optparse::make_option(
          c("-g", "--gc-correct"),
          default = TRUE,
          help = "Path to the GC% file. Default is NULL, i.e. do NOT perform GC correction"
        ),
        optparse::make_option(
          c("-m", "--high-mappability"),
          type = "character",
          help = "Path to the mappability file. Default is NULL, i.e. do NOT exclude fragments from low-mappability regions"
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
          default = "encode.blacklist.hs37-1kg",
          help = "BED files defining regions to be excluded from the analysis, separated by colon. Default is the ENCODE Blacklist: https://www.nature.com/articles/s41598-019-45839-z, which is included in this R package"
        ),
        optparse::make_option(
          c("--exclude-soft-clipping"),
          action = "store_true",
          default = FALSE,
          help = "Exclude fragments with leading soft-clipping from the analysis"
        ),
        optparse::make_option(
          c("-w", "--window-size"),
          default = 200L,
          help = "Size of the sliding window. Default is 200"
        ),
        optparse::make_option(c("-s", "--step-size"), default = 20L, help = "Step size of the sliding window. Default is 20"),
        # optparse::make_option(
        #   c("--cpois"),
        #   action = "store_true",
        #   default = FALSE,
        #   help = "Use continuous Poisson model to call hotspots"
        # ),
        optparse::make_option(c("--fdr"), default = 0.2, help = "FDR cut-off value used in hotspot calling. Default is 0.2"),
        optparse::make_option(c("--pval"), default = 1e-5, help = "Threshold for p-values to call hotspots. Default is 1e-5"),
        optparse::make_option(
          c("--merge-distance"),
          action = "store",
          type = "integer",
          default = NULL,
          help = "During hotspot calling, two hotspots with distance smaller than this will be merged together. If not specified, the sliding window size will by used, i.e. 200bp by default"
        ),
        optparse::make_option(c("--verbose"), default = FALSE, action = "store_true")
      )
    )

    script_args <-
      optparse::parse_args(parser,
                           args = args,
                           convert_hyphens_to_underscores = TRUE)

    library(tidyverse)
    library(magrittr)

    # Process arguments
    if ("input" %in% names(script_args)) {
      script_args$input <-
        str_split(script_args$input, pattern = ":")[[1]]
    }

    if (!("chrom" %in% names(script_args)))
      script_args$chrom <- NULL


    c("chrom", "exclude_chrom", "exclude_region") %>%
      walk(function(arg) {
        if (!is.null(script_args[[arg]])) {
          script_args[[arg]] <<-
            str_split(script_args[[arg]], pattern = ":")[[1]]
        }
      })

    c("min_mapq",
      "min_fraglen",
      "max_fraglen",
      "window_size",
      "step_size") %>%
      walk(function(arg) {
        if (!is_null(script_args))
          script_args[[arg]] <<- as.numeric(script_args[[arg]])
      })

    if (is_null(script_args$merge_distance))
      script_args$merge_distance <- script_args$window_size

    if (script_args$window_size %% script_args$step_size != 0)
      stop("window_size must be multiples of step_size")

    # genome must be hs37-1kg
    if (script_args$genome != "hs37-1kg") {
      stop("Currently, only genome hs37-1kg is supported")
    }

    return(list(subcommand, script_args))
  }
}


# Make adjustment to IFS interval start/end, and write to disk
write_ifs_as_bedgraph <- function(ifs, script_args, comments) {
  # In IFS, each row corresponds to a 200bp region, and they are overlapping.
  # This is not a bedGraph track.
  #
  # We're more interested in a bedGraph-style non-overlapping track
  #
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
  offset <- as.integer((ws_ratio - 1) / 2 * script_args$step_size)

  GenomicRanges::start(ifs) <- GenomicRanges::start(ifs) + offset
  GenomicRanges::width(ifs) <- script_args$step_size

  # # Rearrage orders
  # df <- GenomicRanges::mcols(ifs) %>%
  #   as_tibble() %>%
  #   relocate(c(z_score, score), .after = end) %>% select(-gc, -mappability)
  # GenomicRanges::mcols(ifs) <- df

  bedtorch::write_bed(ifs,
                      file_path = script_args$output_ifs,
                      comments = comments)
}


log_mem <- function(label = "Unknown") {
  if (requireNamespace("lobstr")) {
    mem <- as.numeric(lobstr::mem_used()) / 1024 ** 2
    logging::logdebug(str_interp("[${label}] Memory used: ${mem} MB"))
  }
}


# example
# script_args <- list(
#   input = "sandbox/xh_intermediate/batch2_res_merged_sort.chr21.bed.gz",
#   input_ifs = "sandbox/xh_intermediate/batch2_res_merged_sort.chr21.ifs.bed.gz",
#   input_hotspot = "sandbox/xh_intermediate/batch2_res_merged_sort.chr21.hotspot.bed.gz",
#   genome = "hs37-1kg",
#   gc_correct = TRUE,
#   high_mappability = "sandbox/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
#   chrom = NULL,
#   exclude_chrom = NULL,
#   min_mapq = 30L,
#   min_fraglen = 50L,
#   max_fraglen = 1000L,
#   exclude_region = "encode.blacklist.hs37-1kg",
#   exclude_soft_clipping = FALSE,
#   window_size = 200L,
#   step_size = 20L,
#   cpois = FALSE,
#   fdr = 0.01,
#   local_pval = 1e-5,
#   merge_distance = 200L,
#   verbose = TRUE
# )


# 1. Read fragments
# 2. Calculate raw IFS scores (without GC)
# 3. Calculate GC content for each fragment
subcommand_ifs <- function(script_args) {
  logging::loginfo("Input: fragment data")

  frag <- script_args$input %>% map(function(input_file) {
    logging::loginfo(str_interp("Loading fragments: ${input_file} ..."))
    if (!is_null(script_args$chrom)) {
      frag <-
        read_fragments(
          input_file,
          range = setdiff(script_args$chrom, script_args$exclude_chrom),
          genome = script_args$genome
        )
    } else {
      frag <- read_fragments(input_file, genome = script_args$genome)
      if (!is_null(script_args$exclude_chrom)) {
        frag <- frag[!seqnames(frag) %in% script_args$exclude_chrom]
      }
    }

    frag
  }) %>% do.call(c, args = .)

  logging::loginfo("Fragments summary:")
  print(frag)
  log_mem("Done loading fragments")

  logging::loginfo("Calculating IFS scores ...")
  ifs <- ifs_score(
    frag,
    window_size = script_args$window_size,
    step_size = script_args$step_size,
    gc_correct = script_args$gc_correct,
    blacklist_region = script_args$exclude_region,
    high_mappability_region = script_args$high_mappability,
    min_mapq = script_args$min_mapq,
    min_fraglen = script_args$min_fraglen,
    max_fraglen = script_args$max_fraglen,
    exclude_soft_clipping = script_args$exclude_soft_clipping
  )

  rm(frag)

  logging::loginfo("Raw IFS summary:")
  print(ifs)
  log_mem("Done calculating raw IFS scores")

  if (script_args$gc_correct) {
    ifs <- calc_gc(ifs)
    log_mem("Done calculating GC contents")
  }

  bedtorch::write_bed(ifs, file_path = script_args$output_ifs)
}


subcommand_peak <- function(script_args) {
  # Load IFS score from input file
  logging::loginfo(str_interp("Loading raw IFS scores: ${script_args$input} ..."))
  ifs <-
    bedtorch::read_bed(script_args$input, genome = script_args$genome)

  logging::loginfo("Raw IFS summary:")
  print(ifs)

  if (script_args$gc_correct) {
    logging::loginfo("Performing GC correction ...")
    ifs <- gc_correct(ifs, span = 0.75)
    log_mem("Done GC correction")
  }

  logging::loginfo("Calculating z-scores ...")
  ifs <- calc_ifs_z_score(ifs)

  logging::loginfo("Calculating global p-values ...")
  ifs <- calc_pois_pval(ifs)
  log_mem("Done calculating global p-values")

  logging::loginfo("Calculating local p-values ...")
  ifs <-
    calc_pois_pval_local(
      ifs,
      window_size = script_args$window_size,
      step_size = script_args$step_size,
      local_layout = list(`5k` = 5e3L, `10k` = 10e3L)
    )
  log_mem("Done calculating local p-values")

  logging::loginfo("IFS summary:")
  print(ifs)

  write_ifs_as_bedgraph(ifs, script_args, comments)
}

subcommand_hotspot <- function(script_args) {
  hotspot_standard <-
    call_hotspot(
      ifs,
      fdr_cutoff = script_args$fdr,
      pval_cutoff = script_args$pval,
      local_pval_cutoff = script_args$pval,
      method = script_args$hotspot_method,
      merge_distance = script_args$merge_distance
    )
  # hotspot_standard <-
  #   call_hotspot(
  #     ifs,
  #     fdr_cutoff = script_args$fdr,
  #     pval_cutoff = script_args$local_pval,
  #     local_pval_cutoff = script_args$local_pval,
  #     merge_distance = script_args$merge_distance
  #   )

  if (is.null(hotspot_standard)) {
    logging::loginfo("Called 0 hotspots")
    # Write an empty file anyway. This is useful when you want to use snakemake
    # and cragr together
    system(str_interp("touch ${script_args$output_hotspot}"))
  } else {
    n_hotspot <- length(hotspot_standard)
    logging::loginfo(str_interp("Called ${n_hotspot} hotspots"))
    logging::loginfo("Writing results to disk ...")
    bedtorch::write_bed(hotspot_standard,
                        file_path = script_args$output_hotspot,
                        comments = comments)
  }
}

# Main ----
parse_script_args_result <- parse_script_args()
subcommand <- parse_script_args_result[[1]]
script_args <- parse_script_args_result[[2]]
rm(parse_script_args_result)

# Build comment lines
comments <- c(
  paste0("cragr version: ", as.character(packageVersion("cragr"))),
  paste0("bedtorch version: ", as.character(packageVersion("bedtorch"))),
  # All items in script_args
  names(script_args) %>% purrr::map_chr(function(name) {
    v <- script_args[[name]]
    v_str <- if (!is.null(v))
      v %>% purrr::map_chr(as.character) %>% paste(collapse = ":")
    else
      ""
    paste0(name, "=", v_str)
  })
)

if (script_args$verbose)
  logging::setLevel("DEBUG")

logging::loginfo(str_interp("Argument summary:"))
comments %>% purrr::walk(function(v)
  logging::loginfo(v))

library(cragr)

if (subcommand == "ifs") {
  subcommand_ifs(script_args)
} else if (subcommand == "peak") {
  subcommand_peak(script_args)
} else if (subcommand == "hotspot") {
  subcommand_hotspot(script_args)
}
