stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

library(rlang)

ifs_parser <- optparse::OptionParser(
  option_list = list(
    optparse::make_option(opt_str = c("-i", "--input"),
                          help = "Path to the input fragment file. The file should be in bgzip-compressed BED format, alongside with the .tbi index file."), 
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to the output file."),
    optparse::make_option(c("--genome"), type = "character", help = "Which reference genome the input fragment file is based on. Should be either GRCh37 or GRCh38."),
    optparse::make_option(
      c("-g", "--gc-correct"),
      default = FALSE,
      action = "store_true",
      help = "Perform GC correction."
    ),
    optparse::make_option(
      c("--gc-correct-method"),
      default = "standard",
      help = "GC correction method. Should be either standard or caret [standard]."
    ),
    optparse::make_option(c("--gc-correct-n"),
                          default = 1e6L,
                          help = "Maximal number of data points for GC correction model training [1000000]."),
    optparse::make_option(
      c("-m", "--high-mappability"),
      type = "character",
      help = "Path to the mappability file, which should be in BED format [NULL]."
    ),
    optparse::make_option(
      c("--chrom"),
      type = "character",
      help = "Perform the analysis for the specified chromosome."
    ),
    optparse::make_option(c("--min-mapq"),
                          default = 30L,
                          help = "Minimal MAPQ for fragments included in the analysis."),
    optparse::make_option(c("--min-fraglen"),
                          default = 50L,
                          help = "Minimal length for fragments included in the analysis."),
    optparse::make_option(c("--max-fraglen"),
                          default = 1000L,
                          help = "Maximal length for fragments included in the analysis."),
    optparse::make_option(
      c("--exclude-region"),
      type = "character",
      default = NULL,
      # "encode.blacklist.hs37-1kg",
      help = "BED files defining regions to be excluded from the analysis, separated by colon. Default is the ENCODE Blacklist: https://www.nature.com/articles/s41598-019-45839-z, which is included in this R package"
    ),
    optparse::make_option(
      c("--exclude-soft-clipping"),
      action = "store_true",
      default = FALSE,
      help = "Exclude fragments with leading soft-clipping from the analysis"
    ),
    optparse::make_option(c("-w", "--window-size"),
                          default = 200L,
                          help = "Size of the sliding window [200]"),
    optparse::make_option(c("-s", "--step-size"), default = 20L, help = "Step size of the sliding window [20]"),
    optparse::make_option(c(
      "-t", "--thread", default = 1L, help = "Number of threads [1]"
    )),
    optparse::make_option(c("--verbose"), default = FALSE, action = "store_true")
  )
)


peak_parser <- optparse::OptionParser(
  option_list = list(
    optparse::make_option(c("-i", "--input"), help = "Path to the input file. If there are multiple input files, they should be separated by colons"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output file"),
    optparse::make_option(c("--genome"), type = "character", help = "Genome of the input"),
    optparse::make_option(
      c("-g", "--gc-correct"),
      default = FALSE,
      action = "store_true",
      help = "Whether to perform GC correction"
    ),
    optparse::make_option(
      c("--gc-correct-method"),
      default = "standard",
      help = "Methods used in GC correction. Should be either standard or caret [standard]"
    ),
    optparse::make_option(
      c("--gc-correct-n"),
      default = 1e6L,
      help = "Maximal sample size for GC correction model training [1e6L]"
    ),
    optparse::make_option(
      c("--chrom"),
      type = "character",
      default = NULL,
      help = "Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used"
    ),
    optparse::make_option(
      c("-w", "--window-size"),
      default = 200L,
      help = "Size of the sliding window [200]"
    ),
    optparse::make_option(c("-s", "--step-size"), default = 20L, help = "Step size of the sliding window [20]"),
    optparse::make_option(c("-t", "--thread", default = 1L, help = "Number of threads [1]")),
    optparse::make_option(c("--verbose"), default = FALSE, action = "store_true")
  )
)


signal_parser <- optparse::OptionParser(
  option_list = list(
    optparse::make_option(c("-i", "--input"), help = "Path to the input file. If there are multiple input files, they should be separated by colons"),
    optparse::make_option(c("--hotspot"), help = "Path to the hotspot file"),
    optparse::make_option(c("-o", "--output"), type = "character", help = "Path to output file"),
    optparse::make_option(c("--genome"), type = "character", help = "Genome of the input"),
    optparse::make_option(
      c("-g", "--gc-correct"),
      default = FALSE,
      action = "store_true",
      help = "Whether to perform GC correction"
    ),
    optparse::make_option(
      c("--gc-correct-method"),
      default = "standard",
      help = "Methods used in GC correction. Should be either standard or caret [standard]"
    ),
    optparse::make_option(
      c("--gc-correct-n"),
      default = 1e6L,
      help = "Maximal sample size for GC correction model training [1e6L]"
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
      default = NULL, # "encode.blacklist.hs37-1kg",
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
      help = "Size of the sliding window [200]"
    ),
    # optparse::make_option(c("-s", "--step-size"), default = 20L, help = "Step size of the sliding window [20]"),
    optparse::make_option(c("-t", "--thread", default = 1L, help = "Number of threads [1]")),
    optparse::make_option(c("--verbose"), default = FALSE, action = "store_true")
  )
)



parse_script_args <- function() {
  if (interactive()) {
    script_args <- get0("script_args")
    subcommand <- get0("subcommand")
    
    if (is_null(script_args) || is_null(subcommand))
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
    
    if (!subcommand %in% c("ifs", "peak", "signal"))
      stop("Subcommand should be one of the following: ifs, peak, hotspot, signal")
    
    if (is_true(subcommand == "ifs")) {
      parser <- ifs_parser
    } else if (is_true(subcommand == "peak")) {
      parser <- peak_parser
    } else if (is_true(subcommand == "signal")) {
      parser <- signal_parser
    } else {
      stop("Subcommand should be one of the following: ifs, peak, hotspot, signal")
    }
    
    args <- args[-1]
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
    
    c("chrom", "exclude_region") %>%
      walk(function(arg) {
        if (!is_null(script_args[[arg]])) {
          script_args[[arg]] <<-
            str_split(script_args[[arg]], pattern = "[:,]")[[1]]
        }
      })
    
    c("min_mapq",
      "min_fraglen",
      "max_fraglen",
      "window_size",
      "step_size") %>%
      walk(function(arg) {
        if (!is_null(script_args) && arg %in% names(script_args))
          script_args[[arg]] <<- as.numeric(script_args[[arg]])
      })
    
    if (is_true(script_args$window_size %% script_args$step_size != 0))
      stop("window_size must be multiples of step_size")
    
    # genome must be hs37-1kg
    assertthat::assert_that(
      is_scalar_character(script_args$genome) &&
        script_args$genome %in% c("hs37-1kg", "GRCh37", "GRCh38")
    )
    # if (script_args$genome != "hs37-1kg") {
    #   stop("Currently, only genome hs37-1kg is supported")
    # }
    
    # assertthat::assert_that(script_args$hotspot_method %in% c("pois", "nb"))
    
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
                      file_path = script_args$output,
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
#   verbose = TRUE
# )


# 1. Read fragments
# 2. Calculate raw IFS scores (without GC)
# 3. Calculate GC content for each fragment
subcommand_ifs <- function(script_args) {
  logging::loginfo("Input: fragment data")
  
  # Make sure the genome is available
  bsgenome <- switch(
    script_args$genome,
    "GRCh37" = "BSgenome.Hsapiens.1000genomes.hs37d5",
    "hs37-1kg" = "BSgenome.Hsapiens.1000genomes.hs37d5",
    "GRCh38" = "BSgenome.Hsapiens.NCBI.GRCh38",
    stop(paste0("Invalid genome: ", genome_name))
  )
  
  assertthat::assert_that(requireNamespace(bsgenome), msg = str_interp("${bsgenome} is required"))
  
  frag <- script_args$input %>% map(function(input_file) {
    logging::loginfo(str_interp("Loading fragments: ${input_file} ..."))
    if (!is_null(script_args$chrom)) {
      frag <-
        read_fragments(
          input_file,
          range = script_args$chrom,
          genome = script_args$genome
        )
    } else {
      frag <- read_fragments(input_file, genome = script_args$genome)
    }
    
    frag
  }) %>% do.call(c, args = .)
  
  frag <- sort(frag)
  
  logging::loginfo("Fragments summary:")
  print(frag)
  log_mem("Done loading fragments")
  
  logging::loginfo("Calculating IFS scores ...")
  result <- ifs_score(
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
  
  ifs <- result$ifs
  avg_len <- result$avg_len
  frag_cnt <- result$frag_cnt
  
  comments <-
    c(comments,
      str_interp("avg_len=${avg_len}"),
      str_interp("frag_cnt=${frag_cnt}"))
  
  logging::loginfo("Raw IFS summary:")
  print(ifs)
  log_mem("Done calculating raw IFS scores")
  
  # if (script_args$gc_correct) {
  # Always calculate GC
  if (TRUE) {
    ifs <- calc_gc(ifs)
    log_mem("Done calculating GC contents")
  }
  
  bedtorch::write_bed(ifs, file_path = script_args$output, comments = comments)
}


subcommand_peak <- function(script_args) {
  # Load IFS score from input file
  logging::loginfo(str_interp("Loading raw IFS scores: ${script_args$input} ..."))
  ifs <-
    bedtorch::read_bed(
      script_args$input,
      genome = script_args$genome,
      col.names = c("chrom", "start", "end", "score", "cov", "fraglen", "gc")
    )
  
  logging::loginfo("Raw IFS summary:")
  print(ifs)
  
  if (script_args$gc_correct) {
    logging::loginfo("Performing GC correction ...")
    ifs <-
      gc_correct(
        ifs,
        span = 0.75,
        method = script_args$gc_correct_method,
        max_training_dataset = script_args$gc_correct_n,
        thread = script_args$thread, 
      )
    log_mem("Done GC correction")
  } else {
    # No GC correction, just placeholder
    ifs$score_pre_gc <- ifs$score
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
      local_layout = list(`50k` = 50e3L) #list(`5k` = 5e3L, `10k` = 10e3L, `25k` = 25e3L, `50k` = 50e3L)
    )
  log_mem("Done calculating local p-values")
  
  logging::loginfo("IFS summary:")
  print(ifs)
  
  write_ifs_as_bedgraph(ifs, script_args, comments)
}


# Perform signal-level analysis
subcommand_signal <- function(script_args) {
  # Make sure the genome is available
  bsgenome <- switch(
    script_args$genome,
    "GRCh37" = "BSgenome.Hsapiens.1000genomes.hs37d5",
    "hs37-1kg" = "BSgenome.Hsapiens.1000genomes.hs37d5",
    "GRCh38" = "BSgenome.Hsapiens.NCBI.GRCh38",
    stop(paste0("Invalid genome: ", genome_name))
  )
  assertthat::assert_that(requireNamespace(bsgenome), msg = str_interp("${bsgenome} is required"))
  assertthat::assert_that(
    rlang::is_scalar_character(script_args$chrom),
    msg = "chrom should be the name of a single chromosome"
  )
  
  logging::loginfo("Loading fragment files ...")
  # Calculate IFS scores as usual, which is used to get the GC-correction model
  frag <- script_args$input %>% map(function(input_file) {
    if (!is_null(script_args$chrom)) {
      frag <-
        read_fragments(
          input_file,
          range = script_args$chrom,
          genome = script_args$genome
        )
    } else {
      frag <- read_fragments(input_file, genome = script_args$genome)
    }
    
    frag
  }) %>% do.call(c, args = .)
  frag <- sort(frag)
  
  # IFS from all fragments. Necessary for GC and z-score transformation
  result <- ifs_score(
    frag,
    window_size = script_args$window_size,
    step_size = script_args$window_size,
    gc_correct = script_args$gc_correct,
    blacklist_region = script_args$exclude_region,
    high_mappability_region = script_args$high_mappability,
    min_mapq = script_args$min_mapq,
    min_fraglen = script_args$min_fraglen,
    max_fraglen = script_args$max_fraglen,
    exclude_soft_clipping = script_args$exclude_soft_clipping
  )
  ifs <- result$ifs
  ifs <- calc_gc(ifs)
  ifs$score_pre_gc <- ifs$score
  
  if (script_args$gc_correct) {
    logging::loginfo("Perform GC correction ...")
    
    gc_result <-
      gc_correct(
        ifs,
        span = 0.75,
        method = script_args$gc_correct_method,
        max_training_dataset = script_args$gc_correct_n,
        thread = script_args$thread,
        return_model = TRUE
      )
    gc_model <- gc_result$model
    ifs <- gc_result$ifs
  } else {
    gc_model <- NULL
  }
  score_mean <- mean(ifs$score, na.rm = TRUE)
  score_sd <- sd(ifs$score, na.rm = TRUE)
  
  hotspot <- bedtorch::read_bed(script_args$hotspot, genome = script_args$genome)
  hotspot <- GenomicRanges::resize(hotspot, width = script_args$window_size, fix = "center")
  result <- ifs_score(
    frag,
    interval = hotspot,
    window_size = NULL,
    step_size = NULL,
    gc_correct = script_args$gc_correct,
    blacklist_region = script_args$exclude_region,
    high_mappability_region = script_args$high_mappability,
    min_mapq = script_args$min_mapq,
    min_fraglen = script_args$min_fraglen,
    max_fraglen = script_args$max_fraglen,
    exclude_soft_clipping = script_args$exclude_soft_clipping
  )
  ifs2 <- result$ifs
  ifs2 <- calc_gc(ifs2)
  ifs2$score_pre_gc <- ifs2$score
  
  if (script_args$gc_correct) {
    na_idx <- is.na(ifs2$gc)
    pred <-
      predict(gc_model, newdata = data.frame(gc = ifs2$gc[!na_idx]))
    ifs2$score <- NA
    ifs2$score[!na_idx] <-
      pmax(0, ifs2$score_pre_gc[!na_idx] - pred + mean(ifs$score_pre_gc, na.rm = TRUE))
  }
  
  ifs2$z_score <- (ifs2$score - score_mean) / score_sd
  bedtorch::write_bed(ifs2, file_path = script_args$output, comments = comments)
}

# Main ----
if (interactive()) {
  subcommand <- "ifs"
  
  # example
  script_args <- list(
    input = "sandbox/frag/inhouse_breast_healthy_v2.hg19.frag.bed.gz",
    output = "sandbox/inhouse_breast_healthy_v2.ifs.raw.chr17.bed.gz",
    genome = "GRCh37",
    gc_correct = FALSE,
    gc_correct_method = "standard",
    high_mappability = NULL,
    chrom = "17",
    min_mapq = 30L,
    min_fraglen = 50L,
    max_fraglen = 1000L,
    exclude_region = "sandbox/data/ENCODE.blacklist.hg19.ENCFF001TDO.bed",
    exclude_soft_clipping = TRUE,
    window_size = 200L,
    step_size = 20L,
    verbose = TRUE
  )
} else {
  parse_script_args_result <- parse_script_args()
  subcommand <- parse_script_args_result[[1]]
  script_args <- parse_script_args_result[[2]]
  rm(parse_script_args_result)
}


# Build comment lines
comments <- c(
  paste0("cragr_version=", as.character(packageVersion("cragr"))),
  paste0("bedtorch_version=", as.character(packageVersion("bedtorch"))),
  paste0("timestamp=", lubridate::now() %>% format("%Y-%m-%dT%H:%M:%S%z")),
  # All items in script_args
  names(script_args) %>% purrr::map_chr(function(name) {
    v <- script_args[[name]]
    v_str <- if (!is_null(v))
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
  # } else if (subcommand == "hotspot") {
  #   subcommand_hotspot(script_args)
} else if (subcommand == "signal") {
  subcommand_signal(script_args)
} else {
  stop("Invalid subcommand")
}
