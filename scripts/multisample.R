# Given a list of hotspots, for multiple samples, output the normalized IFS scores within the hotspot regions

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

if (interactive()) {
  if (is.null(get0("script_args")) || is.null(get0("subcommand")))
    stop_quietly()
} else {
  # The first argument is subcommand. Should be multisample
  args <- commandArgs(trailingOnly = TRUE)
  subcommand <- args[1]
  args <- args[-1]

  assertthat::are_equal(subcommand, "multisample")


  # Run in CLI script mode
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("-s", "--hotspot"), help = "Path to hotspot files. For multiple files, separate with colons. For example: hs1.bed.gz:hs2.bed.gz"),
      optparse::make_option(c("--samples"), help = "A .tsv file containing sample information"),
      optparse::make_option(c("-o", "--output"), help = "Output file path")
    )
  )
  script_args <-
    optparse::parse_args(
      parser,
      args = args,
      convert_hyphens_to_underscores = TRUE
    )

  # Process arguments
  script_args$hotspot <- str_split(script_args$hotspot, ":")[[1]]
}

print(script_args)
library(tidyverse)

load_hotspot <- function(hotspot_file) {
  cragr::check_binaries("bedtools", stop_on_fail = TRUE)

  hotspot <- hotspot_file %>% map(cragr::load_bed) %>% data.table::rbindlist()
  hotspot <- hotspot[, 1:3]
  data.table::setkey(hotspot, "chrom", "start", "end")

  merged <-
    bedr::bedr(method = "merge",
               input = list(i = hotspot),
               check.chr = FALSE,
               check.zero.based = FALSE,
               check.valid = FALSE,
               check.sort = FALSE,
               check.merge = FALSE,
               capture.output = "disk")
  data.table::setDT(merged)
  merged[, chrom := factor(chrom, levels = levels(hotspot$chrom))]
  data.table::setkey(merged, "chrom", "start", "end")
  merged
}


intersect_sample <- function(hotspot, sample_info) {
  ifs_list <- 1:nrow(sample_info) %>%
    map(function(row_id) {
      fields <-
        colnames(data.table::fread(file = sample_info$file[row_id], nrows = 100))
      # Remove the leading #-symbol from the 1st column name
      fields[1] <- str_replace(fields[1], "^#[ ]*", "")

      ifs <- bedr::bedr(method = "intersect",
                        input = list(a = sample_info$file[row_id], b = hotspot),
                        check.chr = FALSE,
                        check.zero.based = FALSE,
                        check.valid = FALSE,
                        check.sort = FALSE,
                        check.merge = FALSE,
                        capture.output = "disk")
      data.table::setDT(ifs)
      data.table::setnames(ifs, new = fields)
      ifs <- ifs[, 1:4]
      ifs[, chrom := factor(chrom, levels = levels(hotspot$chrom))]
      data.table::setkey(ifs, "chrom", "start", "end")
      ifs[, name := sample_info$name[row_id]]
      ifs
    })

  browser()
  ifs <- data.table::rbindlist(ifs_list) %>%
    pivot_wider(names_from = "name", values_from = "z_score") %>%
    data.table::as.data.table()
  data.table::setkey(ifs, "chrom", "start", "end")
  ifs
}

hotspot <- load_hotspot(script_args$hotspot)

sample_info <- fread(script_args$samples, col.names = c("name", "file"), header = FALSE)
result <- intersect_sample(hotspot = hotspot, sample_info = sample_info)

data.table::fwrite(result, quote = FALSE, sep = "\t", file = script_args$output)
