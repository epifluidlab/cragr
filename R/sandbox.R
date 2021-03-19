# #!/usr/bin/env Rscript
#
# suppressMessages(library(readr))
# suppressMessages(library(dplyr))
# suppressMessages(library(purrr))
# suppressMessages(library(stringr))
#
# # while(nrow(df <- data.table::fread('file:///dev/stdin')) > 0) {
# #
# # }
#
# lineno <- 0
#
# # f <- file("stdin")#, open = "rt")
# # open(f)
# # while (nrow(df <- tryCatch(read.delim(f, nrows = 150), error = function(e) NULL)) %||% 0 > 0) {
# #   # print(df)
# #   # lineno <- nrow(df)
# # }
# #
#
# while (TRUE) {
#   # lines <- readLines(con = f, n = -1)
#   # if (length(lines) == 0) {
#   #   break
#   # }
#
#   # message(str_interp("Chunk: ${length(lines)}"))
#   # lineno <<- lineno + length(lines)
#
#   # df <- data.table::fread(text = lines)
#   #
#   df <- suppressWarnings(data.table::fread("file:///dev/stdin"))
#   # # # df <- tryCatch(read.delim("/dev/stdin", nrows = 1500, comment.char = "#"), error = function(e) NULL)
#   # #
#   if (is.null(df) || nrow(df) == 0) {
#     print(df)
#     break
#   }
#   #
#   # colnames(df)[2] <- "start"
#   #
#   start_pos <- df$start[1]
#   end_pos <- tail(df$start, n = 1)
#   message(str_interp("Chunk: ${nrow(df)}"))
#   message(str_interp("${start_pos} -> ${end_pos}"))
#   # # print(df)
#   # # tryCatch(print(df), error = function(e) NULL)
#   lineno <<- lineno + nrow(df)
# }
#
# # close(f)
#
# message(str_interp("Total #: ${lineno}"))
#
# # while(nrow(df <- data.table::fread("file:///dev/stdin")) > 0) {
# #
# # }
# #
# # f <- file("stdin")
# # open(f)
# # while (length(lines <- readLines(f, n = -1)) > 0) {
# #   # conn <- textConnection(lines)
# #   # df <- read_tsv(file = conn, col_names = FALSE, col_types = cols())
# #   # df <- read.table(file = conn, sep = "\t") %>% as_tibble()
# #   # write(nrow(df), stderr())
# #   # print(df)
# #   # write(line, stderr())
# #   # process line
# # }


#
# > tmp_maps1
# V1    start
# 1:  0        0
# 2:  0        1
# 3:  0        2
# 4:  0        3
# 5:  0        4
# ---
#   51304562:  0 51304561
# 51304563:  0 51304562
# 51304564:  0 51304563
# 51304565:  0 51304564
# 51304566:  0 51304565

library(tidyverse)

process_mappability <- function(chrom) {
  logging::loginfo(paste0("Loading chr", chrom))
  file_path <- str_interp("../CRAG/Basic_info/mappability/chr${chrom}_m.mat")
  mat_obj <- R.matlab::readMat(file_path)
  mat_obj <- mat_obj$map.b
  colnames(mat_obj) <- "score"
  maps <- as.data.frame(mat_obj)
  rm(mat_obj)
  data.table::setDT(maps)

  maps[, start := 1:.N - 1]
  # Determine boundaries where the score changes
  maps[, boundary := data.table::shift(score, n = 1) != score]
  maps[c(1, length(score)), boundary := TRUE]
  # Clean the dataset
  chrom0 <- chrom
  maps <- maps[boundary == TRUE][, end := data.table::shift(start, n = -1)][1:(length(score) - 1), .(chrom = chrom0, start, end, score)][score > 0]

  # Make 20bp window
  ws <- 20L
  maps[, `:=`(wid1 = start %/% ws, wid2 = (end - 1) %/% ws)]

  maps <- maps[, .(
    chrom = chrom[1],
    start2 = {
      if (wid1 != wid2) {
        v <- wid1:wid2 * ws
        v[1] <- start[1]
        v
      } else {
        start
      }
    },
    end2 = {
      if (wid1 != wid2) {
        v <- (wid1:wid2 + 1) * ws
        v[length(v)] <- end[length(end)]
        v
      } else {
        end
      }
    },
    score = score),
    by = c("wid1", "wid2", "start")][, .(chrom, start = start2, end = end2, score, wid1, wid2)]

  maps[, `:=`(wid1 = NULL, wid2 = NULL, wid = start %/% ws)]
  maps <- maps[, .(
    chrom = chrom[1],
    start = as.integer(wid * ws),
    end = as.integer((wid + 1) * ws),
    score = sum((end - start) * score) / ws
  ), by = wid][, wid := NULL]
  data.table::fwrite(
    maps,
    file = str_interp("wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.chr${chrom}.bedGraph.gz"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE
  )
}

for (chrom in as.character(22:1)) {
  process_mappability(chrom)
}
