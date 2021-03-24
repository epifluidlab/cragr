
args <- commandArgs(trailingOnly = TRUE)
ws <- as.integer(args[1])
chroms <- args[-1]
logging::loginfo(paste0("Process chromosomes: ", paste0("chr", chroms, collapse = ", ")))
logging::loginfo(paste0("Window size: ", ws))

library(tidyverse)

process_mappability <- function(chrom) {
  logging::loginfo(paste0("Loading from Matlab: chr", chrom))
  file_path <- str_interp("/ocean/projects/mcb190124p/clawdeen/CRAG_V3/Basic_info/mappability/chr${chrom}_m.mat")
  mat_obj <- R.matlab::readMat(file_path)
  mat_obj <- mat_obj$map.b
  colnames(mat_obj) <- "score"
  maps <- as.data.frame(mat_obj)
  rm(mat_obj)
  data.table::setDT(maps)

  logging::loginfo("Determin boundaries ...")
  maps[, start := 1:.N - 1]
  # Determine boundaries where the score changes
  maps[, boundary := data.table::shift(score, n = 1) != score]
  maps[c(1, length(score)), boundary := TRUE]
  # Clean the dataset
  chrom0 <- chrom
  maps <- maps[boundary == TRUE][, end := data.table::shift(start, n = -1)][1:(length(score) - 1), .(chrom = chrom0, start, end, score)][score > 0]

  logging::loginfo("Obatining window averages ...")
  # Make 20bp window
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
    file = str_interp("wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.${ws}bp.chr${chrom}.bedGraph.gz"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE
  )
}

for (chrom in chroms) {
  process_mappability(chrom)
}
