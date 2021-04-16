# XH's script generates chr_*.mat files, which are single-column matrix
# containing original IFS scores for each single-bp position.
# This script transforms it into 20bp windows.
rm(list = ls())
res_ifs <- data.table::fread("~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_s_non_chr21.txt.gz")
res_ifs[, id := ((1:.N) - 1) %/% 20]
res_ifs <- res_ifs[, .(score = sum(V1)), by = id]
res_ifs[, `:=`(chrom = "21", start = as.integer(id * 20L))][, end := start + 20]
tmp <- rollsum(res_ifs$score, k = 10, na_pad = TRUE, na_rm = TRUE, align = "left")
res_ifs[, score := tmp][, end := start + 200]
res_ifs <- res_ifs[, .(chrom, start, end, name = ".", score)]
res_ifs <- res_ifs[!is.na(score) & score > 0]
write_bed(res_ifs, file_path = "~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_s_non_chr21.bed.gz", create_index = TRUE)

