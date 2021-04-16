# Xionghui's ma is an interlaced matrix file.
# The chunk order: chr1/group1 => chr1/group2 => chr2/group1 => chr2/group2 => ...
# Split it into two sorted data frames, one for group1 and the other for group2

batch2_s_peak_origin <- data.table::fread("~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_s_peak_origin.txt")
batch2_s_peak_origin <- 
  batch2_s_peak_origin[, .(chrom = as.character(V1), start = V2, end = V3)][, id := 1:.N]

# Find the boundaries
boundaries <-
  batch2_s_peak_origin[, .SD][, `:=`(delta = start - data.table::shift(start, 1)), by = chrom][delta < 0]

gp1 <- (batch2_s_peak_origin[
  , { 
    chrom0 <- chrom
    bid <- boundaries[chrom == chrom0, id]
    .SD[1:(bid - id[1])]
    }, by = chrom])
data.table::setkey(gp1, id)

gp2 <- (batch2_s_peak_origin[
  , { 
    chrom0 <- chrom
    bid <- boundaries[chrom == chrom0, id]
    .SD[(bid - id[1] + 1):.N]
  }, by = chrom])
data.table::setkey(gp2, id)

ma <- data.table::fread("~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_screen_matrix.csv")[, id := 1:.N]
data.table::setkey(ma, id)

ma[gp1][, id := NULL] %>% relocate(c("chrom", "start", "end"), .before = V1) %>%
  write_bed(file_path = "~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_screen_matrix.gp1.tsv.gz", create_index = TRUE)

ma[gp2][, id := NULL] %>% relocate(c("chrom", "start", "end"), .before = V1) %>%
  write_bed(file_path = "~/SynologyDrive/Research/CCHMC/cofrag/results/20210322-week-12/cragr/batch2_screen_matrix.gp2.tsv.gz", create_index = TRUE)

gp1[chrom == "21"]
