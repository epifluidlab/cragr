# Modify XH's code to output intermediate files. Used for troubleshooting: why my results are different from XH's?

rm(list = ls())
library(tidyverse)

devtools::load_all()

## Load XH region
region_chr21 <- data.table::fread("sandbox/xh_intermediate/region_batch2_s_res_chr21.txt.gz")
data.table::setnames(region_chr21, new = c("start", "end", "score", "gc", "maps", "p-val", "score0"))
region_chr21[, start := start - 1L]
data.table::setkey(region_chr21, "start", "end")
region_chr21[start >= 10991360]

## Load XH's raw IFS
xh_raw_ifs <- data.table::fread("sandbox/xh_intermediate/batch2_s_res_chr21_raw.txt.gz")[, start := 1:.N - 1][, group_id := start %/% 20L]
xh_raw_ifs[start >= 10991360][1:200, sum(V1)]

## Load XH's flag_bin
flag_bin <- data.table::fread("sandbox/xh_intermediate/flag_batch2_s_res_chr21.txt.gz")
flag_bin <- flag_bin[, {
  start = as.integer((1:.N - 1)* 20L)
  end = start + 200L
  list(start = start, end = end, flag = V1)
}]
data.table::setkey(flag_bin, "start", "end")
region_chr21[, flag := flag_bin$flag]


xh_raw_ifs_v1 <- xh_raw_ifs[, .(start = start[1], score = sum(V1)), by = group_id]
xh_raw_ifs_v1[, `:=`(start = as.integer(group_id * 20), end = as.integer(group_id * 20 + 20))]
xh_raw_ifs_v1[, score := cragr::rollsum(score, k = 10, na_pad = TRUE)]
xh_raw_ifs_v1[, end := start + 200L]
data.table::setkey(xh_raw_ifs_v1, "start", "end")

xh_raw_ifs_v1[start >= 10991360]

## Test if XH's raw IFS score matches XH region
merged <- xh_raw_ifs_v1[region_chr21, nomatch=0]
merged[abs(score - score0) > 1e-5]  # Empty, means YES, match


## Load my region
my_region_chr21 <- cragr::load_bed("sandbox/xh_intermediate/test_ifs_chr21.bedGraph.gz")
# data.table::setnames(my_region_chr21, 5:12,
my_region_chr21[, `:=`(start = start - 90L)][, end := start + 200L]
data.table::setkey(my_region_chr21, "chrom", "start", "end")

my_region_chr21[start >= 10991360]
region_chr21[start >= 10991360]

## Merge
merged_region <- region_chr21[my_region_chr21, nomatch = 0]
merged_region %>% str()
merged_region

## From XH raw IFS score, mimic:
ifs_new <- region_chr21[flag==1, .(chrom = factor("21"), start, end, score, gc = gc / 100, maps)]
data.table::setkey(ifs_new, "chrom", "start", "end")

## Compare XH's gc/maps with mine
my_region_chr21[, 1:9][ifs_new][abs(gc-i.gc)>1e-5]
my_region_chr21[, 1:9][ifs_new][abs(mappability - maps)>1e-5]

## Mappability is



ifs_new <- my_region_chr21[, .(chrom = factor("21"), start, end)]
data.table::setkey(ifs_new, "start", "end")
ifs_new <- xh_raw_ifs_v1[ifs_new][, .(chrom, start, end, score)]
data.table::setkey(ifs_new, "chrom", "start", "end")

gc <- "sandbox/xh_intermediate/test_gc_chr21.bedGrap.gz"
gc_dt <- cragr::load_bed("test_gc_chr21.bedGrap.gz")
mappability <- "/jet/home/haizizh/data/shared/genomes/hg19/wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz"
ifs_new <- cragr::.ifs_score_feed_raw(tmp2, chrom=tmp2$chrom[1], gc = gc, mappability_region = mappability)

xh_raw_ifs_v1[, .(chrom = "21", start, end, score)]
