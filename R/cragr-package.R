#' @import here
#' @import dplyr
#' @import purrr
#' @import stringr
#' @importFrom magrittr `%<>%`
#' @importFrom data.table setnames setkey fread fwrite shift setDT data.table rbindlist as.data.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps countOverlaps mcols `mcols<-`
#' @importFrom GenomicRanges pintersect ranges `ranges<-` width `width<-` start `start<-` end `end<-` strand `strand<-`
#' @importFrom GenomeInfoDb seqnames `seqnames<-` seqlengths keepSeqlevels seqinfo `seqinfo<-` seqlevels `seqlevels<-`
#' @importFrom S4Vectors from to queryHits subjectHits
NULL

## usethis namespace: start
#' @useDynLib cragr, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
