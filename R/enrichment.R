#' @export
enrichment_analysis <- function(hotspot, feature, half_width = 1000L, flip_rev = TRUE) {
  mcols(hotspot) <- NULL
  mcols(feature) <- NULL
  common_seqlevels <- intersect(seqlevels(feature), seqlevels(hotspot))
  feature <- keepSeqlevels(feature, common_seqlevels, pruning.mode = "coarse")
  hotspot <- keepSeqlevels(hotspot, common_seqlevels, pruning.mode = "coarse")

  # Expand feature from mid point
  feature_mid_point <- as.integer(round(start(feature) + (width(feature) - 1) / 2))
  # Remove potential negative start positions
  positive_feature_span <- feature_mid_point - half_width > 0
  feature <- feature[positive_feature_span]
  feature_mid_point <- feature_mid_point[positive_feature_span]

  ranges(feature) <-
    IRanges::IRanges(start = feature_mid_point - half_width,
                     end = feature_mid_point + half_width)
  hits <- findOverlaps(hotspot, feature)

  # Build the matching data frame between overlapping hotspot and feature
  matched_gr <- hotspot[queryHits(hits)]
  mcols(matched_gr) <- NULL
  matched_gr$origin <- as.integer(start(feature[subjectHits(hits)]) + half_width)

  strand(matched_gr) <- strand(feature[subjectHits(hits)])

  # # Some matched hotspots are large, or far away from the feature origin, thus
  # # go beyond the [origin - half_window, origin + half_window] interval, so we
  # # need to clip them.
  # idx <- start(matched_gr) - matched_gr$origin < -half_width
  # start(matched_gr[idx]) <- matched_gr$origin[idx] - half_width
  # idx <- end(matched_gr) - matched_gr$origin > half_width
  # end(matched_gr[idx]) <- matched_gr$origin[idx] + half_width

  # Build a pseudo BED using relative coordinates, so that we can count the "coverage"
  pseudo_origin <- 10e6L + 1L
  seqnames(matched_gr) <- seqnames(matched_gr)[1]
  # flip - strand
  flip_flag <- strand(matched_gr) == "-"
  ranges(matched_gr) <-
    IRanges::IRanges(
      start = ifelse(flip_flag & flip_rev,
                     matched_gr$origin - end(matched_gr) + pseudo_origin,
                     start(matched_gr) - matched_gr$origin + pseudo_origin),
      width = width(matched_gr)
    )

  scaffold <- GenomicRanges::GRanges(
    seqnames = seqnames(matched_gr)[1],
    ranges = IRanges::IRanges(start = (pseudo_origin - half_width):(pseudo_origin + half_width),
                              width = 1)
  )

  data.table::data.table(
    pos = seq(-half_width, half_width),
    freq = countOverlaps(scaffold, matched_gr) / length(matched_gr)
  )
}
#
# enrichment_tts <- enrichment_analysis(hotspot = hotspot, feature = tts)
# enrichment_tts %>% ggplot(aes(x=offset,y=count))+geom_line()
#
# tss <- bedtorch::read_bed("enrichment/Homo_sapiens.b37.75.TSS.commonChr.sort.bed")
# enrichment_tss <- enrichment_analysis(hotspot = hotspot, feature = tss)
# enrichment_tss %>% ggplot(aes(x=offset,y=count))+geom_line()
#
# hotspot <- bedtorch::read_bed("sandbox/EE86217.hotspot.bed.gz", genome = "hs37-1kg")



