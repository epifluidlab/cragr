# Determine the anchor point
get_anchor <- function(feature, anchor = c("center", "start", "end")) {
  anchor = match.arg(anchor)
  switch(
    anchor,
    center = as.integer(round(start(feature) + (width(
      feature
    ) - 1) / 2)),
    start = start(feature),
    end = end(feature),
    stop(paste0("Invalid anchor: ", anchor))
  )
}


#' @export
enrichment_analysis <-
  function(hotspot,
           feature,
           anchor = c("center", "start", "end"),
           half_width = 1000L,
           flip_rev = TRUE) {
    anchor <- match.arg(anchor)
    mcols(hotspot) <- NULL
    mcols(feature) <- NULL

    # Only process the common part
    common_seqlevels <- intersect(unique(seqnames(hotspot)),
                                  unique(seqnames(feature)))
    # common_seqlevels <-
    #   intersect(seqlevels(feature), seqlevels(hotspot))
    feature <-
      keepSeqlevels(feature, common_seqlevels, pruning.mode = "coarse")
    hotspot <-
      keepSeqlevels(hotspot, common_seqlevels, pruning.mode = "coarse")

    # Determine the anchor point
    anchor_pos <- get_anchor(feature, anchor = anchor)

    # Adjusted feature: anchor - hw -> anchor + hw
    ranges(feature) <-
      IRanges::IRanges(start = pmax(1, anchor_pos - half_width),
                       end = anchor_pos + half_width)

    hits <- findOverlaps(hotspot, feature)

    # Build the matching data frame between overlapping hotspot and feature
    matched_gr <- hotspot[queryHits(hits)]
    matched_gr$anchor <- anchor_pos[subjectHits(hits)]
    feature_strand <- strand(feature[subjectHits(hits)])

    pseudo_origin <- 100e6L
    ranges(matched_gr) <-
      IRanges::IRanges(
        start = ifelse(
          feature_strand == "-" & flip_rev,
          matched_gr$anchor - end(matched_gr),
          start(matched_gr) - matched_gr$anchor
        ) + pseudo_origin,
        width = width(matched_gr)
      )
    seqnames(matched_gr) <- seqnames(matched_gr)[1]

    scaffold <-
      GenomicRanges::GRanges(seqnames = seqnames(matched_gr)[1],
                             ranges = IRanges::IRanges(
                               start = (pseudo_origin - half_width):(pseudo_origin + half_width),
                               width = 1
                             ))

    data.table::data.table(
      pos = seq(-half_width, half_width),
      freq = countOverlaps(scaffold, matched_gr) / length(feature)
    )
  }


#' @export
signal_level_analysis <- function(hotspot, signal, half_width = 1000L) {
  mcols(hotspot) <- NULL
  common_seqlevels <- intersect(seqlevels(signal), seqlevels(hotspot))
  signal <- keepSeqlevels(signal, common_seqlevels, pruning.mode = "coarse")
  hotspot <- keepSeqlevels(hotspot, common_seqlevels, pruning.mode = "coarse")

  hotspot <- GenomicRanges::resize(hotspot, width = 2 * half_width + 1, fix = "center")
  hits <- GenomicRanges::findOverlaps(signal, hotspot)

  matched_gr <- signal[queryHits(hits)]
  matched_gr$origin <- start(hotspot[subjectHits(hits)]) + half_width
  offset <- 1e6L
  ranges(matched_gr) <-
    IRanges::IRanges(
      start = start(matched_gr) - matched_gr$origin + offset,
      width = width(matched_gr)
    )

  seq(offset - half_width, offset + half_width) %>%
    map_dfr(function(x) {
      gr <- matched_gr[start(matched_gr) <= x & end(matched_gr) >= x]
      if (length(gr) > 0)
        value <- mean(gr$score)
      else
        value <- NA

      tibble(offset = x - offset, value = value)
    })
}

