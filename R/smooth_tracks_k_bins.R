smooth_tracks_k_bins <- function(Tracks, smoothing.window = 5) {
    if (smoothing.window < 1) stop("smoothing.window must be >= 1")
    if (smoothing.window %% 2 == 0) stop("smoothing.window must be an odd integer")
    k <- (smoothing.window - 1) %/% 2

    Tracks.smoothed <- Tracks

    for (chrom.ind in seq_along(Tracks)) {
        track.df <- Tracks[[chrom.ind]]
        n <- nrow(track.df)
        if (n == 0) next

        # Replace NAs with 0 before smoothing
        rec <- track.df[, "records"]
        nuc <- track.df[, "nucleotides"]
        rec[is.na(rec)] <- 0
        nuc[is.na(nuc)] <- 0

        for (i in seq_len(n)) {
            left <- max(1L, i - k)
            right <- min(n, i + k)
            Tracks.smoothed[[chrom.ind]][i, "records"] <- mean(rec[left:right], na.rm = TRUE)
            Tracks.smoothed[[chrom.ind]][i, "nucleotides"] <- mean(nuc[left:right], na.rm = TRUE)
        }
    }
    return(Tracks.smoothed)
}