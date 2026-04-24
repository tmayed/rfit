#' @title segment_data
#' @description Split a sample into multiple regimes using an exact
#'   one-dimensional k-medians segmentation.
#' @param data A numeric vector containing the sample to be segmented.
#' @param k Integer. The number of segments to force. If `NULL` (default), 
#'   the optimal number of segments is selected using the specified criterion.
#' @param criterion Character. One of `"gap"`, `"BIC"`, or `"elbow"`. 
#'   `"elbow"` is usually best for finding a "natural" number of regimes (2-4).
#' @param threshold Numeric. Improvement threshold for the `"elbow"` method. 
#'   A new segment is only added if it improves the objective by more than 
#'   this ratio (0 to 1). Defaults to 0.05 (5%).
#' @param penalty Numeric. Penalty multiplier for BIC. Defaults to 1.0.
#' @param min_segment Minimum number of observations required in each regime.
#' @param max_segments Maximum number of contiguous segments to consider.
#' @param n_boot Number of reference datasets used by the gap statistic.
#' @export
segment_data <- function(data, k = NULL, criterion = c("gap", "BIC", "elbow"), 
                        threshold = 0.05, penalty = 1.0, 
                        min_segment = 5L, max_segments = NULL, n_boot = 10L) {
  criterion <- match.arg(criterion)

  # --- 1. Data Pre-processing ---
  x <- sort(data[is.finite(data)])
  n <- length(x)
  if (n == 0) stop("`data` must contain at least one finite value.")

  # --- 2. Constraint Validation ---
  min_segment <- as.integer(min_segment[1])
  max_allowed_segments <- n %/% min_segment
  if (!is.null(k)) {
    k <- as.integer(k[1])
    max_segments <- k
  } else if (is.null(max_segments)) {
    max_segments <- min(10L, max_allowed_segments)
  }
  max_segments <- min(as.integer(max_segments), max_allowed_segments)

  # --- 3. Core Algorithm: Dynamic Programming Segmentation ---
  compute_segmentation <- function(x_sorted, compute_prev = TRUE) {
    x_sum <- c(0, cumsum(x_sorted))
    seg_cost <- matrix(Inf, nrow = n, ncol = n)
    for (start_idx in seq_len(n - min_segment + 1L)) {
      end_indices <- (start_idx + min_segment - 1L):n
      med_indices <- (start_idx + end_indices) %/% 2L
      med_vals <- x_sorted[med_indices]
      left_sums <- x_sum[med_indices + 1L] - x_sum[start_idx]
      right_sums <- x_sum[end_indices + 1L] - x_sum[med_indices + 1L]
      n_left <- med_indices - start_idx + 1L
      n_right <- end_indices - med_indices
      seg_cost[start_idx, end_indices] <- (med_vals * n_left - left_sums) + (right_sums - med_vals * n_right)
    }
    dp <- matrix(Inf, nrow = max_segments, ncol = n)
    prev <- if (compute_prev) matrix(NA_integer_, nrow = max_segments, ncol = n) else NULL
    dp[1L, min_segment:n] <- seg_cost[1L, min_segment:n]
    if (max_segments >= 2L) {
      for (curr_k in 2L:max_segments) {
        min_end <- curr_k * min_segment
        if (min_end > n) next
        for (end_idx in min_end:n) {
          split_indices <- ((curr_k - 1L) * min_segment):(end_idx - min_segment)
          candidates <- dp[curr_k - 1L, split_indices] + seg_cost[split_indices + 1L, end_idx]
          best_idx <- which.min(candidates)
          dp[curr_k, end_idx] <- candidates[best_idx]
          if (compute_prev) prev[curr_k, end_idx] <- split_indices[best_idx]
        }
      }
    }
    list(dp = dp, prev = prev)
  }

  observed <- compute_segmentation(x, compute_prev = TRUE)
  valid_segments <- which(is.finite(observed$dp[, n]))
  obs_costs <- observed$dp[valid_segments, n]
  
  segment_count <- NULL
  method_str <- ""

  if (!is.null(k)) {
    segment_count <- k
    method_str <- sprintf("Forced segmentation into k=%d segments", k)
  } else if (criterion == "elbow") {
    # Elbow method: stop when improvement < threshold
    segment_count <- 1
    if (length(obs_costs) > 1) {
      for (i in 2:length(obs_costs)) {
        improvement <- (obs_costs[i-1] - obs_costs[i]) / obs_costs[i-1]
        if (improvement >= threshold) {
          segment_count <- i
        } else {
          break # Stop adding segments
        }
      }
    }
    method_str <- sprintf("Exact k-medians selected by elbow (threshold=%.2f)", threshold)
  } else if (criterion == "BIC") {
    bics <- n * log(obs_costs / n) + penalty * valid_segments * log(n)
    segment_count <- valid_segments[which.min(bics)]
    method_str <- sprintf("Exact k-medians selected by BIC (penalty=%.1f)", penalty)
  } else {
    # Gap Statistic
    cost_floor <- max(obs_costs[1L], 1) * .Machine$double.eps
    log_w <- log(pmax(obs_costs, cost_floor))
    ref_log_w <- matrix(NA_real_, nrow = n_boot, ncol = length(valid_segments))
    x_range <- range(x)
    for (b in seq_len(n_boot)) {
      ref_x <- sort(stats::runif(n, min = x_range[1L], max = x_range[2L]))
      ref_fit <- compute_segmentation(ref_x, compute_prev = FALSE)
      ref_log_w[b, ] <- log(pmax(ref_fit$dp[valid_segments, n], cost_floor))
    }
    gap <- colMeans(ref_log_w) - log_w
    gap_se <- apply(ref_log_w, 2L, stats::sd) * sqrt(1 + 1 / n_boot)
    selected_idx <- length(valid_segments)
    if (length(valid_segments) > 1L) {
      for (idx in seq_len(length(valid_segments) - 1L)) {
        if (gap[idx] >= gap[idx + 1L] - gap_se[idx + 1L]) {
          selected_idx <- idx; break
        }
      }
    }
    segment_count <- valid_segments[selected_idx]
    method_str <- "Exact k-medians selected by gap statistic"
  }

  # --- 5. Backtracking and Formatting ---
  break_indices <- integer(0); end_idx <- n
  if (segment_count >= 2L) {
    for (curr_k in segment_count:2L) {
      split_idx <- observed$prev[curr_k, end_idx]
      break_indices <- c(split_idx, break_indices)
      end_idx <- split_idx
    }
  }
  splits <- numeric(length(break_indices))
  if (length(break_indices) > 0L) {
    for (idx in seq_along(break_indices)) {
      i <- break_indices[idx]; splits[idx] <- (x[i] + x[i + 1L]) / 2
    }
  }

  structure(
    list(
      splits = splits, n_splits = length(splits), n_segments = segment_count,
      x1 = if (length(splits) >= 1L) splits[1L] else NA_real_,
      x2 = if (length(splits) >= 2L) splits[2L] else NA_real_,
      indices = break_indices, counts = diff(c(0L, break_indices, n)),
      objective = observed$dp[segment_count, n],
      method = method_str
    ),
    class = "segment_data"
  )
}
