#' @title segment_data
#' @description Split a sample into multiple regimes using an exact
#'   one-dimensional k-medians segmentation.
#' @param data A numeric vector containing the sample to be segmented.
#' @param k Integer. The number of segments to force. If `NULL` (default), 
#'   the optimal number of segments is selected using the gap statistic.
#' @param min_segment Minimum number of observations required in each regime.
#'   Ensures statistical significance for each segment. Defaults to 5.
#' @param max_segments Maximum number of contiguous segments to consider when 
#'   `k` is `NULL`. Defaults to `min(8, floor(length(data) / min_segment))`.
#' @param n_boot Number of reference datasets used by the gap statistic
#'   to determine the optimal number of segments. Defaults to 10.
#' @return A list of class `segment_data` containing:
#'   \item{splits}{Numeric vector of split locations (midpoints between points).}
#'   \item{n_splits}{Total number of splits found.}
#'   \item{n_segments}{Total number of segments.}
#'   \item{x1, x2}{Compatibility fields for the first and second splits.}
#'   \item{indices}{Indices in the sorted data where splits occur.}
#'   \item{counts}{Number of observations in each segment.}
#'   \item{objective}{Total absolute deviation (sum of distances to medians).}
#'   \item{gap}{Gap statistic value for the selected model.}
#'   \item{gap_se}{Standard error of the gap statistic.}
#'   \item{method}{String describing the selection criteria used.}
#' @details
#' The function identifies optimal contiguous segments in one-dimensional data.
#' For a fixed number of segments $k$, it minimizes the total absolute deviation:
#' $L = \sum_{j=1}^k \sum_{i \in S_j} |x_i - \text{median}(S_j)|$.
#' This is solved exactly using dynamic programming in $O(k \cdot n^2)$ time.
#' If $k$ is not provided, the optimal number of segments is selected using the 
#' Gap Statistic (Tibshirani et al., 2001).
#' @export
segment_data <- function(data, k = NULL, min_segment = 5L, max_segments = NULL, n_boot = 10L) {
  # --- 1. Data Pre-processing ---
  if (!is.numeric(data)) {
    stop("`data` must be numeric.")
  }

  # Use only finite values and sort for contiguous segmentation
  x <- sort(data[is.finite(data)])
  n <- length(x)
  
  if (n == 0) {
    stop("`data` must contain at least one finite value.")
  }

  # --- 2. Constraint Validation ---
  min_segment <- as.integer(min_segment[1])
  if (is.na(min_segment) || min_segment < 2L) {
    stop("`min_segment` must be an integer >= 2.")
  }

  # Determine search range for number of segments
  max_allowed_segments <- n %/% min_segment
  
  if (!is.null(k)) {
    k <- as.integer(k[1])
    if (is.na(k) || k < 1L) {
      stop("`k` must be an integer >= 1.")
    }
    if (k > max_allowed_segments) {
      stop(sprintf("Requested k=%d segments, but only %d possible with min_segment=%d.", 
                   k, max_allowed_segments, min_segment))
    }
    max_segments <- k
  } else {
    if (is.null(max_segments)) {
      max_segments <- min(8L, max_allowed_segments)
    } else {
      max_segments <- as.integer(max_segments[1])
    }
  }

  if (is.na(max_segments) || max_segments < 1L) {
    stop("`max_segments` must be an integer >= 1.")
  }
  max_segments <- min(max_segments, max_allowed_segments)
  
  n_boot <- as.integer(n_boot[1])
  if (is.na(n_boot) || n_boot < 1L) {
    stop("`n_boot` must be an integer >= 1.")
  }

  # --- 3. Core Algorithm: Dynamic Programming Segmentation ---
  # Internal function to solve the k-medians problem for all k up to max_segments
  compute_segmentation <- function(x_sorted, compute_prev = TRUE) {
    # cumsum allows O(1) calculation of sum(x[i:j])
    x_sum <- c(0, cumsum(x_sorted))
    
    # Pre-calculate costs for all valid [start, end] pairs.
    seg_cost <- matrix(Inf, nrow = n, ncol = n)

    for (start_idx in seq_len(n - min_segment + 1L)) {
      end_indices <- (start_idx + min_segment - 1L):n
      med_indices <- (start_idx + end_indices) %/% 2L
      
      med_vals <- x_sorted[med_indices]
      left_sums <- x_sum[med_indices + 1L] - x_sum[start_idx]
      right_sums <- x_sum[end_indices + 1L] - x_sum[med_indices + 1L]
      
      n_left <- med_indices - start_idx + 1L
      n_right <- end_indices - med_indices
      
      # Vectorized cost calculation for all end_indices for this start_idx
      seg_cost[start_idx, end_indices] <- (med_vals * n_left - left_sums) + 
                                          (right_sums - med_vals * n_right)
    }

    # dp[k, i] = min cost to segment first i points into k segments
    dp <- matrix(Inf, nrow = max_segments, ncol = n)
    # prev stores the optimal split point for backtracking
    prev <- if (compute_prev) matrix(NA_integer_, nrow = max_segments, ncol = n) else NULL

    # Base case: k = 1
    dp[1L, min_segment:n] <- seg_cost[1L, min_segment:n]

    # Recursive step: dp[k, i] = min_{j} (dp[k-1, j] + cost(j+1, i))
    if (max_segments >= 2L) {
      for (curr_k in 2L:max_segments) {
        min_end <- curr_k * min_segment
        if (min_end > n) next
        
        for (end_idx in min_end:n) {
          split_min <- (curr_k - 1L) * min_segment
          split_max <- end_idx - min_segment
          
          split_indices <- split_min:split_max
          candidates <- dp[curr_k - 1L, split_indices] + seg_cost[split_indices + 1L, end_idx]
          
          best_idx <- which.min(candidates)
          dp[curr_k, end_idx] <- candidates[best_idx]
          if (compute_prev) {
            prev[curr_k, end_idx] <- split_indices[best_idx]
          }
        }
      }
    }

    list(dp = dp, prev = prev)
  }

  # Run segmentation on the actual observed data
  observed <- compute_segmentation(x, compute_prev = TRUE)
  valid_segments <- which(is.finite(observed$dp[, n]))
  if (length(valid_segments) == 0L) {
    stop("Failed to compute a valid segmentation.")
  }

  # Extract log-objective values for observed data
  obs_costs <- observed$dp[valid_segments, n]
  cost_floor <- max(obs_costs[1L], 1) * .Machine$double.eps
  log_w <- log(pmax(obs_costs, cost_floor))

  # Special case: All data points are identical
  x_range <- range(x)
  if (x_range[1L] == x_range[2L]) {
    return(structure(
      list(
        splits = numeric(0),
        n_splits = 0L,
        n_segments = 1L,
        x1 = NA_real_,
        x2 = NA_real_,
        indices = integer(0),
        counts = n,
        objective = 0,
        objective_name = "total_absolute_deviation_to_segment_medians",
        gap = 0,
        gap_se = 0,
        method = "All observations are identical, no splits possible."
      ),
      class = "segment_data"
    ))
  }

  # --- 4. Selection Logic ---
  if (!is.null(k)) {
    # Forced number of segments
    segment_count <- k
    selected_idx <- which(valid_segments == k)
    if (length(selected_idx) == 0) {
      stop(sprintf("Could not find a valid segmentation for k=%d.", k))
    }
    
    # Compute gap statistic for reference even if forced
    ref_log_w <- matrix(NA_real_, nrow = n_boot, ncol = length(valid_segments))
    for (b in seq_len(n_boot)) {
      ref_x <- sort(stats::runif(n, min = x_range[1L], max = x_range[2L]))
      ref_fit <- compute_segmentation(ref_x, compute_prev = FALSE)
      ref_costs <- ref_fit$dp[valid_segments, n]
      ref_log_w[b, ] <- log(pmax(ref_costs, cost_floor))
    }
    gap <- colMeans(ref_log_w) - log_w
    gap_sd <- apply(ref_log_w, 2L, stats::sd)
    gap_se <- gap_sd * sqrt(1 + 1 / n_boot)
    method_str <- sprintf("Forced segmentation into k=%d segments", k)
  } else {
    # Automatic selection using Gap Statistic
    ref_log_w <- matrix(NA_real_, nrow = n_boot, ncol = length(valid_segments))
    for (b in seq_len(n_boot)) {
      ref_x <- sort(stats::runif(n, min = x_range[1L], max = x_range[2L]))
      ref_fit <- compute_segmentation(ref_x, compute_prev = FALSE)
      ref_costs <- ref_fit$dp[valid_segments, n]
      ref_log_w[b, ] <- log(pmax(ref_costs, cost_floor))
    }

    gap <- colMeans(ref_log_w) - log_w
    gap_sd <- apply(ref_log_w, 2L, stats::sd)
    gap_se <- gap_sd * sqrt(1 + 1 / n_boot)

    # Pick the smallest k such that Gap(k) >= Gap(k+1) - SE(k+1)
    selected_idx <- length(valid_segments)
    if (length(valid_segments) > 1L) {
      for (idx in seq_len(length(valid_segments) - 1L)) {
        if (gap[idx] >= gap[idx + 1L] - gap_se[idx + 1L]) {
          selected_idx <- idx
          break
        }
      }
    }
    segment_count <- valid_segments[selected_idx]
    method_str <- "Exact k-medians selected by gap statistic"
  }

  # --- 5. Backtracking and Formatting Results ---
  break_indices <- integer(0)
  end_idx <- n
  if (segment_count >= 2L) {
    # Trace back through the prev matrix to find optimal split points
    for (curr_k in segment_count:2L) {
      split_idx <- observed$prev[curr_k, end_idx]
      break_indices <- c(split_idx, break_indices)
      end_idx <- split_idx
    }
  }

  splits <- numeric(length(break_indices))
  if (length(break_indices) > 0L) {
    for (idx in seq_along(break_indices)) {
      i <- break_indices[idx]
      splits[idx] <- (x[i] + x[i + 1L]) / 2
    }
  }

  counts <- diff(c(0L, break_indices, n))

  structure(
    list(
      splits = splits,
      n_splits = length(splits),
      n_segments = segment_count,
      x1 = if (length(splits) >= 1L) splits[1L] else NA_real_,
      x2 = if (length(splits) >= 2L) splits[2L] else NA_real_,
      indices = break_indices,
      counts = counts,
      objective = observed$dp[segment_count, n],
      objective_name = "total_absolute_deviation_to_segment_medians",
      gap = gap[selected_idx],
      gap_se = gap_se[selected_idx],
      method = method_str
    ),
    class = "segment_data"
  )
}
