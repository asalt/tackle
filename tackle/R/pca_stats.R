# Heteroscedastic one-way inference for selected PCA score coordinates.

pca_euclidean_r2 <- function(Y, group) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  group <- as.character(group)
  if (nrow(Y) == 0 || nrow(Y) != length(group) || any(!is.finite(Y))) {
    return(NA_real_)
  }
  grand <- colMeans(Y)
  ss_total <- sum(rowSums((Y - rep(grand, each = nrow(Y)))^2))
  if (!is.finite(ss_total) || ss_total <= 0) return(NA_real_)
  ss_within <- sum(vapply(
    split(seq_len(nrow(Y)), group),
    function(idx) {
      subset <- Y[idx, , drop = FALSE]
      centroid <- colMeans(subset)
      sum(rowSums((subset - rep(centroid, each = nrow(subset)))^2))
    },
    numeric(1)
  ))
  max(0, min(1, 1 - ss_within / ss_total))
}

pca_pairwise_centroid_geometry <- function(Y, group) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  group <- as.character(group)
  empty <- list(
    geometry_group_a = NA_character_,
    geometry_group_b = NA_character_,
    centroid_distance = NA_real_,
    rms_radius_a = NA_real_,
    rms_radius_b = NA_real_,
    pooled_rms_radius = NA_real_,
    standardized_separation = NA_real_
  )
  levels <- unique(group)
  if (nrow(Y) == 0 || nrow(Y) != length(group) || any(!is.finite(Y)) ||
      length(levels) != 2) {
    return(empty)
  }

  left <- Y[group == levels[[1]], , drop = FALSE]
  right <- Y[group == levels[[2]], , drop = FALSE]
  if (nrow(left) == 0 || nrow(right) == 0) return(empty)
  left_centroid <- colMeans(left)
  right_centroid <- colMeans(right)
  centroid_distance <- sqrt(sum((left_centroid - right_centroid)^2))
  radius_a <- sqrt(mean(rowSums((left - rep(left_centroid, each = nrow(left)))^2)))
  radius_b <- sqrt(mean(rowSums((right - rep(right_centroid, each = nrow(right)))^2)))
  pooled_radius <- sqrt((radius_a^2 + radius_b^2) / 2)
  standardized <- if (is.finite(pooled_radius) && pooled_radius > 0) {
    centroid_distance / pooled_radius
  } else {
    NA_real_
  }
  list(
    geometry_group_a = levels[[1]],
    geometry_group_b = levels[[2]],
    centroid_distance = centroid_distance,
    rms_radius_a = radius_a,
    rms_radius_b = radius_b,
    pooled_rms_radius = pooled_radius,
    standardized_separation = standardized
  )
}

.pca_matrix_rank <- function(value) {
  singular_values <- svd(value, nu = 0, nv = 0)$d
  if (length(singular_values) == 0 || max(singular_values) == 0) return(0L)
  tolerance <- max(dim(value)) * max(singular_values) * .Machine$double.eps
  sum(singular_values > tolerance)
}

.pca_block_diag <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1))
  result <- matrix(0, nrow = sum(sizes), ncol = sum(sizes))
  start <- 1L
  for (index in seq_along(blocks)) {
    stop_at <- start + sizes[[index]] - 1L
    result[start:stop_at, start:stop_at] <- blocks[[index]]
    start <- stop_at + 1L
  }
  result
}

.pca_failed_welch_james <- function(status, message) {
  list(
    statistic = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    p_value = NA_real_,
    status = status,
    message = message
  )
}

pca_welch_james <- function(Y, group, condition_limit = 1e12) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  group <- as.character(group)
  if (nrow(Y) != length(group) || ncol(Y) < 1 || any(!is.finite(Y))) {
    return(.pca_failed_welch_james(
      "invalid_values",
      "The selected PCA score matrix contains non-finite values."
    ))
  }

  levels <- unique(group)
  if (length(levels) < 2) {
    return(.pca_failed_welch_james(
      "insufficient_groups", "At least two nonempty groups are required."
    ))
  }
  grouped <- lapply(levels, function(level) Y[group == level, , drop = FALSE])
  group_sizes <- vapply(grouped, nrow, integer(1))
  if (any(group_sizes < 2)) {
    return(.pca_failed_welch_james(
      "insufficient_group_size",
      "Every group needs at least two samples for a covariance estimate."
    ))
  }

  dimensions <- ncol(Y)
  means <- unlist(lapply(grouped, colMeans), use.names = FALSE)
  covariance_blocks <- Map(
    function(values, size) {
      covariance <- as.matrix(stats::cov(values))
      if (!identical(dim(covariance), c(dimensions, dimensions))) {
        covariance <- matrix(covariance, nrow = dimensions, ncol = dimensions)
      }
      covariance / size
    },
    grouped,
    group_sizes
  )
  sigma <- .pca_block_diag(covariance_blocks)

  contrasts <- matrix(0, nrow = length(levels) - 1L, ncol = length(levels))
  contrasts[, 1] <- 1
  for (row in seq_len(nrow(contrasts))) contrasts[row, row + 1L] <- -1
  contrast <- kronecker(contrasts, diag(dimensions))
  numerator_df <- .pca_matrix_rank(contrast)
  contrast_covariance <- contrast %*% sigma %*% t(contrast)
  contrast_covariance <- (contrast_covariance + t(contrast_covariance)) / 2
  covariance_rank <- .pca_matrix_rank(contrast_covariance)
  if (covariance_rank < numerator_df) {
    return(.pca_failed_welch_james(
      "singular_covariance",
      paste0(
        "The covariance contrast has rank ", covariance_rank,
        ", below the required ", numerator_df, "."
      )
    ))
  }
  condition <- kappa(contrast_covariance, exact = TRUE)
  if (!is.finite(condition) || condition > condition_limit) {
    return(.pca_failed_welch_james(
      "ill_conditioned_covariance",
      paste0("The covariance contrast condition number is ", signif(condition, 3), ".")
    ))
  }

  inverse_covariance <- solve(contrast_covariance)
  hypothesis <- contrast %*% means
  twj <- as.numeric(t(hypothesis) %*% inverse_covariance %*% hypothesis)
  projection <- sigma %*% t(contrast) %*% inverse_covariance %*% contrast
  correction_sum <- 0
  for (group_index in seq_along(group_sizes)) {
    selector <- matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
    start <- (group_index - 1L) * dimensions + 1L
    stop_at <- start + dimensions - 1L
    selector[start:stop_at, start:stop_at] <- diag(dimensions)
    component <- projection %*% selector
    correction_sum <- correction_sum + (
      sum(diag(component))^2 + sum(diag(component %*% component))
    ) / (group_sizes[[group_index]] - 1)
  }
  correction <- correction_sum / 2
  if (!is.finite(correction) || correction <= 0) {
    return(.pca_failed_welch_james(
      "invalid_df_correction",
      "The Welch-James approximate-df correction is not positive and finite."
    ))
  }

  denominator_df <- numerator_df * (numerator_df + 2) / (3 * correction)
  divisor <- numerator_df + 2 * correction - 6 * correction / (numerator_df + 2)
  statistic <- twj / divisor
  p_value <- stats::pf(
    statistic, df1 = numerator_df, df2 = denominator_df, lower.tail = FALSE
  )
  if (any(!is.finite(c(statistic, denominator_df, p_value)))) {
    return(.pca_failed_welch_james(
      "nonfinite_result", "The Welch-James approximation produced a non-finite result."
    ))
  }
  list(
    statistic = statistic,
    numerator_df = as.numeric(numerator_df),
    denominator_df = denominator_df,
    p_value = p_value,
    status = "ok",
    message = ""
  )
}

.pca_adjust_pvalues <- function(values, method) {
  result <- rep(NA_real_, length(values))
  keep <- is.finite(values)
  if (!any(keep)) return(result)
  normalized <- tolower(method)
  if (normalized == "bh") normalized <- "BH"
  result[keep] <- stats::p.adjust(values[keep], method = normalized)
  result
}

.pca_clean_group <- function(values) {
  result <- trimws(as.character(values))
  missing <- is.na(values) | tolower(result) %in% c("", "na", "nan", "<na>", "null")
  result[missing] <- NA_character_
  result
}

.pca_group_sizes <- function(group) {
  counts <- table(factor(group, levels = unique(group)))
  paste(paste0(names(counts), "=", as.integer(counts)), collapse = ";")
}

.pca_empty_pairwise <- function() {
  data.frame(
    group_field = character(), scope = character(),
    pcs = character(), n_pcs = integer(), explained_variance_pct = numeric(),
    group_a = character(), group_b = character(), n_a = integer(), n_b = integer(),
    centroid_distance = numeric(), rms_radius_a = numeric(), rms_radius_b = numeric(),
    pooled_rms_radius = numeric(), standardized_separation = numeric(),
    r2 = numeric(), welch_james_f = numeric(), numerator_df = numeric(),
    denominator_df = numeric(), p_value = numeric(), p_adjust_method = character(),
    p_adj = numeric(), p_adj_all_scopes = numeric(), status = character(),
    message = character(), method = character(), stringsAsFactors = FALSE
  )
}

.pca_scope_selection <- function(scope) {
  if (is.null(scope$selection)) return("explicit")
  value <- as.character(unlist(scope$selection, use.names = FALSE))
  if (length(value) == 0 || is.na(value[[1]]) || !nzchar(value[[1]])) {
    return("explicit")
  }
  value[[1]]
}

.pca_select_scope_data <- function(scores, clean_group, pcs) {
  selected <- as.matrix(scores[, pcs, drop = FALSE])
  storage.mode(selected) <- "double"
  keep <- !is.na(clean_group) & apply(selected, 1, function(row) all(is.finite(row)))
  selected <- selected[keep, , drop = FALSE]
  group <- clean_group[keep]
  list(selected = selected, group = group, levels = unique(group))
}

.pca_largest_common_estimable_leading_pcs <- function(scores, clean_group, pcs) {
  for (component_count in rev(seq_along(pcs))) {
    candidate <- pcs[seq_len(component_count)]
    selected_data <- .pca_select_scope_data(scores, clean_group, candidate)
    omnibus <- pca_welch_james(selected_data$selected, selected_data$group)
    if (!identical(omnibus$status, "ok")) next

    pairwise_ok <- TRUE
    if (length(selected_data$levels) > 2) {
      pairs <- utils::combn(selected_data$levels, 2, simplify = FALSE)
      for (pair in pairs) {
        pair_keep <- selected_data$group %in% pair
        pair_test <- pca_welch_james(
          selected_data$selected[pair_keep, , drop = FALSE],
          selected_data$group[pair_keep]
        )
        if (!identical(pair_test$status, "ok")) {
          pairwise_ok <- FALSE
          break
        }
      }
    }
    if (pairwise_ok) return(candidate)
  }
  NULL
}

pca_analyze_separation <- function(
  scores, metadata, group_fields, scopes, p_adjust_method = "holm"
) {
  scores <- as.data.frame(scores, check.names = FALSE)
  metadata <- as.data.frame(metadata, check.names = FALSE, stringsAsFactors = FALSE)
  metadata <- metadata[rownames(scores), , drop = FALSE]
  component_variances <- vapply(
    scores,
    function(values) stats::var(as.numeric(values), na.rm = TRUE),
    numeric(1)
  )
  total_variance <- sum(component_variances[is.finite(component_variances)])
  explained_variance_pct <- function(pcs) {
    selected_variances <- component_variances[pcs]
    if (!is.finite(total_variance) || total_variance <= 0 ||
        any(!is.finite(selected_variances))) {
      return(NA_real_)
    }
    sum(selected_variances) / total_variance * 100
  }
  omnibus_rows <- list()
  pairwise_rows <- list()

  for (field in as.character(unlist(group_fields, use.names = FALSE))) {
    if (!field %in% colnames(metadata)) {
      stop("PCA test metadata field was not found: ", field)
    }
    clean_group <- .pca_clean_group(metadata[[field]])
    field_omnibus_start <- length(omnibus_rows)
    selections <- vapply(scopes, .pca_scope_selection, character(1))
    ordered_scopes <- c(
      scopes[selections != "leading_estimable"],
      scopes[selections == "leading_estimable"]
    )
    for (scope in ordered_scopes) {
      scope_name <- as.character(unlist(scope$name, use.names = FALSE)[[1]])
      pcs <- as.character(unlist(scope$pcs, use.names = FALSE))
      selection <- .pca_scope_selection(scope)
      if (identical(selection, "leading_estimable")) {
        pcs <- .pca_largest_common_estimable_leading_pcs(scores, clean_group, pcs)
        if (is.null(pcs)) next
        scope_name <- "leading_estimable_pcs"
        pcs_text <- paste(pcs, collapse = ",")
        omnibus_indexes <- if (length(omnibus_rows) > field_omnibus_start) {
          seq.int(field_omnibus_start + 1L, length(omnibus_rows))
        } else {
          integer()
        }
        if (length(omnibus_indexes) > 0) {
          duplicate_indexes <- omnibus_indexes[vapply(
            omnibus_rows[omnibus_indexes],
            function(row) identical(row$pcs[[1]], pcs_text),
            logical(1)
          )]
          if (length(duplicate_indexes) > 0) {
            next
          }
        }
      }

      selected_data <- .pca_select_scope_data(scores, clean_group, pcs)
      selected <- selected_data$selected
      group <- selected_data$group
      levels <- selected_data$levels
      test <- pca_welch_james(selected, group)
      geometry <- pca_pairwise_centroid_geometry(selected, group)
      omnibus_rows[[length(omnibus_rows) + 1L]] <- data.frame(
        group_field = field,
        scope = scope_name,
        pcs = paste(pcs, collapse = ","),
        n_pcs = length(pcs),
        explained_variance_pct = explained_variance_pct(pcs),
        geometry_group_a = geometry$geometry_group_a,
        geometry_group_b = geometry$geometry_group_b,
        centroid_distance = geometry$centroid_distance,
        rms_radius_a = geometry$rms_radius_a,
        rms_radius_b = geometry$rms_radius_b,
        pooled_rms_radius = geometry$pooled_rms_radius,
        standardized_separation = geometry$standardized_separation,
        n_samples = nrow(selected),
        n_groups = length(levels),
        group_sizes = .pca_group_sizes(group),
        r2 = pca_euclidean_r2(selected, group),
        welch_james_f = test$statistic,
        numerator_df = test$numerator_df,
        denominator_df = test$denominator_df,
        p_value = test$p_value,
        p_adjust_method = p_adjust_method,
        p_adj = NA_real_,
        status = test$status,
        message = test$message,
        method = "Johansen Welch-James ADF",
        stringsAsFactors = FALSE
      )

      if (length(levels) <= 2) next
      combinations <- utils::combn(levels, 2, simplify = FALSE)
      for (pair in combinations) {
        pair_keep <- group %in% pair
        pair_values <- selected[pair_keep, , drop = FALSE]
        pair_group <- group[pair_keep]
        pair_test <- pca_welch_james(pair_values, pair_group)
        pair_geometry <- pca_pairwise_centroid_geometry(pair_values, pair_group)
        counts <- table(factor(pair_group, levels = pair))
        pairwise_rows[[length(pairwise_rows) + 1L]] <- data.frame(
          group_field = field,
          scope = scope_name,
          pcs = paste(pcs, collapse = ","),
          n_pcs = length(pcs),
          explained_variance_pct = explained_variance_pct(pcs),
          group_a = pair[[1]],
          group_b = pair[[2]],
          n_a = as.integer(counts[[1]]),
          n_b = as.integer(counts[[2]]),
          centroid_distance = pair_geometry$centroid_distance,
          rms_radius_a = pair_geometry$rms_radius_a,
          rms_radius_b = pair_geometry$rms_radius_b,
          pooled_rms_radius = pair_geometry$pooled_rms_radius,
          standardized_separation = pair_geometry$standardized_separation,
          r2 = pca_euclidean_r2(pair_values, pair_group),
          welch_james_f = pair_test$statistic,
          numerator_df = pair_test$numerator_df,
          denominator_df = pair_test$denominator_df,
          p_value = pair_test$p_value,
          p_adjust_method = p_adjust_method,
          p_adj = NA_real_,
          p_adj_all_scopes = NA_real_,
          status = pair_test$status,
          message = pair_test$message,
          method = "Johansen Welch-James ADF",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  omnibus <- if (length(omnibus_rows)) do.call(rbind, omnibus_rows) else data.frame()
  pairwise <- if (length(pairwise_rows)) do.call(rbind, pairwise_rows) else .pca_empty_pairwise()
  for (field in unique(omnibus$group_field)) {
    idx <- which(omnibus$group_field == field)
    omnibus$p_adj[idx] <- .pca_adjust_pvalues(omnibus$p_value[idx], p_adjust_method)
  }
  if (nrow(pairwise)) {
    families <- interaction(pairwise$group_field, pairwise$scope, drop = TRUE)
    for (family in unique(families)) {
      idx <- which(families == family)
      pairwise$p_adj[idx] <- .pca_adjust_pvalues(pairwise$p_value[idx], p_adjust_method)
    }
    for (field in unique(pairwise$group_field)) {
      idx <- which(pairwise$group_field == field)
      pairwise$p_adj_all_scopes[idx] <- .pca_adjust_pvalues(
        pairwise$p_value[idx], p_adjust_method
      )
    }
  }
  list(omnibus = omnibus, pairwise = pairwise)
}

.pca_adjustment_caption_label <- function(method) {
  key <- tolower(gsub("-", "_", trimws(as.character(method))))
  labels <- c(
    none = "p",
    raw = "p",
    holm = "Holm-adjusted p",
    hochberg = "Hochberg-adjusted p",
    bonferroni = "Bonferroni-adjusted p",
    bh = "BH-adjusted p",
    fdr_bh = "BH-adjusted p"
  )
  if (key %in% names(labels)) unname(labels[[key]]) else paste0(method, "-adjusted p")
}

pca_format_test_caption <- function(rows) {
  if (nrow(rows) == 0) return("")
  lines <- vapply(seq_len(nrow(rows)), function(index) {
    row <- rows[index, , drop = FALSE]
    prefix <- if (is.finite(row$r2)) {
      paste0(row$group_field, ": R²=", sprintf("%.3f", row$r2))
    } else {
      paste0(row$group_field, ": R²=NA")
    }
    line <- if (row$status == "ok") {
      p_label <- .pca_adjustment_caption_label(row$p_adjust_method[[1]])
      paste0(
        prefix, "; WJ F*(", sprintf("%.0f", row$numerator_df), ", ",
        sprintf("%.2f", row$denominator_df), ")=", sprintf("%.2f", row$welch_james_f),
        "; ", p_label, "=", format.pval(row$p_adj, digits = 3, eps = 1e-99)
      )
    } else {
      paste0(prefix, "; WJ not estimable (", row$status, ")")
    }
    if ("standardized_separation" %in% colnames(row) &&
        is.finite(row$standardized_separation[[1]])) {
      line <- paste0(
        line,
        "\nCentroid distance = ", sprintf("%.2f", row$centroid_distance),
        "\nRMS radii (", row$geometry_group_a, ", ", row$geometry_group_b, ") = ",
        sprintf("%.2f", row$rms_radius_a), ", ", sprintf("%.2f", row$rms_radius_b),
        "\nStandardized separation = ",
        sprintf("%.2f", row$standardized_separation)
      )
    }
    line
  }, character(1))
  paste(lines, collapse = "\n")
}
