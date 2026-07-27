pca_caption_text_grob <- function(grob) {
  if (inherits(grob, "text")) return(grob)
  children <- grob$children
  if (is.null(children) || length(children) == 0) return(NULL)
  for (child in children) {
    found <- pca_caption_text_grob(child)
    if (!is.null(found)) return(found)
  }
  NULL
}

pca_caption_metrics <- function(plot, fig_width) {
  table <- ggplot2::ggplotGrob(plot)
  caption_index <- which(table$layout$name == "caption")
  if (length(caption_index) == 0) return(NULL)
  caption_index <- caption_index[[1]]
  caption_layout <- table$layout[caption_index, , drop = FALSE]
  caption_columns <- seq.int(caption_layout$l, caption_layout$r)
  outside_columns <- setdiff(seq_along(table$widths), caption_columns)
  outside_width <- if (length(outside_columns)) {
    grid::convertWidth(
      sum(table$widths[outside_columns]),
      "inches",
      valueOnly = TRUE
    )
  } else {
    0
  }

  caption_grob <- table$grobs[[caption_index]]
  horizontal_padding <- 0
  if (!is.null(caption_grob$widths) && length(caption_grob$widths) >= 2) {
    edge_widths <- caption_grob$widths[c(1, length(caption_grob$widths))]
    horizontal_padding <- grid::convertWidth(
      sum(edge_widths),
      "inches",
      valueOnly = TRUE
    )
  }
  text_grob <- pca_caption_text_grob(caption_grob)
  if (is.null(text_grob)) return(NULL)

  list(
    table = table,
    caption_grob = caption_grob,
    text_gp = text_grob$gp,
    available_width = max(
      0,
      as.numeric(fig_width) - outside_width - horizontal_padding
    ),
    height = grid::convertHeight(
      grid::grobHeight(caption_grob),
      "inches",
      valueOnly = TRUE
    )
  )
}

pca_caption_line_width <- function(text, gp) {
  if (!nzchar(text)) return(0)
  grid::convertWidth(
    grid::grobWidth(grid::textGrob(text, gp = gp)),
    "inches",
    valueOnly = TRUE
  )
}

pca_plot_width_metrics <- function(plot) {
  table <- ggplot2::ggplotGrob(plot)
  fixed_width <- grid::convertWidth(
    sum(table$widths),
    "inches",
    valueOnly = TRUE
  )
  guide_indices <- grep("^guide-box", table$layout$name)
  guide_widths <- vapply(
    guide_indices,
    function(index) {
      guide <- table$grobs[[index]]
      if (inherits(guide, "zeroGrob")) return(0)
      grid::convertWidth(
        grid::grobWidth(guide),
        "inches",
        valueOnly = TRUE
      )
    },
    numeric(1)
  )
  list(
    fixed_width = as.numeric(fixed_width),
    legend_width = if (length(guide_widths)) max(guide_widths) else 0
  )
}

pca_split_caption_token <- function(token, max_width, gp) {
  if (pca_caption_line_width(token, gp) <= max_width) return(token)
  characters <- strsplit(token, "", fixed = TRUE)[[1]]
  chunks <- character()
  current <- ""
  for (character in characters) {
    candidate <- paste0(current, character)
    if (nzchar(current) &&
        pca_caption_line_width(candidate, gp) > max_width) {
      chunks <- c(chunks, current)
      current <- character
    } else {
      current <- candidate
    }
  }
  c(chunks, current)
}

pca_wrap_caption_line <- function(line, max_width, gp) {
  line <- trimws(line)
  if (!nzchar(line)) return("")
  words <- strsplit(line, "[[:space:]]+")[[1]]
  wrapped <- character()
  current <- ""
  for (word in words) {
    if (pca_caption_line_width(word, gp) > max_width) {
      if (nzchar(current)) {
        wrapped <- c(wrapped, current)
        current <- ""
      }
      chunks <- pca_split_caption_token(word, max_width, gp)
      if (length(chunks) > 1) wrapped <- c(wrapped, chunks[-length(chunks)])
      current <- chunks[[length(chunks)]]
      next
    }
    candidate <- if (nzchar(current)) paste(current, word) else word
    if (pca_caption_line_width(candidate, gp) <= max_width) {
      current <- candidate
    } else {
      wrapped <- c(wrapped, current)
      current <- word
    }
  }
  c(wrapped, current)
}

pca_wrap_caption_text <- function(text, max_width, gp) {
  paragraphs <- strsplit(as.character(text), "\n", fixed = TRUE)[[1]]
  wrapped <- unlist(
    lapply(
      paragraphs,
      pca_wrap_caption_line,
      max_width = max_width,
      gp = gp
    ),
    use.names = FALSE
  )
  paste(wrapped, collapse = "\n")
}

pca_prepare_plot_for_output <- function(
  plot,
  fig_width,
  fig_height,
  expand_height = TRUE,
  expand_width = FALSE,
  minimum_panel_width = 3.75
) {
  grDevices::pdf(NULL, width = fig_width, height = fig_height)
  on.exit(grDevices::dev.off(), add = TRUE)

  width_metrics <- pca_plot_width_metrics(plot)
  resolved_width <- as.numeric(fig_width)
  if (isTRUE(expand_width)) {
    required_width <- width_metrics$fixed_width + as.numeric(minimum_panel_width)
    resolved_width <- max(resolved_width, ceiling(required_width * 4) / 4)
  }

  caption <- plot$labels$caption
  if (is.null(caption) || length(caption) == 0 ||
      is.na(caption[[1]]) || !nzchar(caption[[1]])) {
    return(list(
      plot = plot,
      fig_width = resolved_width,
      fig_height = as.numeric(fig_height),
      fixed_width = width_metrics$fixed_width,
      legend_width = width_metrics$legend_width,
      caption_width = NA_real_,
      original_line_count = 0L,
      wrapped_line_count = 0L,
      max_line_width = NA_real_
    ))
  }

  original_metrics <- pca_caption_metrics(plot, resolved_width)
  if (is.null(original_metrics) || original_metrics$available_width <= 0) {
    return(list(
      plot = plot,
      fig_width = resolved_width,
      fig_height = as.numeric(fig_height),
      fixed_width = width_metrics$fixed_width,
      legend_width = width_metrics$legend_width,
      caption_width = NA_real_,
      original_line_count = length(strsplit(caption[[1]], "\n", fixed = TRUE)[[1]]),
      wrapped_line_count = length(strsplit(caption[[1]], "\n", fixed = TRUE)[[1]]),
      max_line_width = NA_real_
    ))
  }

  wrapped_caption <- pca_wrap_caption_text(
    caption[[1]],
    max_width = original_metrics$available_width,
    gp = original_metrics$text_gp
  )
  wrapped_plot <- plot + ggplot2::labs(caption = wrapped_caption)
  wrapped_metrics <- pca_caption_metrics(wrapped_plot, resolved_width)
  original_lines <- strsplit(caption[[1]], "\n", fixed = TRUE)[[1]]
  wrapped_lines <- strsplit(wrapped_caption, "\n", fixed = TRUE)[[1]]
  line_widths <- vapply(
    wrapped_lines,
    pca_caption_line_width,
    numeric(1),
    gp = wrapped_metrics$text_gp
  )
  extra_height <- if (isTRUE(expand_height)) {
    max(0, wrapped_metrics$height - original_metrics$height)
  } else {
    0
  }

  list(
    plot = wrapped_plot,
    fig_width = resolved_width,
    fig_height = as.numeric(fig_height) + extra_height,
    fixed_width = width_metrics$fixed_width,
    legend_width = width_metrics$legend_width,
    caption_width = wrapped_metrics$available_width,
    original_line_count = length(original_lines),
    wrapped_line_count = length(wrapped_lines),
    max_line_width = if (length(line_widths)) max(line_widths) else 0
  )
}
