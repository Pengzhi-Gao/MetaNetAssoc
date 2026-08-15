`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop("`", name, "` must be a data.frame.", call. = FALSE)
  }
  invisible(TRUE)
}

.assert_columns <- function(x, columns, name) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns)) {
    stop(
      "`", name, "` is missing required column(s): ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.column_or <- function(x, column, default, n = nrow(x)) {
  if (is.null(column)) {
    return(rep(default, n))
  }
  if (!is.character(column) || length(column) != 1L || !column %in% names(x)) {
    stop("Column `", column, "` was not found.", call. = FALSE)
  }
  x[[column]]
}

.clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

.clean_gene <- function(x) {
  x <- .clean_text(x)
  toupper(x)
}

.collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  paste(x, collapse = ";")
}

.normalise_unit_score <- function(x, name, scale = c("auto", "unit", "string")) {
  scale <- match.arg(scale)
  x <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(x))) {
    stop("`", name, "` must contain only finite numeric values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("`", name, "` cannot contain negative values.", call. = FALSE)
  }
  if (scale == "string") {
    x <- x / 1000
  } else if (scale == "auto" && length(x) && max(x) > 1) {
    if (max(x) <= 1000) {
      x <- x / 1000
    } else {
      stop(
        "`", name, "` exceeds 1000. Supply a score scaled to 0-1 or STRING's 0-1000 range.",
        call. = FALSE
      )
    }
  }
  if (any(x > 1)) {
    stop("`", name, "` must be in the 0-1 range after scaling.", call. = FALSE)
  }
  x
}

.split_rows <- function(x, pattern) {
  pieces <- strsplit(as.character(x), split = pattern, perl = TRUE)
  lapply(pieces, function(z) {
    z <- .clean_gene(z)
    unique(z[!is.na(z)])
  })
}

.cartesian <- function(x, y) {
  x$.join_key <- 1L
  y$.join_key <- 1L
  out <- merge(x, y, by = ".join_key", sort = FALSE)
  out$.join_key <- NULL
  out
}

.safe_rank <- function(x) {
  rank(-x, ties.method = "min", na.last = "keep")
}

