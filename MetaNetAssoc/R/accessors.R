#' Extract tables from a MetaNetAssoc result
#'
#' @param x A `MetaNetResult` object.
#' @return A data frame, except `result_unmatched()` which returns a list.
#' @name result_accessors
NULL

#' @rdname result_accessors
#' @export
result_scores <- function(x) {
  if (!is(x, "MetaNetResult")) stop("`x` must be a MetaNetResult.", call. = FALSE)
  x@scores
}

#' @rdname result_accessors
#' @export
result_edges <- function(x) {
  if (!is(x, "MetaNetResult")) stop("`x` must be a MetaNetResult.", call. = FALSE)
  x@edges
}

#' @rdname result_accessors
#' @export
result_qc <- function(x) {
  if (!is(x, "MetaNetResult")) stop("`x` must be a MetaNetResult.", call. = FALSE)
  x@qc
}

#' @rdname result_accessors
#' @export
result_unmatched <- function(x) {
  if (!is(x, "MetaNetResult")) stop("`x` must be a MetaNetResult.", call. = FALSE)
  x@unmatched
}

#' @rdname result_accessors
#' @export
result_term_reduction <- function(x) {
  if (!is(x, "MetaNetResult")) stop("`x` must be a MetaNetResult.", call. = FALSE)
  x@term_reduction
}
