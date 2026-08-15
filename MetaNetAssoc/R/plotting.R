.require_optional <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install `", package, "` to ", purpose, ".", call. = FALSE)
  }
}

.resolve_term <- function(scores, term) {
  matched <- unique(scores$term_id[scores$term_id == term | scores$term_name == term])
  if (!length(matched)) stop("Term was not found: ", term, call. = FALSE)
  if (length(matched) > 1L) stop("Term name is ambiguous; use a term ID.", call. = FALSE)
  matched
}

#' Plot the metabolite ranking for one term
#'
#' @param x A `MetaNetResult` object.
#' @param term A term ID or exact term name.
#' @param metric Numeric score column.
#' @param top_n Maximum metabolites to show.
#' @return A ggplot object.
#' @export
plot_term_ranking <- function(
    x,
    term,
    metric = "normalized_score",
    top_n = 20) {
  .require_optional("ggplot2", "draw score plots")
  scores <- result_scores(x)
  if (!metric %in% names(scores) || !is.numeric(scores[[metric]])) {
    stop("`metric` must name a numeric score column.", call. = FALSE)
  }
  term_id <- .resolve_term(scores, term)
  plot_data <- scores[scores$term_id == term_id, , drop = FALSE]
  plot_data <- plot_data[order(plot_data[[metric]], decreasing = TRUE), , drop = FALSE]
  plot_data <- utils::head(plot_data, top_n)
  plot_data$value <- plot_data[[metric]]
  plot_data$metabolite_label <- ifelse(
    is.na(plot_data$metabolite_name),
    plot_data$metabolite_id,
    plot_data$metabolite_name
  )
  plot_data$metabolite_label <- stats::reorder(plot_data$metabolite_label, plot_data$value)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = metabolite_label, y = value)) +
    ggplot2::geom_col(fill = "#4C78A8", width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = unique(plot_data$term_name),
      x = NULL,
      y = metric
    ) +
    ggplot2::theme_classic(base_size = 12)
}

#' Plot a metabolite-term score heatmap
#'
#' @param x A `MetaNetResult` object.
#' @param metric Numeric score column.
#' @param top_n_metabolites,top_n_terms Optional display limits based on maximum
#'   score across the opposite dimension.
#' @return A ggplot object.
#' @export
plot_score_heatmap <- function(
    x,
    metric = "normalized_score",
    top_n_metabolites = NULL,
    top_n_terms = NULL) {
  .require_optional("ggplot2", "draw score plots")
  plot_data <- result_scores(x)
  if (!metric %in% names(plot_data) || !is.numeric(plot_data[[metric]])) {
    stop("`metric` must name a numeric score column.", call. = FALSE)
  }
  plot_data$value <- plot_data[[metric]]
  plot_data$metabolite_label <- ifelse(
    is.na(plot_data$metabolite_name),
    plot_data$metabolite_id,
    plot_data$metabolite_name
  )
  if (!is.null(top_n_metabolites)) {
    best <- aggregate(value ~ metabolite_label, plot_data, max)
    keep <- utils::head(best$metabolite_label[order(best$value, decreasing = TRUE)], top_n_metabolites)
    plot_data <- plot_data[plot_data$metabolite_label %in% keep, , drop = FALSE]
  }
  if (!is.null(top_n_terms)) {
    best <- aggregate(value ~ term_id, plot_data, max)
    keep <- utils::head(best$term_id[order(best$value, decreasing = TRUE)], top_n_terms)
    plot_data <- plot_data[plot_data$term_id %in% keep, , drop = FALSE]
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = metabolite_label, y = term_name, fill = value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#B2182B") +
    ggplot2::labs(x = NULL, y = NULL, fill = metric) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Build a four-layer association graph
#'
#' @param x A `MetaNetResult` object.
#' @param metabolites,terms Optional IDs or exact names to retain.
#' @param min_edge_weight Minimum path weight.
#' @return An igraph object with typed nodes and edges.
#' @export
association_graph <- function(x, metabolites = NULL, terms = NULL, min_edge_weight = 0) {
  .require_optional("igraph", "build an association graph")
  edges <- result_edges(x)
  if (!is.null(metabolites)) {
    edges <- edges[
      edges$metabolite_id %in% metabolites | edges$metabolite_name %in% metabolites,
      , drop = FALSE
    ]
  }
  if (!is.null(terms)) {
    edges <- edges[edges$term_id %in% terms | edges$term_name %in% terms, , drop = FALSE]
  }
  edges <- edges[edges$edge_weight >= min_edge_weight, , drop = FALSE]
  if (!nrow(edges)) stop("No paths remained for the requested graph.", call. = FALSE)

  met_node <- paste0("met:", edges$metabolite_key)
  relation_node <- paste0("rel:", edges$relation_gene)
  enrich_node <- paste0("enrich:", edges$enrich_gene)
  term_node <- paste0("term:", edges$term_id)
  graph_edges <- unique(rbind(
    data.frame(from = met_node, to = relation_node, weight = edges$mpi_weight, edge_type = "MPI"),
    data.frame(from = relation_node, to = enrich_node, weight = edges$ppi_score, edge_type = edges$path_type),
    data.frame(from = enrich_node, to = term_node, weight = edges$edge_weight, edge_type = "membership")
  ))
  nodes <- unique(rbind(
    data.frame(
      name = met_node,
      label = ifelse(is.na(edges$metabolite_name), edges$metabolite_id, edges$metabolite_name),
      node_type = "metabolite"
    ),
    data.frame(name = relation_node, label = edges$relation_gene, node_type = "relation_gene"),
    data.frame(name = enrich_node, label = edges$enrich_gene, node_type = "term_gene"),
    data.frame(name = term_node, label = edges$term_name, node_type = "term")
  ))
  igraph::graph_from_data_frame(graph_edges, directed = TRUE, vertices = nodes)
}

#' Plot a four-layer association network
#'
#' @inheritParams association_graph
#' @param ... Additional arguments passed to `igraph::plot.igraph`.
#' @return The plotted igraph object, invisibly.
#' @export
plot_association_network <- function(
    x,
    metabolites = NULL,
    terms = NULL,
    min_edge_weight = 0,
    ...) {
  graph <- association_graph(x, metabolites, terms, min_edge_weight)
  type_colors <- c(
    metabolite = "#6A51A3",
    relation_gene = "#3182BD",
    term_gene = "#31A354",
    term = "#E6550D"
  )
  vertex_colors <- unname(type_colors[igraph::V(graph)$node_type])
  igraph::plot.igraph(
    graph,
    vertex.color = vertex_colors,
    vertex.label = igraph::V(graph)$label,
    vertex.size = 18,
    vertex.label.cex = 0.75,
    edge.arrow.size = 0.35,
    ...
  )
  invisible(graph)
}
