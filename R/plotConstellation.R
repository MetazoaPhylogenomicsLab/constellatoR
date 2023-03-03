#' plotConstellation
#'
#' @description Function used to plot the Orthologous Groups (OGs) clustered by their semantic similarity. A polygon comprising all of the OGs in a cluster and the representative member of that cluster is plotted and their functions are displayed.
#'
#' @param clusters data.frame. Clusters generated with the clusterOGs function.
#' @param annotations data.frame. Annotations for each OG. The rownames should be the OGs and it must contain one or more columns with the annotation with their headers.
#' @param term character.
#' @param color character.
#' @param size character.
#'
#' @return A ggplot2 object
#'
#'
#' @examples
#'
#' annot <- merge(sampleAnnotation, sampleWcl, by.x = 0, by.y = 0)
#' row.names(annot) <- annot[, 1]
#' annot <- annot[, -1]
#' groups <- c("PMAX and LLON", "PMAX and LLON", "PMAX and LLON", "PMAX and LLON", "PMAX and LLON", "PMAX, SMED and LLON", "PMAX, SMED and LLON", "PMAX, SMED and LLON", "PMAX and SMED", "PMAX, SMED and LLON", "PMAX, SMED and LLON", "PMAX and LLON")
#' size <- c(2, 5, 2, 2, 2, 4, 3, 3, 2, 3, 3, 3)
#' annot <- data.frame(annot, groups = groups, size = size)
#' plotConstellation(sampleClusters, annot, term = "BP1", color = "groups", size = "size")
#'
#' @export


plotConstellation <- function(clusters, annotations, term = "BP1", color = NULL, size = NULL){
  hulls <- clusters[, 1:4] %>% dplyr::group_by(Cluster_exemplar) %>% dplyr::slice(grDevices::chull(V1, V2))
  hull_col_pal <- stats::setNames(paletteer::paletteer_d("rcartocolor::Prism", n = dplyr::n_distinct(clusters$Cluster_exemplar), type="continuous"), unique(clusters$Cluster_exemplar))
  hull_col <- hull_col_pal[hulls$Cluster_exemplar]
  clusters <- merge(clusters, annotations, by.x = 0, by.y = 0)
  if(!term %in% colnames(clusters))
    stop(paste("The term ", term, " was not found as a column name in the matrix annotations", sep = ""))
  row.names(clusters) <- clusters$Row.names
  clusters <- clusters[, -1]
  p <- internalPlot(clusters, term, hulls, hull_col, color, size)
  return(p)
}

#################################################################################################
# internalPlot: function used to prepare the object for the plot. The main reason for making it an independent function is passing the label argument to the aes in geom_label_repel. Given that aes quotates its arguments it was not straightforward to pass the ont argument to specify the ontology of the labels.
#
# clusters: the object that contains the data including the coordinates and the name of the representative cluster and the terms that better describe it on different ontological domains.
# term: the ontology to use for labeling each representative cluster. It must be one of "BP", "MF" or "CC".
# hulls: coordinates for the polygon that will group the OGs that belong to the same cluster.
# hull_col: colors to be used for the hulls.
# color:
# size:
#################################################################################################

internalPlot <- function(clusters, term, hulls, hull_col, color, size){
  if(!is.null(color)){
    p <- ggplot2::ggplot(clusters, ggplot2::aes(x = V1, y = V2, color = .data[[color]]))
  } else{
    p <- ggplot2::ggplot(clusters, ggplot2::aes(x = V1, y = V2))
  }
  p <- p + ggplot2::geom_polygon(data = hulls, ggplot2::aes(group = Cluster_exemplar), colour = NA, fill = hull_col, alpha = 0.2, show.legend = F) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_segment(data = clusters, ggplot2::aes(x = X.start, y = Y.start, xend = V1, yend = V2, color = Cluster_exemplar), lwd = 0.3, alpha = 0.4) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = paletteer::paletteer_d("ggthemes::excel_Crop", direction = 1, type = "continuous", n = dplyr::n_distinct(clusters$Cluster_exemplar))) +
    ggplot2::scale_x_continuous(name = "MDS1") +
    ggplot2::scale_y_continuous(name = "MDS2")  +
    ggnewscale::new_scale_color()
  if(!is.null(color) & !is.null(size)){
    p <- p + ggplot2::geom_point(alpha = 0.6, ggplot2::aes(size = .data[[size]], color = .data[[color]]))
  } else if(!is.null(color)){
    p <- p + ggplot2::geom_point(alpha = 0.6, ggplot2::aes(color = .data[[color]]))
  } else if(!is.null(size)){
    p <- p + ggplot2::geom_point(alpha = 0.6, ggplot2::aes(size = .data[[size]]))
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.6)
  }
  p <- p + ggplot2::scale_color_manual(values = paletteer::paletteer_d("dutchmasters::pearl_earring", direction = 1, type = "continuous", n = dplyr::n_distinct(clusters[, colnames(clusters) %in% color]))) +
    ggplot2::scale_size(range=c(3, 9)) +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data[[term]], color = NULL, segment.linetype = 2), color="red", data = subset(clusters, row.names(clusters) == Cluster_exemplar), box.padding = grid::unit(1, "lines"), size = 3, max.overlaps = 220, alpha = 0.7, direction = "both", force = 100)
  if("wordcloud" %in% colnames(clusters)){
      p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = wordcloud, color=NULL, segment.linetype = 2), data = subset(clusters, row.names(clusters) == Cluster_exemplar), box.padding = grid::unit(1, "lines"), size = 3, max.overlaps = 220, alpha = 0.7, direction = "both", force = 100)
  }
  return(p)
}
