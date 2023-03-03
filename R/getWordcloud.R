#' getWordcloud
#'
#' @description This function is used to get the wordcloud for a list of OGs of interest.
#'
#' @param GOsfromOGs data.frame. Contains two columns, the header of the first should be GOs and the header of the second should be OGs. The first column should be a list of GOs and the second one a list of the OGs to which each GO belongs. No empty fields are allowed.
#' @param ont character. The ontology of the GO terms to be evaluated. It must include one or more of the following terms: "BP", "MF" or "CC" (such as c("BP", "MF", "CC"). Defaults to "BP".
#' @param clusters data.frame. Obtained with the clusterOGs function.
#' @param toExclude character. Path to the file containing the terms to be excluded. Defaults to NULL.
#'
#' @return data.frame with the wordcloud for each cluster exemplar.
#'
#' @examples
#'
#'
#' wcl <- getWordcloud(sampleGOsfromOGs, clusters = sampleClusters, toExclude = NULL)
#'
#'
#' @export


getWordcloud <- function(GOsfromOGs, ont = "BP", clusters, toExclude = NULL){
  if(!is.character(ont))
    stop("ont must be an object of the class character.")
  if(!ont %in% c("BP", "CC", "MF") | length(ont) > 1)
    stop("ont can only take one of the following values: BP, CC or MF.")
  message("Processing GOsfromOGs:")
  pangenome <- GOsfromOGs
  if(sum(colnames(pangenome) == c("GOs", "OGs")) != 2)
    stop(paste("The column names in the ", GOsfromOGs, " data.frame should be GOs and OGs.", sep = ""))
  if(sum(pangenome$GOs %in% c("", " ")) > 0 | sum(pangenome$OGs %in% c("", " ")) > 0)
    stop(paste("The ", GOsfromOGs, " data.frame cannot contain empty characters or whitespaces.", sep = ""))
  category <- clusterProfiler::go2ont(unique(pangenome$GOs))
  pangenome_subset <- pangenome[pangenome$GOs %in% category$go_id[category$Ontology == ont],]
  clust <- split(rownames(clusters), clusters$Cluster_exemplar)
  clusters_annotations <- lapply(clust, function(x){
    unique(pangenome_subset[pangenome_subset$OGs %in% x, 1])
  })
  clusters_enrichment <- lapply(clusters_annotations, simplifyEnrichment::keyword_enrichment_from_GO)
  clusters_enrichment <- lapply(clusters_enrichment, function(x){
    x[x$padj <= 0.05, c(1, 5)]
  })
  if(!is.null(toExclude)){
    exclude <- unlist(c(simplifyEnrichment:::GO_EXCLUDE_WORDS, utils::read.table(toExclude)))
    clusters_enrichment <- lapply(clusters_enrichment, function(x){
      x$padj <- -log(x$padj)
      x[!x$keyword %in% exclude, ]
    })
  } else{
    clusters_enrichment <- lapply(clusters_enrichment, function(x){
      x$padj <- -log(x$padj)
      x
    })
  }
  clusters_enrichment <- lapply(clusters_enrichment, function(x){
    x[1:10, ]
  })
  clusters_enrichment <- lapply(clusters_enrichment, function(x){
    x[is.na(x$keyword), 1] <- ""
    return(x)
  })
  clusters_enrichment_label <- lapply(clusters_enrichment, function(x){
    ifelse(x$keyword[1] == "", "" , paste(x$keyword, collapse = "\n"))
  })
  clusters_enrichment_label <- do.call(rbind.data.frame, clusters_enrichment_label)
  row.names(clusters_enrichment_label) <- names(clusters_enrichment)
  colnames(clusters_enrichment_label) <- "wordcloud"
  finalPos <- sapply(clusters$Cluster_exemplar, function(x){
    pos <- which(row.names(clusters_enrichment_label) == x)
  })
  clusters_enrichment_label <- clusters_enrichment_label[finalPos, , drop = FALSE]
  row.names(clusters_enrichment_label) <- row.names(clusters)
  return(clusters_enrichment_label)
}
