#' clusterOGs
#'
#' @description Function used to cluster the Orthologous Groups (OGs) using the semantic similarity of their Gene Ontology annotations (GOs).
#'
#' @param simMat matrix. Similarity matrix generated with the getSimMat function.
#' @param clusterits numeric. Maximal number of iterations to be used by apcluster. Defaults to 20000.
#'
#' @return A data.frame with the clusters which can be used to generate the plot.
#'
#'
#' @examples
#'
#' clusters <- clusterOGs(sampleSimMat)
#'
#' @export


clusterOGs <- function(simMat, clusterits = 20000){
  if(!is.matrix(simMat) | nrow(simMat) != ncol(simMat) | sum(row.names(simMat) != colnames(simMat)) > 0)
    stop("simMat does not contain an appropiate matrix, please build one using te getSimMat function.")
  if(!is.numeric(clusterits) | clusterits < 1 | length(clusterits) > 1)
    stop("clusterits should be a single number greater than 0.")
  xx <- stats::cmdscale(as.matrix(stats::as.dist(1 - simMat)), eig = TRUE, k = 2)
  cluster_apc_linexcsim <- apcluster::apcluster(simMat, details = TRUE, maxits = clusterits)
  clust_lines <- data.frame(xx$points, "Cluster_exemplar" = apcluster::labels(cluster_apc_linexcsim, type = "names"))
  pos <- sapply(1:nrow(clust_lines), function(x){
    which(rownames(clust_lines) == clust_lines$Cluster_exemplar[x])
  })
  clust_lines <- clust_lines[pos, ]
  rownames(clust_lines) <- rownames(xx$points)
  clust_lines <- clust_lines[, c(3, 1, 2)]
  colnames(clust_lines) <-c("Cluster_exemplar", "X.start", "Y.start")
  dfx <- cbind.data.frame(xx$points, clust_lines)
  colnames(dfx)[1:2] <- c("V1", "V2")
  return(dfx)
}
