#' getSimMat
#'
#' @description Function to calculate the similarity between OGs by combining the semantic similarity scores between the GOs from each pair of OGs being analyzed. If needed, it can use the getSemSim function to calculate the semantic similarity between GOs and then calculate the similarity between OGs, however if you are going to calculate the similarities between different sets of OGs you might want to use the getSemSim function to obtain the full matrix, store it in an object and then use it to calculate the similarity scores between different sets of OGs to reduce times.
#'
#' @param GOsfromOGs data.frame. Contains two columns, the header of the first should be GOs and the header of the second should be OGs. The first column should be a list of GOs and the second one a list of the OGs to which each GO belongs. No empty fields are allowed.
#' @param selectedOGs character. The OGs that will be evaluated. They should be included in the OGs column in the GOsfromOGs data.frame, but elements not found will be ignored and a warning will list them.
#' @param ont character. The ontology of the GO terms to be evaluated. It can be either "BP", "MF" or "CC". Defaults to "BP".
#' @param measure character. The method used to calculate the similarity. It can be either "Resnik", "Lin", "Rel", "Jiang", "TCSS" or "Wang". Defaults to "Wang".
#' @param combine character. The method used to combine similarity scores of multiple GOs. It can be either "max", "avg", "rcmax" or "BMA". It defaults to "BMA".
#' @param semSim matrix. Similarity matrix calculated using the getSemSim function. If it has not been calculated it should be set to NULL. Defaults to NULL.
#' @param cores numeric. The number of cores to be used. The greater the number of cores the larger the amount of RAM that will be required. It should be a number greater than 0. Defaults to 1.
#'
#' @return A matrix object with the similarities between the OGs being analysed.
#'
#' @examples
#'
#' OGs <- c("OG0009248", "OG0009250", "OG0009251", "OG0009252", "OG0009253", "OG0009254", "OG0009258", "OG0009261", "OG0009267", "OG0009268", "OG0009269", "OG0009271")
#' simMat <- getSimMat(GOsfromOGs = sampleGOsfromOGs, selectedOGs = OGs)
#'
#' # Using the output of getSemSim
#' simMat <- getSimMat(GOsfromOGs = sampleGOsfromOGs, selectedOGs = OGs, semSim = sampleSemSim)
#'
#' @export


getSimMat <- function(GOsfromOGs, selectedOGs, ont = "BP", measure = "Wang", combine = "BMA", semSim = NULL, cores = 1){
  if(!is.character(combine))
    stop("combine must be an object of the class character.")
  if(!combine %in% c("max", "avg", "rcmax", "BMA") | length(measure) > 1)
    stop("combine can only take one of the following values: max, avg, rcmax or BMA.")
  if(is.null(semSim)){
    semSim <- getSemSim(GOsfromOGs, selectedOGs, ont, measure, cores)
    pangenome <- GOsfromOGs
    annotations_ogs <- pangenome[pangenome$OGs %in% selectedOGs, ]
    annotations_ogs <- split(annotations_ogs$GOs, annotations_ogs$OGs)
  } else{
    if(!is.character(selectedOGs))
      stop("selectedOGs must be an object of the character class")
    selectedOGs <- unique(selectedOGs)
    if(!is.character(ont))
      stop("ont must be an object of the class character.")
    if(!ont %in% c("BP", "CC", "MF") | length(ont) > 1)
      stop("ont can only take one of the following values: BP, CC or MF.")
    if(!is.numeric(cores) | cores < 1 | length(cores) > 1)
      stop("cores should be a single number greater than 0.")
    if(cores > 10)
      warning("A value greater than 10 was selected for cores. Make sure that enough RAM is available.")
    message("Processing GOsfromOGs:")
    pangenome <- GOsfromOGs
    if(sum(colnames(pangenome) == c("GOs", "OGs")) != 2 | !is.data.frame(pangenome))
      stop(paste(GOsfromOGs, " should be a data.frame with two columns: GOs and OGs.", sep = ""))
    if(sum(pangenome$GOs %in% c("", " ")) > 0 | sum(pangenome$OGs %in% c("", " ")) > 0)
      stop(paste("The ", GOsfromOGs, " data.frame cannot contain empty characters or whitespaces.", sep = ""))
    if(sum(pangenome$OGs %in% selectedOGs) == 0)
      stop("None of the OGs of interest are included in the full GOs to OGs dataset")
    if(sum(!selectedOGs %in% pangenome$OGs) > 0)
      warning(paste("The following OGs are not found in the GOsfromOGs object:", paste(selectedOGs[!selectedOGs %in% pangenome$OGs], collapse = ", ")))
    annotations_ogs <- pangenome[pangenome$OGs %in% selectedOGs, ]
    annotations_ogs <- split(annotations_ogs$GOs, annotations_ogs$OGs)
    if(!is(semSim, "dsyMatrix"))
      stop("semSim does not contain a matrix, please generate an appropriate matrix using the getSemSim function.")
    if(!Matrix::isSymmetric(semSim))
      stop("semSim does not contain a symmetrical matrix, please generate an appropriate matrix using the getSemSim function.")
    if(sum(rownames(semSim) != colnames(semSim)))
      stop("row and column names for semSim do not match, please generate a appropriate using the getSemSim function.")
  }
  message("Obtaining similarity matrix")
  simMat <- pbmcapply::pbmclapply(1:(length(annotations_ogs) - 1), function(x){
    out <- sapply((x + 1):length(annotations_ogs), function(y){
      GOSemSim::combineScores(semSim[colnames(semSim) %in% annotations_ogs[[x]], rownames(semSim) %in% annotations_ogs[[y]]], combine = combine)
    })
    out <- c(rep(NA, x), out)
  }, mc.cores = cores)
  cat("\n")
  simMat <- do.call("rbind", simMat)
  simMat <- rbind(simMat, rep(NA, ncol(simMat)))
  simMat <- t(simMat)
  rownames(simMat) <- colnames(simMat) <- names(annotations_ogs)
  return(simMat)
}
