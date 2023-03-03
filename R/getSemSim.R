#' getSemSim
#'
#' @description Function to calculate the semantic similarity between all the Gene Ontology annotations (GOs) in a list using the mgoSim function from the GOSemSim package. The obtained matrix can then be used to obtain the similarity matrix for different sets of Orthologous Groups (OGs), thus to reduce computation times you can use this function to obtain the semantic similarity matrix for the full set of GOs to be evaluated and then use it to evaluate different sets of OGs.
#'
#' @param GOsfromOGs data.frame. Contains two columns, the header of the first should be GOs and the header of the second should be OGs. The first column should be a list of GOs and the second one a list of the OGs to which each GO belongs. No empty fields are allowed.
#' @param selectedOGs character. The OGs that will be evaluated. They should be included in the OGs column in the GOsfromOGs data.frame, but elements not found will be ignored and a warning will list them.
#' @param ont character. The ontology of the GO terms to be evaluated. It can be either "BP", "MF" or "CC". Defaults to "BP".
#' @param measure character. The method used to calculate the similarity. It can be either "Resnik", "Lin", "Rel", "Jiang", "TCSS" or "Wang". Defaults to "Wang".
#' @param cores numeric. The number of cores to be used. The greater the number of cores the larger the amount of RAM that will be required. It should be a number greater than 0. Defaults to 1.
#'
#' @return A matrix object with the semantic similarity between all GOs being analyzed.
#'
#' @examples
#'
#' OGs <- c("OG0009248", "OG0009250", "OG0009251", "OG0009252", "OG0009253", "OG0009254", "OG0009258", "OG0009261", "OG0009267", "OG0009268", "OG0009269", "OG0009271")
#'
#' semSim <- getSemSim(GOsfromOGs = sampleGOsfromOGs, selectedOGs = OGs)
#'
#' @export


getSemSim <- function(GOsfromOGs, selectedOGs, ont = "BP", measure = "Wang", cores = 1){
  if(!is.character(selectedOGs))
    stop("selectedOGs must be an object of the class character.")
  selectedOGs <- unique(selectedOGs)
  if(!is.character(ont))
    stop("ont must be an object of the class character.")
  if(!ont %in% c("BP", "CC", "MF") | length(ont) > 1)
    stop("ont can only take one of the following values: BP, CC or MF.")
  if(!is.character(measure))
    stop("measure must be an object of the class character.")
  if(!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang") | length(measure) > 1)
    stop("measure can only take one of the following values: Resnik, Lin, Rel, Jiang, TCSS or Wang.")
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
  go_terms <- unique(unlist(annotations_ogs))
  if(length(go_terms) == 1)
    stop("The dataset comprises only one GO term, try using a larger set of OGs")
  if(length(go_terms) < 5)
    warning("The dataset comprises less than five GOs, try using a larger set of OGs")
  message("Calculating semantic similarities")
  matBP <- GOSemSim::godata(ont = ont)
  semSim <- pbmcapply::pbmclapply(1:length(go_terms), function(x){
    out <- sapply(x:length(go_terms), function(y){
      out <- GOSemSim::mgoSim(as.character(go_terms[x]), as.character(go_terms[y]), semData = matBP, measure = measure, combine = NULL)
    })
    out <- c(rep(NA, (x - 1)), out)
  }, mc.cores = cores)
  cat("\n")
  semSim <- do.call("rbind", semSim)
  semSim <- Matrix::forceSymmetric(semSim, uplo = "U")
  rownames(semSim) <- colnames(semSim) <- go_terms
  return(semSim)
}
