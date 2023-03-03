#' annotateOGs
#'
#' @description Function used to annotate OGs by identifying the most relevant GO terms associated to each OG. It obtains the three GO terms with the lowest level belonging to each GO and retrieves their description.
#'
#' @param GOsfromOGs data.frame. Contains two columns, the header of the first should be GOs and the header of the second should be OGs. The first column should be a list of GOs and the second one a list of the OGs to which each GO belongs. No empty fields are allowed.
#' @param selectedOGs character. The OGs that will be evaluated. They should be included in the OGs column in the GOsfromOGs data.frame, but elements not found will be ignored and a warning will list them.
#' @param ont character. The ontology of the GO terms to be evaluated. It must include one or more of the following terms: "BP", "MF" or "CC" (such as c("BP", "MF", "CC"). Defaults to "BP".
#' @param ignoreLevels numeric. Levels of the GOs that will be ignored during the annotation. The descriptions for the GOs with lower levels are more general. If no filtering is desired set to NULL. Defaults to 0:3.
#' @param cores numeric. The number of cores to be used. The greater the number of cores the larger the amount of RAM that will be required. It should be a number greater than 0. Defaults to 1.
#'
#' @return A matrix with the annotation for all OGs.
#'
#'
#' @examples
#'
#' OGs <- c("OG0009248", "OG0009250", "OG0009251", "OG0009252", "OG0009253", "OG0009254", "OG0009258", "OG0009261", "OG0009267", "OG0009268", "OG0009269", "OG0009271")
#'
#' annotation <- annotateOGs(sampleGOsfromOGs, OGs)
#'
#' @export


annotateOGs <- function(GOsfromOGs, selectedOGs, ont = "BP", ignoreLevels = 0:3, cores = 1){
  if(!is.character(selectedOGs))
    stop("selectedOGs must be an object of the character class")
  if(sum(ont %in% c("BP", "CC", "MF") == FALSE) > 0 | length(ont) > 3)
    stop("ont can only take one or more of the following values: BP, CC or MF")
  if(!is.numeric(cores) | cores < 1 | length(cores) > 1)
    stop("cores should be a single number greater than 0.")
  if(!is.null(ignoreLevels) & (!is.numeric(ignoreLevels) | sum(ignoreLevels < 0) > 0))
    stop("ignoreLevels should be either a string of numbers equal or greater than 0 or NULL.")
  if(cores > 10)
    warning("A value greater than 10 was selected for cores. Make sure that enough RAM is available.")
  message("Processing GOsfromOGs:")
  pangenome <- GOsfromOGs
  if(sum(colnames(pangenome) == c("GOs", "OGs")) != 2)
    stop(paste("The column names in the ", GOsfromOGs, " data.frame should be GOs and OGs.", sep = ""))
  if(sum(pangenome$GOs %in% c("", " ")) > 0 | sum(pangenome$OGs %in% c("", " ")) > 0)
    stop(paste("The ", GOsfromOGs, " data.frame cannot contain empty characters or whitespaces.", sep = ""))
  if(sum(pangenome$OGs %in% selectedOGs) == 0)
    stop("None of the OGs of interest are included in the full GOs to OGs dataset")
  if(sum(!selectedOGs %in% pangenome$OGs) > 0)
    warning(paste("The following OGs are not found in the GOsfromOGs object:", paste(selectedOGs[!selectedOGs %in% pangenome$OGs], collapse = ", ")))
  annotations_ogs <- pangenome[pangenome$OGs %in% selectedOGs, ]
  annotations_ogs <- split(annotations_ogs$GOs, annotations_ogs$OGs)
  pAll <- pbmcapply::pbmclapply(annotations_ogs, function(x){
    if("BP" %in% ont){
      out <- internalAnnotate(x, ont = "BP", ignoreLevels)
      cNames <- c("BP1", "BP2", "BP3")
    }
    if("MF" %in% ont){
      if("BP" %in% ont){
        out <- c(out, internalAnnotate(x, ont = "MF", ignoreLevels))
        cNames <- c(cNames, "MF1", "MF2", "MF3")
      } else{
        out <- internalAnnotate(x, ont = "MF", ignoreLevels)
        cNames <- c("MF1", "MF2", "MF3")
      }
    }
    if("CC" %in% ont){
      if("BP" %in% ont | "MF" %in% ont){
        out <- c(out, internalAnnotate(x, ont = "CC", ignoreLevels))
        cNames <- c(cNames, "CC1", "CC2", "CC3")
      } else{
        out <- internalAnnotate(x, ont = "CC", ignoreLevels)
        cNames <- c("CC1", "CC2", "CC3")
      }
    }
    names(out) <- cNames
    return(out)
  }, mc.cores = cores)
  cat("\n")
  pAll <- do.call("rbind", pAll)
  out <- lapply(selectedOGs, function(x){
    pos <- which(rownames(pAll) == x)
    if(length(pos) != 0){
      out <- pAll[pos, ]
    } else{
      out <- rep(NA, ncol(pAll))
    }
  })
  out <- do.call("rbind", out)
  rownames(out) <- selectedOGs
  return(out)
}

#################################################################################################
# internalAnnotate: function used to annotate OGs by identifying the most relevant GO terms from a group of GOs. It obtains the top three prioritized GO terms by identifying the GOs on the highest level and thus more specific.
#
# GOs: the list of GOs to be evaluated
# ont: the ontology of the GO terms to be evaluated. It can be either "BP", "MF" or "CC".
# ignoreLevels: an integer string with the levels of the GOs that will be ignored during the annotation. The descriptions for the GOs with lower levels are more general. If no filtering is desired set to NULL. It defaults to 0:3.
#################################################################################################

internalAnnotate <- function(GOs, ont, ignoreLevels){
  if(ont == "BP")
    level <- suppressWarnings(tryCatch(GOxploreR::GOTermBPOnLevel(GOs), error = function(GOs){}))
  if(ont == "MF")
    level <- suppressWarnings(tryCatch(GOxploreR::GOTermMFOnLevel(GOs), error = function(GOs){}))
  if(ont == "CC")
    level <- suppressWarnings(tryCatch(GOxploreR::GOTermCCOnLevel(GOs), error = function(GOs){}))
  if(!is.null(ignoreLevels)){
    level <- level[!level$Level %in% ignoreLevels, ]
  }
  if(is.null(nrow(level))){
    terms <- c(NA, NA, NA)
  } else{
    if(nrow(level) == 0){
      terms <- c(NA, NA, NA)
    } else{
      terms <- clusterProfiler::go2term(level$Term[order(level$Level)][1:3])$Term
      terms <- terms[!is.na(terms)]
      if(length(terms) < 3)
        terms <- c(terms, rep(NA, (3 - length(terms))))
    }
  }
  return(terms)
}
