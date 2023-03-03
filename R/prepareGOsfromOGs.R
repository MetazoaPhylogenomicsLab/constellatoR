#' prepareGOsfromOGs
#'
#' @description Function to prepare a data frame linking the Orthologous Groups (OGs) and the Gene Ontology annotations (GOs) associated to them. It requires the results of OrthoFinder to define the OGs and a directory with the functional annotation results from either eggNOG-mapper or goPredSim for each species.
#'
#' @param OGs character. The path to the Orthogroups.tsv output of OrthoFinder found in the Orthogroups folder. The first column is expected to be Orthogroup and the rest the species analysed.
#' @param GOsdir character. The path to the directory that contains the output of the annotation tool. The file names should be the name of the species followed by ".emapper.annotations" if the method used was eggNOG-mapper (EM) and ".goPredSim.annotations" if method used was goPredSim (GPS).
#' @param method character. The methodology used to functionally annotate the genes for each species. EM stands for eggNOG-mapper and GPS for goPredSim. Defaults to EM.
#' @param species character. The species to be assessed. They should match the column names in the OrthoFinder output and the file names produced by the functional annotation software.
#' @param cores numeric. The number of cores to be used. Defaults to 1.
#'
#' @return A data frame object with two columns, the first with the GOs and the second with the corresponding OGs that contain them.
#'
#' @examples
#'
#' OGs <- system.file("extdata", "Orthogroups.tsv", package="constellatoR")
#' GOsdir <- system.file("extdata", package="constellatoR")
#' species <- c("SMED", "PMAX", "LLON")
#'
#' GOsfromOGs <- prepareGOsfromOGs(OGs, GOsdir, species = species)
#'
#' @export


prepareGOsfromOGs <- function(OGs, GOsdir, method = "EM", species, cores = 1){
  if(!is.character(OGs))
    stop("OGs must be an object of the class character.")
  if(!is.character(GOsdir))
    stop("GOsdir must be an object of the class character.")
  if(!file.exists(OGs))
    stop(paste("The file ", OGs, " does not exist.", sep = ""))
  if(!dir.exists(GOsdir))
    stop(paste("The directory ", GOsdir, " does not exist.", sep = ""))
  if(!is.character(species))
    stop("species must be an object of the class character containing the species to be analysed.")
  if(!method %in% c("GPS", "EM") | length(method) != 1)
    stop("method must be either EM for the output of eggNOG-mapper or GPS for the output of goPredSim")
  if(!is.numeric(cores) | cores < 1 | length(cores) > 1)
    stop("cores should be a single number greater than 0.")
  OGs <- data.table::fread(OGs)
  if(colnames(OGs)[1] != "Orthogroup")
    stop(paste("The file ", OGs, " requires the first column name to be Orthogroup.", sep = ""))
  if(sum(!species %in% colnames(OGs) > 0))
    stop(paste("The following species are not found as column names in the output of OrthoFinder: ", paste(species[!species %in% colnames(OGs)], collapse = ", "), sep = ""))
  if(method == "EM"){
    allAnnotFiles <- list.files(GOsdir, pattern = ".emapper.annotations")
    found <- paste(species, ".emapper.annotations", sep = "") %in% allAnnotFiles
    if(sum(found) != length(species))
      stop(paste("The following species are not found in the annotation files: ", paste(species[!found], collapse = ", "), sep = ""))
    allAnnotFiles <- lapply(paste(species, ".emapper.annotations", sep = ""), function(x){
      cat(paste("Reading file ", x, "\n", sep = ""))
      annot <- data.table::fread(cmd = paste("grep -v '^#' ", GOsdir, "/", x, sep = ""), select = c(1, 10))
      separatedGOs <- lapply(1:nrow(annot), function(y){
        allGOs <- strsplit(annot[[2]][y], ",")
        out <- data.frame(annot[[1]][y], allGOs)
        colnames(out) <- c("genes", "GOs")
        return(out)
      })
      separatedGOs <- do.call("rbind", separatedGOs)
    })
    names(allAnnotFiles) <- species
  } else if(method == "GPS"){
    allAnnotFiles <- list.files(GOsdir, pattern = ".goPredSim.annotations")
    found <- paste(species, ".goPredSim.annotations", sep = "") %in% allAnnotFiles
    if(sum(found) != length(species))
      stop(paste("The following species are not found in the annotation files: ", paste(species[!found], collapse = ", "), sep = ""))
    allAnnotFiles <- lapply(paste(species, ".goPredSim.annotations", sep = ""), function(x){
      cat(paste("Reading file ", x, "\n", sep = ""))
      annot <- data.table::fread(cmd = paste("grep -v '^#' ", GOsdir, "/", x, sep = ""), select = c(1, 2))
      colnames(annot) <- c("genes", "GOs")
      return(annot)
    })
    names(allAnnotFiles) <- species
  }
  posSpecies <- sapply(c("Orthogroup", species), function(x){
    which(colnames(OGs) == x)
  })
  OGs <- OGs[, posSpecies]
  OGs <- OGs[rowSums(OGs != "") > 1, ]
  GOsinOGs <- pbmcapply::pbmclapply(1:nrow(OGs), function(x){
    allGOs <- lapply(1:length(species), function(y){
      adjCol <- y + 1
      genes <- unlist(strsplit(as.character(OGs[, adjCol][x]), ", "))
      GOs <- allAnnotFiles[[y]][which(allAnnotFiles[[y]][, 1] %in% genes), 2]
    })
    allGOs <- unique(unlist(allGOs))
    if(length(allGOs) == 0)
      allGOs <- ""
    out <- data.frame(allGOs, OGs[x, 1])
  }, mc.cores = cores)
  GOsinOGs <- do.call("rbind", GOsinOGs)
  GOsinOGs <- GOsinOGs[!GOsinOGs[, 1] %in% c("-", ""), ]
  GOsinOGs <- GOsinOGs[!duplicated(paste(GOsinOGs[, 1], GOsinOGs[, 2])), ]
  colnames(GOsinOGs) <- c("GOs", "OGs")
  return(GOsinOGs)
}
