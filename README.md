# constellatoR <a><img src='https://github.com/MetazoaPhylogenomicsLab/constellatoR/blob/main/inst/img/constellatoR.png' align="right" height="190" /></a>
R package constellatoR to assess functional convergence in Orthologous Groups (OGs) using the similarity scores of their Gene Ontology annotations (GOs). 

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

-----------------------------------------------------------------------   

# Installation

```{r}

# get the development version from GitHub using devtools :
# install.packages("devtools")
devtools::install_github("MetazoaPhylogenomicsLab/constellatoR",build_vignettes = TRUE)

```

# Usage

#### Foreword : 

The main objective of this package is to cluster Orthologous Groups (OGs) using the similarity scores of their Gene Ontology annotations (GOs). Other features might be used, such as genes along with their associated GOs, however keep in consideration that runtimes will increase the larger the number of features that are used. This tutorial describes the full process going from the outputs of OrthoFinder and [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) but you can choose to use different tools to provide the data given that they are in the appropiate format. Using this package you can obtain constellation plots such as this:

<a><img src='https://github.com/MetazoaPhylogenomicsLab/constellatoR/blob/main/inst/img/Lepto_apclusters_q0.png' align="center" /></a>


### 1 - Prepare the data :

To obtain the OGs we recommend using [OrthoFinder](https://github.com/davidemms/OrthoFinder), however any other tool can be used as long as the output format is the same. For the functional annotation we recommend using eggNOG-mapper or [goPredSim](https://github.com/Rostlab/goPredSim). To minimize the runtime we provide a very reduced dataset.

```{r setup}
library(constellatoR)
```

OGs corresponds to the path to Orthogroups.tsv, one of the output files of OrthoFinder. The first column is expected to be Orthogroup and the rest the species analysed.

```{r}
OGs <- system.file("extdata", "Orthogroups.tsv", package="constellatoR")
rawOGs <- read.delim(OGs)
```

GOsdir corresponds to path to the directory that contains the output of the annotation tool. The file names should be the name of the species followed by ".emapper.annotations" if the method used was eggNOG-mapper (EM) and ".goPredSim.annotations" if method used was goPredSim (GPS). For this tutorial we provide the output of EM.

```{r}
GOsdir <- system.file("extdata", package="constellatoR")
```

Next, we indicate the species that we want to analyse, in this case SMED, PMAX and LLON. Keep in mind that they should match the column names in the output of OrthoFinder and also the file names of the output of the functional annotation software. Finally, we run the `GOsfromOGs` function. We don't indicate the `method` given that we example files are the output of EM and `method` defaults to EM. We don't indicate the number of cores in the `cores` parameter either, given that we will only use one and it defaults to that.  

```{r}
species <- c("SMED", "PMAX", "LLON")
GOsfromOGs <- prepareGOsfromOGs(OGs, GOsdir, species = species)
```

### 2 - Obtain the semantic similarity matrix between all the GOs :

Now we will calculate the semantic similarity between all the GOs in a list using the `mgoSim` function from the GOSemSim package. The obtained matrix can then be used to obtain the similarity matrix for different sets of OGs, thus to reduce computation times we suggest that you use this function to obtain the semantic similarity matrix for the full set of GOs to be evaluated and then use it to evaluate different sets of OGs. Also, we provide the matrices for the full set of GOs for each ontology (BP, MF or CC) from [Figshare](https://figshare.com/articles/dataset/constellatoR_Semantic_similarity_matrices_for_all_GO_terms/22220611).

Using the output of the `prepareGOsfromOGs` we can run the `getSemSim` function. If you did not run the `prepareGOsfromOGs` function you can provide a data.frame with two columns. The header of the first should be `GOs` and the header of the second should be `OGs`. The first column should be a list of GOs and the second one a list of the OGs to which each GO belongs. No empty fields are allowed.

We also have to specify the OGs we are interested in in the `selectedOGs` parameter and the ontology (BP, MF or CC) that we want to evaluate in the `ont` parameter. Given that we are interested in the BP ontology, to which it defaults we do not specify it. Also the method used to calculate the similarity can be specified using the `measure` parameter, check the help of the function to see all available options. In this case we select the Wang method to which the function defaults to. Finally, no cores are specified given that we will only use one.

```{r}
OGs <- c("OG0009248", "OG0009250", "OG0009251", "OG0009252", "OG0009253", "OG0009254", "OG0009258", "OG0009261", "OG0009267", "OG0009268", "OG0009269", "OG0009271")
semSim <- getSemSim(GOsfromOGs = GOsfromOGs, selectedOGs = OGs)
```

### 3 - Calculate the similarity between OGs :

Now we can calculate the similarity between OGs by combining the semantic similarity scores between the GOs from each pair of OGs being analyzed. In this example we obtained the similarity matrix in the previous step, but if needed the `getSimMat` function can internally use the `getSemSim` function to calculate the semantic similarity between GOs and then calculate the similarity between OGs. However, if you are going to calculate the similarities between different sets of OGs we recommend using the getSemSim function to obtain the full matrix, store it in an object and then use it to calculate the similarity scores between different sets of OGs to reduce times as we will do in this example.

We will use the same GOsfromOGs object, we will not specify `ont`, `measure` object given that we are providing the matrix previously created with the `getSimMat` function. We can also select a method to combine the similarity scores from all the GOs in a single OG. In this case we won't specify the `combine` parameter either given that we will use the default method BMA. The matrix created with the `getSimMat` function will be specified using the `semSim` parameter and we will use a single core so we won't use the `cores` parameter. Keep in mind that the greater the number of cores, the larger the amount of RAM that will be required.

```{r}
simMat <- getSimMat(GOsfromOGs = GOsfromOGs, selectedOGs = OGs, semSim = semSim)
```

### 4 - Cluster the OGs :

We can now cluster the OGs using the semantic similarity of their GOs. This can be done using the `clusterOGs` function which uses the `apcluster` function from the apcluster package and requires the output of the getSimMat function. We also need to specify the maximum number of iterations that will  be executed using the `clusterits` parameter. We will not specify it given that we will use its default value of 20000.

```{r}
clusters <- clusterOGs(simMat)
```

### 5 - Annotate the OGs :

Next we can annotate the OGs using the `annotateOGs` function. It identifies the most relevant GO terms associated to each OG and obtains the three GO terms with the lowest level belonging to each GO and retrieves their description.

The input requires the same GOsfromOGs object and the list of selectedOGs. Additionally, the ontology should be specified using the `ont` parameter and it should match the ontology used in the previous steps. Given that it defaults to BP again we will not specify it. We can also specify which levels of the GOs will be ignored during the annotation. The descriptions for the GOs with lower levels are more general so we will ignore them. Given that the `ignoreLevels` parameter defaults to 0:3 we will not specify it. Finally, we will use the default value of 1 for the cores parameter.

```{r}
annotation <- annotateOGs(GOsfromOGs, OGs)
```

If needed you can also obtain a word cloud with the most relevant terms for the OGs of interest using the `getWordcloud` function. This function will obtain the word clouds for all the exemplar OGs. This requires the GOsfromOGs object, an ontology (as before, it should match the ones used previously) and the clusters obtained with the `clusterOGs` function. You can also select a file with terms that should be excluded from the annotation providing the path via the `toExclude` parameter. In this case we will use a file (bio-stopwords.txt) which was downloaded from (Genes2WordCloud from the Ma'ayan Laboratory)[https://www.maayanlab.net/G2W/help.php] and modified to include several additional terms. You can provide any list as long as the file contains one term per line. 

```{r}
stopwords <- system.file("extdata/bio-stopwords.txt", package="constellatoR")
wcl <- getWordcloud(sampleGOsfromOGs, clusters = sampleClusters, toExclude = stopwords)
```

### 6 - Plot the constellations :

Finally, with all of the previous information we can plot our constellations. The `plotConstellation` function is used to plot the OGs that have been previously clustered by their semantic similarity. A polygon comprising all of the OGs in a cluster is plotted and the functions of the representative member of that cluster are displayed. You can also add extra information such as the number of genes per cluster or the species that have members of an OG using the size and color of the dots (or stars). While we provide an automatic function to generate this plot, if changes are needed you can adapt the plot function to better suit your needs. 

This function requires the clusters that were generated with the `clusterOGs` function and an annotation for the OGs. In this case we will use the information obtained with the `annotateOGs` and `getWordcloud` functions, but the user can provide any desired annotation as long as the number of rows in the data.frame matches the number of rows. The rownames should be the OGs and it must contain one or more columns with the annotation with their headers. The `term` parameter is used to specify the column that should be used to annotate the OGs, in this case we will use the column named BP1. Additionally, if one of the columns is named `wordcloud`, the function will identify it and add it to the plot. Finally, you can use the `color` and `size` parameters to specify which column should be used for the sizes and colors of the stars in the constellations.

We will join the annotation and wcl objects in a single data.frame called annot. Then we will calculate the species with members in each OG and the number of genes in each OG and append these values as extra columns in the annot object. Finally we will call the `plotConstellation` function.

```{r}
annot <- merge(annotation, wcl, by.x = 0, by.y = 0)
row.names(annot) <- annot[, 1]
annot <- annot[, -1]
groups <- sapply(row.names(annot), function(x){
  paste(colnames(rawOGs)[rawOGs[which(rawOGs$Orthogroup == x), ] != ""][-1], collapse = ", ")
})
size <- sapply(row.names(annot), function(x){
  sum(unlist(strsplit(as.character(rawOGs[which(rawOGs$Orthogroup == x), -1]), ", ")) != "")
})
annot <- data.frame(annot, groups = groups, size = size)
p <- plotConstellation(clusters, annot, term = "BP1", color = "groups", size = "size")
```

We have now obtained our constellation plot. Of course this is a very reduced data set so we only see three clusters, two comprising a single member (with putative functions 'ion transport' and 'organelle organization' respectively) and one with the remaining OGs with putative 'nucleobase-containing compound metabolic process' function. While this is a very simple data set we can see that the OGs in the larger cluster contain some OGs found only in PMAX and LLON, one only in SMED and PMAX and others in all three species. Thus, all three species have genes in OGs potentially involved in nucleobase-containing compound metabolic process even if all of the OGs don't contain genes shared by all three species.

<a><img src='https://github.com/MetazoaPhylogenomicsLab/constellatoR/blob/main/inst/img/Rplot.png' align="center" /></a>   

# Getting help

Need help, identified a bug, or want to see other features implemented? Open an issue here or send an email to: carlos.vargas@ibe.upf-csic.es.

# Citation
