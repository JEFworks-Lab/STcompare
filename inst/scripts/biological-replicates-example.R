library(STcompare)
library(SpatialExperiment)
library(SEraster)
library(viridis)
library(MERINGUE)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(BiocGenerics)
library(BiocParallel)
devtools::load_all()

s2r2_url <- "https://zenodo.org/records/10724029/files/STalign_S2R2.csv.gz?download=1"
s2r3_url <- "https://zenodo.org/records/10724029/files/STalign_S2R3_to_S2R2.csv.gz?download=1"

# read the data from the url 
read_csv_gz_url <- function(x, ...) {
  u  <- url(x, open = "rb")
  gz <- gzcon(u, text = TRUE)
  on.exit({
    try(close(gz), silent = TRUE)
    try(close(u),  silent = TRUE)
  }, add = TRUE)

  read.csv(gz, ...)
}

target <- read_csv_gz_url(s2r2_url)
source <- read_csv_gz_url(s2r3_url)

# TODO: what are the column names 
head(target)
head(source)


# convert target to spatial experiments

# 1) target ####################################################################
dim(target)
target[1:5, 1:10]

# extract the cell positions (x, y)
pos_target <- target[, c('x', 'y')]
rownames(pos_target) <- target$X
head(pos_target)

# after the thi4thrd column, the remaining columns are the genes
# making this the counts matrix for spatial expression 
gene_target <- target[, 4:dim(target)[2]]
rownames(gene_target) <- target$X
gene_target[1:5, 1:5]

# for spatial experiments, convert from a dataframe to a dgCMatrix
class(gene_target)
target_sparse <- as(t(gene_target), "dgCMatrix")

# format into SpatialExperiment
spe_target <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = target_sparse),
  spatialCoords = as.matrix(pos_target)
)

# 2) source  ###################################################################

dim(source)
source[1:5, 1:10]

# extract the cell positions (x, y) after alignment 
pos_source <- source[, c('STalign_x', 'STalign_y')]
rownames(pos_source) <- source$X
colnames(pos_source) <- c("x", "y")
head(pos_source)

# the genes in source start at the 8th column 
gene_source <- source[, 8:dim(source)[2]]
rownames(gene_source) <- source$X
gene_source[1:5, 1:5]

# need to convert from data.frame to dgCMatrix
class(gene_source)
source_sparse <- as(t(gene_source), "dgCMatrix")

spe_source <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as(t(gene_source), "dgCMatrix")),
  spatialCoords = as.matrix(pos_source)
)

identical(rownames(spe_target), rownames(spe_source))

spe_list <- list(target = spe_target, source = spe_source)

res = 200 
output <- SEraster::rasterizeGeneExpression(
  spe_list, 
  resolution = res, 
  BPPARAM = BiocParallel::SerialParam()
)

moransI <- function(SE_temp, filterDist = resolution){
  
  # Get neighbor-relationships
  w <- MERINGUE::getSpatialNeighbors(spatialCoords(SE_temp), filterDist = filterDist)
  par(mfrow=c(1,1))
  plotNetwork(spatialCoords(SE_temp), w)
  
  # Identify significantly spatially auto-correlated genes
  start_time <- Sys.time()
  I <- MERINGUE::getSpatialPatterns(as.matrix(assay(SE_temp)), w)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  results.filter <- MERINGUE::filterSpatialPatterns(mat = assay(SE_temp),
                                          I = I,
                                          w = w,
                                          adjustPv = TRUE,
                                          alpha = 0.05,
                                          minPercentCells = 0.01,
                                          verbose = TRUE, details = TRUE)
  
  return(results.filter)
}

start_time <- Sys.time()
moransI_merfish_target <- moransI(output$target, filterDist = res)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 46.42798 secs

start_time <- Sys.time()
moransI_merfish_source <- moransI(output$source, filterDist = res)
end_time <- Sys.time()
print(end_time - start_time) #Time difference of 2.000802 hours

svg_merfish_target <- moransI_merfish_target %>% filter(minPercentCells > 0.01) %>% rownames()
svg_merfish_source <- moransI_merfish_source %>% filter(minPercentCells > 0.01) %>% rownames()

svg_int <- intersect(svg_merfish_target, svg_merfish_source)
svg_uni <- union(svg_merfish_target, svg_merfish_source)
non_svg <- rownames(output$target)[!(rownames(output$target) %in% svg_int)]

length(svg_int)
length(non_svg)


set.seed(0)
start_time <- Sys.time()
merfishCorrelation <- spatialCorrelationGeneExpIterPermutations(
  output, 
  nThreads = 22 
)
end_time <- Sys.time()
print(end_time - start_time)

