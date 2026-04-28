
# From the merfish merfish tutorial 

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


# 1) target ####################################################################
dim(target)

# TODO: 
target[1:5, 1:10]

# extract the cell positions (x, y)
pos_target <- target[, c('x', 'y')]

# TODO: 
rownames(pos_target) <- target$X

# TODO:
# after the third column, the remaining columns are the genes
# making this the counts matrix for spatial expression 
gene_target <- target[, 4:dim(target)[2]]
# gene_target <- target[, 5:dim(target)[2]]
rownames(gene_target) <- target$X

# remove blanks 
gene_target <- gene_target[, !grepl("^Blank\\.", colnames(gene_target))]


# for spatial experiments, convert from a dataframe to a dgCMatrix
class(gene_target)
target_sparse <- as(t(gene_target), "dgCMatrix")
class(target_sparse)

# format into SpatialExperiment

# TODO: cpm normalize? 
spe_target <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = target_sparse),
  spatialCoords = as.matrix(pos_target)
)

# 2) source  ###################################################################

dim(source)

source[1:5, 1:10]

# TODO:
pos_source <- source[, c('STalign_x', 'STalign_y')]
rownames(pos_source) <- source$X
colnames(pos_source) <- c("x", "y")
gene_source <- source[, 8:dim(source)[2]]
# gene_source <- source[, 9:dim(source)[2]]

rownames(gene_source) <- source$X

# need to convert from data.frame to dgCMatrix
class(gene_source)

# remove blanks 
gene_source <- gene_source[, !grepl("^Blank\\.", colnames(gene_source))]

source_sparse <- as(t(gene_source), "dgCMatrix")

spe_source <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = as(t(gene_source), "dgCMatrix")),
  spatialCoords = as.matrix(pos_source)
)

spe_list <- list(target = spe_target, source = spe_source)

save(spe_list, file = "test/time_benchmarks/inputs/spe_list.Rdata")
