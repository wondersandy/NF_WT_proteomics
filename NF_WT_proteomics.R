## Import NF-L KO (Beatrix) dataset
library(BiocManager)
remove.packages(c("bayestestR", "beadarray", "boot", "broom", "class", "clubSandwich", "dbplyr", "emmeans", "etm", "fda", "ff", "future",
                "ggeffects", "ggfortify", "graphlayouts", "HDF5Array", "imp4p", "insight", "KernSmooth", "MASS", "mirt", "MuMIn", "nlme", "nnet", 
                "pkgbuild", "PKI", "plotmo", "plotrix", "pmml", "polylabelr", "ProjectTemplate", "proxy", "purrr", "qdap", "qdapTools", "raster", 
                "RcppArmadillo", "RCurl", "rex", "sf", "shinyFiles", "spatial", "stopwords", "svUnit", "systemfonts", "tibble", "timetk",
                "tinytex", "TSP", "withr", "WrightMap", "xml2", "xslt"))
remove.packages(c("xml2", "WrightMap", "withr", "timetk", "tibble", "svUnit", "spatial", "rex", "qdap", "proxy", "pmml", "pkgbuild", "nnet", "nlme", 
                  "mirt", "MASS", "KernSmooth", "insight", "graphlayouts", "ggfortify", "ggeffects", "fda", "class", "broom", "boot", "bayestestR"), 
                  lib = "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
install.packages(c("bayestestR", "beadarray", "boot", "broom", "class", "clubSandwich", "dbplyr", "emmeans", "etm", "fda", "ff", "future",
                   "ggeffects", "ggfortify", "graphlayouts", "HDF5Array", "imp4p", "insight", "KernSmooth", "MASS", "mirt", "MuMIn", "nlme", "nnet", 
                   "pkgbuild", "PKI", "plotmo", "plotrix", "pmml", "polylabelr", "ProjectTemplate", "proxy", "purrr", "qdap", "qdapTools", "raster", 
                   "RcppArmadillo", "RCurl", "rex", "sf", "shinyFiles", "spatial", "stopwords", "svUnit", "systemfonts", "tibble", "timetk",
                   "tinytex", "TSP", "withr", "WrightMap", "xml2", "xslt"))
                # Installing packages into ‘/Users/sandip_home/Library/R/3.6/library’
                #(as ‘lib’ is unspecified)
                #Warning in install.packages :
                #  packages ‘beadarray’, ‘HDF5Array’ are not available (for R version 3.6.1)
install(c("beadarray", "HDF5Array"))
install()
remove.packages(c("beadarray", "clubSandwich", "dbplyr", "emmeans", "etm", "ff", "future", "HDF5Array", "imp4p", "mirt", "MuMIn", "PKI", "plotmo",
  "plotrix", "polylabelr", "ProjectTemplate", "proxy", "purrr", "qdapTools", "raster", "RcppArmadillo", "RCurl", "rex", "sf", "shinyFiles", "stopwords", 
  "svUnit", "systemfonts", "tinytex", "TSP", "withr", "xslt"))
install(c("beadarray", "clubSandwich", "dbplyr", "emmeans", "etm", "ff", "future", "HDF5Array", "imp4p", "mirt", "MuMIn", "PKI", "plotmo",
          "plotrix", "polylabelr", "ProjectTemplate", "proxy", "purrr", "qdapTools", "raster", "RcppArmadillo", "RCurl", "rex", "sf", "shinyFiles", "stopwords", 
          "svUnit", "systemfonts", "tinytex", "TSP", "withr", "xslt")) # Installation path not writeable, unable to update packages:

install.packages(c("beadarray", "clubSandwich", "dbplyr", "emmeans", "etm", "ff", "future", "HDF5Array", "imp4p", "mirt", "MuMIn", "PKI", "plotmo",
                   "plotrix", "polylabelr", "ProjectTemplate", "proxy", "purrr", "qdapTools", "raster", "RcppArmadillo", "RCurl", "rex", "sf", "shinyFiles", "stopwords", 
                   "svUnit", "systemfonts", "tinytex", "TSP", "withr", "xslt"))
install() # # Installation path not writeable, unable to update packages:

### FINALLY, I manually made myself (sandip_home) the owner of the folder/subfolders at "/Library/Frameworks/R.framework/Versions/3.6/Resources/library"
### and gave myself read and write permissions which resolved the issue. 
install()

library(limma)
#remove.packages(c("XLConnect", "XLConnectjars"), lib =  "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
#remove.packages("XLConnectjars", lib =  "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
#library(devtools)
#install_version("XLConnectJars", version = "0.2-12", repos = "http://cran.us.r-project.org")# removed later
#install_version("XLConnect", version = "0.2-12", repos = "http://cran.us.r-project.org") # updated later
options(java.parameters = "-Xmx8000m")
library(XLConnect)
library(Cairo) # nicer graphics, anti-aliased, etc.
library(NMF) # this package has a great annotated heatmap function - aheatmap
library(igraph)
library(gplots)

library(DEP)
library(RColorBrewer)
library(qdap)
library(DESeq2)
library(vsn)

nfl.all.btx <- loadWorkbook("~/Bioinformatics/David Proteomics/Beatrix_NFL KO_061019/ResultsSummary_10June2019_EK.xlsx") # btx ~ Beatrix
nfl.all.btx <- readWorksheet(nfl.all.btx, sheet = 2, header = T)
nfl.all.btx[1:5,1:15]
dim(nfl.all.btx) # 4739   33

## Import INA KO (Beatrix) dataset
ina.all.btx <- loadWorkbook("~/Bioinformatics/David Proteomics/INA KO 041519/INA KO_DataSummary_Beatrix_15Apr2019.xlsx") # btx ~ Beatrix
ina.all.btx <- readWorksheet(ina.all.btx, sheet = 1, header = T)
nfl.all.btx[1:5,1:15]
dim(ina.all.btx) # 3823   31

## Import TKO (Neubert) dataset
tko.all.nbt <- loadWorkbook("~/Bioinformatics/David Proteomics/Data files as of 040620 by David/IHL TKO_Neubert.xlsx") # nbt ~ Neubert
tko.all.nbt <- readWorksheet(tko.all.nbt, sheet = 1, header = T)
tko.all.nbt[1:5,1:9]  ## Make sure that the intensities are log10LFQ as opposed log2LFQ in Beatrix datasets
dim(tko.all.nbt) # 3953   37


## Removing unwanted columns from each data set before combining them
nfl.LFQ.btx <- nfl.all.btx[,c(1:12,22,23)]
head(nfl.LFQ.btx)
row.names(nfl.LFQ.btx) <- nfl.all.btx$Gene # Duplicate warning
sum(duplicated(nfl.all.btx$Gene)) # 27
sum(is.na(nfl.all.btx$Gene)) # 27
sum(nfl.all.btx$Gene == "") # NA

nfl.all.btx[duplicated(nfl.all.btx$Gene) | duplicated(nfl.all.btx$Gene, fromLast = T), c(1:12,23)]

nfl.LFQ.btx_uq <- nfl.all.btx_uq[!is.na(nfl.all.btx_uq$Gene), c(1:12,23)]
dim(nfl.LFQ.btx_uq) # 4711   13

sum(duplicated(nfl.LFQ.btx_uq$Gene)) # 1
nfl.LFQ.btx_uq <- nfl.LFQ.btx_uq[!duplicated(nfl.LFQ.btx_uq$Gene), ]
dim(nfl.LFQ.btx_uq) # 4711   13

head(nfl.LFQ.btx_uq)
colnames(nfl.LFQ.btx_uq)[7:12] <- paste0("WT",".", paste0(rep("LKO",6),1:6))

colnames(nfl.LFQ.btx_uq)

colnames(nfl.LFQ.btx_uq)[1:12] <- paste0(colnames(nfl.LFQ.btx_uq)[1:12], ".btx")
head(nfl.LFQ.btx_uq)

row.names(nfl.LFQ.btx_uq) <- nfl.LFQ.btx_uq$Gene; nfl.LFQ.btx_uq <- nfl.LFQ.btx_uq[,-13]
head(nfl.LFQ.btx_uq)
sapply(nfl.LFQ.btx_uq, class) # ALL character
nfl.LFQ.btx_uq[] <- as.data.frame(sapply(nfl.LFQ.btx_uq, as.numeric))
head(nfl.LFQ.btx_uq)

#nfl.LFQ.btx_uq <- nfl.LFQ.btx_uq[,-13]
#head(nfl.LFQ.btx_uq)

#library(mice)
#md.pattern(nfl.LFQ.btx_uq)
#md.pattern(nfl.LFQ.btx)
## Imputation of missing values
#library(impute)
#library(preprocessCore)

#nfl.LFQ.btx.impute <- impute.knn(as.matrix(nfl.LFQ.btx_uq))

## Maks summarized Experiment object
## create datTrait
datTraits <- loadWorkbook("datTraits_NF.xlsx")
datTraits <- readWorksheet(datTraits, sheet = 1, header = T)
datTraits

row.names(datTraits) <- datTraits$Sample
datTraits
nfl.btx_se <- SummarizedExperiment(assays = nfl.LFQ.btx_uq, 
                                   rowData = row.names(nfl.LFQ.btx_uq), 
                                   colData = datTraits[1:12,])

## Filter on missing values
# Plot a barplot of the protein identification overlap between samples
plot_frequency(nfl.btx_se)
dev.print(pdf, file="RPlot_ProtIdentOverlap_nfl.btx.pdf", height=6, width=8)


# Filter for proteins that are identified in 90% of samples i.e. 11 out of 12 replicates of at least one condition
nfl.btx_se_filt <- filter_missval(nfl.btx_se, thr = 1) 
# Error: 'name' and/or 'ID' columns are not present in 'nfl.btx_se'
nfl.btx_se2 <- nfl.btx_se
fData_nfl.btx <- data.frame("Gene.names"=rownames(nfl.btx_se))
mcols(nfl.btx_se2) <- DataFrame(mcols(nfl.btx_se2), fData_nfl.btx)
mcols(nfl.btx_se2)

head(assays(nfl.btx_se2))


data_filt2
dim(data_filt2) # 2863   12
plot_frequency(data_filt2)

############################################
ina.LFQ.btx <- ina.all.btx[,c(1:12,24,25)]
head(ina.LFQ.btx)
row.names(ina.LFQ.btx) <- ina.all.btx$Gene # Duplicate warning
sum(duplicated(ina.all.btx$Gene)) # 20
sum(is.na(ina.all.btx$Gene)) # 21
sum(ina.all.btx$Gene == "") # NA

ina.all.btx[is.na(ina.all.btx$Gene), c(1:12,24,25)]
# All these genes are NAs and corresponding Uniprots doesn't annotate for any gene (rather pseudogenes etc.)
ina.all.btx[duplicated(ina.all.btx$Gene) | duplicated(ina.all.btx$Gene, fromLast = T), c(1:12,24,25)]

ina.LFQ.btx_uq <- ina.all.btx[!is.na(ina.all.btx$Gene), c(1:12,25)]
head(ina.LFQ.btx_uq)
dim(nfl.LFQ.btx_uq) # 4712   13

sum(duplicated(ina.LFQ.btx_uq$Gene)) # 0
#nfl.LFQ.btx_uq <- nfl.all.LFQ_uq[!duplicated(nfl.LFQ.btx_uq$Gene), ]
#dim(nfl.LFQ.btx_uq) # 4711   13

#head(nfl.LFQ.btx)

row.names(ina.LFQ.btx_uq) <- ina.LFQ.btx_uq$Gene; ina.LFQ.btx_uq <- ina.LFQ.btx_uq[,-13]
head(ina.LFQ.btx_uq)
dim(ina.LFQ.btx_uq) # 3802   12

#nfl.LFQ.btx <- nfl.LFQ.btx[,-13]
#head(nfl.LFQ.btx)

head(ina.LFQ.btx)
colnames(ina.LFQ.btx_uq)[1:6] <- paste0("I", colnames(ina.LFQ.btx_uq)[1:6])
colnames(ina.LFQ.btx_uq)

colnames(ina.LFQ.btx_uq)[7:12] <- paste0("WT",".", paste0(rep("IKO",6),1:6))
head(ina.LFQ.btx_uq)

colnames(ina.LFQ.btx_uq)[1:12] <- paste0(colnames(ina.LFQ.btx_uq)[1:12], ".btx")
head(ina.LFQ.btx_uq)



tko.LFQ.nbt <- tko.all.nbt[,c(1:8,35,37)]
head(tko.LFQ.nbt)
#row.names(nfl.LFQ.btx) <- nfl.all.btx$Uniprot 

## Convert log10 values to log2 values
tko.LFQ.nbt[,-c(9,10)] <- as.data.frame(lapply(tko.LFQ.nbt[,-c(9,10)], function(v)
  ifelse(!is.na(v), log2(10^v), v)))

head(tko.LFQ.nbt)
tail(tko.LFQ.nbt)

colnames(tko.LFQ.nbt)[1:8] <- paste0(rep(c("WT","TKO"),each=4), rep(1:4,2))
head(tko.LFQ.nbt)

colnames(tko.LFQ.nbt)[1:4] <- paste0("WT",".", paste0(rep("TKO",4),1:4))
head(tko.LFQ.nbt)

colnames(tko.LFQ.nbt)[10] <- "Gene"
colnames(tko.LFQ.nbt)[1:9] <- paste0(colnames(tko.LFQ.nbt)[1:9], ".nbt")
head(tko.LFQ.nbt)

tko.LFQ.nbt_uq <- tko.LFQ.nbt[,-9]
sum(duplicated(tko.LFQ.nbt_uq$Gene)) # 176
sum(is.na(tko.LFQ.nbt_uq$Gene)) # 107

# Let's first remove all genes which doesn't code for any protein i.e. remove NAs
tko.LFQ.nbt_uq <- tko.LFQ.nbt_uq[!is.na(tko.LFQ.nbt_uq$Gene), ]
head(tko.LFQ.nbt_uq)
sum(duplicated(tko.LFQ.nbt_uq$Gene)) # 70

tko.LFQ.nbt_uq[duplicated(tko.LFQ.nbt_uq$Gene) | duplicated(tko.LFQ.nbt_uq$Gene, fromLast = T), ]

tko.LFQ.nbt_uq[duplicated(tko.LFQ.nbt_uq$Gene), ]
tko.LFQ.nbt_uq[duplicated(tko.LFQ.nbt_uq$Gene, fromLast = T), ]

#tko.LFQ.nbt_uq[sort(tko.LFQ.nbt_uq$Gene[duplicated(tko.LFQ.nbt_uq$Gene)]),]
#tko.LFQ.nbt_uq[match(tko.LFQ.nbt_uq$Gene[duplicated(tko.LFQ.nbt_uq$Gene, fromLast = T)], tko.LFQ.nbt_uq$Gene[duplicated(tko.LFQ.nbt_uq$Gene)]),]

#tko.LFQ.nbt_uq[duplicated(sort(tko.LFQ.nbt_uq$Gene)), ]
#tko.LFQ.nbt_uq[duplicated(sort(tko.LFQ.nbt_uq$Gene), fromLast = T), ]

sort(tko.LFQ.nbt_uq$Gene[duplicated(tko.LFQ.nbt_uq$Gene)])
sort(tko.LFQ.nbt_uq$Gene[duplicated(tko.LFQ.nbt_uq$Gene, fromLast = T)])

## Keep Genes with higher absolute protein abundance
library(data.table)
#as.data.table(a)[, .SD[which.max(abs(value))], by=id]
setDT(df)[, .SD[which.max(val2)], by = id]
top.test.tr.combd_uniq <- setDT(top.test.tr.combd)[, .SD[which.min(adj.P.Val)], by = ENTREZID]
library(dplyr)
df %>% 
  group_by(id) %>% 
  filter(val2==max(val2))

x %>% 
  group_by(var) %>% 
  summarise_each(funs(max))

#tko.LFQ.nbt_uq2 <- setDT(tko.LFQ.nbt_uq)[, .SD[which.max(rowSums(tko.LFQ.nbt_uq[,-9]))], by=Gene] # OR use this https://stackoverflow.com/a/25962964
#tko.LFQ.nbt_uq2 <- tko.LFQ.nbt_uq %>%
#  group_by(Gene) %>%
#  filter(max(rowSums(tko.LFQ.nbt_uq[,-9])))
tko.LFQ.nbt_uq2 <- tko.LFQ.nbt_uq %>%
  group_by(Gene) %>%
  summarize_each(funs(max))


sum(duplicated(tko.LFQ.nbt_uq2$Gene)) # 0
dim(tko.LFQ.nbt_uq2) # 3776    9
head(tko.LFQ.nbt_uq2)
head(tko.LFQ.nbt_uq[order(tko.LFQ.nbt_uq$Gene), ])
tko.LFQ.nbt_uq2[tko.LFQ.nbt_uq2$Gene == "Slc3a2", ]
tko.LFQ.nbt_uq[tko.LFQ.nbt_uq$Gene == "Slc3a2", ]

tko.LFQ.nbt_uq2 <- as.data.frame(tko.LFQ.nbt_uq2)
head(tko.LFQ.nbt_uq2)
row.names(tko.LFQ.nbt_uq2) <- tko.LFQ.nbt_uq2$Gene; tko.LFQ.nbt_uq2 <- tko.LFQ.nbt_uq2[, -1]
head(tko.LFQ.nbt_uq2)
tail(tko.LFQ.nbt_uq2)
dim(tko.LFQ.nbt_uq2) # 3776    8

## Merging 3 dataframes with Gene column
head(nfl.LFQ.btx)
head(ina.LFQ.btx)
head(tko.LFQ.nbt)
dim(nfl.LFQ.btx) # 4739   14
dim(ina.LFQ.btx) # 3823   14
dim(tko.LFQ.nbt) # 3953   10
###
head(nfl.LFQ.btx_uq)
head(ina.LFQ.btx_uq)
head(tko.LFQ.nbt_uq)
dim(nfl.LFQ.btx) # 4739   14
dim(ina.LFQ.btx) # 3823   14
dim(tko.LFQ.nbt) # 3953   10

df.merged <- merge(nfl.LFQ.btx, ina.LFQ.btx, by="Gene", all = T)
head(df.merged)
tail(df.merged)
dim(df.merged) # 5438   27

df.merged <- merge(df.merged, tko.LFQ.nbt, by="Gene", all = T)
head(df.merged)

run.seq <- function(x) as.numeric(ave(paste(x), x, FUN = seq_along))

L <- list(nfl.LFQ.btx, ina.LFQ.btx, tko.LFQ.nbt)
L2 <- lapply(L, function(x) cbind(x, run.seq = run.seq(x$Gene)))

df.merged <- Reduce(function(...) merge(..., all = TRUE), L2)[-2]
head(df.merged)
tail(df.merged)
dim(df.merged) # 8800   35









































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































