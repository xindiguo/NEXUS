library(GEOquery)
library(affy)
library(gcrma)
library(hgu133plus2hsrefseqprobe)
library(hgu133plus2hsrefseqcdf)
library(hgu133plus2hsrefseq.db)
library(mouse4302mmrefseqprobe)
library(mouse4302mmrefseqcdf)
library(mouse4302mmrefseq.db)
library(data.table)
library(limma)
source("../bin/common.R")

geo_id <- "GSE41747"
geo_study <- getGEO(geo_id)

getGEOSuppFiles(geo_id)
untar(paste0(geo_id,"_RAW.tar"), exdir="data")
cels <- list.files("./data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels <- list.files("data/", pattern = "CEL")

geo_study <- getGEO(geo_id)
phenotype_data1 <- pData(geo_study[[1]])
phenotype_data2 <- pData(geo_study[[2]])

humanData <- list.files(path = "./data/", pattern = "(jan)|(batch)")
mouseData <- list.files(path = "./data/", pattern = "Hyb")

# Human data
############
cdf <- "hgu133plus2hsrefseq"
#Read in the raw data from specified dir of CEL files
raw.data.human <- ReadAffy(verbose=TRUE, celfile.path="./data/", filenames = humanData, cdfname=cdf)
data.human <- gcrma(raw.data.human)
#Get the important stuff out of the data - the expression estimates for each array
gcrma.human <- exprs(data.human)

#Remove control probes: probenames starting with "AFFX"
control <- grep("^AFFX",rownames(gcrma.human),value=TRUE)
gcrma.human <- gcrma.human[!rownames(gcrma.human) %in% control,]

probes.human<- row.names(gcrma.human)
symbol.human <- unlist(mget(probes.human,hgu133plus2hsrefseqSYMBOL))
gcrma.human <- cbind(symbol.human,gcrma.human)

df.h <- as.data.frame(gcrma.human)
df.h <- data.frame(lapply(df.h, as.character), stringsAsFactors=FALSE)
df.h <- df.h[,-1]
df.h <- data.frame(lapply(df.h,as.numeric))
df.h$symbol <- symbol.human

rm(gcrma.human)

df.h <- df.h[complete.cases(df.h),]
df.h <- getMaxIQR(df.h,df.h$symbol)
colnames(df.h) <- sub("_.+","",colnames(df.h))

getMaxIQR <- function(df,refCol){  
  result <- do.call(rbind,lapply(unique(refCol),FUN =function(x){
    temp <- df[df$symbol == x,]
    temp$symbol <- NULL
    if(dim(temp)[1] > 1){
      result <- apply(temp, 1, function(y){
        return(IQR(y))
      })
      final <- temp[which.max(result),]
      row.names(final) <- x
      return(final)
    }else{
      final <- temp[1,]
      row.names(final) <- x
      return(final)
    }
  }))
  return(result)
}

# Mouse data
############
cdf <- "mouse4302mmrefseq"
#Read in the raw data from specified dir of CEL files
raw.data.m <- ReadAffy(verbose=TRUE, celfile.path="./data/", filenames = mouseData, cdfname=cdf)
data.m <- gcrma(raw.data.m)
#Get the important stuff out of the data - the expression estimates for each array
gcrma.m <- exprs(data.m)

#Remove control probes: probenames starting with "AFFX"
control <- grep("^AFFX",rownames(gcrma.m),value=TRUE)
gcrma.m <- gcrma.m[!rownames(gcrma.m) %in% control,]

probes.m<- row.names(gcrma.m)
symbol.m <- unlist(mget(probes.m,mouse4302mmrefseqSYMBOL))
gcrma.m <- cbind(symbol.m,gcrma.m)

df.m <- as.data.frame(gcrma.m)
df.m <- data.frame(lapply(df.m, as.character), stringsAsFactors=FALSE)
df.m <- df.m[,-1]
df.m <- data.frame(lapply(df.m,as.numeric))
df.m$symbol <- symbol.m

rm(gcrma.m)

df.m <- df.m[complete.cases(df.m),]
df.m <- getMaxIQR(df.m,df.m$symbol)
colnames(df.m) <- sub("_.+","",colnames(df.m))


# Combines and normalize 2 rna-seq data with different arrays
combineAndNormalizeData <- function(x,y,method = NULL){
  combinedData <- merge(x,y,by = "row.names")
  rownames(combinedData) <- combinedData$Row.names
  combinedData$Row.names <- NULL
  normalizedData <- normalizeBetweenArrays(as.matrix(combinedData),method = method)
  return(normalizedData)
}

### merge human and mouse data
row.names(df.m) <- toupper(row.names(df.m))
merged_data <- combineAndNormalizeData(df.h,df.m)
write.table(merged_data,"GSE41747_expVal.tsv",sep="\t",quote=FALSE)

# phenotype data
pD1 <- phenotype_data1[,c("title","geo_accession","organism_ch1","characteristics_ch1")]
pD1 <- rename(pD1,c("organism_ch1" = "organsim", "characteristics_ch1" = "tissue"))
pD1$genotype_variation <- NA
pD2 <- phenotype_data2[,c("title","geo_accession","organism_ch1","characteristics_ch1","characteristics_ch1.1")]
pD2 <- rename(pD2,c("organism_ch1" = "organsim", "characteristics_ch1"="genotype_variation","characteristics_ch1.1" = "tissue"))

pD <- rbind(pD1,pD2)

pD$tissue <- sub(".+\\: ","",pD$tissue)
pD$genotype_variation <- sub(".+\\: ","",pD$genotype_variation)
pD <- data.frame(lapply(pD, as.character), stringsAsFactors=FALSE)

pD$title <- NULL
write.table(pD,"GSE41747_phenotype_data.tsv",sep="\t",row.names = FALSE)

############
# GSVA
############
expr <- as.matrix(merged_data)
row.names(pD) <- pD$geo_accession
pD$geo_accession <- NULL

# immune signature
gsva.imm <- gsvaMatrix(expr,geneSet.bindea,"GSE41747","immSig")
gsvaHeatMap(mat = gsva.imm,fileName = "GSE41747_immSig_GSVA_heatmap.png",df = pD)

# hallmarks
gsva.hallmarks <- gsvaMatrix(expr,geneSet.hallmarks,"GSE41747","hallmarks")
gsvaHeatMap(mat = gsva.hallmarks,fileName = "GSE41747_hallmarks_GSVA_heatmap.png",df = pD)


#############
# ESTIMATE
#############
library(ggplot2)
est_score <- estimateMatrix("GSE41747_expVal.tsv", "GSE41747")
write.table(est_score, "./estimate/GSE41747_estimate_score.tsv",quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

est_score <- data.frame(t(est_score))
est_score <- merge(est_score,pD,by="row.names")

ggplot(est_score, aes(x=StromalScore, y=TumorPurity, color=tissue)) + geom_point(size=2)
ggsave("./estimate/stromalScore_tumorPurity_by_tissue.png", width=8, height=4)

ggplot(est_score, aes(x=ImmuneScore, y=TumorPurity, color=tissue)) + geom_point(size=2)
ggsave("./estimate/immuneScore_tumorPurity_by_tissue.png", width=8, height=4)

ggplot(est_score, aes(x=ESTIMATEScore, y=TumorPurity, color=tissue)) + geom_point(size=2)
ggsave("./estimate/ESTIMIATEScore_tumorPurity_by_tissue.png", width=8, height=4)

