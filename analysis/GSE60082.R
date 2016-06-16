library(GEOquery)
library(plyr)
library(data.table)
library(biomaRt)
library(ggplot2)
source("../bin/useful_functions.R")
source("../bin/common.R")

# Get the matrix
gse <- getGEO("GSE60082",GSEMatrix = F)
IDs <- Table(GSMList(gse)[[1]])
IDs$VALUE <- NULL

for (gsm in GSMList(gse)) {
  temp <- Table(gsm)
  colnames(temp)[2] <- gsm@header$geo_accession
  IDs <- merge(IDs, temp, by = "ID_REF",all=T)
  rm(temp)
}

# Convert transcript cluster id to gene symbol
temp <- fread("HuGene-1_0-st-v1.na35.hg19.transcript.csv")

mart.hs <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
probes <- temp[match(IDs$ID_REF,temp$transcript_cluster_id),]
probes <- as.data.frame(probes)
probes <- probes[,c("transcript_cluster_id", "probeset_id")]
genes <- getBM(attributes = c("affy_hugene_1_0_st_v1", "hgnc_symbol"), filters = "affy_hugene_1_0_st_v1", values = probes$probeset_id, mart = mart.hs)

probes <- merge(probes,genes, by.x = "probeset_id", by.y= "affy_hugene_1_0_st_v1")
probes$probeset_id <- NULL

# merge
IDs <- merge(IDs, probes, by.x = "ID_REF", by.y = "transcript_cluster_id")
IDs$ID_REF <- NULL

IDs <- IDs[IDs$hgnc_symbol != "",]
IDs <- rename(IDs, c("hgnc_symbol"= "symbol"))

IDs <- getMaxIQR(IDs, IDs$symbol)

write.table(IDs, "GSE60082_expVal.tsv",sep="\t",quote = FALSE)

# phenotype data
gse <- getGEO("GSE60082")[[1]]
pD <- get_formatted_phenoData(gse)
colnames(pD) <- c("sample_type","host_organism","nf1_associated","original_tumor","passage","geo_accession")
pD$passage <- NULL
pD$original_tumor <- NULL
write.table(pD, "GSE60082_phenotype_data.tsv", row.names = FALSE,sep="\t")

##############
# ESTIMATE
##############
est_score <- estimateMatrix("GSE60082_expVal.tsv","GSE60082")
write.table(est_score, "./estimate/GSE60082_estimate_score.tsv",quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

est_score <- data.frame(t(est_score))
est_score <- merge(est_score,pD,by.x="row.names",by.y="geo_accession")
ggplot(est_score, aes(x=StromalScore, y=TumorPurity, color=sample_type)) + geom_point(size=2)
ggsave("./estimate/stromalScore_tumorPurity_by_sample_type.png", width=8, height=4)

ggplot(est_score, aes(x=ImmuneScore, y=TumorPurity, color=sample_type)) + geom_point(size=2)
ggsave("./estimate/immuneScore_tumorPurity_by_sample_type.png", width=8, height=4)

ggplot(est_score, aes(x=ESTIMATEScore, y=TumorPurity, color=sample_type)) + geom_point(size=2)
ggsave("./estimate/ESTIMIATEScore_tumorPurity_by_sample_type.png", width=8, height=4)

###############
# GSVA
###############
expr <- as.matrix(IDs)
pD$geo_accession <- NULL

# immune signature
gsva.imm <- gsvaMatrix(expr,geneSet.bindea,"GSE60082","immSig")
gsvaHeatMap(mat = gsva.imm,fileName = "GSE60082_immSig_GSVA_heatmap.png",df = pD)

# hallmarks
gsva.hallmarks <- gsvaMatrix(expr,geneSet.hallmarks,"GSE60082","hallmarks")
gsvaHeatMap(mat = gsva.hallmarks,fileName = "GSE60082_hallmarks_GSVA_heatmap.png",df = pD)
