options(timeout = 30000)
library(biomaRt)

testdat <- read.table("/Users/TEHAL/OneDrive - A STAR/Misc/MetS-PH/snp4627.txt")
dim(testdat); head(testdat)

#Select mart
ensembl <- useEnsembl("snp",dataset = "hsapiens_snp",GRCh = "37")


snp_pos_all <- getBM(attributes=c("refsnp_id",
                              "chr_name",
                              "chrom_start",
                              "chrom_end",
                              "ensembl_gene_name"
                              ),
                 filters ="snp_filter", values =testdat$V1[1:1000], 
                 mart = ensembl, uniqueRows=TRUE)

dim(snp_pos_all)
head(snp_pos_all)

K <- 1000

for(i in 1:3) {
  
  snp_pos <- getBM(attributes=c("refsnp_id",
                   "chr_name",
                   "chrom_start",
                   "chrom_end",
                   "ensembl_gene_name"),
      filters ="snp_filter", values =testdat$V1[(1+K):(1000+K)], 
      mart = ensembl, uniqueRows=TRUE)

    snp_pos_all <- rbind(snp_pos_all, snp_pos)
    K <- K+1000  
    print(i)
}

dim(snp_pos_all); head(snp_pos_all)

snp_pos <- getBM(attributes=c("refsnp_id",
                              "chr_name",
                              "chrom_start",
                              "chrom_end",
                              "ensembl_gene_name"),
                 filters ="snp_filter", values =testdat$V1[4001:4627], 
                 mart = ensembl, uniqueRows=TRUE)

snp_pos_all <- rbind(snp_pos_all, snp_pos)
dim(snp_pos_all); head(snp_pos_all)

length(unique(snp_pos_all$refsnp_id))
length(unique(testdat$V1))

testdat$V1[!testdat$V1 %in% snp_pos_all$refsnp_id]

snp_pos <- getBM(attributes=c("refsnp_id",
                              "chr_name",
                              "chrom_start",
                              "chrom_end"),
                 filters ="snp_filter", values =c("rs142662971", "rs4072418",
                                                  "rs58010971", "rs62410026"), 
                 mart = ensembl, uniqueRows=TRUE)
snp_pos

table(snp_pos_all$chr_name)

snp_pos_filter <- snp_pos_all[-grep("HG|HSC", snp_pos_all$chr_name),]
dim(snp_pos_filter); head(snp_pos_filter)
length(unique(snp_pos_filter$refsnp_id))

table(snp_pos_filter$chr_name)
head(snp_pos_filter)

setwd("/Users/TEHAL/OneDrive - A STAR/Misc/MetS-PH/")
write.table(snp_pos_all, "MetS_SNP4627_Chr_Position.txt", row.names = F, 
            sep = "\t")

#snp_pos_all <- fread("MetS_SNPs_Chr_Position.txt", header = T)
#dim(snp_pos_all); head(snp_pos_all)

# to get gene symbol
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = "37")

gene_id <- getBM(attributes=c('ensembl_gene_id',
                   'hgnc_symbol',
                   'chromosome_name',
                   'start_position',
                   'end_position'), 
      filters = 'ensembl_gene_id', 
      values =snp_pos_all$ensembl_gene_name, 
      mart = ensembl)

dim(gene_id)
head(gene_id)
length(unique(gene_id$ensembl_gene_id))
length(unique(snp_pos_all$ensembl_gene_name))

library(dplyr)
gene_id_uniq <- gene_id |> 
  distinct(ensembl_gene_id, hgnc_symbol, start_position, .keep_all = TRUE)
dim(gene_id_uniq)
head(gene_id_uniq)

length(unique(gene_id_uniq$ensembl_gene_id))

write.table(gene_id_uniq, "GeneSymbol_mapped_to_1998ENSEMBLID_HG19_snp4627.txt", row.names = F,
            sep = "\t")

# merge snp and gene symbol

dat_merge <- merge(snp_pos_all, gene_id_uniq, by.x = "ensembl_gene_name",
                   by.y = "ensembl_gene_id", all.x = T)
dim(dat_merge); head(dat_merge)
tail(dat_merge)

length(which(!is.na(dat_merge$hgnc_symbol)))
length(dat_merge$refsnp_id[!is.na(dat_merge$hgnc_symbol)])
length(unique(dat_merge$refsnp_id[!is.na(dat_merge$hgnc_symbol)]))
length(unique(dat_merge$refsnp_id))

tmp <- as.data.frame(cbind(rep(NA, 4),  c("rs142662971", "rs4072418",
  "rs58010971", "rs62410026"), matrix(rep(NA, 28), nrow = 4)))
colnames(tmp) <- colnames(dat_merge)

dat_merge <- rbind(dat_merge, tmp)
length(unique(dat_merge$refsnp_id))

write.table(dat_merge, "All_MetS_snp4627_Annotated.txt",
            row.names = F, sep = "\t")

dat_merge_filter <- dat_merge[!is.na(dat_merge$hgnc_symbol),]
dim(dat_merge_filter); head(dat_merge_filter)

write.table(dat_merge_filter, "MetS_snp4627_with_Gene_Symbol.txt", row.names = F,
            sep = "\t")
