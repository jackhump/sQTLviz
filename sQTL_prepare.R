library(dplyr)
library(data.table)
#library(vcfR)
library(leafcutter)
library(stringr)
library(readr)

setwd("/Users/Jack/Documents/SQTL_LeafViz/")

permutation_res <- "data/permutations.all.CMC.txt.gz.0.05.bh.txt"

# Yang now wants all associations at P < 1e-5
#permutation_res <- "data/permutations.all.CMC.p0.00001.txt"

# other files
VCF = "data/genotypes_MAF1.vcf.gz"
clusters_table <- "data/CMC_perind_numers.counts_renamed.gz" # now just the counts - don't need to strip ratios

#####
# 35 SNPs from Nalls et al 
## present the top SNP x cluster associations
######

GWAS_SNPs <- "data/PD_GWAS_SNPs.txt"
GWAS <- read.table(GWAS_SNPs, header=TRUE, stringsAsFactors = FALSE)

gwas_snps <- GWAS$SNP.orig

GWAS_SNP_associations <- "data/PD_SNPs_associations.all.CMC.txt"
gwas_assoc <- read.table(GWAS_SNP_associations,header=FALSE, stringsAsFactors = FALSE)

gwas_assoc <- gwas_assoc %>%
              mutate( SNP = str_split_fixed(V2, "\\.", 3)[,1] ) %>%
              filter( SNP %in% gwas_snps ) %>%
              rename( "cluster" = V1, "SNP_full" = V2, "P" = V4 )

top_gwas_assoc <- group_by(gwas_assoc, SNP) %>%
                  summarise( P = min(P) ) %>%
                  left_join( gwas_assoc, by = c("SNP", "P" ) )
# for finding in VCF
gwas_full_snps <- top_gwas_assoc$SNP_full
gwas_clusters <- str_split_fixed(top_gwas_assoc$cluster, "\\:", 4)[,4]

#### Yang's Top Two Hits from the TWAS - MAPT and MTOR
# zgrep for them in all.associations
Yang_SNPs <- "data/Yang_chosen_SNPs_associations.all.CMC.txt"
Yang_assoc <- read.table(Yang_SNPs, header=FALSE, stringsAsFactors = FALSE)

Yang_assoc <- Yang_assoc %>%
  mutate( SNP = str_split_fixed(V2, "\\.", 3)[,1] ) %>%
  rename( "cluster" = V1, "SNP_full" = V2, "P" = V4 )

top_yang_assoc <- group_by(Yang_assoc, SNP) %>%
  summarise( P = min(P) ) %>%
  left_join( Yang_assoc, by = c("SNP", "P" ) )

yang_full_snps <- top_yang_assoc$SNP_full
yang_clusters <- str_split_fixed(top_yang_assoc$cluster, "\\:", 4)[,4]

# PRE-REQUISITES:

# annotation
annotation_code <- "/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  as.data.frame(fread(paste("zless",exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}

## CHECK VCFTOOLS IS INSTALLED
if( !file.exists(system("which vcftools", intern=TRUE) )){
  stop("vcftools is not in your $PATH - have you installed it?")
}
## BEDTOOLS
if( !file.exists(system("which bedtools", intern=TRUE) )){
  stop("bedtools is not in your $PATH - have you installed it?")
}

# START

# read in clusters
#clusters <- as.data.frame(fread(paste("zless", clusters_table), sep = " ", row.names = 1,  header = TRUE, stringsAsFactors = FALSE) )

print("reading in clusters")
clusters <- read.table(clusters_table, header=TRUE)

# find and harmonise sample names
samples <- names(clusters)[2:ncol(clusters)]

#convert sample names 
samples <- gsub( "merged_", "", gsub( "_RNA_PFC", "", samples ) ) 
names(clusters)[2:ncol(clusters)] <- samples

# write out samples
samples_file <- "used_samples.txt"
write.table(samples, samples_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

# read in junction x snp results
print("reading in results")
# q < 0.05
res <- as.data.frame(fread(permutation_res, header =TRUE), stringsAsFactors = FALSE)
genotypes <- unique(res$dummy2)

# p < 1e-5
#res <- as.data.frame(fread(permutation_res, header =FALSE), stringsAsFactors = FALSE)
#genotypes <- unique(res$V6)

# add in Nalls' GWAS SNPs and Yang's preferred SNPs
genotypes <- c(genotypes, gwas_full_snps,yang_full_snps)

genotypes_file <- "data/sig_snps_plus_GWAS_Yang_snps.txt"
write.table(genotypes, file = genotypes_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

######
# PREPARE GENOTYPES
######

# use vcftools to filter out just the snps and samples required
print("filtering VCF")
cmd <- paste( "vcftools --gzvcf", VCF,  "--snps", genotypes_file, "--keep", samples_file, "--recode --stdout" )
vcf <- fread(cmd)

# round genotypes to 0,1,2
vcf <- vcf[ complete.cases(vcf), ]
vcf[,10:ncol(vcf)] <- as.data.frame(apply( vcf[,10:ncol(vcf)], MAR = c(1,2), FUN = round ))

get_vcf_meta <- function(vcf){
  # take a VCF and split out the meta data 
  # to use for querying
  vcf_meta <-  as.data.frame(str_split_fixed(vcf$ID, "\\.", 3), stringsAsFactors = FALSE)
  vcf_meta$SNP_pos <- paste(vcf_meta$V2, vcf_meta$V3, sep = ":")
  vcf_meta <- data.frame( SNP = vcf_meta$V1, SNP_pos = vcf_meta$SNP_pos, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
  return(vcf_meta)
}
vcf_meta <- get_vcf_meta(vcf)

##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- str_split_fixed(res[,1], ":",4)[,4]
# add Nalls' GWAS SNP clusters and Yang's too!
sigClusters <- c(sigClusters, gwas_clusters, yang_clusters)

introns <- get_intron_meta(row.names(clusters) )
keepClusters <- match(introns$clu,sigClusters)

# remove non-significant (or non-GWAS SNP-associated) clusters
introns <- introns[ !is.na(keepClusters),]
clusters <- clusters[ !is.na(keepClusters),]

# rearrange sample columns in clusters so they match the VCF
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))
#row.names(clusters) <- clusters$chrom; clusters$chrom <- NULL

# for each cluster work out mean proportion of each junction
# remove junctions < 1% contribution

juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount) 
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate( 
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters) ) %>%
  split( .$clu ) %>%
  purrr::map_df( juncProp ) %>%
  mutate( clu = as.character(.$clu))

introns_to_plot <- introns_to_plot[ splitClusters$prop >= 0.01,]
clusters <- clusters[ splitClusters$prop >= 0.01,]
introns <- introns[ splitClusters$prop >= 0.01,]

####################
# ANNOTATE JUNCTIONS
####################

# functions now in separate script
source("sQTLviz/sQTL_annotation_functions.R")

intersects <- intersect_introns(introns)
threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique( introns$clu ) 
#uniqueClusters <- uniqueClusters[7010:7020]
#annotatedClusters <- lapply( seq_along(uniqueClusters),
#                             FUN = function(i) annotate_single_cluster(introns, clu = uniqueClusters[i], cluIndex = i) )
#annotatedClusters <- do.call(rbind, annotatedClusters)

# replace with purrr?
annotatedClusters <- purrr::map_df( seq_along(uniqueClusters),
                                    ~annotate_single_cluster( introns, clu = uniqueClusters[.], cluIndex = .  ) )

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."


#################
# PREPARE RESULTS
#################

get_snp_meta <- function(snps){
  snp_meta <-  as.data.frame(str_split_fixed(snps, "\\.", 3), stringsAsFactors = FALSE)
  colnames(snp_meta) <- c("snp_ID", "snp_chr", "snp_pos")
  # a few snps don't have any IDs - just the coordinates
  noID <- snp_meta[ grepl("\\.", snp_meta$snp_pos),]
  noID <- select(noID, snp_pos, snp_ID, snp_chr)
  names(noID) <- c("snp_ID", "snp_chr","snp_pos")
  # put back in
  snp_meta[ grepl("\\.", snp_meta$snp_pos),] <- noID
  
  snp_meta$snp_pos <- as.numeric(snp_meta$snp_pos)
  return(snp_meta)
}

## the results table should consist of a list of clusters with their most significant SNP
# bind intron_meta, intron, snp, p value, snp_meta
#sigJunctions <- cbind( get_intron_meta( res[,1]), res[, c(1,6,11)], get_snp_meta(res[,6]))

# use associations with q < 0.05 cut-off. 
# Bind together metadata with original results
sigJunctions <- cbind( get_intron_meta( res[,1]), 
                       res[, c(1,2,3)],
                       get_snp_meta(res[,2]))

# create identical table for GWAS junctions
NallsJunctions <- cbind( get_intron_meta(top_gwas_assoc$cluster),
                         top_gwas_assoc$cluster,
                         top_gwas_assoc$SNP_full,
                         top_gwas_assoc$P,
                         get_snp_meta(top_gwas_assoc$SNP_full))

YangJunctions <- cbind( get_intron_meta(top_yang_assoc$cluster),
                        top_yang_assoc$cluster,
                        top_yang_assoc$SNP_full,
                        top_yang_assoc$P,
                        get_snp_meta(top_yang_assoc$SNP_full))
# bind together - all will be accessible in same table
names(NallsJunctions) <- names(sigJunctions)
names(YangJunctions) <- names(sigJunctions)
sigJunctions <- rbind(sigJunctions, NallsJunctions)
sigJunctions <- rbind(sigJunctions, YangJunctions)
# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

#names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>% 
    dplyr::summarise( chr = first(chr),
                      start = min(start),
                      end = max(end),
                      snp = first(snp_ID),
                      snp_chr = first(snp_chr),
                      pos = first(snp_pos),
                      FDR = first(bpval) ) %>%
    dplyr::arrange(FDR)

# for Nalls and Yang SNPs don't thin at all - we want those SNPs!

####
## PREPARE FOR SHINY
####

code <- "test"
annotation_code <- "gencode_hg19"

resultsByCluster$gene <- annotatedClusters$gene[ match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

# fix coords without "chr"
if( all( !grepl("chr", sample(resultsByCluster$chr, 100)) ) ){
  resultsByCluster$chr <- paste0("chr", resultsByCluster$chr)
}

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr,":", resultsByCluster$start,"-",resultsByCluster$end)

resultsToPlot <- as.data.frame( select( resultsByCluster,
                           SNP = snp,
                           SNP_pos, 
                           gene = gene,
                           cluster_pos,
                           q = FDR
) )

row.names(resultsToPlot) <- resultsByCluster$clu
resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

# create separate table for Nalls results
GWAS_metadata <- read.table("data/PD_GWAS_SNPs.txt", header=TRUE, stringsAsFactors = FALSE)

GWASresults <- select(NallsJunctions,
                      clu,
                      chr,
                      start,
                      end,
                      SNP = snp_ID,
                      snp_chr,
                      pos = snp_pos,
                      q = bpval # is this really FDR or uncorrected P?
                      ) %>%
    mutate( 
      SNP_pos = paste0(snp_chr, ":", pos),
      cluster_pos = paste0("chr", chr, ":", start, "-", end),
      gene = annotatedClusters$gene[ match(clu, annotatedClusters$clusterID)],
      "GWAS P" = GWAS_metadata$P.joint[ match( SNP, GWAS_metadata$SNP.orig)],
      q <- signif(q,  digits = 3)
      ) %>%
  select( "GWAS P", SNP, SNP_pos, gene, cluster_pos, q) %>%
 # arrange("GWAS p") %>%
  as.data.frame(stringsAsFactors = FALSE)
row.names(GWASresults) <- NallsJunctions$clu

YangResults <- select(YangJunctions,
                       clu,
                       chr,
                       start,
                       end,
                       SNP = snp_ID,
                       snp_chr,
                       pos = snp_pos,
                       q = bpval # is this really FDR or uncorrected P?
) %>%
  mutate( 
    SNP_pos = paste0(snp_chr, ":", pos),
    cluster_pos = paste0("chr", chr, ":", start, "-", end),
    gene = annotatedClusters$gene[ match(clu, annotatedClusters$clusterID)],
    q <- signif(q,  digits = 3)
  ) %>%
  select( SNP, SNP_pos, gene, cluster_pos, q) %>%
  # arrange("GWAS p") %>%
  as.data.frame(stringsAsFactors = FALSE)
row.names(YangResults) <- YangJunctions$clu

# get the Betas and per-junction q values
# add in full permutation results to get Beta for each junction
permutation_full_res <- "data/permutations.all.CMC.txt.gz"
perm_full <- read_delim( permutation_full_res,
                         col_names = c("clusterID", "V2","V3","V4","V5","SNP","V7","V8","Beta","V10","FDR"),
                         delim = " "
)
perm_clean <- select(perm_full, clusterID, SNP, Beta, FDR)
perm_clean <- cbind( perm_clean,
                     get_intron_meta(perm_full$clusterID),
                     get_snp_meta(perm_full$SNP) )


# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  rename(snp_ID = SNP) %>%
  left_join( perm_clean, 
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join( annotatedClusters, 
             by = c("clu" = "clusterID", "coord" )
  ) %>%
  select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))



save.image("data/all_data.Rdata")
print("saving objects")
save( annotatedClusters, # every junction needed
      sigJunctions, # every junction x SNP interaction
      resultsToPlot, #significant clusters and the most significant SNP
      GWASresults, # associations with SNPs from a PD GWAS
      YangResults, # associations with Yang's TWAS hit SNPs
      clusters, # junction counts for each sample
      vcf,# the genotypes of each sample
      vcf_meta, # the vcf metadata
      introns_to_plot, # all the intron positions
      #counts, 
      #meta, 
      exons_table, # the annotation
      junctionTable, # the junctions to display for each cluster
      #pca, 
      #intron_summary, 
      #cluster_summary, 
      #introns_to_plot,
      #cluster_ids,
      #sample_table,
      annotation_code,
      code,
      file = paste0( "sQTLviz/sQTL_results.Rdata")
)

# to do - cut down size of exon table to increase speed of querying
#allGenes <- unique(resultsToPlot$gene)
#exons_table_cut <- exons_table[ exons_table$gene_name %in% allGenes ,]

