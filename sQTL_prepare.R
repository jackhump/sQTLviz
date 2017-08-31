library(dplyr)
library(data.table)
#library(vcfR)
library(leafcutter)
library(stringr)

setwd("/Users/Jack/Documents/SQTL_LeafViz/")

permutation_res <- "data/permutations.all.CMC.txt.gz.0.05.bh.txt"
VCF = "data/genotypes_MAF1.vcf.gz"
clusters_table <- "data/CMC_perind_numers.counts_renamed.gz" # now just the counts - don't need to strip ratios

# PRE-REQUISITES:

# annotation
annotation_code <- "/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  #read_table(exon_file)
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
res <- as.data.frame(fread(permutation_res, header =TRUE), stringsAsFactors = FALSE)

genotypes <- unique(res$dummy2)
genotypes_file <- "sig_snps.txt"
write.table(genotypes, genotypes_file, col.names = FALSE, row.names = FALSE, quote = FALSE)

# PREPARE GENOTYPES

# use vcftools to filter out just the snps and samples required
print("filtering VCF")
cmd <- paste( "vcftools --gzvcf", VCF,  "--snps", genotypes_file, "--keep", samples_file, "--recode --stdout" )
vcf <- fread(cmd)

# round genotypes to 0,1,2
vcf <- vcf[ complete.cases(vcf), ]
vcf[,10:ncol(vcf)] <- as.data.frame(apply( vcf[,10:ncol(vcf)], MAR = c(1,2), FUN = round ))

# PREPARE CLUSTERS

sigClusters <- str_split_fixed(res[,1], ":",4)[,4]
introns <- get_intron_meta(row.names(clusters) )
keepClusters <- match(introns$clu,sigClusters)

# remove non-significant clusters
introns <- introns[ !is.na(keepClusters),]
clusters <- clusters[ !is.na(keepClusters),]

#row.names(clusters) <- clusters$chrom; clusters$chrom <- NULL

## if counts are fractions then split off just the numerator
### TODO - make this optional
# now using CMC_perind_numers.counts.gz - no need for this
# clusters[,2:ncol(clusters)] <- as.data.frame( apply( clusters[,2:ncol(clusters)], MAR = c(1,2), FUN = function(x){
#   num <- str_split_fixed(x, "/",2)[,1]
#   return( as.numeric(as.character(num) ) )
# }) )


# PREPARE RESULTS
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
sigJunctions <- cbind( get_intron_meta( res$pid), res, get_snp_meta(res$dummy2))

resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>% 
    dplyr::summarise( chr = first(chr),
                      start = min(start),
                      end = max(end),
                      snp = first(snp_ID),
                      snp_chr = first(snp_chr),
                      pos = first(snp_pos),
                      FDR = first(bpval) ) %>%
    dplyr::arrange(FDR)


# for testing only - generate UCSC link
## genome.ucsc.edu/cgi-bin/hgTracks?&org=human&db=hg19&position=chr5:96075822-96076970&hgFind.matches=rs7724759
# pmin and pmax are vectorised 
# starts <- pmin( resultsByCluster$start, resultsByCluster$pos - 100)
# ends <- pmax( resultsByCluster$end, resultsByCluster$pos + 100)
# URLS <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?&org=human&db=hg19&position=chr", 
#                 resultsByCluster$chr,":",starts,"-", ends,"&hgFind.matches=", resultsByCluster$snp  )

# ANNOTATE JUNCTIONS
intersect_introns <- function(introns){
  all.introns <- introns
  all.introns$start <- as.numeric(all.introns$start)
  all.introns$end <- as.numeric(all.introns$end)
  
  
  # for each splice site write out a bed file  
  all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID = clu)
  
  all.fiveprime <- data.frame( chr = all.introns$chr,
                               start = all.introns$start,
                               end = all.introns$start + 1,
                               clusterID = all.introns$clu)
  all.threeprime <- data.frame( chr = all.introns$chr,
                                start = all.introns$end,
                                end = all.introns$end + 1,
                                clusterID = all.introns$clu)
  all.file <- "all_junctions.bed"
  all.fiveprime.file <- "all_fiveprime.bed"
  all.threeprime.file <- "all_threeprime.bed"
  
  write.table( all.junctions, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  
  print( "BedTools intersect junctions with list of known splice sites")
  # first match junctions 
  all.introns.cmd <- paste0("bedtools intersect -a ", all.file, " -b ", all_introns, " -wa -wb -loj -f 1" )
  
  all.introns_intersect <- fread(all.introns.cmd)
  
  # intersect with bedtools to find the annotations of each splice site
  threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ",threeprime_file, " -wa -wb -loj -f 1" )
  
  threeprime_intersect <- fread(threeprime.cmd)
  
  fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", fiveprime_file, " -wa -wb -loj -f 1" )
  
  fiveprime_intersect <- fread(fiveprime.cmd)
  
  # remove temporary files
  rm.cmd <- paste("rm ", all.file, all.fiveprime.file, all.threeprime.file) 
  system(rm.cmd)
  
  return( list(threeprime_intersect, fiveprime_intersect,all.introns_intersect))
}

intersects <- intersect_introns(introns)
threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]
print("Annotating junctions")
  
annotate_single_cluster <- function(introns, clu, cluIndex){
    #print(clu)
    # for each intron in the cluster, check for coverage of both
    # output a vector of string descriptions 
    cluster <- introns[ introns$clu == clu , ]
    cluster$start <- as.integer(cluster$start)
    cluster$end <- as.integer(cluster$end)
    # subset intersects by clusterID (V4)
    fprimeClu <- fiveprime_intersect[ V4 == clu,]
    tprimeClu <- threeprime_intersect[ V4 == clu,]
    bothSSClu <- all.introns_intersect[ V4 == clu,]
    # for each intron in the cluster:
    #   create vector of overlapping splice sites, indexed by the row of the intersect
    # five prime splice sites
    fprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which( names(cluster) == "start" )
      fprimeClu[   
        V1 == x[chr] & 
          V2 == as.numeric( x[start] ),]
    } )
    # three prime splice sites
    tprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      end <- which( names(cluster) == "end" )
      tprimeClu[   
        V1 == x[chr] & 
        V2 == as.numeric( x[end] ),]
    } )
    
    # both splice sites
    bothSS <-  apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which(names(cluster) == "start")
      end <- which( names(cluster) == "end" )
      
      bothSSClu[   
        V6 == as.numeric( x[start] ) &
        V7 == as.numeric( x[end] ) ,]
    } )
    
    # find gene and ensemblID by the most represented gene among all the splice sites
    cluster_genes <- names(sort(table(do.call( what = rbind, tprime )$V8), decreasing = TRUE ))
    
    cluster_gene <- cluster_genes[ cluster_genes != "." ][1]
    # if no cluster gene found then leave as "."
    if( length(cluster_gene) == 0){
      cluster_gene == "."
    }
    # do the same for EnsemblID
    cluster_ensemblIDs <- names(sort(table(do.call( what = rbind, tprime )$V9), decreasing = TRUE ))
    cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
    if( length( cluster_ensemblID ) == 0 ){
      cluster_ensemblID == "."
    }
    
    verdict <- c()
    coord <- c()
    gene <- c()
    ensemblID <- c()
    transcripts <- list() 
    
    for( intron in 1:nrow(cluster) ){
      coord[intron] <- paste0(cluster[intron,]$chr,":", cluster[intron,]$start,"-", cluster[intron,]$end )
      
      gene[intron] <- cluster_gene
      ensemblID[intron] <- cluster_ensemblID
      
      # for each intron create vector of all transcripts that contain both splice sites
      #transcripts[[intron]] <- unique( intersect( tprime[[intron]]$V10, fprime[[intron]]$V10 ) ) 
      
      verdict[intron] <- "error"
      if( # if neither are annotated
        all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 == "." )
      ){ verdict[intron] <- "cryptic_unanchored"
      }
      if( # if only one is annotated
        all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 != "." )
      ){ verdict[intron] <- "cryptic_threeprime"
      }
      if(
        all( tprime[[intron]]$V5 != ".") & all( fprime[[intron]]$V5 == "." )
      ){ verdict[intron] <- "cryptic_fiveprime"
      }
      if( # if both splice sites are annotated
        all( tprime[[intron]]$V5 != "." ) & all( fprime[[intron]]$V5 != "." )
      ){ 
        # test if the splice sites are paired in a known intron
        if( all(bothSS[[intron]]$V5 != ".") ){
          verdict[intron] <- "annotated"
        }else{ # both are annotated but never in the same junction
          verdict[intron] <- "novel annotated pair"
        }
      }

    }
    if( cluIndex %% 500 == 0 ){
      print( paste("processed", cluIndex, "clusters" ))
    }
    #print(clu)
    
    return(
      data.frame(
        clusterID = clu,
        coord = coord,
        gene = gene,
        ensemblID = ensemblID,
        verdict = verdict,
        stringsAsFactors = FALSE)
    )
}

uniqueClusters <- unique( introns$clu ) 
#uniqueClusters <- uniqueClusters[7010:7020]
annotatedClusters <- lapply( 1:length(uniqueClusters), FUN = function(i) annotate_single_cluster(introns, clu = uniqueClusters[i], cluIndex = i) )
annotatedClusters <- do.call(rbind, annotatedClusters)

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."


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
                           SNP_pos,                           gene = gene,
                           cluster_pos,
                           q = FDR
) )
row.names(resultsToPlot) <- resultsByCluster$clu
resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

save.image("all_data.Rdata")
print("saving objects")
save( annotatedClusters, # every junction needed
      resultsToPlot, #significant clusters and the most significant SNP
      clusters, # junction counts for each sample
      vcf, # the genotypes of each sample
      #counts, 
      #meta, 
      exons_table, # the annotation
      #pca, 
      #intron_summary, 
      #cluster_summary, 
      #introns_to_plot,
      #cluster_ids,
      #sample_table,
      annotation_code,
      code,
      file = paste0( "sQTL_results.Rdata")
)

# to do - cut down size of exon table to increase speed of querying
#allGenes <- unique(resultsToPlot$gene)
#exons_table_cut <- exons_table[ exons_table$gene_name %in% allGenes ,]

