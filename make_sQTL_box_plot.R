#  sel <- 1
#  cluster_to_plot = row.names(resultsToPlot)[sel]
#  # get lowest p value junction in cluster to start
#  all_junctions = sigJunctions
#  junction_to_plot= all_junctions[all_junctions$clu==cluster_to_plot,]
#  junction_to_plot = junction_to_plot[ which( junction_to_plot$bpval == min(junction_to_plot$bpval) ), ]$pid
# #                    main_title = "test"
#                    vcf=vcf
#                     vcf_meta=vcf_meta
#                     exons_table = exons_table
#                     counts = clusters
#                     introns = annotatedClusters
#                     cluster_ids = annotatedClusters$clusterID
#                     snp_pos = resultsToPlot[sel,]$SNP_pos
#                     snp = resultsToPlot[sel,]$SNP
# # 
# #                    
# 


#' Make genotype x junction count box plots
#'
#' @import ggplot2
#' @import ggbeeswarm
#' @export
make_sQTL_box_plot <- function(
  cluster_to_plot,
  all_junctions = NA,
  junction_to_plot = NA,
  main_title = NA,
  exons_table = NULL,
  vcf = NULL,
  vcf_meta = NULL,
  cluster_ids = NULL,
  counts = NULL,
  introns = NULL,
  snp_pos=NA,
  junctionTable = NA,
  snp = snp ){
  
  #print("HELLO JACK")
  #print(junction_to_plot)
  
  #print(snp)
  
  stopifnot( snp %in% vcf_meta$SNP ) 
  
  # sometimes chr is missing
  
  if( !grepl("^chr", junction_to_plot)){
    junction_to_plot <- paste0("chr", junction_to_plot)
  }
  
  # subset VCF and get genotype groups
  vcfIndex <- which( vcf_meta$SNP == snp)
  VCF <- vcf[vcfIndex,10:ncol(vcf)]
  #table(t(VCF) )
  #for testing!
  VCF_meta <- vcf_meta[vcfIndex,]
  meta <- as.data.frame(t(VCF))
  meta$group=as.factor(meta$V1)
  group_names <- c(0,1,2)
  names(group_names) <- c(
    paste0( VCF_meta$REF, "/", VCF_meta$REF),
    paste0( VCF_meta$REF, "/", VCF_meta$ALT),
    paste0( VCF_meta$ALT, "/", VCF_meta$ALT)
  )
  
  #print(group_names)
  y <- t(counts[ cluster_ids==cluster_to_plot, ])
  
  # for each sample divide each junction count by the total for that sample
  normalisedCounts <- as.data.frame(sweep(y, 1, rowSums(y), "/"))
  
  genotypes <- as.data.frame(t(VCF))
  
  normalisedCounts$genotypeCode <- genotypes$V1[ match( row.names(normalisedCounts), row.names(genotypes))] 
  normalisedCounts <- normalisedCounts[ complete.cases(normalisedCounts),]
  
  normalisedCounts$genotype <- names(group_names)[ match(normalisedCounts$genotypeCode, group_names)]
  
  toPlot <- select( normalisedCounts, 
                    junction = junction_to_plot,
                    geno =  "genotype")
  toPlot$geno <- factor(toPlot$geno, levels = (names(group_names))) # this was reversed
  
  # get junction information for title
  junc <- mutate(junctionTable,
                 j = paste0( gsub("-",":", coord),":", clu)
                 ) %>%
    filter( j == junction_to_plot )
  
  values <- paste( signif(as.numeric(junc$Beta),3)," q =", signif(as.numeric(junc$q),3) )
  
  plot <- ggplot( data = toPlot, aes(x = geno, y = junction, group = geno ) ) + 
    geom_boxplot(outlier.colour = NA, fill = "orange") +
    geom_quasirandom(size = 0.8) + #coord_flip() +
    theme_classic() +
    theme(axis.title.y = element_text( angle = 90 ) ) +
    ylab("contribution to cluster") +
    xlab("") +
    labs(title = junc$coord,
         subtitle = bquote(beta==.(values)) )
      
  return(plot)
}
  