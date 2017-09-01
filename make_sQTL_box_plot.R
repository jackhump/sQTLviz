sel <- 2
make_sQTL_box_plot(cluster_to_plot = row.names(resultsToPlot)[sel],
                   all_junctions = sigJunctions,
                   junction_to_plot = 
                   main_title = "test",
                   vcf=vcf,
                   vcf_meta=vcf_meta,
                   exons_table = exons_table,
                   counts = clusters,
                   introns = annotatedClusters,
                   cluster_ids = annotatedClusters$clusterID,
                   snp_pos = resultsToPlot[sel,]$SNP_pos,
                   snp = resultsToPlot[sel,]$SNP )




#' Make genotype x junction count box plots
#'
#' @import ggplot2
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
  snp = snp ){
  
  
}
  