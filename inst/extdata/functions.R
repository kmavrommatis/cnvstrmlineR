
#' runQC
#' @description
#' runs the CNAqc tool
#' and stores rds data with the CNV object including cNAqc
#'
#' @returns an updated cnvobj which includes CNAqc object
#' @param cnvobj (cnvstrmlinr) the object with the CNV segments
#' @param muts (data.frame) a dataframe with mutations produced by the cnvstrmliner::loadMutationsFromVCF
#' @param n_bootstrap (numeric) the number of iterations to run to get the peaks QC
#' @param min_absolute_karyotype_mutations (integer) How many mutations to require in order to analyze a segment
#' @param rds_cache (character) the name of the file to store the cnvobj including the CNAqc
#' @param qc_output (character) the name of the file to store a summary of the QC
#'

runQC <- function( cnvobj, muts,
                 n_bootstrap=4,
                 min_absolute_karyotype_mutations=100,
                 rds_cache,
                 qc_output){

  logger::log_info("Calculating fragmentation statistics for {cnvobj@sample}")
  frag_summary=lapply(c(5,20,10), function(x){
    cnvobj=cnvstrmliner::detectOverfragmentation(cnvobj, minimum_segments_for_testing = x)
    list( min.segments=x,
          overfragmented.arms.count=cnvobj@quality$fragmentation$table %>% as.data.frame() %>%
            dplyr::filter( significant==TRUE) %>% dplyr::pull(significant)%>% length(),
          overfragmented.arms=cnvobj@quality$fragmentation$table %>% as.data.frame() %>%
            dplyr::filter( significant==TRUE) %>% dplyr::mutate(s=paste(seqnames,arm))%>% dplyr::pull(s) %>% paste(., collapse = ":")
    )
  }) %>% bind_rows()
  logger::log_info("Computing QC metrics for {cnvobj@sample}")
  if( ! is.na(cnvobj@purity)){
    cnvobj=cnvstrmliner::detectOverfragmentation(cnvobj, minimum_segments_for_testing = 20)
    cnvobj=cnvstrmliner::applyQCMetrics( cnvobj, muts, n_bootstrap=n_bootstrap, min_absolute_karyotype_mutations = min_absolute_karyotype_mutations)
  }else{
    logger::log_warn("Computing QC metrics for {cnvobj@sample} is not possible because there isn't a valid purity computed.")
  }

  fn= rds_cache
  logger::log_info("Saving objects with QC metrics in file {fn}")
  saveRDS(cnvobj, fn )


  if(!is.na(cnvobj@purity)){
    failed_chr=cnvobj@quality$cnaqcobj$scores_chr %>%subset(is.na(.)) %>% names()
  }
  logger::log_info("Building the summary")
  # #print( "Summary information")
  # logger::log_info("n_bootstrap {n_bootstrap}")
  # logger::log_info("min_absolute_karyotype {min_absolute_karyotype_mutations}")
  qc=  list(
    sample=cnvobj@sample,
    method=cnvobj@method,
    purity=cnvobj@purity,
    ploidy=cnvobj@effective_ploidy,
    n_bootstrap=n_bootstrap,
    min_absolute_karyotype_mutations=min_absolute_karyotype_mutations,
    score=ifelse(!is.na( cnvobj@purity), cnvobj@quality$cnaqcobj$peaks_analysis$score, NA),
    peaksQC=ifelse(!is.na( cnvobj@purity),cnvobj@quality$cnaqcobj$peaks_analysis$QC,NA),
    karyotypes.failed=ifelse(!is.na( cnvobj@purity),cnvobj@quality$cnaqcobj$peaks_analysis$matches %>%
                               dplyr::filter( QC == 'FAIL')%>%dplyr::pull( karyotype)%>% unique()%>%sort()%>% paste(., collapse = ','), NA),
    chromosomes.failed.count=ifelse(!is.na( cnvobj@purity),length(failed_chr),NA),
    chromosomes.failed=ifelse(!is.na( cnvobj@purity),paste(failed_chr , collapse=","),NA)
  ) %>% bind_cols()

  fn= qc_output
  logger::log_info("Saving summary of QC metrics for in file {fn}")
  write.table( cbind(qc, frag_summary), file=fn, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  return(cnvobj)
}
