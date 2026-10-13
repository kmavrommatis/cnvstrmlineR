#' Get the  seqinfo for a genome
#'
#' Because the seqinfo functions in DBGenomeInfo require internet connection
#' we override them with a new function that reads this information from the package contents:

#'
#' @param genome (character) name of the genome ('hg38','hg19')
#' @docType methods
#' @return a seqinfo object
#'
getSeqInfo=function(genome){
  # seqinfo=list()
  # seqinfo[[genome]]= chrom_size[[ genome ]]
  # seqinfo[[genome]] <- seqinfo[[genome]] %>%
  #   magrittr::set_colnames(c('seqnames', 'seqlengths')) %>%
  #   dplyr::mutate( isCircular=FALSE,
  #                  genome=genome)
  #
  # seqinfo[[genome]]=GenomeInfoDb::Seqinfo(seqnames=seqinfo[[genome]]$seqnames,
  #                               seqlengths = seqinfo[[genome]]$seqlengths,
  #                               isCircular = seqinfo[[genome]]$isCircular,
  #                               genome=seqinfo[[genome]]$genome)
  #
  #
  # seqinfo=usethis::use_data(seqinfo)
  cnvstrmliner::seqinfo[[ genome ]]
}


