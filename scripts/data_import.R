library(GenomicFeatures)

#' Import CLIPper peaks
#'
#' @param path Path to the CLIPper BED files
#' @param snames String vector of sample names. All files must be named
#'   snames_clipper_peaks.bed
#'
#' @return List GRanges objects. The names are snames and each GRange is a peak.
#' @export
#'
#' @examples
read_clipper <- function(path = NA, snames = NA) {
  # The score is the p-value
  clipper <- list()
  for (sample in snames){
    clipper[[sample]] <- import(file.path(path, paste0(sample, 
                                                       "_deduplicated_clipper_peaks.bed")))
    clipper[[sample]]$gene_id <- str_split(clipper[[sample]]$name, "_", 
                                           simplify = TRUE)[,1]
  }
  clipper
}


#' Prepare gtf annotation
#'
#' Split GTF annotation into intron, exon, 3'UTR and 5'UTR. The exonic parts 
#' that overlap with UTRs are counted as UTRs! Intronic parts that overlap 
#' with exons or UTRs are not considered introns.
#'
#' @param gtf GRanges object
#'
#' @return GRangesList with exon, intron, 3'UTR and 5'UTR annotation
#' @export
#' 
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#'
#' @examples
prepare_anno <- function(gtf){
  exon <- gtf[gtf$type == "exon"] %>% unique
  gene <- gtf[gtf$type == "gene"] %>% unique
  five_utr <- gtf[gtf$type == "five_prime_utr"] %>% unique
  three_utr <- gtf[gtf$type == "three_prime_utr"] %>% unique
  
  ## We remove all 3' and 5' UTR regions that overlap with any exons
  exon_utr <- GenomicRanges::setdiff(exon, three_utr)
  exon_unique <- GenomicRanges::setdiff(exon_utr, five_utr) %>% unique
  ## copy metadata information from original exons
  olaps <- findOverlaps(exon_unique, exon)
  idx <- which(!duplicated(queryHits(olaps))) ## take the first match
  mcols(exon_unique) <- mcols(exon[subjectHits(olaps)[idx]])
  anno <- GRangesList(gene = gene, exon = exon_unique, three_prime_utr = three_utr, 
             five_prime_utr = five_utr)
  ## intron annotation
  txdb <- makeTxDbFromGRanges(gtf)
  introns <- unlist(intronsByTranscript(txdb))
  ## match intron metadata
  ## add one nucleotide upstream of the exon to allow overlap with intron
  olaps <- findOverlaps(introns, resize(exon, width(exon) +1, fix="end"))
  idx <- which(!duplicated(queryHits(olaps))) ## take the first match
  mcols(introns) <- mcols(exon[subjectHits(olaps)[idx]])
  
  ## remove the intronic parts that overlap with exons from other transcripts
  anno[["intron"]] <- GenomicRanges::setdiff(introns, c(anno[["exon"]], 
                                         anno[["three_prime_utr"]], 
                                         anno[["five_prime_utr"]]))
  olaps <- findOverlaps(anno[["intron"]], introns)
  idx <- which(!duplicated(queryHits(olaps))) ## take the first match
  mcols(anno[["intron"]]) <- mcols(introns[subjectHits(olaps)[idx]])
  anno
}


#' Annotate peaks
#'
#' @param peaks GRanges object with peaks
#' @param genes GRanges object with GTF gene annotations
#'
#' @return GRanges object with new metadata columns 'gene_id', 'gene_name' and
#'   'gene_biotype'
#' @export
#'
#' @examples
add_gene_annotation <- function(peaks, genes){
  m <- match(peaks$gene_id, genes$gene_id)
  res <- peaks[, "score"][!is.na(m)]
  mcols(res) <- cbind(mcols(res), mcols(genes[m[!is.na(m)], 
                                              c("gene_id", "gene_name", 
                                                "gene_biotype")]))
  res
}


