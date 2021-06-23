## Given a list of regions and an IGV session, write a batch script that takes screenshots of the regions.

## Customize IGV manually and load the batch scrip: Tools --> Load Batch Script

library(here)
library(rtracklayer)

gdf <- read.table(here("analysis", "clipper_analysis", "gene_cluster_count_merge_peaks.txt"), 
                  header = TRUE)
out_dir <- here("IGV", "RBDm_6M_specific_peaks")
gtf_file <- here("reference", "Homo_sapiens.GRCh38.98.sorted.gtf")
gtf <- import(gtf_file)

rg <- list("6M" =  gdf %>% dplyr::filter(nclus_WT == 0 & nclus_RBDm == 0) %>% 
             pull("gene_id") %>% 
             sapply(., function(x) gtf[gtf$gene_id == x & gtf$type == "gene"]), 
           RBDm = gdf %>% dplyr::filter(nclus_WT == 0 & nclus_6M == 0) %>% 
             pull("gene_id") %>% 
             sapply(., function(x) gtf[gtf$gene_id == x & gtf$type == "gene"]))


snapshot_gene <- function(peaks, prefix) {
  s <- min(start(peaks))
  e <- max(end(peaks))
  
  if(e-s < 1000){ ## at least 1000bp are shown
    border <- floor(1000-(e-s)/2)
  }
  else{
    border <- 150
  }
  s <- s - border
  e <- e + border
  
  cat(paste0("goto chr", unique(seqnames(peaks)), ":", s, "-", e), 
      fill = TRUE)
  cat(paste0("snapshot ", prefix, "_",  unique(peaks$gene_name), "_", unique(seqnames(peaks)), "-", s, "-", e, ".png"), 
      fill = TRUE)
}

## Write the batch script to file
sink(file.path(out_dir, "IGV_screenshot.bat"))
cat(paste0("snapshotDirectory ", "/Volumes/kathi/Manu_TDP_CLIP/IGV/RBDm_6M_specific_peaks"), 
    fill = TRUE)
for (g in names(rg)) {
  for(i in names(rg[[g]])){
    snapshot_gene(rg[[g]][[i]], paste0(g, "_specific"))
  }
}
sink()
