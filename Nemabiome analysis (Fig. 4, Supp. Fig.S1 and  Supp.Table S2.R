#Nemabiome analysis (Fig. 4, Supp. Fig.S1 and  Supp.Table S2)

library("knitr")
library("BiocStyle")
require('ggplot2')
require('gridExtra')
require('dada2')
require('phyloseq')
require(dplyr)
library(data.table)
require(stringi)
require(stringr)
library(readr)
library(DECIPHER)
library(ShortRead)
library(Biostrings)
library(Hmisc)
theme_set(theme_bw())

set.seed(106)
path <- "~/Documents/CHIRON/WP3/Sainfoin/"
setwd(path)

###########----------------------------- Read in dada2 output with BAND_SIZE parameters modified

for(mxee in c(11,'2.05.0')){
  for(truncL in c(200, 217)){
    for(BS in c(-1,16,32)){
      
      ### Read in seqtab
      seqtab_nochim = read.table(file = paste0('./dada2_R_mxee',mxee,'_trunc',truncL,'_BS',BS,'/output.tsv'),
                                 header=T)
      
      dim(seqtab_nochim)
      #[1] 189 105
      
      colnames(seqtab_nochim) = sapply(stringr::str_split(colnames(seqtab_nochim),'_'),
                                       function(x) x[1])
      
      dna <- DNAStringSet(getSequences(seqtab_nochim$OTUID))
      seqtab_nochim$OTUID = paste0('ASV_',1:nrow(seqtab_nochim))
      rownames(seqtab_nochim) = seqtab_nochim$OTUID
      names(dna) = seqtab_nochim$OTUID
      seqtab_nochim$OTUID = NULL
      
      ## samp_dat
      samples = unique(colnames(seqtab_nochim))
      samp_data <-  data.frame(
        row.names = samples,
        sample = samples
      )
      
      for(taxmeth in c('idtaxa','assigntaxa')){
        
        if(taxmeth=='idtaxa'){
          train <- readDNAStringSet("~/Documents/CYATHOMIX/WP3/NMB/db/idtaxa_03022022.fasta") 
          tax <- read_tsv("~/Documents/CYATHOMIX/WP3/NMB/db/idtaxa_03022022.tax") 
          
          trainingSet <- LearnTaxa(train, names(train), tax)
          #dna <- DNAStringSet(getSequences(seqtab_nochim))
          
          ids <- IdTaxa(dna,
                        trainingSet,
                        strand = "both",
                        threshold = 50,
                        bootstraps = 100,
                        processors = NULL,
                        verbose = TRUE,
                        type = "extended")
          
          ranks <- c("root","domain", "phylum", "class", "order", "family", "subfamily", "genus", "species") # ranks of interest
          # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
          taxidITS <- t(sapply(ids, function(x) {
            m <- match(ranks, x$rank)
            taxa <- x$taxon[m]
            taxa[startsWith(taxa, "unclassified_")] <- NA
            taxa
          }))
          colnames(taxidITS) <- ranks; rownames(taxidITS) <- names(dna) #getSequences(seqtab_nochim)
          
          # Phyloseq obj
          assign(paste0('ps_',taxmeth,'_mxee',mxee,'_BS',BS,'_Tl',truncL), phyloseq(
            otu_table(seqtab_nochim, taxa_are_rows = TRUE),
            tax_table(taxidITS),
            sample_data(samp_data) #,dna
          ))
        }else{
          
          ###------ Use assignTaxonomy as in Poissant et al. 2020
          idsITS <- assignTaxonomy(
            dna, #paste0("./dada2_output_",bc,"_mxe",mxee,"_exported/dna-sequences.fasta"),
            "~/Documents/CYATHOMIX/WP3/NMB/db/idtaxa_03022022.fasta",
            minBoot = 0,
            tryRC = TRUE,
            outputBootstraps = TRUE,
            taxLevels = c("domain", "phylum", "class", "order", "family", "subfamily", "genus", "species"),
            multithread = FALSE,
            verbose = FALSE
          )
          taxidITS = idsITS
          ranks <- c("root","domain", "phylum", "class", "order", "family", "subfamily", "genus", "species") # ranks of interest
          colnames(taxidITS$tax) <- ranks
          rownames(taxidITS$tax) = names(rownames(taxidITS$tax))
          #taxidITS = taxidITS$tax[which(taxidITS$boot[,9]>=50),]
          
          # Phyloseq obj
          assign(paste0('ps_',taxmeth,'_mxee',mxee,'_BS',BS,'_Tl',truncL), phyloseq(
            otu_table(seqtab_nochim, taxa_are_rows = TRUE),
            tax_table(taxidITS$tax),
            sample_data(samp_data) #,dna
          ))
          assign(paste0('boot_',taxmeth,'_mxee',mxee,'_BS',BS,'_Tl',truncL),taxidITS$boot)
        }
        
      }
      rm(taxidITS,samples,seqtab_nochim)
    }
  }
}

taxGlomRank='species'

##############---- AssignTAXA
## mxee (Maximal Error Rate) 1,1
require(viridis)
`ps_assigntaxa_mxee11_BS-1_Tl200`
length(which(`boot_assigntaxa_mxee11_BS-1_Tl200`[,9]>50))
#[1] 111
length(which(`boot_assigntaxa_mxee11_BS-1_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS-1_Tl200`)
#[1] 0.8161765
`mock_assigntaxa_mxee11_BS-1_Tl200` = subset_samples(`ps_assigntaxa_mxee11_BS-1_Tl200`,
                                           sample %in% c('ITSM1','ITSM2','ITSM3',
                                                            'ITSM4','ITSM5'))
`mock_assigntaxa_mxee11_BS-1_Tl200.tf` = filter_taxa(`mock_assigntaxa_mxee11_BS-1_Tl200`,
                                                     function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_assigntaxa_mxee11_BS-1_Tl200.tf`,taxrank = taxGlomRank)

mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)

o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS-1_Tl200') )
rm(mock)

###
`ps_assigntaxa_mxee11_BS-1_Tl217`
length(which(`boot_assigntaxa_mxee11_BS-1_Tl217`[,9]>50))
#[1] 88
length(which(`boot_assigntaxa_mxee11_BS-1_Tl217`[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS-1_Tl217`)
#[1] 0.8301887

`mock_assigntaxa_mxee11_BS-1_Tl217` = subset_samples(`ps_assigntaxa_mxee11_BS-1_Tl217`,
                                                     sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                   'ITSM4','ITSM5'))
`mock_assigntaxa_mxee11_BS-1_Tl217.tf` = filter_taxa(`mock_assigntaxa_mxee11_BS-1_Tl217`,
                                                     function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_assigntaxa_mxee11_BS-1_Tl217.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)

o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS-1_Tl217') )
rm(mock)

###
ps_assigntaxa_mxee11_BS16_Tl200
length(which(`boot_assigntaxa_mxee11_BS16_Tl200`[,9]>50))
#[1] 112
length(which(`boot_assigntaxa_mxee11_BS16_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS16_Tl200`)
#[1] 0.8235294
mock_assigntaxa_mxee11_BS16_Tl200 = subset_samples(ps_assigntaxa_mxee11_BS16_Tl200,
                                                     sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                   'ITSM4','ITSM5'))
mock_assigntaxa_mxee11_BS16_Tl200.tf = filter_taxa(mock_assigntaxa_mxee11_BS16_Tl200,
                                                     function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee11_BS16_Tl200.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)


o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS16_Tl200') )
rm(mock)

##
ps_assigntaxa_mxee11_BS16_Tl217
length(which(`boot_assigntaxa_mxee11_BS16_Tl217`[,9]>50))
#[1] 89
length(which(`boot_assigntaxa_mxee11_BS16_Tl217`[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS16_Tl217`)
#[1] 0.8396226
mock_assigntaxa_mxee11_BS16_Tl217 = subset_samples(ps_assigntaxa_mxee11_BS16_Tl217,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_assigntaxa_mxee11_BS16_Tl217.tf = filter_taxa(mock_assigntaxa_mxee11_BS16_Tl217,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee11_BS16_Tl217.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS16_Tl217') )
rm(mock)

###
ps_assigntaxa_mxee11_BS32_Tl200
length(which(`boot_assigntaxa_mxee11_BS32_Tl200`[,9]>50))
#[1] 110
length(which(`boot_assigntaxa_mxee11_BS32_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS32_Tl200`)
#[1] 0.8088235
mock_assigntaxa_mxee11_BS32_Tl200 = subset_samples(ps_assigntaxa_mxee11_BS32_Tl200,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_assigntaxa_mxee11_BS32_Tl200.tf = filter_taxa(mock_assigntaxa_mxee11_BS32_Tl200,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee11_BS32_Tl200.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS32_Tl200') )
rm(mock)

###
ps_assigntaxa_mxee11_BS32_Tl217
length(which(boot_assigntaxa_mxee11_BS32_Tl217[,9]>50))
#[1] 89
length(which(boot_assigntaxa_mxee11_BS32_Tl217[,9]>50))/nrow(`boot_assigntaxa_mxee11_BS32_Tl217`)
#[1] 0.8396226
mock_assigntaxa_mxee11_BS32_Tl217 = subset_samples(ps_assigntaxa_mxee11_BS32_Tl217,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_assigntaxa_mxee11_BS32_Tl217.tf = filter_taxa(mock_assigntaxa_mxee11_BS32_Tl217,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee11_BS32_Tl217.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee11_BS32_Tl217') )
rm(mock)


##------- mxee 2,5
`ps_assigntaxa_mxee2.05.0_BS-1_Tl200`
length(which(`boot_assigntaxa_mxee2.05.0_BS-1_Tl200`[,9]>50))
#[1] 117
length(which(`boot_assigntaxa_mxee2.05.0_BS-1_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS-1_Tl200`)
#[1] 0.8731343
`mock_assigntaxa_mxee2.05.0BS-1_Tl200` = subset_samples(`ps_assigntaxa_mxee2.05.0BS-1_Tl200`,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
`mock_assigntaxa_mxee2.05.0BS-1_Tl200.tf` = filter_taxa(`mock_assigntaxa_mxee2.05.0BS-1_Tl200`,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_assigntaxa_mxee2.05.0BS-1_Tl200.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS-1_Tl200') )
rm(mock)

###
`ps_assigntaxa_mxee2.05.0_BS-1_Tl217`
length(which(`boot_assigntaxa_mxee2.05.0_BS-1_Tl217`[,9]>50))
#[1] 97
length(which(`boot_assigntaxa_mxee2.05.0_BS-1_Tl217`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS-1_Tl217`)
#[1] 0.8818182
`mock_assigntaxa_mxee2.05.0BS-1_Tl217` = subset_samples(`ps_assigntaxa_mxee2.05.0_BS-1_Tl217`,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
`mock_assigntaxa_mxee2.05.0BS-1_Tl217.tf` = filter_taxa(`mock_assigntaxa_mxee2.05.0BS-1_Tl217`,
                                                        function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_assigntaxa_mxee2.05.0BS-1_Tl217.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS-1_Tl217') )
rm(mock)

###
ps_assigntaxa_mxee2.05.0_BS16_Tl200
length(which(`boot_assigntaxa_mxee2.05.0_BS16_Tl200`[,9]>50))
#[1] 118
length(which(`boot_assigntaxa_mxee2.05.0_BS16_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS16_Tl200`)
#[1] 0.880597
`mock_assigntaxa_mxee2.05.0BS16_Tl217` = subset_samples(ps_assigntaxa_mxee2.05.0_BS16_Tl200,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
`mock_assigntaxa_mxee2.05.0BS16_Tl217.tf` = filter_taxa(`mock_assigntaxa_mxee2.05.0BS16_Tl217`,
                                                        function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_assigntaxa_mxee2.05.0BS16_Tl217.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS16_Tl200') )
rm(mock)

###
ps_assigntaxa_mxee2.05.0_BS16_Tl217
length(which(`boot_assigntaxa_mxee2.05.0_BS16_Tl217`[,9]>50))
#[1] 96
length(which(`boot_assigntaxa_mxee2.05.0_BS16_Tl217`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS16_Tl217`)
#[1] 0.8727273
mock_assigntaxa_mxee2.05.0BS16_Tl217 = subset_samples(ps_assigntaxa_mxee2.05.0_BS16_Tl217,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
mock_assigntaxa_mxee2.05.0BS16_Tl217.tf = filter_taxa(mock_assigntaxa_mxee2.05.0BS16_Tl217,
                                                        function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee2.05.0BS16_Tl217.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS16_Tl217') )
rm(mock)

###
ps_assigntaxa_mxee2.05.0_BS32_Tl200
length(which(`boot_assigntaxa_mxee2.05.0_BS32_Tl200`[,9]>50))
#[1] 118
length(which(`boot_assigntaxa_mxee2.05.0_BS32_Tl200`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS32_Tl200`)
#[1] 0.880597
mock_assigntaxa_mxee2.05.0BS32_Tl200 = subset_samples(ps_assigntaxa_mxee2.05.0_BS32_Tl200,
                                                      sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                    'ITSM4','ITSM5'))
mock_assigntaxa_mxee2.05.0BS32_Tl200.tf = filter_taxa(mock_assigntaxa_mxee2.05.0BS32_Tl200,
                                                      function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee2.05.0BS32_Tl200.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS32_Tl200') )
rm(mock)

###
ps_assigntaxa_mxee2.05.0_BS32_Tl217
length(which(`boot_assigntaxa_mxee2.05.0_BS32_Tl217`[,9]>50))
#[1] 96
length(which(`boot_assigntaxa_mxee2.05.0_BS32_Tl217`[,9]>50))/nrow(`boot_assigntaxa_mxee2.05.0_BS32_Tl217`)
#[1] 0.8727273
mock_assigntaxa_mxee2.05.0BS32_Tl217 = subset_samples(ps_assigntaxa_mxee2.05.0_BS32_Tl217,
                                           sample %in% c('ITSM1','ITSM2','ITSM3',
                                                            'ITSM4','ITSM5'))
mock_assigntaxa_mxee2.05.0BS32_Tl217.tf = filter_taxa(mock_assigntaxa_mxee2.05.0BS32_Tl217,
                                                      function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_assigntaxa_mxee2.05.0BS32_Tl217.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_assigntaxa_mxee2.05.0BS32_Tl200') )
rm(mock)

##############---- IdTAXA
## mxee 1,1
require(viridis)
`ps_idtaxa_mxee11_BS-1_Tl200`
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS-1_Tl200`)[,9]))))
nona
#[1] 102
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS-1_Tl200`)[,9])
#[1] 0.75
`mock_idtaxa_mxee11_BS-1_Tl200` = subset_samples(`ps_idtaxa_mxee11_BS-1_Tl200`,
                                                     sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                   'ITSM4','ITSM5'))
`mock_idtaxa_mxee11_BS-1_Tl200.tf` = filter_taxa(`mock_idtaxa_mxee11_BS-1_Tl200`,
                                                     function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_idtaxa_mxee11_BS-1_Tl200.tf`,taxrank = taxGlomRank)

mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)

o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS-1_Tl200') )
rm(mock)

###
`ps_idtaxa_mxee11_BS-1_Tl217`
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS-1_Tl217`)[,9]))))
nona
#[1] 80
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS-1_Tl217`)[,9])
#[1] 0.754717

`mock_idtaxa_mxee11_BS-1_Tl217` = subset_samples(`ps_idtaxa_mxee11_BS-1_Tl217`,
                                                     sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                   'ITSM4','ITSM5'))
`mock_idtaxa_mxee11_BS-1_Tl217.tf` = filter_taxa(`mock_idtaxa_mxee11_BS-1_Tl217`,
                                                     function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_idtaxa_mxee11_BS-1_Tl217.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)

o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS-1_Tl217') )
rm(mock)

###
ps_idtaxa_mxee11_BS16_Tl200
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS16_Tl200`)[,9]))))
nona
#[1] 105
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS16_Tl200`)[,9])
#[1] 0.7720588
mock_idtaxa_mxee11_BS16_Tl200 = subset_samples(ps_idtaxa_mxee11_BS16_Tl200,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_idtaxa_mxee11_BS16_Tl200.tf = filter_taxa(mock_idtaxa_mxee11_BS16_Tl200,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee11_BS16_Tl200.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)


o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS16_Tl200') )
rm(mock)

##
ps_idtaxa_mxee11_BS16_Tl217
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS16_Tl217`)[,9]))))
nona
#[1] 81
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS16_Tl217`)[,9])
#[1] 0.7641509

mock_idtaxa_mxee11_BS16_Tl217 = subset_samples(ps_idtaxa_mxee11_BS16_Tl217,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_idtaxa_mxee11_BS16_Tl217.tf = filter_taxa(mock_idtaxa_mxee11_BS16_Tl217,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee11_BS16_Tl217.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS16_Tl217') )
rm(mock)

###
ps_idtaxa_mxee11_BS32_Tl200
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS32_Tl200`)[,9]))))
nona
#[1] 103
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS32_Tl200`)[,9])
#[1] 0.7573529
mock_idtaxa_mxee11_BS32_Tl200 = subset_samples(ps_idtaxa_mxee11_BS32_Tl200,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_idtaxa_mxee11_BS32_Tl200.tf = filter_taxa(mock_idtaxa_mxee11_BS32_Tl200,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee11_BS32_Tl200.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS32_Tl200') )
rm(mock)

###
ps_idtaxa_mxee11_BS32_Tl217
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee11_BS32_Tl217`)[,9]))))
nona
#[1] 83
nona/nrow(tax_table(`ps_idtaxa_mxee11_BS32_Tl217`)[,9])
#[1] 0.7830189
mock_idtaxa_mxee11_BS32_Tl217 = subset_samples(ps_idtaxa_mxee11_BS32_Tl217,
                                                   sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                 'ITSM4','ITSM5'))
mock_idtaxa_mxee11_BS32_Tl217.tf = filter_taxa(mock_idtaxa_mxee11_BS32_Tl217,
                                                   function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee11_BS32_Tl217.tf,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee11_BS32_Tl217') )
rm(mock)


##------- mxee 2,5
`ps_idtaxa_mxee2.05.0_BS-1_Tl200`
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS-1_Tl200`)[,9]))))
nona
#[1] 102
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS-1_Tl200`)[,9])
#[1] 0.761194
`mock_idtaxa_mxee2.05.0BS-1_Tl200` = subset_samples(`ps_idtaxa_mxee2.05.0_BS-1_Tl200`,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
`mock_idtaxa_mxee2.05.0BS-1_Tl200.tf` = filter_taxa(`mock_idtaxa_mxee2.05.0BS-1_Tl200`,
                                                        function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_idtaxa_mxee2.05.0BS-1_Tl200.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS-1_Tl200') )
rm(mock)

###
`ps_idtaxa_mxee2.05.0_BS-1_Tl217`
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS-1_Tl217`)[,9]))))
nona
#[1] 92
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS-1_Tl217`)[,9])
#[1] 0.8363636
`mock_idtaxa_mxee2.05.0BS-1_Tl217` = subset_samples(`ps_idtaxa_mxee2.05.0_BS-1_Tl217`,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
`mock_idtaxa_mxee2.05.0BS-1_Tl217.tf` = filter_taxa(`mock_idtaxa_mxee2.05.0BS-1_Tl217`,
                                                        function(x) sum(x) > 50, TRUE)

mock = tax_glom(`mock_idtaxa_mxee2.05.0BS-1_Tl217.tf`,taxrank = taxGlomRank)

mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS-1_Tl217') )

#(## Check Cyathostomum %
`mock_idtaxa_mxee2.05.0BS-1_Tl217.tf` = filter_taxa(`mock_idtaxa_mxee2.05.0BS-1_Tl217`,
                                                    function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_idtaxa_mxee2.05.0BS-1_Tl217.tf`,taxrank = taxGlomRank)
mock = transform_sample_counts(mock,function(x) x/sum(x))
# OTU Table:          [5 taxa and 5 samples]
# taxa are rows
#             ITSM1      ITSM2     ITSM3     ITSM4      ITSM5
# ASV_3  0.03384718 0.05555556 0.1486340 0.0000000 0.00000000
# ASV_6  0.65750670 0.48653199 0.5740083 0.3768393 0.07227755
# ASV_13 0.11863271 0.21604938 0.2773576 0.6231607 0.92772245
# ASV_21 0.16253351 0.15937149 0.0000000 0.0000000 0.00000000
# ASV_58 0.02747989 0.08249158 0.0000000 0.0000000 0.00000000
# tax_table(mock)[,9]
# Taxonomy Table:     [5 taxa by 1 taxonomic ranks]:
#   species                 
# ASV_3  "Cylicocyclus_nassatus" 
# ASV_6  "Cyathostomum_pateratum"
# ASV_13 "Cyathostomum_catinatum"
# ASV_21 "Coronocyclus_labiatus" 
# ASV_58 "Cylicocyclus_insigne" )
rm(mock)

###
ps_idtaxa_mxee2.05.0_BS16_Tl200
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS16_Tl200`)[,9]))))
nona
#[1] 108
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS16_Tl200`)[,9])
#[1] 0.8059701
`mock_idtaxa_mxee2.05.0BS16_Tl217` = subset_samples(ps_idtaxa_mxee2.05.0_BS16_Tl200,
                                                        sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                      'ITSM4','ITSM5'))
`mock_idtaxa_mxee2.05.0BS16_Tl217.tf` = filter_taxa(`mock_idtaxa_mxee2.05.0BS16_Tl217`,
                                                        function(x) sum(x) > 50, TRUE)
mock = tax_glom(`mock_idtaxa_mxee2.05.0BS16_Tl217.tf`,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS16_Tl200') )
rm(mock)

###
ps_idtaxa_mxee2.05.0_BS16_Tl217
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS16_Tl217`)[,9]))))
nona
#[1] 90
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS16_Tl217`)[,9])
#[1] 0.8181818
mock_idtaxa_mxee2.05.0BS16_Tl217 = subset_samples(ps_idtaxa_mxee2.05.0_BS16_Tl217,
                                                      sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                    'ITSM4','ITSM5'))
mock_idtaxa_mxee2.05.0BS16_Tl217.tf = filter_taxa(mock_idtaxa_mxee2.05.0BS16_Tl217,
                                                      function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee2.05.0BS16_Tl217.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS16_Tl217') )
rm(mock)

###
ps_idtaxa_mxee2.05.0_BS32_Tl200
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS32_Tl200`)[,9]))))
nona
#[1] 106
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS32_Tl200`)[,9])
#[1] 0.7910448
mock_idtaxa_mxee2.05.0BS32_Tl200 = subset_samples(ps_idtaxa_mxee2.05.0_BS32_Tl200,
                                                      sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                    'ITSM4','ITSM5'))
mock_idtaxa_mxee2.05.0BS32_Tl200.tf = filter_taxa(mock_idtaxa_mxee2.05.0BS32_Tl200,
                                                      function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee2.05.0BS32_Tl200.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS32_Tl200') )
rm(mock)

###
ps_idtaxa_mxee2.05.0_BS32_Tl217
nona = length(which(!(is.na(tax_table(`ps_idtaxa_mxee2.05.0_BS32_Tl217`)[,9]))))
nona
#[1] 89
nona/nrow(tax_table(`ps_idtaxa_mxee2.05.0_BS32_Tl217`)[,9])
#[1] 0.8090909
mock_idtaxa_mxee2.05.0BS32_Tl217 = subset_samples(ps_idtaxa_mxee2.05.0_BS32_Tl217,
                                                      sample %in% c('ITSM1','ITSM2','ITSM3',
                                                                    'ITSM4','ITSM5'))
mock_idtaxa_mxee2.05.0BS32_Tl217.tf = filter_taxa(mock_idtaxa_mxee2.05.0BS32_Tl217,
                                                      function(x) sum(x) > 50, TRUE)
mock = tax_glom(mock_idtaxa_mxee2.05.0BS32_Tl217.tf ,taxrank = taxGlomRank)
mock <- filter_taxa(transform_sample_counts(mock, 
                                            function(x) ifelse(x>0,1,0)),
                    function(x) sum(x) > 0, TRUE)
o1 = order(otu_table(mock)[,1],decreasing =T)
print(pheatmap::pheatmap(otu_table(mock)[o1,],cluster_rows = F, 
                         cluster_cols = T,fontsize_col = 6,#cellheight = 14,
                         color = viridis_pal(option='D')(2),
                         legend_breaks = c(0,1),legend_labels = c('Absence','Presence'),
                         labels_row = tax_table(mock)[o1,ncol(tax_table(mock))],
                         border_color = "grey60",main = 'ps_idtaxa_mxee2.05.0BS32_Tl217') )
rm(mock)

rm(list=ls())

########################################################## SAINFOIN study - 
### paired mxee 2,5 - 217 bp
### pseudo-pooled, consensus chimera
### BS 32
#ps_idtaxa_mxee2.05.0_BS32_Tl217
BS = -1
truncL = 217
mxee = '2.05.0'
taxmeth='idtaxa'
### Sequencing stats
stats = read.table(file = paste0('./dada2_R_mxee',mxee,'_trunc',truncL,'_BS',BS,'/track.tsv'),header=T)

plot(stats$non.chimeric[order(stats$non.chimeric)])
abline(h = 50)

summary(stats$input[stats$non.chimeric>50])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13079   16474   17808   25924   28376   74679

summary(stats$non.chimeric[stats$non.chimeric>50])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2269    3894    4816    6768    8120   25434

### Import dada
seqtab_nochim = read.table(file = paste0('./dada2_R_mxee',mxee,'_trunc',truncL,'_BS',BS,'/output.tsv'),
                           header=T)

dim(seqtab_nochim)
#[1] 110 39

colnames(seqtab_nochim) = sapply(stringr::str_split(colnames(seqtab_nochim),'_'),function(x) x[1])

###--------- IDTAXA for txid

dna <- DNAStringSet(getSequences(seqtab_nochim$OTUID))
seqtab_nochim$OTUID = paste0('ASV_',1:nrow(seqtab_nochim))
rownames(seqtab_nochim) = seqtab_nochim$OTUID
names(dna) = seqtab_nochim$OTUID

seqtab_nochim$OTUID=NULL
## samp_dat
samples = unique(colnames(seqtab_nochim))
samp_data <-  data.frame(
  row.names = samples,
  sample = samples
)


train <- readDNAStringSet("~/Documents/CYATHOMIX/WP3/NMB/db/idtaxa_03022022.fasta") 
tax <- read_tsv("~/Documents/CYATHOMIX/WP3/NMB/db/idtaxa_03022022.tax") 
trainingSet <- LearnTaxa(train, names(train), tax)

ids <- IdTaxa(dna,
              trainingSet,
              strand = "both",
              threshold = 50,
              bootstraps = 100,
              processors = NULL,
              verbose = TRUE,
              type = "extended")

ranks <- c("root","domain", "phylum", "class", "order", "family", "subfamily", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxidITS <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxidITS) <- ranks; rownames(taxidITS) <- names(dna)

# Phyloseq obj
assign(paste0('ps_',taxmeth,'_mxee',mxee,'_BS',BS,'_Tl',truncL), phyloseq(
  otu_table(seqtab_nochim, taxa_are_rows = TRUE),
  tax_table(taxidITS),
  sample_data(samp_data)
))

`ps_idtaxa_mxee2.05.0_BS-1_Tl217`
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 110 taxa and 38 samples ]
# sample_data() Sample Data:       [ 38 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 110 taxa by 9 taxonomic ranks ]

### Remove rare counts (contaminants) 
### this is to remove cases where low counts match a single ASV
ps = filter_taxa(`ps_idtaxa_mxee2.05.0_BS-1_Tl217`,function(x) sum(x) > 50, TRUE)
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 95 taxa and 38 samples ]
# sample_data() Sample Data:       [ 38 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 95 taxa by 9 taxonomic ranks ]

### Discard 4 samples with low counts following contaminant removal
colSums(otu_table(ps))[which(colSums(otu_table(ps)) < 50)]
# ITSH2O ITSJ28.22 
#      0        12  

ps2 = prune_samples(sample_sums(ps) >= 50, ps)
ps2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 95 taxa and 36 samples ]
# sample_data() Sample Data:       [ 36 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 95 taxa by 9 taxonomic ranks ]

length(which(is.na(tax_table(ps2)[,ncol(tax_table(ps2))])))
#[1] 18
length(which(is.na(tax_table(ps2)[,ncol(tax_table(ps2))-1])))
#[1] 4

### Relative abundance of undetermined species
undet = which(is.na(tax_table(ps2)[,ncol(tax_table(ps2))-1]))
100*sum(otu_table(ps2)[undet,])/sum(otu_table(ps2))
#[1] 0.6073756

### Contribution of each genus to missing species
undet_gen = which(is.na(tax_table(ps2)[,ncol(tax_table(ps2))]) & !(is.na(tax_table(ps2)[,ncol(tax_table(ps2))-1])))
undet_gen_asv = data.frame(tax_table(ps2)[undet_gen,])
undet_cylicocyclus = undet_gen_asv[undet_gen_asv$genus=='Cylicocyclus',]
undet_cyathostomum = undet_gen_asv[undet_gen_asv$genus=='Cyathostomum',]
undet_cylicostephanus = undet_gen_asv[undet_gen_asv$genus=='Cylicostephanus',]

100*sum(otu_table(ps2)[match(rownames(undet_cylicocyclus),rownames(otu_table(ps2)))])/sum(otu_table(ps2))
#[1] 2.749217
100*sum(otu_table(ps2)[match(rownames(undet_cyathostomum),rownames(otu_table(ps2)))])/sum(otu_table(ps2))
#[1] 5.786095
100*sum(otu_table(ps2)[match(rownames(undet_cylicostephanus),rownames(otu_table(ps2)))])/sum(otu_table(ps2))
#[1] 0.09246246

###########------------------------------------------- Taxonomic + phylo agglomeration
# Taxonomic glom
# How many genera are present after filtering?
taxGlomRank = "species"
length(which(!is.na(get_taxa_unique(ps2, taxonomic.rank = taxGlomRank))))
## [1] 13
get_taxa_unique(ps2, taxonomic.rank = taxGlomRank)
# [1] "Cylicostephanus_minutus"       "Cylicocyclus_nassatus"         "Cyathostomum_pateratum"        NA                             
# [5] "Coronocyclus_coronatus"        "Cylicostephanus_longibursatus" "Cylicocyclus_ashworthi"        "Cyathostomum_catinatum"       
# [9] "Cylicocyclus_leptostomus"      "Coronocyclus_labiatus"         "Cylicostephanus_goldi"         "Cylicocyclus_insigne"         
# [13] "Cylicostephanus_calicatus"     "Cylicostephanus_bidentatus"   

pss = tax_glom(ps2,taxrank = taxGlomRank)
pss
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 13 taxa and 46 samples ]
# sample_data() Sample Data:       [ 46 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 13 taxa by 9 taxonomic ranks ]

length(which(is.na(tax_table(pss)[,9])))
#[1] 0
length(which(is.na(tax_table(pss)[,8])))
#[1] 0

rowSums(otu_table(pss))
pss.t = transform_sample_counts(pss,function(x) x / sum(x))
tot = data.frame(rel.abundance = rowSums(otu_table(pss.t)))
tot$species = tax_table(pss)[match(rownames(tot),rownames(tax_table(pss))),9]

tot[order(tot[,1],decreasing = T),]
#        rel.abundance                       species
# ASV_1    13.82566810       Cylicostephanus_minutus
# ASV_7     5.47439508        Cylicocyclus_ashworthi
# ASV_3     5.05943218         Cylicocyclus_nassatus
# ASV_10    3.21142270 Cylicostephanus_longibursatus
# ASV_13    2.78190253        Cyathostomum_catinatum
# ASV_6     2.52859764        Cyathostomum_pateratum
# ASV_9     1.15391125        Coronocyclus_coronatus
# ASV_21    0.75984728         Coronocyclus_labiatus
# ASV_25    0.38851636      Cylicocyclus_leptostomus
# ASV_58    0.29109349          Cylicocyclus_insigne
# ASV_41    0.28105594     Cylicostephanus_calicatus
# ASV_43    0.21118802         Cylicostephanus_goldi
# ASV_78    0.03296941    Cylicostephanus_bidentatus

## Define color palette
pal = qualpalr::qualpal(n=nrow(tax_table(pss)),
                        colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))
df = sample_data(pss.t)
df$sample.id = gsub('ITSM1','Mock-6sp-raw',df$sample) ## rename mock community 
df$sample.id = gsub('ITSM2','Mock-6sp-equimolar',df$sample.id) ## rename mock community 
df$sample.id = gsub('ITSM4','Mock-2sp-1/2',df$sample.id) ## rename mock community 
df$sample.id = gsub('ITSM5','Mock-2sp-3/4',df$sample.id) ## rename mock community 
df$sample.id = gsub('ITSJ7','D0',df$sample.id) ## recode from the start of the sainfoin regime
df$sample.id = gsub('ITSJ28','D21',df$sample.id) ## recode from the start of the sainfoin regime
df$sample.id = gsub('\\.','\\.H',df$sample.id) ## recode Horse name
sample_data(pss.t) = df

pss.t.f = subset_samples(pss.t,!(sample.id %in% c('ITSM3'))) ## remove contaminated mock

pdf(file = './supplementary_Figure1.pdf',width = 14, height = 8)
plot_bar(pss.t.f, x = 'sample.id', fill="species") +
  scale_fill_manual(values = pal$hex) +
  theme(legend.position = 'bottom',text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank())+
  guides(fill = guide_legend(nrow = 3))
dev.off()

###------------------------------------------------ Sainfoin experiment metadata
sainf.design = data.frame(Horse=c(1, 2, 3, 4, 5,
                                  6, 7,  8, 9, 10,
                                  11, 12,13, 14,15,
                                  16, 17, 19, 21,22),
                          GP = c('Sainfoin','Sainfoin','Control',
                                 'Sainfoin','Control','Sainfoin','Control',
                                 'Sainfoin','Sainfoin','Sainfoin','Sainfoin',
                                 'Control','Control','Control','Sainfoin',
                                 'Control','Control','Control','Sainfoin','Control'))

sample_data(pss)$Horse = factor(as.integer(sapply(stringr::str_split(sample_data(pss)$sample,'\\.'),
                                                  function(x) x[2])))
sample_data(pss)$Day = sapply(stringr::str_split(sample_data(pss)$sample,'\\.'),
                              function(x) gsub('ITSJ','',x[1]))
sample_data(pss)$Day = factor(sample_data(pss)$Day,levels=c(7,28)) ## coded from the start of the transition period; xp starts at 7day post-transition
sample_data(pss)$GP = sainf.design$GP[match(sample_data(pss)$Horse,sainf.design$Horse)]

###------------------------------------------------ Remove mock communities
psITSsainf = subset_samples(pss,GP %in% c('Sainfoin','Control') & Day %in% c(7,28))
psITSsainf
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 13 taxa and 31 samples ]
# sample_data() Sample Data:       [ 31 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 13 taxa by 9 taxonomic ranks ]

get_taxa_unique(psITSsainf, taxonomic.rank = taxGlomRank)
# [1] "Cylicostephanus_minutus"       "Cylicocyclus_nassatus"         "Coronocyclus_coronatus"       
# [4] "Cylicocyclus_leptostomus"      "Cyathostomum_pateratum"        "Cylicocyclus_ashworthi"       
# [7] "Cylicostephanus_longibursatus" "Cyathostomum_catinatum"        "Coronocyclus_labiatus"        
# [10] "Cylicocyclus_insigne"          "Cylicostephanus_goldi"         "Cylicostephanus_calicatus"    
# [13] "Cylicostephanus_bidentatus"   

table(sample_data(psITSsainf)$GP, sample_data(psITSsainf)$Day)
#          7 28
# Control  8  5
# Sainfoin 9  9

### Most abundant taxa on day 0 (7 days after the transition period)
ma0 = sort(100*rowSums(otu_table(psITSsainf)[,grep('J7',colnames(otu_table(psITSsainf)))]),
           decreasing = T)/nrow(otu_table(psITSsainf)[,grep('J7',colnames(otu_table(psITSsainf)))])
ma0 = data.frame(ma0)
ma0$sp = tax_table(pss)[match(rownames(ma0),rownames(tax_table(pss))),9]
ma0
#                ma0                       species
# ASV_1  602207.6923       Cylicostephanus_minutus
# ASV_7  120446.1538        Cylicocyclus_ashworthi
# ASV_3   99015.3846         Cylicocyclus_nassatus
# ASV_10  64300.0000 Cylicostephanus_longibursatus
# ASV_9   29576.9231        Coronocyclus_coronatus
# ASV_6   15269.2308        Cyathostomum_pateratum
# ASV_13  11400.0000        Cyathostomum_catinatum
# ASV_41   6115.3846     Cylicostephanus_calicatus
# ASV_25   5653.8462      Cylicocyclus_leptostomus
# ASV_21   4800.0000         Coronocyclus_labiatus
# ASV_43   4100.0000         Cylicostephanus_goldi
# ASV_58   1730.7692          Cylicocyclus_insigne
# ASV_78    923.0769    Cylicostephanus_bidentatus

###----- Alpha sainfoin to be tested within each day
alpha_div <- estimate_richness(psITSsainf, split = TRUE, measure = "Shannon")
alpha_div$sample <- rownames(alpha_div) %>%  as.factor()

alphaplot <- sample_data(psITSsainf) %>%
  unclass() %>%
  data.frame() %>%
  left_join(alpha_div, by = "sample") %>%
  reshape2::melt(measure.vars = "Shannon",
                 variable.name = "diversity_measure",
                 value.name = "alpha_diversity")
ggplot(alphaplot) +
  geom_boxplot(aes(x = GP, y = alpha_diversity, fill = GP), alpha = .4) +
  facet_wrap(~ Day) +
  scale_y_continuous(limits = c(0, 2.2), breaks = seq(0, 3, .2)) +
  labs(x = "Type", y = "Shannon Diversity", color = "Pipeline") +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 

ggplot(alphaplot) +
  geom_point(aes(x = Day, y = alpha_diversity, col = GP,group = Horse), alpha = .4) +
  geom_line(aes(x = Day, y = alpha_diversity, col = GP,group = Horse), alpha = .4) +
  scale_y_continuous(limits = c(0, 2.2), breaks = seq(0, 3, .2)) +
  labs(x = "Type", y = "Shannon Diversity", color = "Pipeline") +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 

t.test(alpha_diversity ~ GP, 
       data = alphaplot[alphaplot$Day==7,])
# data:  alpha_diversity by GP
# t = -0.74824, df = 11.101, p-value = 0.4699
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7123224  0.3505815
# sample estimates:
#   mean in group Control mean in group Sainfoin 
# 1.090634               1.271504

t.test(alpha_diversity ~ GP, 
       data = alphaplot[alphaplot$Day==28,])
# data:  alpha_diversity by GP
# t = -0.66654, df = 6.9188, p-value = 0.5267
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.7237944  0.4060651
# sample estimates:
#   mean in group Control mean in group Sainfoin 
# 1.318642               1.477507 

### NMDS Sainfoin
pslog <- transform_sample_counts(psITSsainf, function(x) log(1 + x))
psbin  = transform_sample_counts(psITSsainf, function(x) ifelse(x>0,1,0))

######------- PCoA Bray
out.bra.log <- ordinate(pslog, method = "PCoA", distance = "bray")
evals <- out.bra.log$values$Eigenvalues

## Method
p = plot_ordination(pslog, out.bra.log,
                    color = "Day", shape  = "GP") +
  #labs(col = "Pipeline") +
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) +
  ggtitle(paste0('PCoA - Bray - ITS')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16))
print(p)

######------- PCoA Jaccard
out.jac.log <- ordinate(psbin, method = "PCoA", distance = "jaccard")
evals <- out.jac.log$values$Eigenvalues

## Method
p = plot_ordination(psbin, out.jac.log,
                    color = "Day", shape  = "GP") +
  #labs(col = "Pipeline") +
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) +
  ggtitle(paste0('PCoA - Jaccard - ITS')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16))
print(p)

######------- PERMANOVA - Bray
library(vegan)
meta.sainf = data.frame(sample_data(pslog))
pslog.sp = as.matrix(otu_table(pslog))

pslog.sp7 = as.matrix(otu_table(pslog))[,meta.sainf$Day==7]
pslog.sp28 = as.matrix(otu_table(pslog))[,meta.sainf$Day==28]

## Day 7
ps7 = subset_samples(pslog,Day==7)
out.bra.log7 <- ordinate(ps7, method = "PCoA", distance = "bray")
evals7 <- out.bra.log7$values$Eigenvalues

## Method
p7 = plot_ordination(ps7, out.bra.log7, 
                     color = "GP", shape  = "GP") +
  #labs(col = "Pipeline") + 
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) + 
  ggtitle(paste0('PCoA - Bray - Day7')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 
p7

ps7.fr = filter_taxa(ps7, function(x) sum(x) > 0, TRUE)
ps7.fr
adonis(vegdist(t(as.matrix(otu_table(ps7.fr))), method="bray") ~ GP,
       data = meta.sainf[meta.sainf$Day==7,])
# GP         1   0.04080 0.040801 0.93678 0.05878  0.472
# Residuals 15   0.65331 0.043554         0.94122       
# Total     16   0.69411                  1.00000 

## Day 28
ps28 = subset_samples(pslog,Day==28)
out.bra.log28 <- ordinate(ps28, method = "PCoA", distance = "bray")
evals28 <- out.bra.log28$values$Eigenvalues

## Method
p28 = plot_ordination(ps28, out.bra.log28, 
                      color = "GP", shape  = "GP") +
  #labs(col = "Pipeline") + 
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) + 
  ggtitle(paste0('PCoA - Bray - Day 28')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 
p28

ps28.fr = filter_taxa(ps28, function(x) sum(x) > 0, TRUE)
ps28.fr
adonis(vegdist(t(as.matrix(otu_table(ps28.fr))), method="bray") ~ GP,
       data = meta.sainf[meta.sainf$Day==28,])
# GP         1   0.04156 0.041561  1.1302 0.08608  0.331
# Residuals 12   0.44126 0.036772         0.91392       
# Total     13   0.48282                  1.00000

source('~/Documents/Scripts/multiplot.R')
multiplot(p7,p28,cols=2)

######------- PERMANOVA - Jaccard
meta.sainf = data.frame(sample_data(psbin))
psbin.sp = as.matrix(otu_table(psbin))

psbin.sp7 = as.matrix(otu_table(psbin))[,meta.sainf$Day==7]
psbin.sp28 = as.matrix(otu_table(psbin))[,meta.sainf$Day==28]

## Day 7
ps7 = subset_samples(psbin,Day==7)
out.bra.log7 <- ordinate(ps7, method = "PCoA", distance = "jaccard")
evals7 <- out.bra.log7$values$Eigenvalues

## Method
p7 = plot_ordination(ps7, out.bra.log7, 
                     color = "GP", shape  = "GP") +
  #labs(col = "Pipeline") + 
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) + 
  ggtitle(paste0('PCoA - jaccard - Day7')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 
p7

adonis(vegdist(t(as.matrix(otu_table(ps7))), method="jaccard") ~ GP,
       data = meta.sainf[meta.sainf$Day==7,])
# GP         1   0.09747 0.097468  1.5277 0.09244  0.204
# Residuals 15   0.95698 0.063799         0.90756       
# Total     16   1.05445                  1.00000

## Day 28
ps28 = subset_samples(psbin,Day==28)
out.bra.log28 <- ordinate(ps28, method = "PCoA", distance = "jaccard")
evals28 <- out.bra.log28$values$Eigenvalues

## Method
p28 = plot_ordination(ps28, out.bra.log28, 
                      color = "GP", shape  = "GP") +
  #labs(col = "Pipeline") + 
  geom_point(size = 6, alpha = .4) +
  #coord_fixed(sqrt(evals[2] / evals[1])) + 
  ggtitle(paste0('PCoA - jaccard - Day 28')) + theme(legend.position = 'bottom') +
  theme(legend.position = 'bottom', text = element_text(size = 16)) 
p28

adonis(vegdist(t(as.matrix(otu_table(ps28))), method="jaccard") ~ GP,
       data = meta.sainf[meta.sainf$Day==28,])
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# GP         1   0.01617 0.016174 0.37128 0.03001   0.72
# Residuals 12   0.52276 0.043563         0.96999       
# Total     13   0.53893                  1.00000     

multiplot(p7,p28,cols=2)

#########----- Samples with complete data throughout the experiment

#####--- Differential abundance in species using horses present at day7 and day28

Hok = names(table(sample_data(psITSsainf)$Horse))[table(sample_data(psITSsainf)$Horse)==2]

pss2 = subset_samples(psITSsainf,Horse %in% Hok)
sample_data(pss2)$Day = factor(sample_data(pss2)$Day,levels=c(7,28))
sample_data(pss2)$Group = paste0(sample_data(pss2)$Day,'-',
                                 sample_data(pss2)$Horse)
sample_data(pss2)$Horse = factor(sample_data(pss2)$Horse)
pss2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 13 taxa and 24 samples ]
# sample_data() Sample Data:       [ 24 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 13 taxa by 9 taxonomic ranks ]

table(sample_data(pss2)$GP, sample_data(pss2)$Day)
#          7 28
# Control  4  4
# Sainfoin 8  8

pss2.t = transform_sample_counts(pss2,function(x) x/sum(x))


### Check species prevalence
dfTraj$bin = ifelse(dfTraj$count>0,1,0)

dfprev = data.frame(table(dfTraj$bin,dfTraj$Group, dfTraj$species,dfTraj$Day))

ggplot(dfprev,aes(x = Var3, y = Freq, fill = Var1)) + 
  geom_bar(stat='identity') + coord_flip() + ylab('') +
  facet_wrap(~ Var2 + Var4, scales = 'free_y')

dfpres = dfprev[dfprev$Var1==1,]

dfpres$Freq[dfpres$Var2=='Control']=dfpres$Freq[dfpres$Var2=='Control']/4
dfpres$Freq[dfpres$Var2=='Sainfoin']=dfpres$Freq[dfpres$Var2=='Control']/8

aggregate(Freq ~ Var2 + Var4, FUN = summary, data = dfpres)
#       Var2 Var4 Freq.Min. Freq.1st Qu. Freq.Median Freq.Mean Freq.3rd Qu. Freq.Max.
# 1  Control    7 0.2500000    0.7500000   1.0000000 0.8461538    1.0000000 1.0000000
# 2 Sainfoin    7 0.0312500    0.0937500   0.1250000 0.1057692    0.1250000 0.1250000
# 3  Control   28 0.2500000    0.5000000   1.0000000 0.8076923    1.0000000 1.0000000
# 4 Sainfoin   28 0.0312500    0.0625000   0.1250000 0.1009615    0.1250000 0.1250000

## Species present everywhere (group and day)

### Check species abundance
pss2.toplot = pss2.t
sample_data(pss2.toplot)$Day = as.integer(as.character(sample_data(pss2.toplot)$Day))
sample_data(pss2.toplot)$Day[sample_data(pss2.toplot)$Day==7] = 0 ##relevel relative to the beginning of the sainfoin regime
sample_data(pss2.toplot)$Day[sample_data(pss2.toplot)$Day==28] = 21  ##relevel relative to the beginning of the sainfoin regime
sample_data(pss2.toplot)$Day = factor(sample_data(pss2.toplot)$Day,levels=c(0,21))
tax_table(pss2.toplot)[,9] = gsub('_',' ',tax_table(pss2.toplot)[,9])

p1 = plot_bar(pss2.toplot, x = 'Horse', fill="species") +
  scale_fill_manual(values = pal$hex) +
  facet_wrap(~ paste0(GP, ' - Day ', Day) , ncol = 2, scales = 'free_x') +
  theme_classic() + 
  theme(legend.position='bottom',strip.background = element_blank(), 
        text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14), legend.title = element_blank()) + 
  guides(fill = guide_legend(nrow = 5))


dfTraj0 = reshape2::melt(as.matrix(otu_table(pss2)))
colnames(dfTraj0) = c('otu','sample.id','count')
dfTraj0$Group = sample_data(pss2)$GP[match(dfTraj0$sample.id,sample_data(pss2)$sample)]
dfTraj0$Id = factor(sapply(str_split(dfTraj0$sample.id,'\\.'), function(x) x[2]))
dfTraj0$Day = factor(gsub('ITSJ','',sapply(str_split(dfTraj0$sample.id,'\\.'), function(x) x[1])))
dfTraj0$Day = factor(dfTraj0$Day, levels = c(7, 28))
idx = match(dfTraj0$otu,rownames(tax_table(pss2)))
dfTraj0$species = tax_table(pss2)[idx,9]
dfTraj0$species = factor(dfTraj0$species)

aggregate(count ~ species, 
          FUN = sum, data = dfTraj0)
#                          species count
# 1         Coronocyclus_coronatus  4762
# 2          Coronocyclus_labiatus  2259
# 3         Cyathostomum_catinatum   800
# 4         Cyathostomum_pateratum  2969
# 5         Cylicocyclus_ashworthi 24050
# 6           Cylicocyclus_insigne   576
# 7       Cylicocyclus_leptostomus  2289
# 8          Cylicocyclus_nassatus 17400
# 9     Cylicostephanus_bidentatus   142
# 10     Cylicostephanus_calicatus   781
# 11         Cylicostephanus_goldi   926
# 12 Cylicostephanus_longibursatus 14406
# 13       Cylicostephanus_minutus 58609

####------ Retain top 4 species
abu.thr = 10000 
psCountTraj.f = metagMisc::phyloseq_filter_prevalence(pss2,
                                                      prev.trh = .3, ## 4 horses at least
                                                      abund.trh = abu.thr, abund.type = "total",
                                                      threshold_condition = "AND")  # 13 taxa
psCountTraj.f
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4 taxa and 24 samples ]
# sample_data() Sample Data:       [ 24 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 4 taxa by 9 taxonomic ranks ]

dfTraj = reshape2::melt(as.matrix(otu_table(psCountTraj.f)))
colnames(dfTraj) = c('otu','sample.id','count')
dfTraj$Group = sample_data(psCountTraj.f)$GP[match(dfTraj$sample.id,sample_data(psCountTraj.f)$sample)]

idx = match(dfTraj$otu,rownames(tax_table(psCountTraj.f)))
dfTraj$species = tax_table(psCountTraj.f)[idx,9]
dfTraj$species = factor(dfTraj$species)
dfTraj$Group = factor(dfTraj$Group)
dfTraj$Id = factor(sapply(str_split(dfTraj$sample.id,'\\.'), function(x) x[2]))
dfTraj$Day = factor(gsub('ITSJ','',sapply(str_split(dfTraj$sample.id,'\\.'), function(x) x[1])))
dfTraj$Day = factor(dfTraj$Day, levels = c(7, 28))

plot(density(dfTraj$count))
shapiro.test(dfTraj$count)
shapiro.test(sqrt(sqrt(dfTraj$count)))
# data:  sqrt(sqrt(dfTraj$count))
# W = 0.9618, p-value = 4.535e-05
# dfTrajd = merge(dfTraj21,dfTraj0,by=c('otu','Group','species','Id'))
# dfTrajd$abudif = dfTrajd$count.x - dfTrajd$count.y
# dfTrajd = dfTrajd[dfTrajd$abudif> -10000,]

mTraj = nlme::lme(sqrt(sqrt(count)) ~ species*Group*Day ,
                  random=~ 1|Id,
                  data = dfTraj)
#summary(mTraj)
printCoefmat(summary(mTraj)$tTable)
#                                                              Value Std.Error        DF t-value   p-value    
# (Intercept)                                               5.804075  0.898493 70.000000  6.4598 1.188e-08 ***
# speciesCylicocyclus_nassatus                             -0.814961  1.270661 70.000000 -0.6414   0.52338    
# speciesCylicostephanus_longibursatus                     -2.747597  1.270661 70.000000 -2.1623   0.03401 *  
# speciesCylicostephanus_minutus                            1.502024  1.270661 70.000000  1.1821   0.24117    
# GroupSainfoin                                            -1.321021  1.100425 10.000000 -1.2005   0.25762    
# Day28                                                    -0.959106  1.270661 70.000000 -0.7548   0.45290    
# speciesCylicocyclus_nassatus:GroupSainfoin                1.131211  1.556236 70.000000  0.7269   0.46972    
# speciesCylicostephanus_longibursatus:GroupSainfoin        1.913151  1.556236 70.000000  1.2293   0.22306    
# speciesCylicostephanus_minutus:GroupSainfoin              0.018768  1.556236 70.000000  0.0121   0.99041    
# speciesCylicocyclus_nassatus:Day28                        0.897798  1.796986 70.000000  0.4996   0.61891    
# speciesCylicostephanus_longibursatus:Day28                2.082428  1.796986 70.000000  1.1588   0.25046    
# speciesCylicostephanus_minutus:Day28                     -0.562425  1.796986 70.000000 -0.3130   0.75522    
# GroupSainfoin:Day28                                       0.687347  1.556236 70.000000  0.4417   0.66009    
# speciesCylicocyclus_nassatus:GroupSainfoin:Day28         -0.586206  2.200850 70.000000 -0.2664   0.79075    
# speciesCylicostephanus_longibursatus:GroupSainfoin:Day28 -0.087516  2.200850 70.000000 -0.0398   0.96839    
# speciesCylicostephanus_minutus:GroupSainfoin:Day28        1.173555  2.200850 70.000000  0.5332   0.59556 

dfTraj$species = gsub('_',' ',dfTraj$species)
dfTraj$Day = as.integer(as.character(dfTraj$Day))
dfTraj$Day[dfTraj$Day==7] = 0 ##relevel relative to the beginning of the sainfoin regime
dfTraj$Day[dfTraj$Day==28] = 21  ##relevel relative to the beginning of the sainfoin regime
dfTraj$Day = factor(dfTraj$Day,levels=c(0,21))

p2 = ggplot(dfTraj,aes(x = Day, y = sqrt(sqrt(count)), group = paste0(Day,Group), fill = Group)) +
  #geom_point(size = 3,alpha = .6,position = position_dodge(width = .5))  + 
  geom_boxplot() +
  theme_classic() +
  theme(legend.position='bottom',strip.background = element_blank(), text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  ylab('4th-root transformed counts') +
  scale_fill_manual(values = c('#bdbdbd','#016450')) +
  facet_wrap(~ species, scales='free') 

pdf(file = './Figure_nemabiome.pdf',width = 14, height = 8)
multiplot(p1,p2,cols=2)
dev.off()
