## Duplicate Protein Screening
### This script loops through a directories fasta files (gunzipped) and identifies
### duplicated protein sequences (100% homology) so that they can be filtered
### in later analyses.
setwd("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/")
file_names <- list.files()
file_names <- file_names[grepl(".fasta.gz", file_names)]
length(file_names) == 219

duplicated_proteins <- data.frame()
for(active_file in file_names){
  fasta <- seqinr::read.fasta(file = active_file,
            seqtype = "AA", as.string = TRUE, set.attributes = FALSE, whole.header = TRUE)
  
  
  if(sum(duplicated(fasta)) > 0){
    duplicated_proteins <- rbind.data.frame(duplicated_proteins,
        data.frame(query_seqname = names(fasta)[duplicated(fasta)],
                   pathogen = gsub(".fasta.gz", "", active_file),
                   `Proteome ID` = gsub(".*_UP", "UP", gsub(".fasta.gz", "", active_file))))
  }
  
}
duplicated_proteins$viral_uniprot <- sapply(strsplit(duplicated_proteins$query_seqname, split="\\|"), function(x) paste(x[2]))

write.csv(duplicated_proteins, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/duplicated_proteins.csv",
          row.names = FALSE)
################################################## Human Duplicated Proteins
library(tidyverse)
library(data.table)
source("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Utilities.R")
setwd("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Human_Proteome/")
file_names <- "homo_sapiens_proteome.fasta.gz"
### Removed proteins 
dead_proteins_tr <- fread("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/delac_tr.txt",skip = 26)%>%
  dplyr::rename("uniprot" = 1)
dead_proteins_sp <- fread("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/delac_sp.txt",skip = 26) %>%
  dplyr::rename("uniprot" = 1)
dead_proteins <- rbind.data.frame(dead_proteins_tr, dead_proteins_sp)
rm(dead_proteins_sp, dead_proteins_tr)

fasta <- seqinr::read.fasta(file = file_names,
                            seqtype = "AA", as.string = TRUE, set.attributes = FALSE, whole.header = TRUE)
fasta_df <- data.frame(seq = unlist(fasta)) %>%
  rownames_to_column("query_seqname") %>%
  rowwise() %>%
  mutate(uniprot = id_split(query_seqname)) 

dead_proteins_small <- dead_proteins %>%
  filter(uniprot %in% fasta_df$uniprot)

removed <- fasta_df %>%
  filter(uniprot %in% dead_proteins_small$uniprot)

duplicated <- fasta_df %>%
  filter(!uniprot %in% removed$uniprot) %>%
  filter(duplicated(seq))

duplicated_or_removed_proteins <- unique(c(removed$uniprot, duplicated$uniprot))
write.csv(duplicated_or_removed_proteins, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/human_duplicated_or_removed_proteins.csv",
          row.names = FALSE)

## 12mer length filtering #####################################################
setwd("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/")
file_names <- list.files()
file_names <- file_names[grepl(".fasta.gz", file_names)]
length(file_names) == 219

short_proteins <- data.frame()
for(active_file in file_names){
  fasta <- seqinr::read.fasta(file = active_file,
                              seqtype = "AA", as.string = TRUE, set.attributes = FALSE, whole.header = TRUE)
  if(sum(nchar(fasta) < 12)){
    short_proteins <- rbind.data.frame(short_proteins,
                                            data.frame(query_seqname = names(fasta)[nchar(fasta) < 12],
                                                       ##sequence_length = nchar(fasta)[nchar(fasta) < 12],
                                                       pathogen = gsub(".fasta.gz", "", active_file),
                                                       `Proteome ID` = gsub(".*_UP", "UP", gsub(".fasta.gz", "", active_file))))
  }
}
short_proteins$viral_uniprot <- sapply(strsplit(short_proteins$query_seqname, split="\\|"), function(x) paste(x[2]))

write.csv(short_proteins, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/12mer_length_filtered_proteins.csv",
          row.names = FALSE)

## 18mer length filtering #####################################################
setwd("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/")
file_names <- list.files()
file_names <- file_names[grepl(".fasta.gz", file_names)]
length(file_names) == 219

short_proteins <- data.frame()
for(active_file in file_names){
  fasta <- seqinr::read.fasta(file = active_file,
                              seqtype = "AA", as.string = TRUE, set.attributes = FALSE, whole.header = TRUE)
  if(sum(nchar(fasta) < 18)){
    short_proteins <- rbind.data.frame(short_proteins,
                                       data.frame(query_seqname = names(fasta)[nchar(fasta) < 18],
                                                  ## sequence_length = nchar(fasta)[nchar(fasta) < 18],
                                                  pathogen = gsub(".fasta.gz", "", active_file),
                                                  `Proteome ID` = gsub(".*_UP", "UP", gsub(".fasta.gz", "", active_file))))
  }
}
short_proteins$viral_uniprot <- sapply(strsplit(short_proteins$query_seqname, split="\\|"), function(x) paste(x[2]))

write.csv(short_proteins, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/18mer_length_filtered_proteins.csv",
          row.names = FALSE)

## 8mer length filtering #####################################################
setwd("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/")
file_names <- list.files()
file_names <- file_names[grepl(".fasta.gz", file_names)]
length(file_names) == 219

short_proteins <- data.frame()
for(active_file in file_names){
  fasta <- seqinr::read.fasta(file = active_file,
                              seqtype = "AA", as.string = TRUE, set.attributes = FALSE, whole.header = TRUE)
  if(sum(nchar(fasta) < 8)){
    short_proteins <- rbind.data.frame(short_proteins,
                                       data.frame(query_seqname = names(fasta)[nchar(fasta) < 8],
                                                  ## sequence_length = nchar(fasta)[nchar(fasta) < 18],
                                                  pathogen = gsub(".fasta.gz", "", active_file),
                                                  `Proteome ID` = gsub(".*_UP", "UP", gsub(".fasta.gz", "", active_file))))
  }
}
short_proteins$viral_uniprot <- sapply(strsplit(short_proteins$query_seqname, split="\\|"), function(x) paste(x[2]))
## There are none
write.csv(short_proteins, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/Metadata/8mer_length_filtered_proteins.csv",
          row.names = FALSE)
