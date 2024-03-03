### Permutation Practice ##############
### Molecular Mimicry #################
## Full Protein Scrmable, 100 times Epstein Barr Virus 
library(seqinr)
set.seed(123)

# viral.files <- list.files("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Complete_Final_Set_easy")
# viral.files <-viral.files[grepl("virus", viral.files)]

viral.files <- list.files("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/")
viral.files <-viral.files[grepl("virus|MERS|UP000109776", viral.files)]
files_ran <- paste0(gsub(".*\\/", "", list.dirs("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Permutation_Experiment")), ".fasta.gz")
### 07-06-2023 Let's run the last few missing viruses
viral.files <- viral.files[!viral.files %in% files_ran]

for(active_file in viral.files){
  active_viral_base <-
    seqinr::read.fasta(paste0("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Final_Cohort/", active_file),
               seqtype = "AA")
  
  dir.create(paste0("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Permutation_Experiment/",
                    gsub(".fasta.gz", "", active_file)), showWarnings = FALSE)
  
  for(round_i in 1:30){
    set.seed(round_i)
    active_viral_base_modded <- active_viral_base
    for(active_prot in 1:length(active_viral_base)){
      active_viral_base_modded[active_prot][[1]] <- sample(active_viral_base_modded[active_prot][[1]])
    }
    seqinr::write.fasta(active_viral_base_modded, names = names(active_viral_base_modded),
                file.out = paste0("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Permutation_Experiment/",
                                  gsub(".fasta.gz", "", active_file),
                                  "/Full_Scramble_",  gsub(".fasta.gz", "", active_file), "-",round_i ,".fasta.gz" ))
  }
  
  ## Reverse all the proteins
  active_viral_base_reverse <- active_viral_base
  for(active_prot in 1:length(active_viral_base_reverse)){
    active_viral_base_reverse[active_prot][[1]] <- rev(active_viral_base_modded[active_prot][[1]])
  }
  seqinr::write.fasta(active_viral_base_reverse, names = names(active_viral_base_reverse),
              file.out = paste0("/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Viral_Proteomes/Permutation_Experiment/Reverse/",
                                "Reversed_",  gsub(".fasta.gz", "", active_file), ".fasta.gz" ))
  
}


