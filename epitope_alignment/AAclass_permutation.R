## Script to scramble proteins by AA class
library(stringr)
library(seqinr)

viral.files <- list.files("./Viral_Proteomes/Complete_Final_Set_easy")
viral.files <-viral.files[grepl("virus", viral.files)]

hydrophobic <- rep("1", 5); names(hydrophobic) <- c("A", "V", "L", "I", "M")
hydrophilic <- rep("2", 4); names(hydrophilic) <- c("N", "Q", "A", "T")
positive <- rep("3", 3); names(positive) <- c("K", "R", "H")
negative <- rep("4", 2); names(negative) <- c("D", "E")
aromatic <- rep("5", 3); names(aromatic) <- c("F", "Y", "W")
## Essentially this forms the dictionary, note I'm using numbers because other lettes could complicate the pattern sub latter
replace_with <- c(hydrophobic, hydrophilic, positive, negative, aromatic)

## Loop through every fasta file and read it in, create a directory if need be for the 30 permutations
## Note for the 30 permutations, we keep the active_viral_base unedited so we don't keep rereading it in
for(active_file in viral.files){
  active_viral_base <-
    read.fasta(paste0("./Viral_Proteomes/Complete_Final_Set_easy/", active_file),
               seqtype = "AA")
  dir.create(paste0("./Viral_Proteomes/Permutation_Experiment_AAClass/",
                    gsub(".fasta.gz", "", active_file)), showWarnings = FALSE)
  ## Run the permutation 30 times
  for(round_i in 1:30){
    set.seed(round_i) ## Critical to set seed for reproducibility
    active_viral_base_modded <- active_viral_base ## Create a temporary list where we will be scrambling
    for(active_prot in 1:length(active_viral_base)){ ## Loop through the proteins in the proteome
      protein_sequence <- active_viral_base_modded[active_prot][[1]]
      scramble.seq <- unlist(strsplit(str_replace_all(protein_sequence, replace_with), split = "")) ## Convert to numeric key
      for(i in 1:5){ ## for 1:5 in the AA class naming used above
        pool <- unlist(strsplit(protein_sequence, split = ""))[scramble.seq == i] ## Calculate sampling pool for AA Class
        scramble.seq[scramble.seq == i] <- sample(pool, replace = F) ## Sample the class back in to replace key
      }
      active_viral_base_modded[active_prot][[1]] <- scramble.seq ## Write the scrambled protein back to the object
    }
    ## Save the object, indicating the permutation number in the file name
    write.fasta(active_viral_base_modded, names = names(active_viral_base_modded),
                file.out = paste0("./Viral_Proteomes/Permutation_Experiment_AAClass/",
                                  gsub(".fasta.gz", "", active_file),
                                  "/AAClass_Scramble_",  gsub(".fasta.gz", "", active_file), "-",round_i ,".fasta.gz" ))
  }
}
