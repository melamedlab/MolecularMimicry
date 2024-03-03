### Human Proteome 12mer count

human_proteome <- seqinr::read.fasta(file = "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Human_Proteome/homo_sapiens_proteome.fasta.gz", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
Sys.time()
split <- data.frame()
for(active_prot in 1:length(human_proteome)){
  print(paste0("On protein ", active_prot, " out of 75776"))
  frags <- vector()
  if(nchar(unlist(human_proteome[active_prot])) < 12){next}
  for(active_position in 1:(nchar(unlist(human_proteome[active_prot]))- 11)){
    frags <- c(frags, paste(unlist(strsplit(unlist(human_proteome[active_prot]), split = ""))[active_position:(active_position+11)], collapse = ""))
    names(frags)[length(frags)] <- active_position
  }
  split <- rbind.data.frame(split, data.frame("mer" = frags, protein = names(human_proteome)[active_prot], position = names(frags)))
  browser()
}
## write.csv(split, "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/all_human_12mers.csv", row.names = F)
Sys.time()

################################################################################

human_proteome <- seqinr::read.fasta(file = "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/epitope_alignment/Human_Proteome/homo_sapiens_proteome.fasta.gz", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
mTEC_ready <- mTEC_ready %>%
  filter(!is.na(uniprot_gn_id) & uniprot_gn_id != "")
mTEC_ready_distinct <- mTEC_ready %>%
  distinct(uniprot_gn_id, in_mTEC)

pass1 <- mTEC_ready_distinct$uniprot_gn_id[mTEC_ready_distinct$in_mTEC == "non-mTEC"]
pass1 <- pass1[! pass1 %in% mTEC_ready_distinct$uniprot_gn_id[mTEC_ready_distinct$in_mTEC == "mTEC"]]
non_mTEC_cells <- human_proteome[sapply(strsplit(names(human_proteome), split="\\|"), function(x) paste(x[2])) %in% pass1]
mTEC_cells <- human_proteome[sapply(strsplit(names(human_proteome), split="\\|"), function(x) paste(x[2])) %in% mTEC_ready_distinct$uniprot_gn_id[mTEC_ready_distinct$in_mTEC == "mTEC"]]
sum(mTEC_ready_distinct$uniprot_gn_id[mTEC_ready_distinct$in_mTEC == "non-mTEC"] %in% mTEC_ready_distinct$uniprot_gn_id[mTEC_ready_distinct$in_mTEC == "mTEC"])

write.fasta(non_mTEC_cells, names = names(non_mTEC_cells), "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/mTEC_Human_12mer/non_mTEC.fasta.gz")
write.fasta(mTEC_cells, names = names(mTEC_cells), "/stor/work/Ehrlich_COVID19/Cole/Molecular_Mimicry/Analysis/mTEC_Human_12mer/mTEC.fasta.gz")
