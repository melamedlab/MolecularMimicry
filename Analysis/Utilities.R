## Utilities File
## Contains functions and color palettes

DNA_Diff <- function(DNA1, DNA2){
  DNA1 <- strsplit(DNA1, split = "")
  DNA2 <- strsplit(DNA2, split = "")
  DNAcon <- rbind(unlist(DNA1), unlist(DNA2))
  DNAdiff <- DNAcon[,DNAcon[1,] != DNAcon[2,]]
  ncol(as.data.frame(DNAdiff))
}

run_identifier <- function(x){
  runs <- split(seq_along(x), cumsum(c(0, diff(x) > 1)))
  position <- setNames(unlist(runs, use.names=F),rep(names(runs), lengths(runs)))
  out <- x[position]
  names(out) <- names(position)
  return(names(out))
}

dual_run_identifier <- function(x){
  runs <- split(seq_along(x), cumsum(c(0, diff(x) > 1)))
  position <- setNames(unlist(runs, use.names=F),rep(names(runs), lengths(runs)))
  out <- x[position]
  names(out) <- names(position)
  return(names(out))
}

target_merge <- function(x){
  target_seq <- paste0(x[1], paste(sapply(strsplit(x[-1], split = ""),tail,1), collapse = ""))
  return(target_seq)
}

mismatch_calling <- function(x, y){
  x <- strsplit(x, split = "")
  y <- strsplit(y, split = "")
  xycon <- rbind(unlist(x), unlist(y))
  xyDiff <- xycon[,xycon[1,] != xycon[2,]]
  ncol(as.data.frame(xyDiff))
}
mismatch_calling_fast <- function(x,y){
  return(sum(unlist(strsplit(x, split = "")) != unlist(strsplit(y, split = ""))))
}


id_split <- function(x){
  target_seq <- sapply(strsplit(x, split="\\|"), function(x) paste(x[2]))
  return(target_seq)
}
id_split_fast <- function(x){return(strsplit(x, split="\\|")[[1]][2])}

fast_ec_screen <- function(ec_number.human, ec_number.viral,
                           grepl_ec_on = TRUE){
  tmp.viral <- unlist(strsplit(ec_number.viral, ";"))
  tmp.viral <- tmp.viral[!is.na(tmp.viral)]
  if(length(tmp.viral) == 0){
    return(FALSE)
  }
  tmp.human <- unlist(strsplit(ec_number.human, ";"))
  tmp.human <- tmp.human[!is.na(tmp.human)]
  if(length(tmp.human) == 0){
    return(FALSE)
  }
  if(grepl_ec_on){
    tmp.human.pattern <- paste(paste0("^", gsub(".-.*$", "", tmp.human)), collapse = "|")
    tmp.human.pattern <- gsub("^\\^$", "No pattern", tmp.human.pattern)
    ## No human enzyme should have two or more ec numbers,
    ## but this just defensively allows for multiple
    return((sum(grepl(paste0(tmp.human.pattern, collapse = "|"), tmp.viral)) > 0 ))
  } else {
    return((sum(tmp.human %in% tmp.viral) > 0))
  }
}

## fast_ec_screen("3.1.-.-; 3.6.4.-", "2.7.7.48; 3.4.22.28; 3.6.4.13")

Chronic_Acute_palette <- c("Chronic" = "#BF5700", "Acute" = "#1975EE")

## The following is from 
### library(rcartocolor)
### nColor <- 12
### scales::show_col(carto_pal(nColor, "Safe"))
Family_palette_v2 <- c("Arenaviridae" = "#888888", "Phenuiviridae" = "#332288",
                       "Coronaviridae" = "#CC6677",
                       "Flaviviridae" = "#88CCEE", "Herpesviridae" = "#AA4499",
                       "Paramyxoviridae" = "#117733", "Picornaviridae" = "#999933",
                       "Polyomaviridae" = "#DDCC77",
                       "Poxviridae" = "#44AA99", "Retroviridae" = "#6699CC",
                       "Rhabdoviridae" = "#661100",
                       "Togaviridae" = "#882255")

Baltimore_classification_colors <- c("circular ssRNA-" = "#E69F00", "dsDNA" = "#56B4E9","dsDNA-RT" = "#009E73","dsRNA" = "#F0E442","ssDNA" = "#0072B2","ssRNA-" = "#D55E00","ssRNA+" = "#CC79A7","ssRNA+-RT" = "#999999")

protein_compare <- function(mimicry_df, fasta, protein_ids, kmer_length,
                            feature_name, star_heights = c(9,15,25, 25),
                            bar_height = 27, virus_name){
  library(scales)
  library(ggpubr)
  mimicry_df$viral_uniprot <- sapply(strsplit(mimicry_df$query_seqname, split="\\|"), function(x) paste(x[2]))
  mimicry_df$uniprot <- sapply(strsplit(mimicry_df$target_seqname, split="\\|"), function(x) paste(x[2]))
  PLength <- as.data.frame(nchar(fasta))
  
  data <- mimicry_df %>%
    filter(!viral_uniprot %in% duplicated_proteins$viral_uniprot) %>%
    filter(!uniprot %in% human_duplicated_proteins$uniprot) %>%
    distinct(viral_uniprot, query_seqname, query_start, mismatches) %>%
    group_by(viral_uniprot, query_seqname, query_start) %>%
    slice_min(mismatches) %>%
    group_by(viral_uniprot, query_seqname, mismatches) %>%
    dplyr::summarize(count = n()) %>%
    full_join(PLength %>% rownames_to_column("query_seqname")) 
  data$viral_uniprot <- sapply(strsplit(data$query_seqname, split="\\|"), function(x) paste(x[2]))
  
  data <- data %>%
    mutate(feature = ifelse(viral_uniprot %in% protein_ids,
                            feature_name[1], feature_name[2])) %>%
    dplyr::rename("protein_length" = "nchar(fasta)") %>%
    mutate(count_corr = count / (protein_length -(kmer_length-1)) * 100) %>%
    pivot_wider(names_from = mismatches, values_from = count_corr,
                id_cols = -count) %>% 
    mutate_at(c("0", "1", "2", "3"), ~replace_na(.,0))%>%
    rowwise() %>%
    mutate(`3` = sum(`0`, `1`, `2`, `3`,na.rm = T),
           `2` = sum(`0`, `1`, `2`,na.rm = T),
           `1` = sum(`0`, `1`,na.rm = T)) %>%
    ungroup()
  
  testing_data <- data %>%
    pivot_longer(cols = c("0", "1", "2", "3"), names_to = "mismatches",
                 values_to = "count_corr")
  
  wt <- compare_means(data = testing_data, formula = count_corr ~ feature,
                      group.by = c("mismatches"), method = "wilcox.test") %>%
    mutate(p.adj.star = ifelse(p.adj <= 0.001, "*\n*\n*", ifelse(p.adj <= 0.01, "*\n*",
                                                                 ifelse(p.adj<=0.05, "*", NA))))
  
  mean_measures <- data %>%
    pivot_longer(cols = c("0", "1", "2", "3"), names_to = "mismatches",
                 values_to = "count_corr") %>%
    mutate(mismatches =ifelse(mismatches != "0",  paste0("≤", mismatches), mismatches)) %>%
    mutate(mismatches = factor(mismatches, c("0", "≤1", "≤2", "≤3"))) %>%
    group_by(feature, mismatches) %>%
    dplyr::summarize(mean_value = mean(count_corr), count = n(), sd = sd(count_corr)) %>%
    mutate(value_stderr = sd/sqrt(count))
 
  gg2 <- data %>%
    pivot_longer(cols = c("0", "1", "2", "3"), names_to = "mismatches",
                 values_to = "count_corr") %>%
    mutate(mismatches =ifelse(mismatches != "0",  paste0("≤", mismatches), mismatches)) %>%
    mutate(mismatches = factor(mismatches, c("0", "≤1", "≤2", "≤3"))) %>%
    mutate(count_corr = ifelse(count_corr == 0, 0.1, count_corr)) %>%
    arrange(feature, mismatches) %>% #viral_uniprot
    ggplot(aes(x=as.numeric(interaction(feature,mismatches)), y=count_corr, color = `feature`, group = viral_uniprot)) +
    geom_boxplot(inherit.aes = FALSE, aes(x=as.numeric(interaction(feature,mismatches)),
                                          y=count_corr, color = `feature`, group = as.numeric(interaction(feature,mismatches))),
                 outlier.size = 0, outlier.shape = NA, alpha = 0.3, width = 0.4,position=position_dodge(1)) +
    geom_point(inherit.aes = FALSE,
               aes(x=as.numeric(interaction(feature,mismatches)), y=count_corr, color = `feature`),
               position = position_jitter(width = 0.1, seed = 3922), size = 0.4) +
    geom_line(alpha = 0.05, position = position_jitter(width = 0.1, seed = 3922)) +
    theme_classic() +
    xlab("Number of mismatches") +
    scale_color_manual(values =c("#B77737","#56709E"), name = "Expression") +
    ylab(paste0("Percent ", kmer_length, "mer mimics (%)")) +
    annotate("text", x= c(1.5,3.5,5.5,7.5), y= star_heights,
             label = gsub("\n", "", wt$p.adj.star), size = 6, lineheight = 0.5)+
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) round(10^x, 3))) +
    scale_x_continuous(breaks=c(1.5, 3.5, 5.5, 7.5), labels = c("0", "≤1", "≤2", "≤3"))
  
  return(c(list(gg2), list(wt)))
}


