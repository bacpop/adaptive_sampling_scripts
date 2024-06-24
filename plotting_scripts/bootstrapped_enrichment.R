library(ggplot2)
library(data.table)
library(purrr)
library(ggpubr)
library(ggsci)
library(stringr)
library(dplyr)
library(ggbeeswarm)
library(hash)
library("plotrix")
library(scales)
library(ggh4x)
library(svglite)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

experiment <- "Mixed_V14"

if (experiment == "V14_CPS_graph_remclust2")
{
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V14_minimap2_v_graph_remclust2/"
} else if (experiment == "V14_CPS_graph_fulldb")
{
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V14_minimap2_v_graph_fulldb/"
}else if (experiment == "V12_CPS_graph")
{
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V12_minimap2_v_graph/"
} else if (experiment == "V12_WGS")
{
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V12_WGS_size_selection/"
} else if (experiment == "V12_CPS")
{
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V12_CPS_size_selection/"
} else if (experiment == "V12_WGS_v_CPS"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V12_WGS_v_CPS/"
} else if (experiment == "Minimap2_partialdb_fulldb"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/Minimap2_partialdb_fulldb/"
} else if (experiment == "Graph_partialdb_fulldb"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/Graph_partialdb_fulldb/"
} else if (experiment == "V12_v_V14"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V12_v_V14/"
} else if (experiment == "Mixed_V14"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V14_Mixed_graph_fulldb/"
} else if (experiment == "Mixed_V14_all_CBL"){
  indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/V14_Mixed_graph_fulldb_all_CBL/"
}

bootstrap_files <- Sys.glob(paste(indir,"*_bootstrap.txt", sep = ""))
params <- gsub("_bootstrap.txt", "", gsub(".*/", "", bootstrap_files))
summary_files <- Sys.glob(paste(indir,"*_summary.txt", sep = ""))

#target <- "FM211187.1"
h <- hash()

if (experiment == "V12_CPS")
{
  h[["barcode01"]] <- c("23F")
  h[["barcode02"]] <- c("19A")
  h[["barcode03"]] <- c("19F")
  h[["barcode04"]] <- c("03")
  h[["barcode05"]] <- c("06B")
  h[["barcode06"]] <- c("19F")
  h[["barcode07"]] <- c("23F", "19A")
  h[["barcode08"]] <- c("23F", "19F")
  h[["barcode09"]] <- c("23F", "03")
  h[["barcode10"]] <- c("23F", "06B")
  h[["barcode11"]] <- c("23F", "19F")
  h[["barcode12"]] <- c("23F")
} else if (experiment == "Mixed_V14_all_CBL") {
  target <- c("01","02","03","04","05","06A","06B","06D","07A","07B","07C","07F","08","09A","09L",
              "09N","09V","10A","10B","10C","10F","11A","11B","11C","11D","11F","12A","12B","12F",
              "13","15A","15B","15C","15F","16A","16F","17F","18A","18B","18C","18F","19A","19B",
              "19C","19F","21","22A","23A","24B","24F","25A","25F","27","28A","28F","29","31","32F",
              "33A","33B","33C","33D","34","35A","35B","35C","37","39","40","41A","41F","42","43","44",
              "45","46","47A","47F","48","06E","19AF","35D","20B","24C","06G","06F","11E","06H","22F",
              "14","17A","23B1","23B","24A","32A","33F","35F","36","38","06C","23F","15D","07D","10D","20A","10X")
  h[["barcode03"]] <- c("18C", "18B", "18F", "23F", "23A")
  h[["barcode04"]] <-  c("23F")
  h[["barcode05"]] <- c("19F", "23F")
  h[["barcode06"]] <- c("19A", "23F")
} else {
  # set up dictionary for WGS
  target <- c("23F", "FM211187.1")
  h[["barcode01"]] <- target
  h[["barcode02"]] <- target
  h[["barcode03"]] <- target
  h[["barcode04"]] <- target
  h[["barcode05"]] <- target
  h[["barcode06"]] <- target
  h[["barcode07"]] <- target
  h[["barcode08"]] <- target
  h[["barcode09"]] <- target
  h[["barcode10"]] <- target
  h[["barcode11"]] <- target
  h[["barcode12"]] <- target
  h[["barcode13"]] <- target
  h[["barcode14"]] <- target
  h[["barcode15"]] <- target
  h[["barcode16"]] <- target
  h[["barcode17"]] <- target
  h[["barcode18"]] <- target
  h[["barcode19"]] <- target
  h[["barcode20"]] <- target
  h[["barcode21"]] <- target
  h[["barcode22"]] <- target
  h[["barcode23"]] <- target
  h[["barcode24"]] <- target
}

target.df <- data.frame(Run = c(), Barcode = c(), Alignment = c(), Enrichment = c(), Type = c())

i <- 1
for (i in 1:length(bootstrap_files))
{
  filename <- bootstrap_files[i]
  summary.name <- summary_files[i]
  df <- read.table(filename, sep = "\t", comment.char = "", header = 1)
  summary.df <- read.table(summary.name, sep = "\t", comment.char = "", header = 1)
  summary.df <- subset(summary.df, Statistic == "Enrichment")
  name <- params[i]
  
  # subset data, only get enrichment for desired sequences
  barcodes <- unique(df$Barcode)
  
  for (barcode in barcodes)
  {
    # ignore non-expected barcodes
    if (!is.null(h[[barcode]]))
    {
      subsample <- subset(df, Barcode == barcode & Alignment %in% h[[barcode]])
      subsample.summary <- subset(summary.df, Barcode == barcode & Alignment %in% h[[barcode]])
      
      if (nrow(subsample.summary) != 0)
      {
        # append total stats
        to.append <- data.frame(Run = name, Barcode = subsample$Barcode, Alignment = subsample$Alignment, Enrichment = subsample$Enrichment, Type = "Bootstrap")
        target.df <- rbind(target.df, to.append)
        
        to.append <- data.frame(Run = name, Barcode = subsample.summary$Barcode, Alignment = subsample.summary$Alignment, Enrichment = subsample.summary$Value, Type = "Observed")
        target.df <- rbind(target.df, to.append)
      } 
    }
  }
}

# rename datasets
if (experiment == "V14_CPS_graph_fulldb" | experiment == "V14_CPS_graph_remclust2" | experiment == "V12_CPS_graph" | experiment == "V12_WGS_v_CPS")
{
  if (experiment == "V14_CPS_graph_remclust2" | experiment == "V14_CPS_graph_fulldb")
  {
    target.df <- subset(target.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
    target.df$Contaminant <- NA
    target.df$Contaminant_species <- NA
    target.df$Concentration <- 0
    target.df$Contaminant <- NA
    target.df$Contaminant_species <- NA
    target.df$Concentration <- 0
    target.df$Contaminant[target.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode01"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode02"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode03"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode03"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode04"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode04"] <- 0.5
    target.df$Contaminant[target.df$Barcode == "barcode05"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode06"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode07"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode07"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode08"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode08"] <- 0.5
    target.df$Contaminant[target.df$Barcode == "barcode09"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode09"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode10"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode10"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode11"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode11"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode12"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode12"] <- 0.5
    
    if (experiment == "V14_CPS_graph_remclust2")
    {
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_minimap2_rep2_all"] <- "Minimap2"
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p90_all"] <- "Graph k19 (S=90%)"
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%)"
      
      target.df$Run <- factor(target.df$Run, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    } else if (experiment == "V14_CPS_graph_fulldb")
    {
      target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_minimap2_all"] <- "Minimap2"
      target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%)"
      
      target.df$Run <- factor(target.df$Run, levels = c("Minimap2", "Graph k19 (S=75%)"))
    }
  } else if (experiment == "V12_CPS_graph" | experiment == "V12_WGS_v_CPS")
  {
    # change barcode names
    target.df <- subset(target.df, Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12")
    target.df$Contaminant <- NA
    target.df$Contaminant_species <- NA
    target.df$Concentration <- 0
    target.df$Contaminant <- NA
    target.df$Contaminant_species <- NA
    target.df$Concentration <- 0
    target.df$Contaminant[target.df$Barcode == "barcode13"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode13"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode14"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode14"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode15"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode15"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode16"] <- "S. pneumoniae R6"
    target.df$Concentration[target.df$Barcode == "barcode16"] <- 0.5
    target.df$Contaminant[target.df$Barcode == "barcode17"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode17"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode18"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode18"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode19"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode19"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode20"] <- "S. mitis"
    target.df$Concentration[target.df$Barcode == "barcode20"] <- 0.5
    target.df$Contaminant[target.df$Barcode == "barcode21"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode21"] <- 0.005
    target.df$Contaminant[target.df$Barcode == "barcode22"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode22"] <- 0.01
    target.df$Contaminant[target.df$Barcode == "barcode23"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode23"] <- 0.1
    target.df$Contaminant[target.df$Barcode == "barcode24"] <- "E. coli DH5a"
    target.df$Concentration[target.df$Barcode == "barcode24"] <- 0.5
    
    if (experiment == "V12_CPS_graph")
    {
      target.df$Run[target.df$Run == "23F_CPS_WGS_v_CPS_23F_only_all"] <- "Minimap2"
      target.df$Run[target.df$Run == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50_all"] <- "Graph k31 (S=90%)"
      
      target.df$Run <- factor(target.df$Run, levels = c("Minimap2", "Graph k31 (S=90%)"))
    } else if (experiment == "V12_WGS_v_CPS")
    {
      target.df$Run[target.df$Run == "23F_CPS_WGS_v_CPS_23F_only_all"] <- "CBL"
      target.df$Run[target.df$Run == "23F_WGS_WGS_v_CPS_all"] <- "Whole Genome"
      
      target.df$Run <- factor(target.df$Run, levels = c("Whole Genome", "CBL"))
    }
  }

  #target.df$Concentration <- target.df$Concentration * 100
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
} else if (experiment == "V12_WGS")
{
  target.df <- subset(target.df, Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04")
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant[target.df$Barcode == "barcode05"] <- "M. catarrhalis + H. influenzae"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode06"] <- "M. catarrhalis + H. influenzae"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode07"] <- "M. catarrhalis + H. influenzae"
  target.df$Concentration[target.df$Barcode == "barcode07"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode08"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode08"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode09"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode09"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode10"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode10"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode11"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode11"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode12"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode12"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode13"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode13"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode14"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode14"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode15"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode15"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode16"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode16"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode17"] <- "S. oralis"
  target.df$Concentration[target.df$Barcode == "barcode17"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode18"] <- "S. oralis"
  target.df$Concentration[target.df$Barcode == "barcode18"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode19"] <- "S. oralis"
  target.df$Concentration[target.df$Barcode == "barcode19"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode20"] <- "S. pneumoniae 110.58"
  target.df$Concentration[target.df$Barcode == "barcode20"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode21"] <- "S. pneumoniae 110.58"
  target.df$Concentration[target.df$Barcode == "barcode21"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode22"] <- "S. pneumoniae 110.58"
  target.df$Concentration[target.df$Barcode == "barcode22"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode23"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode23"] <- 0.001
  target.df$Contaminant[target.df$Barcode == "barcode24"] <- "S. pneumoniae 110.58"
  target.df$Concentration[target.df$Barcode == "barcode24"] <- 0.001
  
  target.df$Run[target.df$Run == "23F_WGS_size_selection_all"] <- "Size-selected"
  target.df$Run[target.df$Run == "23F_WGS_wo_size_selection_all"] <- "Unselected"
  #target.df$Concentration <- target.df$Concentration * 100
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a" | target.df$Contaminant == "M. catarrhalis + H. influenzae"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. oralis" | target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae 110.58" | target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Run <- factor(target.df$Run, levels = c("Unselected", "Size-selected"))
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "M. catarrhalis + H. influenzae", "S. oralis", "S. mitis", "S. pneumoniae 110.58", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
} else if (experiment == "V12_CPS")
{
  target.df <- subset(target.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
  
  target.df$Strain <- NA
  target.df$Strain.long <- NA
  target.df$Serotype <- NA
  target.df$Concentration <- 0
  
  target.df$Strain[target.df$Barcode == "barcode01"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode01"] <- "PMEN1-23F"
  target.df$Serotype[target.df$Barcode == "barcode01"] <- "23F"
  target.df$Concentration[target.df$Barcode == "barcode01"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode02"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode02"] <- "PMEN1-19A"
  target.df$Serotype[target.df$Barcode == "barcode02"] <- "19A"
  target.df$Concentration[target.df$Barcode == "barcode02"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode03"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode03"] <- "PMEN1-19F"
  target.df$Serotype[target.df$Barcode == "barcode03"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode04"] <- "PMEN33"
  target.df$Strain.long[target.df$Barcode == "barcode04"] <- "PMEN33"
  target.df$Serotype[target.df$Barcode == "barcode04"] <- "03"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode05"] <- "PMEN2"
  target.df$Strain.long[target.df$Barcode == "barcode05"] <- "PMEN2"
  target.df$Serotype[target.df$Barcode == "barcode05"] <- "06B"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode06"] <- "PMEN14"
  target.df$Strain.long[target.df$Barcode == "barcode06"] <- "PMEN14"
  target.df$Serotype[target.df$Barcode == "barcode06"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- "1.0"
  
  target.df$Strain[target.df$Barcode == "barcode07"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode07"] <- "PMEN1-19A"
  target.df$Serotype[target.df$Barcode == "barcode07"] <- "19A"
  target.df$Concentration[target.df$Barcode == "barcode07"] <- "0.5"
  
  target.df$Strain[target.df$Barcode == "barcode08"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode08"] <- "PMEN1-19F"
  target.df$Serotype[target.df$Barcode == "barcode08"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode08"] <- "0.5"
  
  target.df$Strain[target.df$Barcode == "barcode09"] <- "PMEN33"
  target.df$Strain.long[target.df$Barcode == "barcode09"] <- "PMEN33"
  target.df$Serotype[target.df$Barcode == "barcode09"] <- "03"
  target.df$Concentration[target.df$Barcode == "barcode09"] <- "0.5"
  
  target.df$Strain[target.df$Barcode == "barcode10"] <- "PMEN2"
  target.df$Strain.long[target.df$Barcode == "barcode10"] <- "PMEN2"
  target.df$Serotype[target.df$Barcode == "barcode10"] <- "06B"
  target.df$Concentration[target.df$Barcode == "barcode10"] <- "0.5"
  
  target.df$Strain[target.df$Barcode == "barcode11"] <- "PMEN14"
  target.df$Strain.long[target.df$Barcode == "barcode11"] <- "PMEN14"
  target.df$Serotype[target.df$Barcode == "barcode11"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode11"] <- "0.5"
  
  target.df$Strain[target.df$Barcode == "barcode12"] <- "R6"
  target.df$Strain.long[target.df$Barcode == "barcode12"] <- "R6"
  target.df$Serotype[target.df$Barcode == "barcode12"] <- "NT"
  target.df$Concentration[target.df$Barcode == "barcode12"] <- "0.5"
  
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  #target.df$Concentration <- target.df$Concentration * perc.genome
  
  target.df$Run[target.df$Run == "CPS_w_size_selection_all"] <- "Size-selected"
  target.df$Run[target.df$Run == "CPS_wo_size_selection_all"] <- "Unselected"
  target.df$Run <- factor(target.df$Run, levels = c("Unselected", "Size-selected"))
  target.df$Concentration <- factor(target.df$Concentration, levels = c("0.5", "1.0"))
} else if (experiment == "Minimap2_partialdb_fulldb")
{
  V12.target.df <- subset(target.df, Run == "23F_CPS_WGS_v_CPS_23F_only_all" & Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12")
  V14.target.df <- subset(target.df, Run != "23F_CPS_WGS_v_CPS_23F_only_all" & Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
  
  V12.target.df$Contaminant <- NA
  V12.target.df$Contaminant_species <- NA
  V12.target.df$Concentration <- 0
  V12.target.df$Contaminant <- NA
  V12.target.df$Contaminant_species <- NA
  V12.target.df$Concentration <- 0
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode13"] <- "S. pneumoniae R6"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode13"] <- 0.005
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode14"] <- "S. pneumoniae R6"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode14"] <- 0.01
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode15"] <- "S. pneumoniae R6"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode15"] <- 0.1
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode16"] <- "S. pneumoniae R6"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode16"] <- 0.5
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode17"] <- "S. mitis"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode17"] <- 0.005
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode18"] <- "S. mitis"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode18"] <- 0.01
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode19"] <- "S. mitis"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode19"] <- 0.1
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode20"] <- "S. mitis"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode20"] <- 0.5
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode21"] <- "E. coli DH5a"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode21"] <- 0.005
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode22"] <- "E. coli DH5a"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode22"] <- 0.01
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode23"] <- "E. coli DH5a"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode23"] <- 0.1
  V12.target.df$Contaminant[V12.target.df$Barcode == "barcode24"] <- "E. coli DH5a"
  V12.target.df$Concentration[V12.target.df$Barcode == "barcode24"] <- 0.5
  
  
  V14.target.df$Contaminant <- NA
  V14.target.df$Contaminant_species <- NA
  V14.target.df$Concentration <- 0
  V14.target.df$Contaminant <- NA
  V14.target.df$Contaminant_species <- NA
  V14.target.df$Concentration <- 0
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode01"] <- 0.005
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode02"] <- 0.01
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode03"] <- "S. pneumoniae R6"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode03"] <- 0.1
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode04"] <- "S. pneumoniae R6"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode04"] <- 0.5
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode05"] <- "S. mitis"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode05"] <- 0.005
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode06"] <- "S. mitis"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode06"] <- 0.01
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode07"] <- "S. mitis"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode07"] <- 0.1
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode08"] <- "S. mitis"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode08"] <- 0.5
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode09"] <- "E. coli DH5a"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode09"] <- 0.005
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode10"] <- "E. coli DH5a"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode10"] <- 0.01
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode11"] <- "E. coli DH5a"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode11"] <- 0.1
  V14.target.df$Contaminant[V14.target.df$Barcode == "barcode12"] <- "E. coli DH5a"
  V14.target.df$Concentration[V14.target.df$Barcode == "barcode12"] <- 0.5
  
  target.df <- rbind(V12.target.df, V14.target.df)
  
  target.df$Run[target.df$Run == "23F_CPS_WGS_v_CPS_23F_only_all"] <- "Minimap2 full db (V12)"
  target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_minimap2_all"] <- "Minimap2 full db (V14)"
  target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_minimap2_rep2_all"] <- "Minimap2 partial db (V14)"
  
  #target.df$Concentration <- target.df$Concentration * 100
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
  
} else if (experiment == "Minimap2_partialdb_fulldb")
{
  target.df <- subset(target.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant[target.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode01"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode02"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode03"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode04"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode05"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode06"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode07"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode07"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode08"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode08"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode09"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode09"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode10"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode10"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode11"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode11"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode12"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode12"] <- 0.5
  
  target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%) full db (V14)"
  target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%) partial db (V14)"
  
  #target.df$Concentration <- target.df$Concentration * 100
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
  
} else if (experiment == "V12_v_V14")
{
  
  target.df$Run[target.df$Run == "V12_CPS_23F_only_all"] <- "V12"
  target.df$Run[target.df$Run == "V14_CPS_23F_only_all"] <- "V14"
  
  target.df$Run <- factor(target.df$Run, levels = c("V12", "V14"))
  
  #remove specific barcodes
  target.df <- subset(target.df, (Run == "V14" & (Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")) | 
                        (Run == "V12" & (Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12")))
  
  
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant[target.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode01"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode02"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode03"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode04"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode05"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode06"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode07"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode07"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode08"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode08"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode09"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode09"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode10"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode10"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode11"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode11"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode12"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode12"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode13"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode13"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode14"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode14"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode15"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode15"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode16"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode16"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode17"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode17"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode18"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode18"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode19"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode19"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode20"] <- "S. mitis"
  target.df$Concentration[target.df$Barcode == "barcode20"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode21"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode21"] <- 0.005
  target.df$Contaminant[target.df$Barcode == "barcode22"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode22"] <- 0.01
  target.df$Contaminant[target.df$Barcode == "barcode23"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode23"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode24"] <- "E. coli DH5a"
  target.df$Concentration[target.df$Barcode == "barcode24"] <- 0.5
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
  
} else if (experiment == "Mixed_V14")
{
  target.df$Size.Type <- NA
  
  target.df$Size.Type[target.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Size-selected"
  target.df$Size.Type[target.df$Run == "Mixed_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Unselected"
  
  target.df$Size.Type <- factor(target.df$Size.Type, levels = c("Unselected", "Size-selected"))
  
  #remove specific barcodes
  target.df <- subset(target.df, (Size.Type == "Unselected" & (Barcode != "barcode07" & Barcode != "barcode08" &  Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" &  Barcode != "barcode12" &Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")) | 
                        (Size.Type == "Size-selected" & (Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12" & Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")))
  
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant[target.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode01"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode02"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode03"] <- "PCV-C-1464-1"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode04"] <- "PCV-C-0657-1"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode05"] <- "09B10326"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode06"] <- "PCV-C-0720-1"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.1
  
  # change back for size-selected mixture
  target.df$Contaminant[target.df$Barcode == "barcode03" & target.df$Size.Type == "Size-selected"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode03" & target.df$Size.Type == "Size-selected"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode04" & target.df$Size.Type == "Size-selected"] <- "S. pneumoniae R6"
  target.df$Concentration[target.df$Barcode == "barcode04" & target.df$Size.Type == "Size-selected"] <- 0.5
  
  target.df$Run[target.df$Run == "Mixed_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%)"
  target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%)"
  
  target.df$Run <- factor(target.df$Run, levels = c("Graph k19 (S=75%)"))
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-1464-1"] <- "Sample 1"
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-0657-1"] <- "Sample 2"
  target.df$Contaminant_species[target.df$Contaminant == "09B10326"] <- "Sample 3"
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-0720-1"] <- "Sample 4"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "Pneumo R6"
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Pneumo R6"))
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
} else if (experiment == "Mixed_V14_all_CBL")
{
  target.df$Size.Type <- NA
  
  target.df$Size.Type[target.df$Run == "Mixed_23F_fulldb_V14_400T_graphk19_p75_all_CBL_all"] <- "Unselected"
  
  target.df$Size.Type <- factor(target.df$Size.Type, levels = c("Unselected"))
  
  #remove specific barcodes
  target.df <- subset(target.df, (Size.Type == "Unselected" & (Barcode != "barcode07" & Barcode != "barcode08" &  Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" &  Barcode != "barcode12" &Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")) | 
                        (Size.Type == "Size-selected" & (Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12" & Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")))
  
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Contaminant <- NA
  target.df$Contaminant_species <- NA
  target.df$Concentration <- 0
  target.df$Concentration[target.df$Barcode == "barcode02"] <- 0.5
  target.df$Contaminant[target.df$Barcode == "barcode03"] <- "PCV-C-1464-1"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode04"] <- "PCV-C-0657-1"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode05"] <- "09B10326"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- 0.1
  target.df$Contaminant[target.df$Barcode == "barcode06"] <- "PCV-C-0720-1"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- 0.1
  
  
  target.df$Run[target.df$Run == "Mixed_23F_fulldb_V14_400T_graphk19_p75_all_CBL_all"] <- "Graph k19 (S=75%)"
  
  target.df$Run <- factor(target.df$Run, levels = c("Graph k19 (S=75%)"))
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-1464-1"] <- "Sample 1"
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-0657-1"] <- "Sample 2"
  target.df$Contaminant_species[target.df$Contaminant == "09B10326"] <- "Sample 3"
  target.df$Contaminant_species[target.df$Contaminant == "PCV-C-0720-1"] <- "Sample 4"
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"))
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
}


count.df <- target.df %>%
  group_by(Run, Barcode) %>%
  summarise(Freq = n())

# subset to remove observed
boostrap <- subset(target.df, Type == "Bootstrap")
observed <- subset(target.df, Type == "Observed")

max.min.df <- boostrap %>%
  group_by(Run, Barcode, Type) %>%
  summarise(Max = max(Enrichment), Min = min(Enrichment))

max.min.df <- merge(max.min.df, observed, by=c("Run",  "Barcode"))
max.min.df$in.range <- ifelse((max.min.df$Enrichment <= max.min.df$Max & max.min.df$Enrichment >= max.min.df$Min), TRUE, FALSE)

# plot mean comparison
p <- ggplot(observed, aes(x = Run, y = Enrichment, group = Run, colour = Run)) + theme_light() + xlab("Alignment") + ylab("Enrichment") + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 1, linetype = "dashed") + geom_jitter(size=2, alpha=0.9) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)# + coord_trans(y = 'log10') + scale_y_continuous(breaks = 10^seq(0, 6, by = 1))## + scale_y_continuous(breaks = seq(0,100,5)) 
p
ggsave(file="CPS_v_WGS_v12_mean_comp.svg", plot=p)


# plot standard error
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run, group=interaction(Concentration, Run, Contaminant_species))) + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE)) + scale_y_log10() + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1) + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_V14_minimap2_v_graph_fulldb_concentration.svg", plot=p, height = 10, width = 15)

# just plot observed values
p <- ggplot(observed, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_point(size=3) + geom_line() + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE)) + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_minimap2_v_graph_fulldb_pres_obs.svg", plot=p, height = 6, width = 15)


# just plot observed values with se
boostrap$Enrichment[boostrap$Enrichment == 0.0] <- 0.1
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, fun.min = function(z) {mean(z) - std.error(z)}, fun.max = function(z) {mean(z) + std.error(z)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment), linetype="solid") + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE)) + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_minimap2_v_graph_fulldb_pres_obs.svg", plot=p, height = 6, width = 15)

# plot quantiles
# change 0 values to very low values
boostrap$Enrichment[boostrap$Enrichment == 0.0] <- 0.1
observed$Enrichment[observed$Enrichment == 0.0] <- 0.1

p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment), linetype="solid") + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE)) + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="Graph_part_fulldb_enrichment_comp.svg", plot=p, height = 6, width = 15)

# try boxplot
p <- ggplot(boostrap, aes(x = as.factor(signif(Concentration, 1)), y = Enrichment, colour = Run)) + geom_boxplot() + geom_point(data = observed, aes(x = as.factor(signif(Concentration, 1)), y = Enrichment), size=3) + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_y_log10()# + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE)) + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_minimap2_v_graph_remclust2_concentration.svg", plot=p, height = 10, width = 15)

# plot IQR quantiles
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + facet_grid(.~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = label_log(digits = 1)) + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_V14_minimap2_v_graph_fulldb_quantile.svg", plot=p, height = 12, width = 15)

#split by group for publication with V12 WGS vs. CPS
{
  # change 0 values to very low values
  boostrap$Enrichment <- boostrap$Enrichment + 0.01
  observed$Enrichment <- observed$Enrichment + 0.01
  
  boostrap$Contaminant <- as.character(boostrap$Contaminant)
  observed$Contaminant <- as.character(observed$Contaminant)
  
  boostrap$Contaminant[boostrap$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  boostrap$Contaminant[boostrap$Contaminant == "E. coli DH5a"] <- "E. coli"
  observed$Contaminant[observed$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  observed$Contaminant[observed$Contaminant == "E. coli DH5a"] <- "E. coli"
  
  boostrap$Contaminant <- factor(boostrap$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  observed$Contaminant <- factor(observed$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  boostrap$Run <- factor(boostrap$Run, levels = c("CBL", "Whole Genome"))
  observed$Run <- factor(observed$Run, levels = c("CBL", "Whole Genome"))
  
  p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + facet_grid(~Contaminant, scales = "free_x") + geom_errorbar(stat = "summary", linewidth=1, width=0.35, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=2.5) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 14, face = "italic"), strip.text.y = element_text(size = 14, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + scale_color_npg() + scale_y_log10(limits = c(0.01, NA), breaks = c(1, 10, 100)) + scale_x_log10(labels = label_log(digits = 1)) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + guides(colour=guide_legend(title="Target Type")) + coord_cartesian(ylim=c(1, NA))
  p
  ggsave(file="V12_WGS_CPS_combined_quantile.svg", plot=p, height = 6, width = 12)
  
  operon.df.bootstrap <- subset(boostrap, Run == "CBL")
  operon.df.observed <- subset(observed, Run == "CBL")
  wgs.df.bootstrap <- subset(boostrap, Run == "Whole Genome")
  wgs.df.observed <- subset(observed, Run == "Whole Genome")
  
  # plot IQR quantiles
  p <- ggplot(operon.df.bootstrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = operon.df.observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = operon.df.observed, aes(x = Concentration, y = Enrichment)) + facet_grid(.~Contaminant) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.position = "none") + scale_color_npg() + scale_x_log10(labels = label_log(digits = 1), breaks=c(10^-5,10^-4,10^-3, 10^-2), limits=c(10^-5, 10^-2)) + scale_y_log10() + theme(strip.text.x = element_text(face = "italic")) + coord_cartesian(ylim=c(1, NA))
  p
  ggsave(file="V12_WGS_vs_CPS_CPS_only_quantile.svg", plot=p, height = 6, width = 12)
  
  p <- ggplot(wgs.df.bootstrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = wgs.df.observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = wgs.df.observed, aes(x = Concentration, y = Enrichment)) + facet_grid(.~Contaminant) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.position = "none") + scale_color_manual(values = c("#4DBBD5B2")) + scale_x_log10(labels = label_log(digits = 1), breaks=c(10^-3,10^-2,10^-1, 10^-0), limits=c(10^-3, 10^-0)) + scale_y_log10() + theme(strip.text.x = element_text(face = "italic")) + coord_cartesian(ylim=c(1, NA))
  p
  ggsave(file="V12_WGS_vs_CPS_WGS_only_quantile.svg", plot=p, height = 6, width = 12)
  
  p <- ggplot(wgs.df.bootstrap, aes(x = Concentration, y = Enrichment, colour = Contaminant)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = wgs.df.observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = wgs.df.observed, aes(x = Concentration, y = Enrichment)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.position = "none") + scale_color_npg() + scale_x_log10(labels = scientific_10, breaks=c(10^-3,10^-2,10^-1, 10^0), limits=c(10^-3, 10^0)) + scale_y_log10() + theme(strip.text.x = element_text(face = "italic"))
  p
  ggsave(file="V12_WGS_vs_CPS_WGS_only_quantile_combined.svg", plot=p, height = 6, width = 12)
}

# plot only dilutions for publication for V12 CPS
{
  strain_labeller <- function(variable,value){
    return(cluster.labs[value])
  }
  
  colour_vec <- c("23F" = "#e6194B", "3" = "#911eb4",  "6B" = "#f032e6", "19A" = "#000075", "19F" = "#42d4f4")
  #colour_vec <- c("23F" = "#e4543d", "3" = "#818181",  "6B" = "#818181", "19A" = "#818181", "19F" = "#818181")
  #colour_vec <- c("23F" = "black", "3" = "#e4543d",  "6B" = "#53b9d3", "19A" = "#0a9e87", "19F" = "#aa7ae6")
  #colour_vec <- c("23F" = "#e4543d", "3" = "#2c2c2c",  "6B" = "#555555", "19A" = "#818181", "19F" = "#b0b0b0")
  #colour_vec <- c("23F" = "#e4543d", "3" = "#555555",  "6B" = "#b0b0b0", "19A" = "#818181", "19F" = "#2c2c2c")
  
  # compare dilutions of non-serotype 23F operons, mixed and single
  mixed.iso<- subset(observed, Concentration == 0.5 & Run == "Size-selected")
  mixed.iso$Strain[mixed.iso$Strain.long == "PMEN1-19A"] <- "PMEN1 (A)"
  mixed.iso$Strain[mixed.iso$Strain.long == "PMEN1-19F"] <- "PMEN1 (B)"
  mixed.iso$Alignment[mixed.iso$Alignment == "03"] <- "3"
  mixed.iso$Alignment[mixed.iso$Alignment == "06B"] <- "6B"
  #mixed.iso$Alignment[mixed.iso$Strain == "PMEN14" & mixed.iso$Alignment == "19F"] <- "19F (PMEN14)"
  #mixed.iso$Alignment[mixed.iso$Strain == "PMEN1" & mixed.iso$Alignment == "19F"] <- "19F (PMEN1)"
  #mixed.iso$Serotype[mixed.iso$Strain == "PMEN14" & mixed.iso$Serotype == "19F"] <- "19F (PMEN14)"
  #mixed.iso$Serotype[mixed.iso$Strain == "PMEN1" & mixed.iso$Serotype == "19F"] <- "19F (PMEN1)"
  mixed.iso.bootstrap <- subset(boostrap, Concentration == 0.5 & Run == "Size-selected")
  mixed.iso.bootstrap$Strain[mixed.iso.bootstrap$Strain.long == "PMEN1-19A"] <- "PMEN1 (A)"
  mixed.iso.bootstrap$Strain[mixed.iso.bootstrap$Strain.long == "PMEN1-19F"] <- "PMEN1 (B)"
  mixed.iso.bootstrap$Alignment[mixed.iso.bootstrap$Alignment == "03"] <- "3"
  mixed.iso.bootstrap$Alignment[mixed.iso.bootstrap$Alignment == "06B"] <- "6B"
  #mixed.iso.bootstrap$Alignment[mixed.iso.bootstrap$Strain == "PMEN14" & mixed.iso.bootstrap$Alignment == "19F"] <- "19F (PMEN14)"
  #mixed.iso.bootstrap$Alignment[mixed.iso.bootstrap$Strain == "PMEN1" & mixed.iso.bootstrap$Alignment == "19F"] <- "19F (PMEN1)"
  #mixed.iso.bootstrap$Serotype[mixed.iso.bootstrap$Strain == "PMEN14" & mixed.iso.bootstrap$Serotype == "19F"] <- "19F (PMEN14)"
  #mixed.iso.bootstrap$Serotype[mixed.iso.bootstrap$Strain == "PMEN1" & mixed.iso.bootstrap$Serotype == "19F"] <- "19F (PMEN1)"
  
  # add serotype and strain together
  mixed.iso$id <- paste(mixed.iso$Strain, " (", mixed.iso$Serotype, ")", sep = "")
  mixed.iso.bootstrap$id <- paste(mixed.iso.bootstrap$Strain, " (", mixed.iso.bootstrap$Serotype, ")", sep = "")
  
  mixed.iso$Alignment <- factor(mixed.iso$Alignment, levels = c("3", "6B", "19A", "19F", "23F"))
  mixed.iso.bootstrap$Alignment <- factor(mixed.iso.bootstrap$Alignment, levels = c("3", "6B", "19A", "19F", "23F"))
  #mixed.iso$Strain <- factor(mixed.iso$Strain, levels = c("PMEN1 (A)", "PMEN1 (B)", "PMEN2", "PMEN14", "PMEN33", "R6"))
  #mixed.iso.bootstrap$Strain <- factor(mixed.iso.bootstrap$Strain, levels = c("PMEN1 (A)", "PMEN1 (B)", "PMEN2", "PMEN14", "PMEN33", "R6"))
  
  #cluster.labs <- list("PMEN1 (A)"="GPSC16", "PMEN1 (B)"="GPSC16", "PMEN2"="GPSC23", "PMEN14"="GPSC1", "PMEN33"="GPSC3", "R6"="GPSC622")
  cluster.labs <- list("GPSC16 (A)"="PMEN1 (A)", "GPSC16 (B)"="PMEN1 (B)", "GPSC23"="PMEN2", "GPSC1"="PMEN14", "GPSC3"="PMEN33", "GPSC622"="R6")
  
  cluster.labs <- unlist(cluster.labs)
  mixed.iso.bootstrap$Strain <- names(cluster.labs)[match(mixed.iso.bootstrap$Strain, cluster.labs)]
  mixed.iso$Strain <- names(cluster.labs)[match(mixed.iso$Strain, cluster.labs)]
  mixed.iso$Strain <- factor(mixed.iso$Strain, levels = c("GPSC1", "GPSC3", "GPSC16 (A)", "GPSC16 (B)", "GPSC23", "GPSC622"))
  mixed.iso.bootstrap$Strain <- factor(mixed.iso.bootstrap$Strain, levels = c("GPSC1", "GPSC3", "GPSC16 (A)", "GPSC16 (B)", "GPSC23", "GPSC622"))
  
  cluster.labs <- list("GPSC1"="GPSC1", "GPSC3"="GPSC3","GPSC16 (A)"="GPSC16", "GPSC16 (B)"="GPSC16", "GPSC23"="GPSC23", "GPSC622"="GPSC622")
  
  p <- ggplot(mixed.iso, aes(x = Alignment, y = Enrichment, colour=Alignment)) + geom_errorbar(data = mixed.iso.bootstrap,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + facet_grid(.~Strain, scales = "free_x", labeller = strain_labeller) + theme_light() + xlab("Serotype") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none", strip.placement = "outside") + scale_colour_manual(values = colour_vec)
  p
  ggsave(file="V12_CPS_mixtures_quantiles.svg", plot=p, height = 6, width = 12)
  
  #p <- ggplot(mixed.iso, aes(x = Run, y = Enrichment, colour=Alignment)) + geom_errorbar(data = mixed.iso.bootstrap,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + facet_grid(.~id, switch = "x") + theme_light() + xlab("Serotype-strain of isolate mixed with Spn23F") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14, colour = "black"), strip.background = element_rect(fill = "white"), strip.text.y = element_text(size = 14), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), strip.placement = "outside") + guides(colour=guide_legend(title="Target Serotype")) + scale_color_npg()
  #p
  
}

# plot quantiles for remclust
{
  # change 0 values to very low values
  boostrap$Enrichment <- boostrap$Enrichment + 0.01
  observed$Enrichment <- observed$Enrichment + 0.01
  
  boostrap$Contaminant <- as.character(boostrap$Contaminant)
  observed$Contaminant <- as.character(observed$Contaminant)
  
  boostrap$Contaminant[boostrap$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  boostrap$Contaminant[boostrap$Contaminant == "E. coli DH5a"] <- "E. coli"
  observed$Contaminant[observed$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  observed$Contaminant[observed$Contaminant == "E. coli DH5a"] <- "E. coli"
  
  boostrap$Contaminant <- factor(boostrap$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  observed$Contaminant <- factor(observed$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  
  # plot IQR quantiles
  p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.2, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + facet_grid(Contaminant~Run, scales = "free_y") + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16, face = "italic"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14), legend.position = "none") + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = label_log(digits = 1)) + scale_y_log10()
  p
  ggsave(file="CPS_V14_minimap2_v_graph_remclust2_quantile.svg", plot=p, height = 8, width = 10)
  
  # for presentation
  observed.subset <- subset(observed, Run != "Graph k19 (S=90%)")
  bootstrap.subset <- subset(boostrap, Run != "Graph k19 (S=90%)")
  
  observed.subset$Run <- as.character(observed.subset$Run)
  observed.subset$Run[observed.subset$Run == "Minimap2"] <- "Linear"
  observed.subset$Run[observed.subset$Run == "Graph k19 (S=75%)"] <- "Graph"
  observed.subset$Run <- factor(observed.subset$Run, levels = c("Linear", "Graph"))
  bootstrap.subset$Run <- as.character(bootstrap.subset$Run)
  bootstrap.subset$Run[bootstrap.subset$Run == "Minimap2"] <- "Linear"
  bootstrap.subset$Run[bootstrap.subset$Run == "Graph k19 (S=75%)"] <- "Graph"
  bootstrap.subset$Run <- factor(bootstrap.subset$Run, levels = c("Linear", "Graph"))
  
  p <- ggplot(bootstrap.subset, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.2, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed.subset, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed.subset, aes(x = Concentration, y = Enrichment)) + facet_grid(~Contaminant, scales = "free_y") + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16, face = "italic"), strip.text.x = element_text(size = 16, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = label_log(digits = 1)) + scale_y_log10()
  p
  ggsave(file="CPS_V14_minimap2_v_graph_remclust2_quantile_presentation.svg", plot=p, height = 6, width = 12)
}

#split by group for publication with V12 vs. V14
{
  # change 0 values to very low values
  boostrap$Enrichment <- boostrap$Enrichment + 0.01
  observed$Enrichment <- observed$Enrichment + 0.01
  
  boostrap$Contaminant <- as.character(boostrap$Contaminant)
  observed$Contaminant <- as.character(observed$Contaminant)
  
  boostrap$Contaminant[boostrap$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  boostrap$Contaminant[boostrap$Contaminant == "E. coli DH5a"] <- "E. coli"
  observed$Contaminant[observed$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  observed$Contaminant[observed$Contaminant == "E. coli DH5a"] <- "E. coli"
  
  boostrap$Contaminant <- factor(boostrap$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  observed$Contaminant <- factor(observed$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
  #boostrap$Run <- factor(boostrap$Run, levels = c("CBL", "Whole Genome"))
  #observed$Run <- factor(observed$Run, levels = c("CBL", "Whole Genome"))
  
  p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + facet_grid(Contaminant~Run, scales = "free_x") + geom_errorbar(stat = "summary", linewidth=1, width=0.1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=2.5) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16), legend.position = "none") + scale_color_npg() + scale_y_log10(limits = c(0.01, NA), breaks = c(1, 10, 100)) + scale_x_log10(labels = label_log(digits = 1)) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + guides(colour=guide_legend(title="Chemistry")) + coord_cartesian(ylim=c(1, NA))
  p
  ggsave(file="V12_v_V14_CPS_combined_quantile.svg", plot=p, height = 6, width = 8)
  
  p <- ggplot(observed, aes(x = Run, y = Enrichment, group = Run, colour = Run)) + theme_light() + xlab("Chemistry") + ylab("Enrichment") + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 1, linetype = "dashed") + geom_jitter(size=2, alpha=0.9) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)# + coord_trans(y = 'log10') + scale_y_continuous(breaks = 10^seq(0, 6, by = 1))## + scale_y_continuous(breaks = seq(0,100,5)) 
  p
  ggsave(file="V12_v_V14_CPS_mean_comp.svg", plot=p, height = 6, width = 8)
}

# for real mixed samples:
{
  boostrap$Enrichment <- boostrap$Enrichment + 0.01
  observed$Enrichment <- observed$Enrichment + 0.01
  
  sample.bootstrap <- subset(boostrap, Concentration < 0.001)
  sample.observed <- subset(observed, Concentration < 0.001)
  
  p <- ggplot(sample.bootstrap, aes(x = Contaminant_species, y = Enrichment, colour = Size.Type)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = sample.observed, aes(x = Contaminant_species, y = Enrichment), size=3) + theme_light() + xlab("Mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Library type"), shape=guide_legend(title="Target")) + scale_color_npg() 
  p
  
  observed$Concentration <- signif(observed$Concentration, digits = 1)
  boostrap$Concentration <- signif(boostrap$Concentration, digits = 1)
  boostrap$Concentration[boostrap$Concentration == 0.0008] <- "`8x10`^`-4`"
  boostrap$Concentration[boostrap$Concentration == 0.004] <- "`4x10`^`-3`"
  observed$Concentration[observed$Concentration == 0.0008] <- "`8x10`^`-4`"
  observed$Concentration[observed$Concentration == 0.004] <- "`4x10`^`-3`" 
  observed$Concentration <- factor(observed$Concentration, levels = c("`8x10`^`-4`", "`4x10`^`-3`"))
  boostrap$Concentration <- factor(boostrap$Concentration, levels = c("`8x10`^`-4`", "`4x10^`-3`"))
  
  lower.quartile <- boostrap %>% group_by(Run, Barcode, Size.Type) %>%
    summarize(low.quartile=quantile(Enrichment,probs=0.25))
  upper.quartile <- boostrap %>% group_by(Run, Barcode, Size.Type) %>%
    summarize(up.quartile=quantile(Enrichment,probs=0.75))
  
  new.df <- merge(observed, lower.quartile, by = c("Run", "Barcode", "Size.Type"))
  new.df <- merge(new.df, upper.quartile, by = c("Run", "Barcode", "Size.Type"))
  
  p <- ggplot(new.df, aes(x = Contaminant_species, y = Enrichment, group = Concentration, colour = Size.Type, shape = Concentration)) + geom_point(size=3, position = position_dodge(width = 0.5)) + geom_errorbar(aes(ymin=low.quartile, ymax=up.quartile), linewidth=1, width=0.3, position = position_dodge(width = 0.5)) + theme_light() + xlab("Mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Library type"), shape=guide_legend(title="Target\nproportion")) + scale_color_npg() + scale_shape(labels = function(x) parse(text = x))
  p
  
  ggsave(file="V14_Mixed_graphk19_p75_quantile.svg", plot=p, height = 6, width = 10)
}

# for real mixed samples looking for multiple CBL
{
  sample.bootstrap <- subset(boostrap, Concentration < 0.001)
  sample.observed <- subset(observed, Concentration < 0.001)
  
  p <- ggplot(sample.bootstrap, aes(x = Contaminant_species, y = Enrichment, colour = Size.Type)) + geom_errorbar(stat = "summary", linewidth=1, width=0.3, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = sample.observed, aes(x = Contaminant_species, y = Enrichment), size=3) + theme_light() + xlab("Mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Library type"), shape=guide_legend(title="Target")) + scale_color_npg() 
  p
  
  observed$Concentration <- signif(observed$Concentration, digits = 1)
  boostrap$Concentration <- signif(boostrap$Concentration, digits = 1)
  boostrap$Concentration[boostrap$Concentration == 0.0008] <- "`8x10`^`-4`"
  boostrap$Concentration[boostrap$Concentration == 0.004] <- "`4x10`^`-3`"
  observed$Concentration[observed$Concentration == 0.0008] <- "`8x10`^`-4`"
  observed$Concentration[observed$Concentration == 0.004] <- "`4x10`^`-3`" 
  observed$Concentration <- factor(observed$Concentration, levels = c("`8x10`^`-4`", "`4x10`^`-3`"))
  boostrap$Concentration <- factor(boostrap$Concentration, levels = c("`8x10`^`-4`", "`4x10^`-3`"))
  
  lower.quartile <- boostrap %>% group_by(Run, Barcode, Size.Type, Alignment) %>%
    summarize(low.quartile=quantile(Enrichment,probs=0.25))
  upper.quartile <- boostrap %>% group_by(Run, Barcode, Size.Type, Alignment) %>%
    summarize(up.quartile=quantile(Enrichment,probs=0.75))
  
  new.df <- merge(observed, lower.quartile, by = c("Run", "Barcode", "Size.Type", "Alignment"))
  new.df <- merge(new.df, upper.quartile, by = c("Run", "Barcode", "Size.Type", "Alignment"))
  
  new.df$colour <- FALSE
  new.df$colour[new.df$Alignment == "23F"] <- TRUE
  
  p <- ggplot(new.df, aes(x = Alignment, y = Enrichment, group = Alignment, colour = colour)) + facet_grid(~Contaminant_species, scales = "free_x") + geom_errorbar(aes(ymin=low.quartile, ymax=up.quartile), linewidth=1, width=0.3, position = position_dodge(width = 0.5)) + geom_point(size=3, position = position_dodge(width = 0.5)) + theme_light() + xlab("Serotype detected") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14), legend.position = "none") + guides(colour=guide_legend(title="Library type"), shape=guide_legend(title="Target\nproportion")) + scale_color_npg() + scale_shape(labels = function(x) parse(text = x)) #+ scale_y_log10() # +
  p
  
  ggsave(file="V14_Mixed_graphk19_p75_multiCBL_quantile.svg", plot=p, height = 6, width = 12)
}

# for presentations
# remove k19 90% pid
boostrap <- subset(boostrap, Run != "Graph k19 (S=90%)")
boostrap$Run <- as.character(boostrap$Run)
boostrap$Run[boostrap$Run == "Graph k19 (S=75%)"] <- "Graph"
boostrap$Run <- factor(boostrap$Run, levels = c("Minimap2", "Graph"))
observed <- subset(observed, Run != "Graph k19 (S=90%)")
observed$Run <- as.character(observed$Run)
observed$Run[observed$Run == "Graph k19 (S=75%)"] <- "Graph"
observed$Run <- factor(observed$Run, levels = c("Minimap2", "Graph"))


# plot quantiles
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=2.5) + geom_line(data = observed, aes(x = Concentration, y = Enrichment), linetype="dashed") + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Alignment"), shape=guide_legend(title="Target")) + scale_y_log10() + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_minimap2_v_graph_remclust2_quantiles_for_pres.svg", plot=p, height = 6, width = 15)


# for V12 WGS
# just plot observed values
p <- ggplot(observed, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_point(size=3) + geom_line() + facet_wrap(Contaminant_species~Contaminant, nrow=3) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Experiment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="WGS_23F_obs.svg", plot=p, height = 12, width = 15)

# plot IQR
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.4, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + facet_wrap(Contaminant_species~Contaminant, nrow=3) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=14)) + guides(colour=guide_legend(title="Library type"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = label_log(digits = 1))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="WGS_23F_quantiles.svg", plot=p, height = 12, width = 15)

boostrap <- subset(boostrap, Contaminant == "S. pneumoniae 110.58" | Contaminant == "S. mitis" |  Contaminant == "E. coli DH5a")
observed <- subset(observed, Contaminant == "S. pneumoniae 110.58" | Contaminant == "S. mitis" |  Contaminant == "E. coli DH5a")
p <- ggplot(boostrap, aes(x = Concentration, y = Enrichment, colour = Run)) + geom_errorbar(stat = "summary", linewidth=1, width=0.4, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(data = observed, aes(x = Concentration, y = Enrichment), size=3) + geom_line(data = observed, aes(x = Concentration, y = Enrichment)) + facet_grid(~Contaminant_species) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Experiment"), shape=guide_legend(title="Target")) + scale_color_npg() + scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="WGS_23F_quantiles_for_pres.svg", plot=p, height = 4, width = 16)


# for V12 CPS
# just plot observed values
# single isolate only
# single.iso <- subset(observed, Concentration== "100%")
# p <- ggplot(single.iso, aes(x = Strain, y = Enrichment, colour = Alignment)) + geom_point(size=3) + theme_light() + facet_grid(~Run) + xlab("Strain background") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Target CBL")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
# ggsave(file="CPS_v12_singleiso_obs.svg", plot=p, height = 12, width = 15)
# # plot IQR
# single.iso.bootstrap <- subset(boostrap, Concentration == "100%")
# p <- ggplot(single.iso.bootstrap, aes(x = Strain, y = Enrichment, colour = Alignment)) + geom_errorbar(position=position_dodge(width = 0.9), stat = "summary", linewidth=1, width=0.5, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), data = single.iso, aes(x = Strain, y = Enrichment), size=3) + facet_grid(~Run) + theme_light() + xlab("Strain Background") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Target CBL")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
# ggsave(file="CPS_v12_singleiso_quantiles.svg", plot=p, height = 12, width = 15)
# 
# # mixture isolate only
# mixed.iso <- subset(observed, Concentration== "50%")
# p <- ggplot(mixed.iso, aes(x = Strain, y = Enrichment, colour = Alignment)) + geom_point(size=3) + theme_light() + facet_grid(~Run) + xlab("Strain mixed with PMEN1 (23F)") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Target CBL")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + scale_y_log10()# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
# ggsave(file="CPS_v12_mixediso_obs.svg", plot=p, height = 12, width = 15)
# # plot IQR of strain background
# mixed.iso.23F <- subset(observed, Concentration == "50%" & Alignment == "23F")
# mixed.iso.bootstrap.23F <- subset(boostrap, Concentration == "50%" & Alignment == "23F")
# p <- ggplot(mixed.iso.bootstrap.23F, aes(x = Strain.long, y = Enrichment, colour=Strain)) + geom_errorbar(position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), data = mixed.iso.23F, aes(x = Strain.long, y = Enrichment), size=3) + facet_grid(~Run) + theme_light() + xlab("Strain mixed with PMEN1-23F") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Mixture Strain Background")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
# ggsave(file="CPS_v12_mixediso_23Fonly_quantiles.svg", plot=p, height = 12, width = 15)
# 
# # IQR of serotype
# p <- ggplot(mixed.iso.bootstrap.23F, aes(x = Serotype, y = Enrichment, colour=Strain)) + geom_errorbar(position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), data = mixed.iso.23F, aes(x = Serotype, y = Enrichment), size=3) + facet_grid(~Run) + theme_light() + xlab("Serotype mixed with PMEN1-23F") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Mixture Strain Background")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
# ggsave(file="CPS_v12_mixediso_23Fonly_serotype_quantiles.svg", plot=p, height = 12, width = 15)
# 
# # compare dilutions of non-serotype 23F operons, mixed only
# mixed.iso.non23F <- subset(observed, Alignment != "23F" & Concentration== "50%")
# #mixed.iso.non23F$Strain[mixed.iso.non23F$Strain == "PMEN1-19A"] <- "PMEN1"
# #mixed.iso.non23F$Strain[mixed.iso.non23F$Strain == "PMEN1-19F"] <- "PMEN1"
# #mixed.iso.non23F$Alignment[mixed.iso.non23F$Strain == "PMEN14"] <- "19F (PMEN14)"
# #mixed.iso.non23F$Alignment[mixed.iso.non23F$Strain == "PMEN1" & mixed.iso.non23F$Alignment == "19F"] <- "19F (PMEN1)"
# mixed.iso.bootstrap.non23F <- subset(boostrap,Alignment != "23F" & Concentration== "50%")
# #mixed.iso.bootstrap.non23F$Strain[mixed.iso.bootstrap.non23F$Strain == "PMEN1-19A"] <- "PMEN1"
# #mixed.iso.bootstrap.non23F$Strain[mixed.iso.bootstrap.non23F$Strain == "PMEN1-19F"] <- "PMEN1"
# #mixed.iso.bootstrap.non23F$Alignment[mixed.iso.bootstrap.non23F$Strain == "PMEN14"] <- "19F (PMEN14)"
# #mixed.iso.bootstrap.non23F$Alignment[mixed.iso.bootstrap.non23F$Strain == "PMEN1" & mixed.iso.bootstrap.non23F$Alignment == "19F"] <- "19F (PMEN1)"
# 
# p <- ggplot(mixed.iso.bootstrap.non23F, aes(x = Serotype, y = Enrichment, colour=Strain)) + geom_errorbar(position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), data = mixed.iso.non23F, aes(x = Serotype, y = Enrichment), size=3) + facet_grid(~Run) + theme_light() + xlab("Serotype mixed with PMEN1-23F") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust=1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Mixture Strain Background")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
# p
#ggsave(file="CPS_v12_mixediso_non23F_serotype_quantiles.svg", plot=p, height = 12, width = 15)

# look at single isolate mixtures
single.iso <- subset(observed, Concentration == "100%")
single.iso$Strain[single.iso$Strain == "PMEN1-19A"] <- "PMEN1"
single.iso$Strain[single.iso$Strain == "PMEN1-19F"] <- "PMEN1"
single.iso$Alignment[single.iso$Strain == "PMEN14"] <- "19F (PMEN14)"
single.iso$Alignment[single.iso$Strain == "PMEN1" & single.iso$Alignment == "19F"] <- "19F (PMEN1)"
single.iso.bootstrap <- subset(boostrap, Concentration == "100%")
single.iso.bootstrap$Strain[single.iso.bootstrap$Strain == "PMEN1-19A"] <- "PMEN1"
single.iso.bootstrap$Strain[single.iso.bootstrap$Strain == "PMEN1-19F"] <- "PMEN1"
single.iso.bootstrap$Alignment[single.iso.bootstrap$Strain == "PMEN14"] <- "19F (PMEN14)"
single.iso.bootstrap$Alignment[single.iso.bootstrap$Strain == "PMEN1" & single.iso.bootstrap$Alignment == "19F"] <- "19F (PMEN1)"

single.iso <- subset(single.iso, Run == "Size-selected")
single.iso.bootstrap <- subset(single.iso.bootstrap, Run == "Size-selected")
p <- ggplot(single.iso, aes(x = Serotype, y = Enrichment, colour=Strain)) + geom_errorbar(data = single.iso.bootstrap,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + theme_light() + xlab("% Strain in mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Strain Background"), shape=guide_legend(title="Target")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_v12_singleiso_serotype_quantiles_for_pres.svg", plot=p, height = 6, width = 12)

# compare dilutions of non-serotype 23F operons, mixed and single
mixed.iso.non23F <- subset(observed, Alignment != "23F")
mixed.iso.non23F$Strain[mixed.iso.non23F$Strain == "PMEN1-19A"] <- "PMEN1"
mixed.iso.non23F$Strain[mixed.iso.non23F$Strain == "PMEN1-19F"] <- "PMEN1"
mixed.iso.non23F$Alignment[mixed.iso.non23F$Strain == "PMEN14"] <- "19F (PMEN14)"
mixed.iso.non23F$Alignment[mixed.iso.non23F$Strain == "PMEN1" & mixed.iso.non23F$Alignment == "19F"] <- "19F (PMEN1)"
mixed.iso.bootstrap.non23F <- subset(boostrap,Alignment != "23F")
mixed.iso.bootstrap.non23F$Strain[mixed.iso.bootstrap.non23F$Strain == "PMEN1-19A"] <- "PMEN1"
mixed.iso.bootstrap.non23F$Strain[mixed.iso.bootstrap.non23F$Strain == "PMEN1-19F"] <- "PMEN1"
mixed.iso.bootstrap.non23F$Alignment[mixed.iso.bootstrap.non23F$Strain == "PMEN14"] <- "19F (PMEN14)"
mixed.iso.bootstrap.non23F$Alignment[mixed.iso.bootstrap.non23F$Strain == "PMEN1" & mixed.iso.bootstrap.non23F$Alignment == "19F"] <- "19F (PMEN1)"

p <- ggplot(mixed.iso.non23F, aes(x = Concentration, y = Enrichment, group=Strain.long, colour=Strain)) + geom_errorbar(data = mixed.iso.bootstrap.non23F,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + facet_grid(Run~Alignment) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Strain Background"), shape=guide_legend(title="Target")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_v12_mixediso_non23F_serotype_quantiles.svg", plot=p, height = 12, width = 15)

mixed.iso.non23F <- subset(mixed.iso.non23F, Run == "Size-selected")
mixed.iso.bootstrap.non23F <- subset(mixed.iso.bootstrap.non23F, Run == "Size-selected")
p <- ggplot(mixed.iso.non23F, aes(x = Concentration, y = Enrichment, group=Strain.long, colour=Strain)) + geom_errorbar(data = mixed.iso.bootstrap.non23F,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + facet_grid(~Alignment) + theme_light() + xlab("% Strain in mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 14), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Strain Background"), shape=guide_legend(title="Target")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_v12_mixediso_non23F_serotype_quantiles_for_pres.svg", plot=p, height = 6, width = 12)

# compare dilutions of 23F operons, mixed and single
mixed.iso.23F <- subset(observed, Alignment == "23F")
mixed.iso.23F$Strain[mixed.iso.23F$Strain == "PMEN1-19A"] <- "PMEN1"
mixed.iso.23F$Strain[mixed.iso.23F$Strain == "PMEN1-19F"] <- "PMEN1"
mixed.iso.23F$Serotype[mixed.iso.23F$Strain == "PMEN14"] <- "19F (PMEN14)"
mixed.iso.23F$Serotype[mixed.iso.23F$Strain == "PMEN1" & mixed.iso.23F$Serotype == "19F"] <- "19F (PMEN1)"

# add 100% 23F data to all points
# mixed.iso.23F.temp <- subset(mixed.iso.23F, Serotype != "23F")
# mixed.iso.23F.temp$Concentration <- "100%"
# mixed.iso.23F.ss <- mixed.iso.23F[mixed.iso.23F$Serotype == "23F" & mixed.iso.23F$Run == "Size-selected",]
# mixed.iso.23F.ori <- mixed.iso.23F[mixed.iso.23F$Serotype == "23F" & mixed.iso.23F$Run == "Unselected",]
# for(i in 1:nrow(mixed.iso.23F.temp)) {
#   row <- mixed.iso.23F.temp[i,]
#   # do stuff with row
#   mixed.iso.23F.temp[i,]$Enrichment <- ifelse(row$Run == "Unselected", mixed.iso.23F.ori$Enrichment, mixed.iso.23F.ss$Enrichment)
# }
# mixed.iso.23F <- subset(mixed.iso.23F, Serotype != "23F")
# mixed.iso.23F <- rbind(mixed.iso.23F, mixed.iso.23F.temp)

# repeat for bootstrapped samples
mixed.iso.bootstrap.23F <- subset(boostrap,Alignment == "23F")
mixed.iso.bootstrap.23F$Strain[mixed.iso.bootstrap.23F$Strain == "PMEN1-19A"] <- "PMEN1"
mixed.iso.bootstrap.23F$Strain[mixed.iso.bootstrap.23F$Strain == "PMEN1-19F"] <- "PMEN1"
mixed.iso.bootstrap.23F$Serotype[mixed.iso.bootstrap.23F$Strain == "PMEN14"] <- "19F (PMEN14)"
mixed.iso.bootstrap.23F$Serotype[mixed.iso.bootstrap.23F$Strain == "PMEN1" & mixed.iso.bootstrap.23F$Serotype == "19F"] <- "19F (PMEN1)"

# # add 100% 23F data to all points
# mixed.iso.bootstrap.23F.temp <- subset(mixed.iso.bootstrap.23F, Serotype != "23F")
# mixed.iso.bootstrap.23F.temp$Concentration <- "100%"
# mixed.iso.23F.ss <- mixed.iso.bootstrap.23F[mixed.iso.bootstrap.23F$Serotype == "23F" & mixed.iso.bootstrap.23F$Run == "Size-selected",]
# mixed.iso.23F.ori <- mixed.iso.bootstrap.23F[mixed.iso.bootstrap.23F$Serotype == "23F" & mixed.iso.bootstrap.23F$Run == "Unselected",]
# barcodes <- unique(mixed.iso.bootstrap.23F.temp$Barcode)
# runs <- unique(mixed.iso.bootstrap.23F.temp$Run)
# 
# for(b in barcodes) {
#   for (r in runs)
#   {
#     df <- subset(mixed.iso.bootstrap.23F.temp, Barcode == b & Run == r)
#     if (r == "Unselected")
#     {
#       df$Enrichment <- mixed.iso.23F.ori$Enrichment
#     } else
#     {
#       df$Enrichment <- mixed.iso.23F.ss$Enrichment
#     }
#     mixed.iso.bootstrap.23F.temp[mixed.iso.bootstrap.23F.temp$Barcode == b & mixed.iso.bootstrap.23F.temp$Run == r,] <- df
#   }
# }
# mixed.iso.bootstrap.23F <- subset(mixed.iso.bootstrap.23F, Serotype != "23F")
# mixed.iso.bootstrap.23F <- rbind(mixed.iso.bootstrap.23F, mixed.iso.bootstrap.23F.temp)

mixed.iso.23F$Serotype[mixed.iso.23F$Serotype == "23F"] <- "Undiluted"
mixed.iso.bootstrap.23F$Serotype[mixed.iso.bootstrap.23F$Serotype == "23F"] <- "Undiluted"
mixed.iso.23F$Serotype <- factor(mixed.iso.23F$Serotype, level = c("Undiluted", "03", "06B", "19A", "19F (PMEN1)", "19F (PMEN14)", "NT"))
mixed.iso.bootstrap.23F$Serotype <- factor(mixed.iso.bootstrap.23F$Serotype, level = c("Undiluted", "03", "06B", "19A", "19F (PMEN1)", "19F (PMEN14)", "NT"))
p <- ggplot(mixed.iso.23F, aes(x = Serotype, y = Enrichment, colour=Strain)) + geom_errorbar(data = mixed.iso.bootstrap.23F,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + facet_grid(Run~.) + theme_light() + xlab("Contaminant Serotype") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 14), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Strain Background"), shape=guide_legend(title="Target")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_v12_mixediso_23F_serotype_quantiles.svg", plot=p, height = 12, width = 15)

mixed.iso.23F <- subset(mixed.iso.23F, Run == "Size-selected")
mixed.iso.bootstrap.23F <- subset(mixed.iso.bootstrap.23F, Run == "Size-selected")
p <- ggplot(mixed.iso.23F, aes(x = Concentration, y = Enrichment, group=Strain.long, colour=Strain)) + geom_errorbar(data = mixed.iso.bootstrap.23F,position=position_dodge(width = 0.9), width=0.5, stat = "summary", linewidth=1, fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}) + geom_point(position=position_dodge(width = 0.9), size=3) + geom_line() + facet_grid(~Serotype) + theme_light() + xlab("% PMEN1-23F in mixture") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Mixture Strain Background"), shape=guide_legend(title="Target")) + scale_color_npg() #+ scale_x_log10(labels = function(x) format(x, scientific = TRUE))# + stat_summary(fun.data = mean_se, geom = "errorbar", linewidth=1)# + stat_summary(fun.y="mean", geom="line", linewidth=1, linetype = "dotted", aes(group=interaction(Run, Contaminant_species))) + stat_summary(fun.y="mean", geom="point", aes(group=interaction(Run, Contaminant_species)))  
p
ggsave(file="CPS_v12_mixediso_23F_serotype_quantiles_for_pres.svg", plot=p, height = 6, width = 12)
