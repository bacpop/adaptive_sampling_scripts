library(ggplot2)
library(data.table)
library(purrr)
library(ggpubr)
library(ggsci)
library(stringr)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(hash)
library(ggtext)


# read in files
indir = "./coverage_analysis/CPS_full_genome_alignment/"
prefix <- ""
experiment <- "V12_WGS_v_CPS"
hist_files <- Sys.glob(paste(indir, experiment,"/", prefix, "*.csv", sep = ""))
summary_files <- Sys.glob(paste(indir, experiment,"/*.txt", sep = ""))
histify <- TRUE
write.output <- FALSE
smooth <- TRUE
sample <- FALSE
CPS.only <- FALSE

i <- 1
for (i in 1:length(hist_files))
{
  filename <- hist_files[i]
  #summary.name <- summary_files[i]
  df <- read.table(filename, sep = ",", comment.char = "", header = 0)
  
  if (CPS.only == TRUE)
  {
    temp.df <- df
    temp.df$Position <- seq(1:nrow(temp.df))
    temp.df <- subset(temp.df, Position >= 303559 & Position <= 322212)
    df <- data.frame(V1 = temp.df$V1)
  }
  
  #summary.df <- read.table(summary.name, sep = "\t", comment.char = "", header = 1)
  #summary.df <- subset(summary.df, Statistic == "Enrichment")
  params <- gsub("_hist.csv", "", gsub(".*/", "", filename))
  params <- str_split(params, "_all_")[[1]]
  
  run.name <- params[1]
  details <- str_split(params[2], "_")[[1]]
  barcode <- details[length(details) - 2]
  align <- details[length(details) - 1]
  channel <- details[length(details)]
  
  if (sample == FALSE)
  {
    if (grepl("sample", filename) == TRUE)
    {
      next
    }
  } else
  {
    if (grepl("sample", filename) == FALSE)
    {
      next
    }
  }
  
  if (histify == TRUE)
  {
    df$index <- seq(1, nrow(df))
    
    if (grepl("WGS", run.name) == TRUE)
    {
      if (grepl("CPS_only", run.name) == TRUE | grepl("23F_only", run.name) == TRUE | grepl("23F", run.name) == TRUE)
      {
        type = "CPS"
      }
      else{
        type = "WGS"
      }
    } else {
      type = "CPS"
    }
    
    if (smooth == TRUE)
    {
      bins <- 20000
    } else if (type == "WGS")
    {
      bins <- 20000
    } else {
      bins <- 250
    }
    
    df2 <- aggregate(df, #the data frame
                     by=list(cut(df$index, seq(1,nrow(df),bins))), #the bins (see below)
                     mean) #the aggregating function
    
    # df2 %>%
    #   extract(Group.1, c("start", "end"), "(-?\\d+),(-?\\d+)")
    
    df2$Group.1 <- as.character(df2$Group.1)
    
    df.parsed <- df2 %>% separate(Group.1, c("A", "B"), sep = ",")
    df.parsed$A <- gsub("\\[|\\]", "", df.parsed$A)
    df.parsed$A <- gsub("\\(|\\)", "", df.parsed$A)
    df.parsed$B <- gsub("\\[|\\]", "", df.parsed$B)
    df.parsed$B <- gsub("\\(|\\)", "", df.parsed$B)
    df.parsed$A <- as.numeric(df.parsed$A)
    df.parsed$B <- as.numeric(df.parsed$B)
    
    to.file <- data.frame(chr = align, start = df.parsed$A, end = df.parsed$B, coverage = df.parsed$V1)
    
    if (write.output == TRUE)
    {
      outdir <- paste("./assemblies/Inspector_all_structural_errors/", run.name, "_sup", sep = "")
      outdir <- paste(outdir, "/bc_", barcode, "_ch_", channel, "_al_combined", sep = "")
      
      error <- try(write.table(to.file, file = paste(outdir, "/coverage_ref", ".bed", sep = ""), row.names=FALSE, col.names = FALSE, sep="\t"), silent = TRUE)
      if(inherits(error, "try-error"))
      {
        next
      }
    } else
    {
      if (i == 1)
      {
        target.df <- data.frame(Run = run.name, Barcode = barcode, Alignment = to.file$chr, Start = to.file$start, End=to.file$end,  Coverage = to.file$coverage, Channel = channel)
      } else 
      {
        to.append <- data.frame(Run = run.name, Barcode = barcode, Alignment = to.file$chr, Start = to.file$start, End=to.file$end,  Coverage = to.file$coverage, Channel = channel)
        target.df <- rbind(target.df, to.append)
      }
    }
    
  } else {
    #append total stats
    if (i == 1)
    {
      target.df <- data.frame(Run = run.name, Barcode = barcode, Alignment = align, Position = seq(1:length(df$V1)), Coverage = df$V1, Channel = channel)
    } else 
    {
      to.append <- data.frame(Run = run.name, Barcode = barcode, Alignment = align, Position = seq(1:length(df$V1)), Coverage = df$V1, Channel = channel)
      target.df <- rbind(target.df, to.append)
    }
  }
  
  prev.barcode <- barcode
  prev.channel <- channel
}


for (i in 1:length(summary_files))
{
  filename <- summary_files[i]
  #summary.name <- summary_files[i]
  df <- read.table(filename, sep = "\t", comment.char = "", header = 1)
  #summary.df <- read.table(summary.name, sep = "\t", comment.char = "", header = 1)
  #summary.df <- subset(summary.df, Statistic == "Enrichment")
  params <- gsub("_all_summary.txt", "", gsub(".*/", "", filename))
  params <- gsub("_all_Spn23F_summary.txt", "", gsub(".*/", "", params))

  df$Run <- params

  # append total stats
  if (i == 1)
  {
    summary.df <- df
  } else 
  {
    summary.df <- rbind(summary.df, df)
  }
}
# just take control bases
#summary.df <- subset(summary.df, Channel == "control")

# normalise coverage by total bases
total.summary.df <- summary.df %>%
  group_by(Run, Channel, Barcode) %>%
  summarise(Total.bases = sum(Total_bases) / 1000000000)

target.df <- merge(target.df, total.summary.df, by = c("Run", "Channel", "Barcode"))
target.df$norm.coverage <- target.df$Coverage / target.df$Total.bases

# remove NA barcodes
target.df <- subset(target.df, Barcode != "NA")
target.df <- na.omit(target.df)

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
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "Minimap2"
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "Graph k19 (S=90%)"
      target.df$Run[target.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "Graph k19 (S=75%)"
      
      target.df$Run <- factor(target.df$Run, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    } else if (experiment == "V14_CPS_graph_fulldb")
    {
      target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_minimap2"] <- "Minimap2"
      target.df$Run[target.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "Graph k19 (S=75%)"
      
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
      target.df$Run[target.df$Run == "CPS_experiment3_23F_minimap2_23F_only"] <- "Minimap2"
      target.df$Run[target.df$Run == "CPS_experiment3_23F_graph_k31_90_c50_23F_only"] <- "Graph k31 (S=90%)"
      
      target.df$Run <- factor(target.df$Run, levels = c("Minimap2", "Graph k31 (S=90%)"))
    } else if (experiment == "V12_WGS_v_CPS")
    {
      target.df$Run[target.df$Run == "23F_CPS_WGS_v_CPS_23F_only"] <- "CBL"
      target.df$Run[target.df$Run == "23F_WGS_WGS_v_CPS"] <- "Whole Genome"
      
      target.df$Run <- factor(target.df$Run, levels = c("Whole Genome", "CBL"))
    }
    
  }
  
  #target.df$Concentration <- target.df$Concentration * 100
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  target.df$Concentration <- as.numeric(target.df$Concentration)
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a"] <- "E. coli"
  target.df$Contaminant_species[target.df$Contaminant == "S. mitis"] <- "S. mitis"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  
  target.df$Concentration[target.df$Run != "Whole Genome"] <- target.df$Concentration[target.df$Run != "Whole Genome"] * perc.genome
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
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
  
  target.df$Run[target.df$Run == "23F_WGS_size_selection"] <- "Size-selected"
  target.df$Run[target.df$Run == "23F_WGS_wo_size_selection"] <- "Unselected"
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
  target.df$Concentration[target.df$Barcode == "barcode01"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode02"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode02"] <- "PMEN1-19A"
  target.df$Serotype[target.df$Barcode == "barcode02"] <- "19A"
  target.df$Concentration[target.df$Barcode == "barcode02"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode03"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode03"] <- "PMEN1-19F"
  target.df$Serotype[target.df$Barcode == "barcode03"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode03"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode04"] <- "PMEN33"
  target.df$Strain.long[target.df$Barcode == "barcode04"] <- "PMEN33"
  target.df$Serotype[target.df$Barcode == "barcode04"] <- "03"
  target.df$Concentration[target.df$Barcode == "barcode04"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode05"] <- "PMEN2"
  target.df$Strain.long[target.df$Barcode == "barcode05"] <- "PMEN2"
  target.df$Serotype[target.df$Barcode == "barcode05"] <- "06B"
  target.df$Concentration[target.df$Barcode == "barcode05"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode06"] <- "PMEN14"
  target.df$Strain.long[target.df$Barcode == "barcode06"] <- "PMEN14"
  target.df$Serotype[target.df$Barcode == "barcode06"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode06"] <- "100%"
  
  target.df$Strain[target.df$Barcode == "barcode07"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode07"] <- "PMEN1-19A"
  target.df$Serotype[target.df$Barcode == "barcode07"] <- "19A"
  target.df$Concentration[target.df$Barcode == "barcode07"] <- "50%"
  
  target.df$Strain[target.df$Barcode == "barcode08"] <- "PMEN1"
  target.df$Strain.long[target.df$Barcode == "barcode08"] <- "PMEN1-19F"
  target.df$Serotype[target.df$Barcode == "barcode08"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode08"] <- "50%"
  
  target.df$Strain[target.df$Barcode == "barcode09"] <- "PMEN33"
  target.df$Strain.long[target.df$Barcode == "barcode09"] <- "PMEN33"
  target.df$Serotype[target.df$Barcode == "barcode09"] <- "03"
  target.df$Concentration[target.df$Barcode == "barcode09"] <- "50%"
  
  target.df$Strain[target.df$Barcode == "barcode10"] <- "PMEN2"
  target.df$Strain.long[target.df$Barcode == "barcode10"] <- "PMEN2"
  target.df$Serotype[target.df$Barcode == "barcode10"] <- "06B"
  target.df$Concentration[target.df$Barcode == "barcode10"] <- "50%"
  
  target.df$Strain[target.df$Barcode == "barcode11"] <- "PMEN14"
  target.df$Strain.long[target.df$Barcode == "barcode11"] <- "PMEN14"
  target.df$Serotype[target.df$Barcode == "barcode11"] <- "19F"
  target.df$Concentration[target.df$Barcode == "barcode11"] <- "50%"
  
  target.df$Strain[target.df$Barcode == "barcode12"] <- "R6"
  target.df$Strain.long[target.df$Barcode == "barcode12"] <- "R6"
  target.df$Serotype[target.df$Barcode == "barcode12"] <- "NT"
  target.df$Concentration[target.df$Barcode == "barcode12"] <- "50%"
  
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  #target.df$Concentration <- target.df$Concentration * perc.genome
  
  target.df$Run[target.df$Run == "CPS_w_size_selection"] <- "Size-selected"
  target.df$Run[target.df$Run == "CPS_wo_size_selection"] <- "Unselected"
  target.df$Run <- factor(target.df$Run, levels = c("Unselected", "Size-selected"))
  target.df$Concentration <- factor(target.df$Concentration, levels = c("50%", "100%"))
}

# remove NA barcodes
target.df <- subset(target.df, Barcode != "NA")
target.df <- na.omit(target.df)

# get positions
target.df$Concentration <- signif(target.df$Concentration, 1)
target.df <- target.df[order( target.df[,1], target.df[,2], target.df[,3], target.df[,4], target.df[,5]),]
target.df$Position <- round((target.df$Start + target.df$End) / 2)

# subsample for CPS
if (type == "CPS" & max(target.df$Position) > 1000000)
{
  target.df.ori <- target.df
  target.df <- subset(target.df, Position >= 303559 & Position <= 322212)
  
  # adjust coordinates
  min.coord <- min(target.df$Start)
  target.df$Start <- (target.df$Start - min.coord)
  target.df$End <- (target.df$End - min.coord)
  target.df$Position <- round((target.df$Start + target.df$End) / 2)
  
  # change names from non-streptococcus etc.
  target.df$Contaminant <- as.character(target.df$Contaminant)
  
  target.df$Contaminant[target.df$Contaminant == "S. pneumoniae R6"] <- "S. pneumoniae"
  target.df$Contaminant[target.df$Contaminant == "E. coli DH5a"] <- "E. coli"
  
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli", "S. mitis", "S. pneumoniae"))
}

#plot coverage for CPS operon
if (histify == TRUE)
{
  Spn23F.label.bed <- data.frame(
    chr = "NA",
    Start = c(303559,1487660, 1207639),
    End = c(322212, 1526939, 1288735),
    label = c("CBL", "\u03C6MM1", "ICE*Sp*23FST81")
  )
  Spn23F.label.bed$Start <- Spn23F.label.bed$Start / 1000000
  Spn23F.label.bed$End <- Spn23F.label.bed$End / 1000000
  
  
  if (experiment == "V12_CPS")
  {
    target.df$Run[target.df$Run == "CPS_w_size_selection"] <- "Size-selected"
    target.df$Run[target.df$Run == "CPS_wo_size_selection"] <- "Unselected"
    target.df$Run <- factor(target.df$Run, levels = c("Unselected", "Size-selected"))
    sample.df <- subset(target.df, Barcode == "barcode01")
  } else if (experiment == "V12_WGS_v_CPS")
  {
    sample.df <- subset(target.df, (Run == "Whole Genome" & Concentration == 0.1) | (Run == "CBL" & Concentration == 0.0008))
  }
  
  sample.df <- sample.df[order(sample.df$Start),]
  
  sample.df$Channel[sample.df$Channel == "adaptive"] <- "NAS"
  sample.df$Channel[sample.df$Channel == "control"] <- "Control"
  sample.df$Channel <- factor(sample.df$Channel, levels = c("NAS", "Control"))
  
  if (experiment == "V12_CPS")
  {
    Spn23F.label.bed$norm.coverage <- 600
    Spn23F.label.bed2 <- Spn23F.label.bed
    Spn23F.label.bed$Channel <- "NAS"
    Spn23F.label.bed2$Channel <- "Control"
    Spn23F.label.bed <- rbind(Spn23F.label.bed, Spn23F.label.bed2)
    Spn23F.label.bed$Channel <- factor(Spn23F.label.bed$Channel, levels = c("NAS", "Control"))
    p <- ggplot(sample.df, aes(x = (Start/1000000), y = norm.coverage, colour = Channel)) + facet_grid(.~Run) + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_smooth(span = 0.125, se = FALSE) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + geom_rect(data = Spn23F.label.bed, aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf), colour = NA, fill = "grey", alpha = 0.4) + geom_richtext(data = Spn23F.label.bed, aes(x = Start, label = label), colour="black", hjust = 1, size = 4, fill = NA, label.color = NA,) #+ scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    ggsave(file="CPS_V12_23F_operon_coverage.svg", plot=p, width = 15, height = 6)
  }
  else
  {
    p <- ggplot(sample.df, aes(x = (Start/1000000), y = norm.coverage, colour = Channel)) + facet_grid(Contaminant_species~Run) + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_smooth(span = 0.125, se = FALSE) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    p <- ggplot(sample.df, aes(x = (Start/1000000), y = norm.coverage, colour = Channel)) + facet_grid(Contaminant_species~Run) + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    p <- ggplot(subset(sample.df, Channel == "NAS"), aes(x = (Start/1000000), y = Coverage, colour = Run)) + facet_grid(Contaminant_species~., scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    p <- ggplot(subset(sample.df, Channel == "NAS"), aes(x = (Start/1000000), y = Coverage, colour = Run)) + facet_grid(Contaminant_species~., scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_smooth(span = 0.125, se = FALSE)  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    p <- ggplot(sample.df, aes(x = (Start/1000000), y = Coverage, colour = Channel)) + facet_grid(Contaminant_species~Run, scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    #p <- ggplot(sample.df, aes(x = (Start/1000000), y = norm.coverage, colour = Channel)) + facet_grid(Contaminant_species~Run, scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_line() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    #p
    
    adaptive.channels <- subset(sample.df, Channel == "NAS")
    control.channels <- subset(sample.df, Channel == "Control")
    
    sample.df <- merge(adaptive.channels, control.channels, by = c("Run", "Barcode", "Alignment", "Position", "Contaminant", "Contaminant_species", "Concentration"))
    sample.df$cov.diff <- sample.df$Coverage.x - sample.df$Coverage.y
    sample.df$norm.cov.diff <- sample.df$norm.coverage.x - sample.df$norm.coverage.y
    sample.df$cov.prop <- sample.df$Coverage.x / sample.df$Coverage.y
    sample.df$norm.cov.prop <- sample.df$norm.coverage.x / sample.df$norm.coverage.y
    sample.df[is.na(sample.df)] <- 0
    sample.df <- sample.df[order(sample.df$Position),]
    sample.df$Run <- factor(sample.df$Run, levels = c("CBL", "Whole Genome"))
    
    p <- ggplot(sample.df, aes(x = (Position/1000000), y = norm.cov.prop, colour = Run)) + facet_grid(Contaminant_species~., scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Channel")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    p <- ggplot(sample.df, aes(x = (Position/1000000), y = norm.cov.prop, colour = Run)) + facet_grid(Contaminant_species~., scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("NAS/Control coverage per Gb") + geom_smooth(span = 0.125, se = FALSE) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Target")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    p <- ggplot(sample.df, aes(x = (Position/1000000), y = norm.cov.diff, colour = Run)) + facet_grid(Contaminant_species~., scales = "free_y") + theme_light() + xlab("Position in Spn23F genome (Mb)") + ylab("NAS-Control coverage difference per Gb") + geom_smooth(span = 0.125, se = FALSE) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12, face="italic"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + guides(colour=guide_legend(title="Target")) + scale_color_npg() + scale_x_continuous(breaks = seq(0,10,0.2)) + annotate("rect", xmin = 303559/1000000, xmax = 322212/1000000, ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.6) + geom_hline(yintercept = 0, linetype="dashed", colour = "black", alpha=0.4)# + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
    p
    
    ggsave(file="V12_WGS_v_CPS_operon_normalised_coverage_0.1_dilution.svg", plot=p, width = 12, height = 6)
  }
}

# plot coverage
adaptive.channels <- subset(target.df, Channel == "adaptive")
#p <- ggplot(adaptive.channels, aes(x = Position, y = norm.coverage, group = Run, colour = Run)) + facet_grid(Concentration~Contaminant_species, scales = "free_y") + theme_light() + xlab("Position (bp)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + scale_color_npg() 
#p
#ggsave(file="CPS_V14_minimap2_v_graph_remclust2_adaptive_norm_coverage.svg", plot=p, width = 12, height = 8)

control.channels <- subset(target.df, Channel == "control")
#p <- ggplot(control.channels, aes(x = Position, y = norm.coverage, group = Run, colour = Run)) + facet_grid(Concentration~Contaminant_species, scales = "free_y") + theme_light() + xlab("Position (bp)") + ylab("Coverage per Gb") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + scale_color_npg() 
#p
#ggsave(file="CPS_V14_minimap2_v_graph_remclust2_control_norm_coverage.svg", plot=p, width = 12, height = 8)

# calculate difference between adaptive and control
merged <- merge(adaptive.channels, control.channels, by = c("Run", "Barcode", "Alignment", "Position", "Contaminant", "Contaminant_species", "Concentration"))
merged$cov.diff <- merged$Coverage.x - merged$Coverage.y
merged$norm.cov.diff <- merged$norm.coverage.x - merged$norm.coverage.y
merged$cov.prop <- merged$Coverage.x / merged$Coverage.y
merged$norm.cov.prop <- merged$norm.coverage.x / merged$norm.coverage.y

merged[is.na(merged)] <- 0
merged <- merged[order( merged[,1], merged[,2], merged[,3], merged[,4]),]

merged$Concentration <- signif(merged$Concentration, 1)
merged$Concentration[merged$Concentration == 0.00004] <- "`4x10`^`-5`"
merged$Concentration[merged$Concentration == 0.00008] <- "`8x10`^`-5`"
merged$Concentration[merged$Concentration == 0.0008] <- "`8x10`^`-4`"
merged$Concentration[merged$Concentration == 0.004] <- "`4x10`^`-3`"

merged$Contaminant <- as.character(merged$Contaminant)
merged$Contaminant[merged$Contaminant == "S. pneumoniae R6"] <- "italic(`S. pneumoniae`)"
merged$Contaminant[merged$Contaminant == "S. mitis"] <- "italic(`S. mitis`)"
merged$Contaminant[merged$Contaminant == "E. coli DH5a"] <- "italic(`E. coli`)"
merged$Contaminant <- factor(merged$Contaminant, levels = c("italic(`E. coli`)", "italic(`S. mitis`)", "italic(`S. pneumoniae`)"))
merged$Concentration <- factor(merged$Concentration, levels = c("`4x10`^`-5`", "`8x10`^`-5`", "`8x10`^`-4`", "`4x10`^`-3`"))


p <- ggplot(merged, aes(x = Position, y = norm.cov.diff, group = Run, colour = Run)) + facet_grid(Concentration~Contaminant, scales = "free_y", labeller = label_parsed) + theme_light() + xlab("Position (bp)") + ylab("NAS-Control coverage difference per Gb") + stat_smooth(geom='line', alpha=0.7, se=FALSE, linewidth=1, span=0.1) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 16)) + guides(colour=guide_legend(title="Aligner")) + scale_color_npg() 
p

#p <- ggplot(merged, aes(x = Position, y = norm.cov.prop, group = Run, colour = Run)) + facet_grid(Concentration~Contaminant, scales = "free_y", labeller = label_parsed) + theme_light() + xlab("Position (bp)") + ylab("NAS-Control coverage difference per Gb") + stat_smooth(geom='line', alpha=0.7, se=FALSE, linewidth=1, span=0.1) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 16)) + guides(colour=guide_legend(title="Aligner")) + scale_color_npg() 
#p
ggsave(file="CPS_V14_minimap2_v_graph_fulldb_norm_diff_coverage.svg", plot=p, width = 12, height = 6)

# for presentations
merged <- subset(merged, Run != "Graph k19 (S=90%)")
merged$Run <- as.character(merged$Run)
merged$Run[merged$Run == "Graph k19 (S=75%)"] <- "Graph"
merged$Run <- factor(merged$Run, levels = c("Minimap2", "Graph"))

p <- ggplot(merged, aes(x = Position, y = norm.cov.diff, group = Run, colour = Run)) + facet_grid(Concentration~Contaminant_species, scales = "free_y") + theme_light() + xlab("Position (bp)") + ylab("Adaptive-Control normalised coverage difference") + geom_path() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() 
p
ggsave(file="CPS_V14_minimap2_v_graph_remclust2_norm_diff_coverage_for_pres.svg", plot=p, width = 12, height = 8)
