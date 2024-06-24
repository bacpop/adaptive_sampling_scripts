library(ggplot2)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggtext)
library(dplyr)
library(ggbeeswarm)

parse_filename <- function(indir)
{
  experiment <- str_split(indir, "/")[[1]]
  run <- experiment[length(experiment) - 2]
  experiment <- experiment[length(experiment) - 1]
  
  params <- str_split(experiment, "_")[[1]]
  
  barcode <- params[2]
  channel <- if(params[4] == "adaptive") "Adaptive" else "Control"
  min.kmer <- as.numeric(gsub("[^0-9.-]", "", params[7]))
  min.percentage <- as.numeric(gsub("[^0-9.-]", "", params[8]))
  
  array <- c(barcode, channel, min.kmer, min.percentage, run)
  
  array 
}

# read in files
#indir = "./pneumokitty/remclust2/"
#experiment = "V14_CPS_graph_remclust2"

indir = "./pneumokitty/Mixed_V14/"
experiment = "Mixed_V14"

analyse.dirs <- Sys.glob(paste(indir, "*/*/pneumo_capsular_typing", sep = ""))

top.n.rows <- 200

i <- 1
for (i in 1:length(analyse.dirs)){
  indir <- analyse.dirs[i]
  #iterate over bed files, identify number of barcodes
  all.data <- Sys.glob(paste(indir, "*/*_alldata.csv", sep = ""))
  #result.data <- Sys.glob(paste(indir, "*/*_result_data.csv", sep = ""))
  mixed.data <- Sys.glob(paste(indir, "*/*_mixed_serotypes.csv", sep = ""))
  
  mixed.present <- if(length(mixed.data) == 0) FALSE else TRUE
  valid <- if(length(all.data) == 0) FALSE else TRUE
  
  if (valid == FALSE)
  {
    next
  }
  
  all.data.df <- read.table(all.data, sep = ",", comment.char = "", header = 1)
  #result.data.df <- read.table(result.data, sep = ",", comment.char = "", header = 1)
  
  # if (mixed.present == TRUE)
  # {
  #   mixed.data.df <- read.table(mixed.data, sep = ",", comment.char = "", header = 1)
  # }
  
  parse.file <- parse_filename(indir)
  
  top.rows <- top_n(all.data.df, top.n.rows, percent)
  top.rows$Barcode <- parse.file[1]
  top.rows$Channel <- parse.file[2]
  top.rows$Min_kmer <- parse.file[3]
  top.rows$min.percentage <- parse.file[4]
  top.rows$Run <- parse.file[5]
  
  if (i == 1)
  {
    total.df <- top.rows
  } else
  {
    total.df <- rbind(total.df, top.rows)
  }
}


# rename datasets
if (experiment == "V14_CPS_graph_fulldb" | experiment == "V14_CPS_graph_remclust2" | experiment == "V12_CPS_graph" | experiment == "V12_WGS_v_CPS")
{
  if (experiment == "V14_CPS_graph_remclust2" | experiment == "V14_CPS_graph_fulldb")
  {
    total.df <- subset(total.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
    total.df$Contaminant <- NA
    total.df$Contaminant_species <- NA
    total.df$Concentration <- 0
    total.df$Contaminant <- NA
    total.df$Contaminant_species <- NA
    total.df$Concentration <- 0
    total.df$Contaminant[total.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode01"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode02"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode03"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode03"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode04"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode04"] <- 0.5
    total.df$Contaminant[total.df$Barcode == "barcode05"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode05"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode06"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode06"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode07"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode07"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode08"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode08"] <- 0.5
    total.df$Contaminant[total.df$Barcode == "barcode09"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode09"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode10"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode10"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode11"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode11"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode12"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode12"] <- 0.5
    
    if (experiment == "V14_CPS_graph_remclust2")
    {
      total.df$Run[total.df$Run == "CPS_23F_remclust2_V14_400T_minimap2_rep2_aligned"] <- "Minimap2"
      total.df$Run[total.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p90_aligned"] <- "Graph k19 (S=90%)"
      total.df$Run[total.df$Run == "CPS_23F_remclust2_V14_400T_graphk19_p75_aligned"] <- "Graph k19 (S=75%)"
      
      total.df$Run <- factor(total.df$Run, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    } else if (experiment == "V14_CPS_graph_fulldb")
    {
      total.df$Run[total.df$Run == "CPS_23F_fulldb_V14_400T_minimap2_all"] <- "Minimap2"
      total.df$Run[total.df$Run == "CPS_23F_fulldb_V14_400T_graphk19_p75_all"] <- "Graph k19 (S=75%)"
      
      total.df$Run <- factor(total.df$Run, levels = c("Minimap2", "Graph k19 (S=75%)"))
    }
  } else if (experiment == "V12_CPS_graph" | experiment == "V12_WGS_v_CPS")
  {
    # change barcode names
    total.df <- subset(total.df, Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12")
    total.df$Contaminant <- NA
    total.df$Contaminant_species <- NA
    total.df$Concentration <- 0
    total.df$Contaminant <- NA
    total.df$Contaminant_species <- NA
    total.df$Concentration <- 0
    total.df$Contaminant[total.df$Barcode == "barcode13"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode13"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode14"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode14"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode15"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode15"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode16"] <- "S. pneumoniae R6"
    total.df$Concentration[total.df$Barcode == "barcode16"] <- 0.5
    total.df$Contaminant[total.df$Barcode == "barcode17"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode17"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode18"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode18"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode19"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode19"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode20"] <- "S. mitis"
    total.df$Concentration[total.df$Barcode == "barcode20"] <- 0.5
    total.df$Contaminant[total.df$Barcode == "barcode21"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode21"] <- 0.005
    total.df$Contaminant[total.df$Barcode == "barcode22"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode22"] <- 0.01
    total.df$Contaminant[total.df$Barcode == "barcode23"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode23"] <- 0.1
    total.df$Contaminant[total.df$Barcode == "barcode24"] <- "E. coli DH5a"
    total.df$Concentration[total.df$Barcode == "barcode24"] <- 0.5
    
    if (experiment == "V12_CPS_graph")
    {
      total.df$Run[total.df$Run == "23F_CPS_WGS_v_CPS_23F_only_all"] <- "Minimap2"
      total.df$Run[total.df$Run == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50_all"] <- "Graph k31 (S=90%)"
      
      total.df$Run <- factor(total.df$Run, levels = c("Minimap2", "Graph k31 (S=90%)"))
    } else if (experiment == "V12_WGS_v_CPS")
    {
      total.df$Run[total.df$Run == "23F_CPS_WGS_v_CPS_23F_only_all"] <- "CBL"
      total.df$Run[total.df$Run == "23F_WGS_WGS_v_CPS_all"] <- "Whole Genome"
      
      total.df$Run <- factor(total.df$Run, levels = c("Whole Genome", "CBL"))
    }
  }
  
  #total.df$Concentration <- total.df$Concentration * 100
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  total.df$Contaminant_species[total.df$Contaminant == "E. coli DH5a"] <- "Non-Streptococcus"
  total.df$Contaminant_species[total.df$Contaminant == "S. mitis"] <- "Streptococcus"
  total.df$Contaminant_species[total.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  total.df$Concentration[total.df$Run != "Whole Genome"] <- total.df$Concentration[total.df$Run != "Whole Genome"] * perc.genome
  
  total.df$Contaminant <- factor(total.df$Contaminant, levels = c("E. coli DH5a", "S. mitis", "S. pneumoniae R6"))
  total.df$Contaminant_species <- factor(total.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
} else if (experiment == "Mixed_V14")
{
  
  #remove specific barcodes
  total.df <- subset(total.df, Barcode != "barcode07" & Barcode != "barcode08" &  Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" &  Barcode != "barcode12" &Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
  
  total.df$Contaminant <- NA
  total.df$Contaminant_species <- NA
  total.df$Concentration <- 0
  total.df$Contaminant <- NA
  total.df$Contaminant_species <- NA
  total.df$Concentration <- 0
  total.df$Contaminant[total.df$Barcode == "barcode01"] <- "S. pneumoniae R6"
  total.df$Concentration[total.df$Barcode == "barcode01"] <- 0.1
  total.df$Contaminant[total.df$Barcode == "barcode02"] <- "S. pneumoniae R6"
  total.df$Concentration[total.df$Barcode == "barcode02"] <- 0.5
  total.df$Contaminant[total.df$Barcode == "barcode03"] <- "PCV-C-1464-1"
  total.df$Concentration[total.df$Barcode == "barcode03"] <- 0.1
  total.df$Contaminant[total.df$Barcode == "barcode04"] <- "PCV-C-0657-1"
  total.df$Concentration[total.df$Barcode == "barcode04"] <- 0.1
  total.df$Contaminant[total.df$Barcode == "barcode05"] <- "09B10326"
  total.df$Concentration[total.df$Barcode == "barcode05"] <- 0.1
  total.df$Contaminant[total.df$Barcode == "barcode06"] <- "PCV-C-0720-1"
  total.df$Concentration[total.df$Barcode == "barcode06"] <- 0.1
  
  total.df$Run[total.df$Run == "Mixed_23F_fulldb_V14_400T_graphk19_p75_aligned"] <- "Graph k19 (S=75%)"
  
  total.df$Run <- factor(total.df$Run, levels = c("Graph k19 (S=75%)"))
  
  operon.length <- 18654
  genome.length <- 2221315
  perc.genome <- operon.length / genome.length
  
  total.df$Contaminant_species[total.df$Contaminant == "PCV-C-1464-1"] <- "`Sample 1`"
  total.df$Contaminant_species[total.df$Contaminant == "PCV-C-0657-1"] <- "`Sample 2`"
  total.df$Contaminant_species[total.df$Contaminant == "09B10326"] <- "`Sample 3`"
  total.df$Contaminant_species[total.df$Contaminant == "PCV-C-0720-1"] <- "`Sample 4`"
  total.df$Contaminant_species[total.df$Contaminant == "S. pneumoniae R6"] <- "`Pneumo R6`"
  total.df$Contaminant_species <- factor(total.df$Contaminant_species, levels = c("`Sample 1`", "`Sample 2`", "`Sample 3`", "`Sample 4`", "`Pneumo R6`"))
  
  total.df$Concentration[total.df$Run != "Whole Genome"] <- total.df$Concentration[total.df$Run != "Whole Genome"] * perc.genome
}

total.df$highlight <- ifelse(total.df$Serotype == "23F", "23F", "other")
#total.df.subset <- subset(total.df, Channel == "Adaptive")
total.df.subset <- total.df
#total.df.subset$highlight <- factor(total.df.subset$highlight, levels = c("highlight", "normal"))
total.df.subset$Concentration <- signif(total.df.subset$Concentration, 1)
total.df.subset$Concentration[total.df.subset$Concentration == 0.00004] <- "`4x10`^`-5`"
total.df.subset$Concentration[total.df.subset$Concentration == 0.00008] <- "`8x10`^`-5`"
total.df.subset$Concentration[total.df.subset$Concentration == 0.0008] <- "`8x10`^`-4`"
total.df.subset$Concentration[total.df.subset$Concentration == 0.004] <- "`4x10`^`-3`"

# remove any derivative serotypes
total.df.subset <- total.df.subset[!grepl("-", total.df.subset$Serotype),]

total.df.subset <- total.df.subset[order(total.df.subset$highlight, decreasing = TRUE),]
colour_vec <- c("23F" = '#FF6666', "other" = "#717171")

if (experiment != "Mixed_V14")
{
  total.df.subset$Contaminant <- as.character(total.df.subset$Contaminant)
  total.df.subset$Contaminant[total.df.subset$Contaminant == "S. pneumoniae R6"] <- "italic(`S. pneumoniae`)"
  total.df.subset$Contaminant[total.df.subset$Contaminant == "S. mitis"] <- "italic(`S. mitis`)"
  total.df.subset$Contaminant[total.df.subset$Contaminant == "E. coli DH5a"] <- "italic(`E. coli`)"
  total.df.subset$Contaminant <- factor(total.df.subset$Contaminant, levels = c("italic(`E. coli`)", "italic(`S. mitis`)", "italic(`S. pneumoniae`)"))
  total.df.subset$Concentration <- factor(total.df.subset$Concentration, levels = c("`4x10`^`-5`", "`8x10`^`-5`", "`8x10`^`-4`", "`4x10`^`-3`"))
  
  p <- ggplot(total.df.subset, aes(x = Run, y = percent, colour = factor(highlight), shape=Channel, group=Channel)) + geom_point(alpha=1, size=2, position = position_dodge(width = 0.5)) + theme_light() + facet_grid(Contaminant~Concentration, scales = "free_y", labeller = label_parsed) + xlab("Alignment method") + ylab("Proportion of reference matched (%)") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + geom_hline(yintercept = 70, linetype="dashed", colour = "grey") + guides(colour=guide_legend(title="Serotype"), shape=guide_legend(title="Channel")) + scale_colour_manual(values = colour_vec)
  p
  ggsave(file="CPS_remclust2_pnuemoKITy.svg", plot=p, height = 8, width = 12)
} else
{
  total.df.subset$Concentration <- factor(total.df.subset$Concentration, levels = c("`8x10`^`-4`", "`4x10`^`-3`"))
  predicted.serotypes <- subset(total.df.subset, percent >= 70 & Channel == "Adaptive" & Barcode != "barcode01" & Barcode != "barcode02")
  p <- ggplot(total.df.subset, aes(x = Concentration, y = percent, colour = factor(highlight), shape=Channel, group=Channel)) + geom_point(alpha=1, size=2, position = position_dodge(width = 0.5)) + theme_light() + facet_wrap(.~Contaminant_species, scales = "free_x", labeller = label_parsed, ncol=5) + xlab("Proportion of total DNA targeted for enrichment") + ylab("Proportion of reference matched (%)") + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + geom_hline(yintercept = 70, linetype="dashed", colour = "grey") + guides(colour=guide_legend(title="Serotype"), shape=guide_legend(title="Channel")) + scale_colour_manual(values = colour_vec) + scale_x_discrete(labels=parse(text=levels(total.df.subset$Concentration)))
  p
  ggsave(file="Mixed_V14_pnuemoKITy.svg", plot=p, height = 8, width = 12)

  #for presentation
  total.df.subset.pres <- subset(total.df.subset, Concentration !="`4x10`^`-3`" & Contaminant_species != "`Pneumo R6`")
  total.df.subset.pres$Channel <- as.character(total.df.subset.pres$Channel)
  total.df.subset.pres$Channel[total.df.subset.pres$Channel == "Adaptive"] <- "NAS"
  total.df.subset.pres$Channel <- factor(total.df.subset.pres$Channel, levels = c("NAS", "Control"))
  p <- ggplot(total.df.subset.pres, aes(x = Channel, y = percent, colour = factor(highlight), group=Channel)) + geom_point(alpha=1, size=4, position = position_dodge(width = 0.5)) + theme_light() + facet_wrap(.~Contaminant_species, scales = "free_x", labeller = label_parsed, ncol=4) + xlab("Channel") + ylab("Proportion of reference matched (%)") + theme(axis.text.x = element_text(size = 16, angle = 45, hjust=1), axis.text.y = element_text(size = 14), axis.title=element_text(size=20,face="bold"), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14, face = "italic"), legend.title=element_text(size=18,face="bold"), legend.text=element_text(size=16)) + geom_hline(yintercept = 70, linetype="dashed", colour = "grey") + guides(colour=guide_legend(title="Serotype")) + scale_colour_manual(values = colour_vec) + scale_y_continuous(limits = c(50, 100))
  p
  ggsave(file="Mixed_V14_pnuemoKITy_for_pres.svg", plot=p, height = 8, width = 12)
}



