library(ggplot2)
library(purrr)
library(ggpubr)
library(ggsci)
library(stringr)
library(dplyr)
library(ggbeeswarm)
library(hash)

indir = "./alignment_params/alignment_only/"
length_files <- Sys.glob(paste(indir,"*.txt", sep = ""))
params <- gsub("_summary.txt", "", gsub(".*/", "", length_files))
target <- c("23F", "FM211187.1")

target.df <- data.frame(Run = c(), Statistic = c(), Channel = c(), Barcode = c(), pid = c(), Alignment = c(), Value = c())


h <- hash()
{
  # set up barcodes
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


for (i in 1:length(length_files))
{
  filename <- length_files[i]
  df <- read.table(filename, sep = "\t", comment.char = "", header = 1)
  size_selection <- "Original"
  pid <- 0
  
  if (!grepl("_wo_", params[i]))
  {
    size_selection <- "Size selection"
  }
  
  names = str_split(params[i], "_p")[[1]]
  perc_id <- as.integer(names[[2]])
  
  # subset data, only get enrichment for desired sequences
  barcodes <- unique(df$Barcode)
  
  for (barcode in barcodes)
  {
    # ignore non-expected barcodes
    if (!is.null(h[[barcode]]))
    {
      subsample <- subset(df, Barcode == barcode & Alignment %in% h[[barcode]])
      
      if (nrow(subsample) == 0)
      {
        stats <- c("Reads_mapped", "Bases_mapped","Reads_mapped", "Bases_mapped", "Reads_mapped", "Bases_mapped","Reads_mapped", "Bases_mapped", "Enrichment")
        channels <- c("Target", "Target", "Target", "Target", "Non-target", "Non-target", "Non-target", "Non-target", NA)
        to.append <- data.frame(Run = size_selection, pid = perc_id, Statistic = stats, Channel = channels, Barcode = barcode, Alignment = NA, Value =0)
      } else {
        to.append <- data.frame(Run = size_selection, pid = perc_id, Statistic = subsample$Statistic, Channel = subsample$Channel, Barcode = subsample$Barcode, Alignment = subsample$Alignment, Value = subsample$Value)
      }
      
      target.df <- rbind(target.df, to.append)
      
      # append total stats
      subsample <- subset(df, Barcode == barcode & Alignment == "Total")
      to.append <- data.frame(Run = size_selection, pid = perc_id, Statistic = subsample$Statistic, Channel = subsample$Channel, Barcode = subsample$Barcode, Alignment = subsample$Alignment, Value = subsample$Value)
      target.df <- rbind(target.df, to.append)
    }
  }
}

# change barcode names
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
  
  target.df$Concentration <- target.df$Concentration * 100
  
  target.df$Contaminant_species[target.df$Contaminant == "E. coli DH5a" | target.df$Contaminant == "M. catarrhalis + H. influenzae"] <- "Non-Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. oralis" | target.df$Contaminant == "S. mitis"] <- "Streptococcus"
  target.df$Contaminant_species[target.df$Contaminant == "S. pneumoniae 110.58" | target.df$Contaminant == "S. pneumoniae R6"] <- "Streptococcus pneumoniae"
  
  target.df$Run <- factor(target.df$Run, levels = c("Original", "Size selection"))
  target.df$Contaminant <- factor(target.df$Contaminant, levels = c("E. coli DH5a", "M. catarrhalis + H. influenzae", "S. oralis", "S. mitis", "S. pneumoniae 110.58", "S. pneumoniae R6"))
  target.df$Contaminant_species <- factor(target.df$Contaminant_species, levels = c("Non-Streptococcus", "Streptococcus", "Streptococcus pneumoniae"))
}

sampled.data <- subset(target.df, Statistic == "Enrichment")
sampled.data <- sampled.data[!is.infinite(sampled.data$Value),]
sampled.data <- subset(sampled.data, Run == "Original")

min.max.df <- sampled.data %>%
  group_by(pid) %>%
  summarise(min = min(Value), max = max(Value))

p <- ggplot(sampled.data, aes(x = pid, y = Value, group = interaction(Run, Contaminant), colour = Contaminant)) + facet_grid(Contaminant_species~Concentration) + theme_light() + xlab("Minimum alignment identity (%)") + ylab("Enrichment") + geom_hline(yintercept = 1, linetype = "dashed") + geom_line() + geom_point() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), strip.text.x = element_text(size = 12)) + scale_color_npg() + scale_x_continuous(breaks = seq(0,100,4)) + labs(color = "Contaminant(s) in library") + scale_y_log10()
p
ggsave(file="alignment_param_comp.svg", plot=p, width=12, height=8)

p <- ggplot(sampled.data, aes(x = as.factor(pid), y = Value, colour = "red")) + theme_light() + xlab("Minimum alignment identity (%)") + ylab("Enrichment") + geom_boxplot() + geom_hline(yintercept = 1, linetype = "dashed") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.position = "none", strip.text.x = element_text(size = 12)) + scale_color_npg() + scale_y_log10() #+ stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.2) + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 4, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE)
p
ggsave(file="alignment_param_stats.svg", plot=p)
