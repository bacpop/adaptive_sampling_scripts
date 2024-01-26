library(hash)
library(openxlsx)
library(tidyverse)
library(stringr)
library(data.table)

indir.top = "./assemblies/Inspector_all_structural_errors_lq/"

bed.dirs <- Sys.glob(paste(indir.top, "*", sep = ""))

parse_filename <- function(file, indir)
{
  prefix <- str_remove(file, indir)
  
  params <- str_split(prefix, "/")[[1]]
  
  bc.details <- str_split(params[2], "_")[[1]]
  bed.details <- str_split(params[3], "_")[[1]]
  
  barcode <- bc.details[2]
  channel <- bc.details[4]
  type <- bed.details[1]
  
  array <- c(barcode, channel, type)
  
  array 
}

i <- 1
for (i in 1:length(bed.dirs))
{
  indir <- bed.dirs[i]
  bed_files <- Sys.glob(paste(indir, "*/*/summary_statistics", sep = ""))
  
  Experiment <- str_split(indir, "/")[[1]]
  Experiment <- Experiment[length(Experiment)]
  
  Experiment <- str_remove(Experiment, "_sup")
  
  # determine if WGS, CPS-multi or just 23F
  if (grepl("WGS", indir) == TRUE)
  {
    if (grepl("CPS_only", indir) == TRUE)
    {
      type = "CPS"
    }
    else{
      type = "WGS"
    }
  } else if (grepl("CPS", indir) == TRUE & grepl("size_selection", indir) == TRUE)
  {
    type = "multiCPS"
  } else {
    type = "CPS"
  }
  
  if (type == "WGS")
  {
    Alignment = "Spn23F"
  } else
  {
    Alignment = "23F CBL"
  }
  
  j <- 1
  for (j in 1:length(bed_files))
  {
    file <- bed_files[j]
    parse.file <- parse_filename(file, indir)
    curr.barcode <- parse.file[1]
    curr.channel <- parse.file[2]
    
    summary.df <- fread(file, sep = "\t", header = FALSE)
    colnames(summary.df) <- c("Statistic", "Value")
    
    
    summary.df$Barcode <- curr.barcode
    summary.df$Channel <- curr.channel
    summary.df$Experiment <- Experiment
    summary.df$Target <- Alignment
    
    
    # Experiment datasets
    if (Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75" | Experiment == "CPS_23F_fulldb_V14_400T_minimap2" | Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75" | Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90" | Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2")
    {
      summary.df.temp <- subset(summary.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
      summary.df.temp$Mixture <- NA
      summary.df.temp$Concentration <- 0
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode01"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode01"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode02"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode02"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode03"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode03"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode04"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode04"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode05"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode05"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode06"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode06"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode07"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode07"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode08"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode08"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode09"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode09"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode10"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode10"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode11"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode11"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode12"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode12"] <- 0.5
      
    } else if (Experiment == "23F_WGS_WGS_v_CPS" | Experiment == "23F_CPS_WGS_v_CPS_23F_only")
    {
      # change barcode Experiments
      summary.df.temp <- subset(summary.df, Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04" & Barcode != "barcode05" & Barcode != "barcode06" & Barcode != "barcode07" & Barcode != "barcode08" & Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" & Barcode != "barcode12")
      summary.df.temp$Mixture <- NA
      summary.df.temp$Concentration <- 0
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode13"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode13"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode14"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode14"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode15"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode15"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode16"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode16"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode17"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode17"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode18"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode18"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode19"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode19"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode20"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode20"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode21"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode21"] <- 0.005
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode22"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode22"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode23"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode23"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode24"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode24"] <- 0.5
      
      
    } else if (Experiment == "23F_WGS_w_size_selection" | Experiment == "23F_WGS_wo_size_selection")
    {
      summary.df.temp <- subset(summary.df, Barcode != "barcode01" & Barcode != "barcode02" & Barcode != "barcode03" & Barcode != "barcode04")
      summary.df.temp$Mixture <- NA
      summary.df.temp$Concentration <- 0
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode05"] <- "M. catarrhalis + H. influenzae  + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode05"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode06"] <- "M. catarrhalis + H. influenzae + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode06"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode07"] <- "M. catarrhalis + H. influenzae + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode07"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode08"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode08"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode09"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode09"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode10"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode10"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode11"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode11"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode12"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode12"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode13"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode13"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode14"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode14"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode15"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode15"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode16"] <- "S. mitis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode16"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode17"] <- "S. oralis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode17"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode18"] <- "S. oralis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode18"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode19"] <- "S. oralis + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode19"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode20"] <- "S. pneumoniae 110.58 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode20"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode21"] <- "S. pneumoniae 110.58 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode21"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode22"] <- "S. pneumoniae 110.58 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode22"] <- 0.01
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode23"] <- "E. coli DH5a + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode23"] <- 0.001
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode24"] <- "S. pneumoniae 110.58 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode24"] <- 0.001
      
      
      #summary.df.temp$Concentration <- summary.df.temp$Concentration * 100
      
      
    } else if (Experiment == "CPS_w_size_selection" | Experiment == "CPS_wo_size_selection")
    {
      h <- hash()
      
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
      
      
      summary.df.temp <- subset(summary.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
      
      barcodes <- unique(summary.df.temp$Barcode)
      
      summary.df.temp.temp <- data.frame(matrix(nrow = 0, ncol = ncol(summary.df.temp)))
      names(summary.df.temp.temp) <- names(summary.df.temp)
      
      
      # for (barcode in barcodes)
      # {
      #   subsample.summary <- subset(summary.df.temp, Barcode == barcode & Alignment %in% h[[barcode]])
      #   
      #   # ignore non-expected barcodes
      #   if (nrow(subsample.summary) != 0)
      #   {
      #     summary.df.temp.temp <- rbind(summary.df.temp.temp, subsample.summary)
      #   }
      # }
      
      summary.df.temp <- summary.df.temp.temp
      
      summary.df.temp$Mixture <- NA
      summary.df.temp$Concentration <- 0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode01"] <- "GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode01"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode02"] <- "GPSC16-19A"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode02"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode03"] <- "GPSC16-19F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode03"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode04"] <- "GPSC3-03"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode04"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode05"] <- "GPSC23-06B"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode05"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode06"] <- "GPSC1-19F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode06"] <- 1.0
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode07"] <- "GPSC16-19A + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode07"] <- 0.5
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode08"] <- "GPSC16-19F + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode08"] <- 0.5
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode09"] <- "GPSC3-03 + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode09"] <- 0.5
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode10"] <- "GPSC23-06B + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode10"] <- 0.5
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode11"] <- "GPSC1-19F + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode11"] <- 0.5
      
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode12"] <- "GPSC622 + GPSC16-23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode12"] <- 0.5
    } else if (Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75")
    {
      #remove specific barcodes
      summary.df.temp <- subset(summary.df, Barcode != "barcode07" & Barcode != "barcode08" &  Barcode != "barcode09" & Barcode != "barcode10" & Barcode != "barcode11" &  Barcode != "barcode12" &Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
      
      summary.df.temp$Mixture <- NA
      summary.df.temp$Concentration <- 0
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode01"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode01"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode02"] <- "S. pneumoniae R6 + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode02"] <- 0.5
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode03"] <- "Sample 1 (PCV-C-1464-1) + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode03"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode04"] <- "Sample 2 (PCV-C-0657-1) + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode04"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode05"] <- "Sample 3 (09B10326) + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode05"] <- 0.1
      summary.df.temp$Mixture[summary.df.temp$Barcode == "barcode06"] <- "Sample 4 (PCV-C-0720-1) + Spn23F"
      summary.df.temp$Concentration[summary.df.temp$Barcode == "barcode06"] <- 0.1
    }
    
    if (i > 1)
    {
      summary.df.final <- rbind(summary.df.final, summary.df.temp)
    } else {
      summary.df.final <- summary.df.temp
    }
  }
}

summary.df.final$Target <- NA
summary.df.final$Aligner <- NA
summary.df.final$Library <- NA
summary.df.final$Chemistry <- NA

# WGS vs. CBL CBL only
summary.df.final$Target[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_23F_only"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_23F_only"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_23F_only"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_23F_only"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_23F_only"] <- "WGS vs. CBL enrichment"

# CBL vs. WGS, WGS only
summary.df.final$Target[summary.df.final$Experiment == "23F_WGS_WGS_v_CPS"] <- "Spn23F Whole Genome"
summary.df.final$Aligner[summary.df.final$Experiment == "23F_WGS_WGS_v_CPS"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "23F_WGS_WGS_v_CPS"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "23F_WGS_WGS_v_CPS"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "23F_WGS_WGS_v_CPS"] <- "WGS vs. CBL enrichment"

# WGS size selection
summary.df.final$Target[summary.df.final$Experiment == "23F_WGS_w_size_selection"] <- "Spn23F Whole Genome"
summary.df.final$Aligner[summary.df.final$Experiment == "23F_WGS_w_size_selection"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "23F_WGS_w_size_selection"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "23F_WGS_w_size_selection"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "23F_WGS_w_size_selection"] <- "Spn23F whole genome enrichment"

# WGS without size selection
summary.df.final$Target[summary.df.final$Experiment == "23F_WGS_wo_size_selection"] <- "Spn23F Whole Genome"
summary.df.final$Aligner[summary.df.final$Experiment == "23F_WGS_wo_size_selection"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "23F_WGS_wo_size_selection"] <- "Unselected"
summary.df.final$Chemistry[summary.df.final$Experiment == "23F_WGS_wo_size_selection"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "23F_WGS_wo_size_selection"] <- "Spn23F whole genome enrichment"

# CBL with size selection
summary.df.final$Target[summary.df.final$Experiment == "CPS_w_size_selection"] <- "Multi-CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_w_size_selection"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "CPS_w_size_selection"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_w_size_selection"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_w_size_selection"] <- "CBL enrichment"

# CBL without size selection
summary.df.final$Target[summary.df.final$Experiment == "CPS_wo_size_selection"] <- "Multi-CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_wo_size_selection"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "CPS_wo_size_selection"] <- "Unselected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_wo_size_selection"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_wo_size_selection"] <- "CBL enrichment"

# CBL full db graph p75
summary.df.final$Target[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "GP (k=19, S=75%, min. read=50 bp)"
summary.df.final$Library[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_graphk19_p75"] <- "Full CBL database"

# CBL full db minimap
summary.df.final$Target[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_minimap2"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_minimap2"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_minimap2"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_minimap2"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_23F_fulldb_V14_400T_minimap2"] <- "Full CBL database"

#CBL partial database minimap2
summary.df.final$Target[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_minimap2_rep2"] <- "Partial CBL database"

#CBL partial database graph s90
summary.df.final$Target[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "GP (k=19, S=90%, min. read=50 bp)"
summary.df.final$Library[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p90"] <- "Partial CBL database"

#CBL partial database graph s75
summary.df.final$Target[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "GP (k=19, S=75%, min. read=50 bp)"
summary.df.final$Library[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "Partial CBL database"

#Mixed culture full database graph s75
summary.df.final$Target[summary.df.final$Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75"] <- "23F CBL"
summary.df.final$Aligner[summary.df.final$Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75"] <- "GP (k=19, S=75%, min. read=50 bp)"
summary.df.final$Library[summary.df.final$Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75"] <- "Unselected"
summary.df.final$Chemistry[summary.df.final$Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75"] <- "V14"
summary.df.final$Experiment[summary.df.final$Experiment == "Mixed_23F_fulldb_V14_400T_graphk19_p75"] <- "Mixed culture Full CBL database"

# # CBL graph vs. minimap
# summary.df.final$Target[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50"] <- "CBL"
# summary.df.final$Aligner[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50"] <- "GP (k=31, S=90%, min. read=50 bp)"
# summary.df.final$Library[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50"] <- "Size-selected"
# summary.df.final$Chemistry[summary.df.final$Experiment == "CPS_23F_remclust2_V14_400T_graphk19_p75"] <- "V14"
# summary.df.final$Experiment[summary.df.final$Experiment == "23F_CPS_WGS_v_CPS_graph_k31_p90_c50"] <- "CBL vs WGS and Minimap2 vs Graph (V12)"

summary.df.final$Alignment[summary.df.final$Alignment == "FM211187.1"] <- "Spn23F"

summary.df.final$Experiment <- factor(summary.df.final$Experiment, levels = c("Spn23F whole genome enrichment", "CBL enrichment", "WGS vs. CBL enrichment", "Partial CBL database", "Full CBL database", "Mixed culture Full CBL database"))
summary.df.final$Library <- factor(summary.df.final$Library, levels = c("Unselected", "Size-selected"))
summary.df.final$Target <- factor(summary.df.final$Target, levels = c("Spn23F Whole Genome", "23F CBL", "Multi-CBL" ))
summary.df.final$Aligner <- factor(summary.df.final$Aligner, levels = c("Minimap2", "GP (k=19, S=75%, min. read=50 bp)", "GP (k=19, S=90%, min. read=50 bp)"))
summary.df.final$Alignment <- factor(summary.df.final$Alignment, levels = c("03", "06B", "19A", "19F", "23F", "Spn23F", "unaligned", "Total"))

summary.df.final$Concentration <- summary.df.final$Concentration * 100
operon.length <- 18654
genome.length <- 2221315
perc.genome <- operon.length / genome.length

summary.df.final$Concentration[summary.df.final$Target == "23F CBL"] <- summary.df.final$Concentration[summary.df.final$Target == "23F CBL"] * perc.genome
summary.df.final$Concentration <- signif(summary.df.final$Concentration, 1)
summary.df.final$Channel[summary.df.final$Channel == "adaptive"] <- "Adaptive"
summary.df.final$Channel[summary.df.final$Channel == "control"] <- "Control"

#reorder
summary.df.final <- summary.df.final[order(summary.df.final$Experiment, summary.df.final$Barcode),]
summary.df.final <- summary.df.final[,c(5, 7, 8, 6, 9, 10, 4, 3, 1, 2)]

# correct
summary.df.final$Statistic[grepl("Reference base with Depth=0 (including Ns)*", summary.df.final$Statistic)] <- "Reference base with Depth=0 (including Ns)"
summary.df.final$Statistic[grepl("Reference base with Depth=1 *", summary.df.final$Statistic)] <- "Reference base with Depth=1"
summary.df.final$Statistic[grepl("Reference base with Depth=2 *", summary.df.final$Statistic)] <- "Reference base with Depth=2"
summary.df.final$Statistic[grepl("Reference base with Depth>2 *", summary.df.final$Statistic)] <- "Reference base with Depth>2"
summary.df.final$Statistic[summary.df.final$Statistic == "Depth in large conigs"] <- "Depth in large contigs"

names(summary.df.final)[names(summary.df.final) == "Concentration"] <- "Concentration_perc"

contiguity <- subset(summary.df.final, Statistic == "Number of contigs" | Statistic == "Number of contigs > 10000 bp" | Statistic == "Number of contigs >1000000 bp" | Statistic == "Total length" | Statistic == "Total length of contigs > 10000 bp" | Statistic == "Total length of contigs >1000000bp" | Statistic == "Longest contig" | Statistic == "Second longest contig length" | Statistic == "N50" | Statistic == "N50 of contigs >1Mbp")
read.alignment <- subset(summary.df.final, Statistic == "Mapping rate /%" | Statistic == "Split-read rate /%" | Statistic == "Depth" | Statistic == "Mapping rate in large contigs /%" | Statistic == "Split-read rate in large contigs /%" | Statistic == "Depth in large contigs" | Statistic == "Structural error" | Statistic == "Expansion" | Statistic == "Collapse", Statistic == "Haplotype switch" | Statistic == "Inversion" | Statistic == "Small-scale assembly error /per Mbp", Statistic == "Total small-scale assembly error" | Statistic == "Base substitution" | Statistic == "Small-scale expansion" | Statistic == "Small-scale collapse" | Statistic == "QV")
ref.alignment <- subset(summary.df.final, Statistic == "Genome Coverage /%" | Statistic == "Reference base with Depth=0 (including Ns)" | Statistic == "Reference base with Depth=1" | Statistic == "Reference base with Depth=2" | Statistic == "Reference base with Depth>2" | Statistic == "Assembly contig mapping ratio (length) /%" | Statistic == "Assembly contig NA50" | Statistic == "Number of assembly collapse" | Statistic == "Number of assembly expansion" | Statistic == "Number of assembly inversion")

wb <- createWorkbook()
addWorksheet(wb, "Contiguity")
addWorksheet(wb, "Read_to_assembly_alignment")
addWorksheet(wb, "Assembly_to_reference_alignment")

writeData(wb, "Contiguity", contiguity, startRow = 1, startCol = 1)
writeData(wb, "Read_to_assembly_alignment", read.alignment, startRow = 1, startCol = 1)
writeData(wb, "Assembly_to_reference_alignment", ref.alignment, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "assembly_metadata.xlsx", overwrite = TRUE)

