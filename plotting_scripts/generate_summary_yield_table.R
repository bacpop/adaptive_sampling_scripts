library(hash)
library(openxlsx)
library(tidyverse)

indir = "./enrichment_analysis/bootstrapped_enrichment_alignment_only/all_summaries_manuscript/"
summary_files <- Sys.glob(paste(indir,"*_summary.txt", sep = ""))


for (i in 1:length(summary_files))
{
  summary.Experiment <- summary_files[i]
  summary.df <- read.table(summary.Experiment, sep = "\t", comment.char = "", header = 1)
  
  Experiment <- gsub("_all_summary.txt", "", gsub(".*/", "", summary.Experiment))
  summary.df$Experiment <- Experiment
  
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
    
    
  } else if (Experiment == "23F_WGS_size_selection" | Experiment == "23F_WGS_wo_size_selection")
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
    
    h[["barcode01"]] <- c("23F", "unaligned", "Total")
    h[["barcode02"]] <- c("19A", "unaligned", "Total")
    h[["barcode03"]] <- c("19F", "unaligned", "Total")
    h[["barcode04"]] <- c("03", "unaligned", "Total")
    h[["barcode05"]] <- c("06B", "unaligned", "Total")
    h[["barcode06"]] <- c("19F", "unaligned", "Total")
    h[["barcode07"]] <- c("23F", "19A", "unaligned", "Total")
    h[["barcode08"]] <- c("23F", "19F", "unaligned", "Total")
    h[["barcode09"]] <- c("23F", "03", "unaligned", "Total")
    h[["barcode10"]] <- c("23F", "06B", "unaligned", "Total")
    h[["barcode11"]] <- c("23F", "19F", "unaligned", "Total")
    h[["barcode12"]] <- c("23F", "unaligned", "Total")
    
    
    summary.df.temp <- subset(summary.df, Barcode != "barcode13" & Barcode != "barcode14" & Barcode != "barcode15" & Barcode != "barcode16" & Barcode != "barcode17" & Barcode != "barcode18" & Barcode != "barcode19" & Barcode != "barcode20" & Barcode != "barcode21" & Barcode != "barcode22" & Barcode != "barcode23" & Barcode != "barcode24")
    
    barcodes <- unique(summary.df.temp$Barcode)
    
    summary.df.temp.temp <- data.frame(matrix(nrow = 0, ncol = ncol(summary.df.temp)))
    names(summary.df.temp.temp) <- names(summary.df.temp)
    
    
    for (barcode in barcodes)
    {
      subsample.summary <- subset(summary.df.temp, Barcode == barcode & Alignment %in% h[[barcode]])
      
      # ignore non-expected barcodes
      if (nrow(subsample.summary) != 0)
      {
        summary.df.temp.temp <- rbind(summary.df.temp.temp, subsample.summary)
      }
    }
    
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

test <- subset(summary.df.final, Experiment == "CPS_wo_size_selection")

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
summary.df.final$Target[summary.df.final$Experiment == "23F_WGS_size_selection"] <- "Spn23F Whole Genome"
summary.df.final$Aligner[summary.df.final$Experiment == "23F_WGS_size_selection"] <- "Minimap2"
summary.df.final$Library[summary.df.final$Experiment == "23F_WGS_size_selection"] <- "Size-selected"
summary.df.final$Chemistry[summary.df.final$Experiment == "23F_WGS_size_selection"] <- "V12"
summary.df.final$Experiment[summary.df.final$Experiment == "23F_WGS_size_selection"] <- "Spn23F whole genome enrichment"

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
summary.df.final$Concentration[summary.df.final$Target == "Multi-CBL"] <- summary.df.final$Concentration[summary.df.final$Target == "23F CBL"] * perc.genome
summary.df.final$Concentration <- signif(summary.df.final$Concentration, 1)
summary.df.final$Channel[summary.df.final$Channel == "Target"] <- "Adaptive"
summary.df.final$Channel[summary.df.final$Channel == "Non-target"] <- "Control"
summary.df.final$Statistic <- gsub("_", " ", summary.df.final$Statistic)

names(summary.df.final)[names(summary.df.final) == "Concentration"] <- "Concentration_perc"

# reorder tables
bases.mapped <- subset(summary.df.final, summary.df.final$Statistic == "Bases mapped")
bases.mapped <- bases.mapped[order(bases.mapped[,6], bases.mapped[,11], bases.mapped[,1],  bases.mapped[,3], bases.mapped[,4], bases.mapped[,2]),]
bases.mapped <- bases.mapped[,c(6, 7, 8, 9, 10, 11, 2, 3, 4, 5)]
bases.mapped <- bases.mapped %>% 
  rename(
    Bases_mapped = Value
  )

bases.total <- subset(summary.df.final, summary.df.final$Statistic == "Bases total")
bases.total <- bases.total[order(bases.total[,6], bases.total[,11], bases.total[,1],  bases.total[,3], bases.total[,4], bases.total[,2]),]
bases.total <- bases.total[,c(6, 7, 8, 9, 10, 11, 2, 3, 4, 5)]
bases.total <- bases.total %>% 
  rename(
    Bases_total = Value
  )

read.mapped <- subset(summary.df.final, summary.df.final$Statistic == "Reads mapped")
read.mapped <- read.mapped[order(read.mapped[,6], read.mapped[,11], read.mapped[,1],  read.mapped[,3], read.mapped[,4], read.mapped[,2]),]
read.mapped <- read.mapped[,c(6, 7, 8, 9, 10, 11, 2, 3, 4, 5)]
read.mapped <- read.mapped %>% 
  rename(
    Reads_mapped = Value
  )

read.total <- subset(summary.df.final, summary.df.final$Statistic == "Reads total")
read.total <- read.total[order(read.total[,6], read.total[,11], read.total[,1],  read.total[,3], read.total[,4], read.total[,2]),]
read.total <- read.total[,c(6, 7, 8, 9, 10, 11, 2, 3, 4, 5)]
read.total <- read.total %>% 
  rename(
    Reads_total = Value
  )

enrichment <- subset(summary.df.final, summary.df.final$Statistic == "Enrichment")
enrichment <- enrichment[order(enrichment[,6], enrichment[,11], enrichment[,1],  enrichment[,3], enrichment[,4], enrichment[,2]),]
enrichment <- enrichment[,c(6, 7, 8, 9, 10, 11, 3, 4, 5)]
enrichment <- enrichment %>% 
  rename(
    Enrichment = Value
  )

wb <- createWorkbook()
addWorksheet(wb, "bases_mapped")
addWorksheet(wb, "bases_total")
addWorksheet(wb, "reads_mapped")
addWorksheet(wb, "reads_total")
addWorksheet(wb, "enrichment")

writeData(wb, "bases_mapped", bases.mapped, startRow = 1, startCol = 1)
writeData(wb, "bases_total", bases.total, startRow = 1, startCol = 1)
writeData(wb, "reads_mapped", read.mapped, startRow = 1, startCol = 1)
writeData(wb, "reads_total", read.total, startRow = 1, startCol = 1)
writeData(wb, "enrichment", enrichment, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "run_metadata.xlsx", overwrite = TRUE)

# write.csv(bases.mapped, "NAS_bases_mapped_summary.csv", row.names=FALSE)
# write.csv(read.mapped, "NAS_reads_mapped_summary.csv", row.names=FALSE)
# write.csv(enrichment, "NAS_enrichment_summary.csv", row.names=FALSE)
