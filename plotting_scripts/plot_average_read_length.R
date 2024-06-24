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
library(openxlsx)

filename <- "Supplementary_File_1_run_metadata.xlsx"

all.mapped.df <- read.xlsx(filename, "read_len_mapped",)
subset.df <- subset(all.mapped.df, Alignment != "Total" & Library == "Size-selected")
subset.df$aligned <- NA
subset.df$aligned[subset.df$Alignment == "unaligned"] <- "Unaligned"
subset.df$aligned[subset.df$Alignment != "unaligned"] <- "Aligned"

unique.exp <- unique(subset.df$Experiment)
subset.df$Experiment <- factor(subset.df$Experiment, levels = c("Spn23F whole genome enrichment", "CBL enrichment", "WGS vs. CBL enrichment", "Partial CBL database", "Full CBL database"))

subset.df$Channel <- as.character(subset.df$Channel)
subset.df$Channel[subset.df$Channel == "Adaptive"] <- "NAS"
subset.df$Channel <- factor(subset.df$Channel, levels = c("NAS", "Control"))

p <- ggplot(subset.df, aes(x = aligned, y = Avg_read_length, group = aligned, colour = aligned)) + theme_light() + xlab("Read Type") + ylab("Average read length (bp)") + facet_grid(Channel~Experiment, scales="free_x") + geom_boxplot()  + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = FALSE, size = 4.5, label.x = 1.0, label.y = 7000)
p
ggsave(file="average_read_length_comparison_size_selected.svg", plot=p, height = 6, width = 15)

subset.df <- subset(all.mapped.df, Alignment != "Total" & Library == "Unselected")
subset.df$aligned <- NA
subset.df$aligned[subset.df$Alignment == "unaligned"] <- "Unaligned"
subset.df$aligned[subset.df$Alignment != "unaligned"] <- "Aligned"

unique.exp <- unique(subset.df$Experiment)
subset.df$Experiment <- factor(subset.df$Experiment, levels = c("Spn23F whole genome enrichment", "CBL enrichment", "Mixed culture Full CBL database"))

subset.df$Channel <- as.character(subset.df$Channel)
subset.df$Channel[subset.df$Channel == "Adaptive"] <- "NAS"
subset.df$Channel <- factor(subset.df$Channel, levels = c("NAS", "Control"))

p <- ggplot(subset.df, aes(x = aligned, y = Avg_read_length, group = aligned, colour = aligned)) + theme_light() + xlab("Read Type") + ylab("Average read length (bp)") + facet_grid(Channel~Experiment, scales="free_x") + geom_boxplot()  + theme(axis.text.x = element_text(size = 14, angle = 45, hjust=1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = FALSE, size = 4, label.x = 1.0, label.y = 4500)
p
ggsave(file="average_read_length_comparison_unselected.svg", plot=p, height = 6, width = 10)

