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

all.mapped.df <- read.xlsx(filename, "bases_mapped",)
all.total.df <- read.xlsx(filename, "bases_total",)

#df$Library <- factor(df$Library, levels = c("Unselected", "Size-selected"))
#df$Mixture <- gsub('  + Spn23F', '', df$Mixture, fixed = TRUE)
#df$Mixture <- gsub(' + Spn23F', '', df$Mixture, fixed = TRUE)
#df$Mixture <- gsub('PMEN1-23F', 'Spn23F', df$Mixture, fixed = TRUE)

experiments <- unique(all.mapped.df$Experiment)
i <- 5

for (i in 1:length(experiments))
{
  e = experiments[i]
  subsample <- subset(all.mapped.df, Experiment == e)
  total.df <- subset(all.total.df, Experiment == e)
  
  subsample$Channel[subsample$Channel == "Adaptive"] <- "NAS"
  subsample$Channel[subsample$Channel == "Control"] <- "Control"
  subsample$Channel <- factor(subsample$Channel, levels = c("NAS", "Control"))
  
  total.df$Channel[total.df$Channel == "Adaptive"] <- "NAS"
  total.df$Channel[total.df$Channel == "Control"] <- "Control"
  total.df$Channel <- factor(total.df$Channel, levels = c("NAS", "Control"))
  
  if (e == "Spn23F whole genome enrichment")
  {
    p <- ggplot(aligned.bases, aes(x = as.factor(Concentration), y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Bases aligning to target (Mb)") + facet_grid(Mixture~Library, scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_WGS_aligned_yield.svg", plot=p, height = 12, width = 12)

    p <- ggplot(total.df, aes(x = as.factor(Concentration), y = Bases_total / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Total bases generated (Mb)") + facet_grid(Mixture~Library, scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_WGS_total_yield.svg", plot=p, height = 12, width = 12)
  }

  if (e == "CBL enrichment")
  {
    aligned.bases <- subset(subsample, Alignment == "23F")
    aligned.bases$Mixture[aligned.bases$Mixture == "Spn23F"] <- "Undiluted"

    p <- ggplot(aligned.bases, aes(x = Mixture, y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Contaminant strain-serotype combination") + ylab("Bases aligning to target (Mb)") + facet_grid(~Library, scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14,angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_CPS_23F_aligned_yield.svg", plot=p, height = 8, width = 12)

    aligned.bases <- subset(subsample, Alignment != "23F")

    p <- ggplot(aligned.bases, aes(x = Mixture, y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Contaminant strain-serotype combination") + ylab("Bases aligning to target (Mb)") + facet_grid(Concentration~Library, scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14,angle = 45, hjust = 1, vjust = 1), strip.text.y = element_text(size = 16), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 16), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_CPS_non23F_aligned_yield.svg", plot=p, height = 8, width = 12)
    
    aligned.bases <- subsample
    aligned.bases <- aligned.bases[grepl(" ", aligned.bases$Mixture),]
    
    p <- ggplot(aligned.bases, aes(x = Channel, y = Bases_mapped / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.5)
    p
    ggsave(file="V12_mixed_serotypes_aligned_yield_mean_comp.svg", plot=p, height = 6, width = 8)
  }
  
  if (e == "WGS vs. CBL enrichment")
  {
    aligned.bases <- subset(subsample, Alignment == "23F" | Alignment == "Spn23F")
    p <- ggplot(aligned.bases, aes(x = Channel, y = Bases_mapped / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(Target~., scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.y = element_text(size = 14), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.5, vjust=1)
    p
    
    ggsave(file="WGS_CPS_yield_comparison_direct_split.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(aligned.bases, aes(x = as.factor(Concentration_perc), y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Bases aligning to target (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_CPS_v_WGS_aligned_yield.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = as.factor(Concentration), y = Bases_total / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Total bases generated (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V12_CPS_v_WGS_total_yield.svg", plot=p, height = 6, width = 12)
  }
  
  if (e == "Partial CBL database")
  {
    aligned.bases <- subset(subsample, Alignment == "23F")
    
    aligned.bases$Aligner[aligned.bases$Aligner == "Minimap2"] <- "Minimap2"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=90%, min. read=50 bp)"] <- "Graph k19 (S=90%)"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=75%, min. read=50 bp)"] <- "Graph k19 (S=75%)"
    aligned.bases$Aligner <- factor(aligned.bases$Aligner, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    
    p <- ggplot(aligned.bases, aes(x = Channel, y = Bases_mapped / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_partialdb_aligned_yield_mean_comp.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = Channel, y = Bases_total / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_partialdb_total_yield_mean_comp.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(aligned.bases, aes(x = as.factor(Concentration), y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Bases aligning to target (Mb)") + facet_grid(Mixture~Aligner, scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_partialdb_aligned_yield.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = as.factor(Concentration), y = Bases_total / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Total bases generated (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_partialdb_total_yield.svg", plot=p, height = 6, width = 12)
  }
  
  if (e == "Full CBL database")
  {
    aligned.bases <- subset(subsample, Alignment == "23F")
    
    aligned.bases$Aligner[aligned.bases$Aligner == "Minimap2"] <- "Minimap2"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=90%, min. read=50 bp)"] <- "Graph k19 (S=90%)"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=75%, min. read=50 bp)"] <- "Graph k19 (S=75%)"
    aligned.bases$Aligner <- factor(aligned.bases$Aligner, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    
    p <- ggplot(aligned.bases, aes(x = Channel, y = Bases_mapped / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_fulldb_aligned_yield_mean_comp.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = Channel, y = Bases_total / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_fulldb_total_yield_mean_comp.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(aligned.bases, aes(x = as.factor(Concentration), y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Bases aligning to target (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_fulldb_aligned_yield.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = as.factor(Concentration), y = Bases_total / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Total bases generated (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_fulldb_total_yield.svg", plot=p, height = 6, width = 12)
  }
  
  if (e == "Mixed culture Full CBL database")
  {
    aligned.bases <- subset(subsample, Alignment == "23F")
    
    aligned.bases$Aligner[aligned.bases$Aligner == "Minimap2"] <- "Minimap2"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=90%, min. read=50 bp)"] <- "Graph k19 (S=90%)"
    aligned.bases$Aligner[aligned.bases$Aligner == "GP (k=19, S=75%, min. read=50 bp)"] <- "Graph k19 (S=75%)"
    aligned.bases$Aligner <- factor(aligned.bases$Aligner, levels = c("Minimap2", "Graph k19 (S=75%)", "Graph k19 (S=90%)"))
    
    p <- ggplot(aligned.bases, aes(x = Channel, y = Bases_mapped / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_Mixed_aligned_yield_mean_comp.svg", plot=p, height = 6, width = 8)
    
    p <- ggplot(total.df, aes(x = Channel, y = Bases_total / 1000000, group = Channel, colour = Channel)) + theme_light() + xlab("Channel") + ylab("Bases aligned to Target (Mb)") + geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) + facet_grid(.~Aligner, scales="free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.position="none") + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + stat_compare_means(paired = TRUE, size = 5.5, label.x = 1.0)
    p
    ggsave(file="V14_Mixed_total_yield_mean_comp.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(aligned.bases, aes(x = as.factor(Concentration), y = Bases_mapped / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Bases aligning to target (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_fulldb_aligned_yield.svg", plot=p, height = 6, width = 12)
    
    p <- ggplot(total.df, aes(x = as.factor(Concentration), y = Bases_total / 1000000, group = Channel, fill = Channel)) + theme_light() + xlab("Proportion of total DNA targeted for enrichment") + ylab("Total bases generated (Mb)") + facet_grid(Mixture~., scales="free_y") + geom_bar(stat="identity", position="dodge")  + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12)) + scale_color_npg()
    p
    ggsave(file="V14_fulldb_total_yield.svg", plot=p, height = 6, width = 12)
  }
}
