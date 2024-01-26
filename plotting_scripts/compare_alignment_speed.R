library(ggplot2)
library(ggsci)
library(dplyr)

indir = "./compare_alignment_speed/ref_CPS_run/"

FP.file <- Sys.glob(paste(indir, "*FP_time.txt", sep = ""))
TP.file <- Sys.glob(paste(indir, "*TP_time.txt", sep = ""))

FP.df <- read.table(FP.file, sep = "\t", comment.char = "", header = TRUE)
TP.df <- read.table(TP.file, sep = "\t", comment.char = "", header = TRUE)

# assign true positive and true negatives
FP.df$Type <- "TN"
FP.df$Type[FP.df$Rejection == 0] <- "FP"
TP.df$Type <- "TP"
TP.df$Type[TP.df$Rejection == 1] <- "FN"


full.df <- rbind(FP.df, TP.df)
full.df$Tool[full.df$Tool == "Mappy"] <- "Minimap2"
full.df$Time <- full.df$Time * 1000
graph.subsample <- subset(full.df, Tool == "Graph")
mappy.subsample <- subset(full.df, Tool == "Minimap2")


stacked.df <- data.frame(Graph_len = graph.subsample$Seq_len, Graph_time = graph.subsample$Time, Graph_type = graph.subsample$Type, Mappy_len = mappy.subsample$Seq_len, Mappy_time = mappy.subsample$Time, Mappy_type = mappy.subsample$Type)
san.check <- stacked.df$Graph_len == stacked.df$Mappy_len
sum(san.check)

stacked.df$time.diff <- stacked.df$Graph_time - stacked.df$Mappy_time
stacked.df$time.ratio <- stacked.df$Graph_time / stacked.df$Mappy_time

p <- ggplot(stacked.df, aes(x = Mappy_time, y = Graph_time)) + theme_light() + xlab("Mappy alignment time (millisecs)") + ylab("Graph psuedoalignment time (millisecs)") + geom_point(alpha=0.3, colour="blue") + geom_abline(slope=1, intercept = 0, linetype="dashed", colour="grey") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() + scale_y_log10() 
p

p <- ggplot(stacked.df, aes(x = Graph_len, y = time.ratio)) + theme_light() + xlab("Read length (bp)") + ylab("Ratio Graph-Mappy alignment time") + geom_point(alpha=0.3, colour="#2471A3") + geom_abline(slope=1, intercept = 0, linetype="dashed", colour="#grey") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() #+ scale_y_log10()
p

p <- ggplot(stacked.df, aes(x = Graph_len, y = time.diff)) + theme_light() + xlab("Read length (bp)") + ylab("Graph-Mappy alignment time difference (millisecs)") + geom_point(alpha=0.3, colour="#2471A3") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() #+ scale_y_log10()
p

# plot time
p <- ggplot(full.df, aes(x = Seq_len, y = Time, group = Tool, colour = Tool)) + theme_light() + xlab("Position (bp)") + ylab("Adaptive-Control normalised coverage diffence") + geom_point() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() 
p

# boxplot
p <- ggplot(full.df, aes(x = Tool, y = Time, group = Tool, colour = Tool)) + theme_light() + xlab("Tool") + ylab("Alignment time (mililsecs)") + geom_boxplot() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() 
p

p <- ggplot(full.df, aes(x = Tool, y = Time, group = Tool, colour = Tool)) + theme_light() + xlab("Tool") + ylab("Alignment time (mililsecs)") + geom_boxplot(outlier.shape = NA) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() + scale_y_continuous(limits = quantile(full.df$Time, c(0.1, 0.9)))
p

# plot histogram
p <- ggplot(full.df, aes(x = Time, group = Tool, fill = Tool)) + theme_light() + xlab("Alignment time (millisecs)") + ylab("Count") + geom_histogram(bins = 100) + facet_grid(Tool~., scales = "free_y") + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), legend.position = "none") + guides(colour=guide_legend(title="Method")) + scale_color_npg() + scale_x_continuous(limits = c(NA, 0.12))
p
ggsave(file="mappy_vs_graph_alignment_speed.svg", plot=p, width = 8, height = 5)

full.df%>%
  group_by(Tool)%>% 
  summarise(Mean=mean(Time), Max=max(Time), Min=min(Time), Median=median(Time), Std=sd(Time))

#plot lengths of each category
p <- ggplot(full.df, aes(x = Type, y = Seq_len, fill = Tool)) + theme_light() + xlab("Read classification") + ylab("Read length (bp)") + geom_violin() + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) + guides(colour=guide_legend(title="Method")) + scale_color_npg() + stat_summary(fun=mean, colour="black", geom="text", position = position_dodge(0.9), show_guide = FALSE, vjust=-1, size = 4, aes( label=round(after_stat(y), digits=1))) + stat_summary(fun.y=mean, position = position_dodge(0.9), colour="black", geom="point", shape=18, size=2, show_guide = FALSE)
p
ggsave(file="mappy_vs_graph_read_classification.svg", plot=p, width = 8, height = 5)

#summarise results
count.df <- full.df %>%
  group_by(Tool, Type) %>%
  summarise(Total = n())

#subset df
unique.tool <- unique(count.df$Tool)
tool <- "Graph"
unstacked.count.df <- data.frame(Tool = character(), FN = numeric(), FP = numeric(), TN = numeric(), TP = numeric())

for (tool in unique.tool)
{
  subsample <- subset(count.df, Tool == tool)
  FN. <- subsample$Total[subsample$Type == "FN"]
  TN. <- subsample$Total[subsample$Type == "TN"]
  TP. <- subsample$Total[subsample$Type == "TP"]
  FP. <- subsample$Total[subsample$Type == "FP"]
  
  to.append <- data.frame(Tool = tool, FN = FN., FP = FP., TN = TN., TP = TP.)
  unstacked.count.df <- rbind(unstacked.count.df, to.append)
}

# negatives are sum of false positives and true negatives
unstacked.count.df$negatives <- unstacked.count.df$FP + unstacked.count.df$TN
# positives are sum of true positives and false negatives
unstacked.count.df$positives <- unstacked.count.df$TP + unstacked.count.df$FN

# calculate rates
unstacked.count.df$FP.rate <- unstacked.count.df$FP / unstacked.count.df$negatives
unstacked.count.df$TP.rate <- unstacked.count.df$TP / unstacked.count.df$positives
unstacked.count.df$FN.rate <- unstacked.count.df$FN / unstacked.count.df$positives
unstacked.count.df$TN.rate <- unstacked.count.df$TN / unstacked.count.df$negatives
unstacked.count.df$f1 = (2 * unstacked.count.df$TP) / ((2 * unstacked.count.df$TP) + unstacked.count.df$FP + unstacked.count.df$FN)

unstacked.count.df$precision <- unstacked.count.df$TP / (unstacked.count.df$TP + unstacked.count.df$FP)
unstacked.count.df$recall <- unstacked.count.df$TP / (unstacked.count.df$TP + unstacked.count.df$FN)

write.csv(unstacked.count.df, "minimap2_graph_simulated_accuracy.csv", row.names=FALSE)

#unstacked.count.df$f1.2 <- 2 * ((unstacked.count.df$precision * unstacked.count.df$recall) / (unstacked.count.df$precision + unstacked.count.df$recall))

# plot accuracy results
accuracy.file <- read.csv("minimap2_graph_simulated_accuracy.csv", header = 1)
accuracy.subset <- data.frame(Tool = accuracy.file$Tool, TP.rate = accuracy.file$TP.rate, TN.rate = accuracy.file$TN.rate)
accuracy.subset <- cbind(accuracy.subset[1], stack(accuracy.subset[2:3]))
accuracy.subset$ind <- as.character(accuracy.subset$ind)
accuracy.subset$Tool <- as.character(accuracy.subset$Tool)

accuracy.subset$ind[accuracy.subset$ind == "TP.rate"] <- "TP rate"
accuracy.subset$ind[accuracy.subset$ind == "TN.rate"] <- "TN rate"
accuracy.subset$ind <- factor(accuracy.subset$ind, levels = c("TP rate", "TN rate"))

accuracy.subset$Tool[accuracy.subset$Tool == "Graph"] <- "Graph k19 (S=75%)"
accuracy.subset$Tool <- factor(accuracy.subset$Tool, levels = c("Minimap2", "Graph k19 (S=75%)"))

p <- ggplot(accuracy.subset, aes(x=Tool,y=values, fill=Tool)) + geom_col(colour = "black") + facet_grid(ind~., switch = "y") + theme_light() + xlab("Tool") + ylab("") + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 18, colour = "black", face = "bold"), strip.background = element_rect(fill = "white"), legend.position = "none", strip.placement = "outside") + scale_color_npg()
p
ggsave(file="nanosim_accuracy_comparison.svg", plot=p, height = 8, width = 12)
