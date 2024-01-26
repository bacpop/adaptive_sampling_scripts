library(ggplot2)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggtext)

roundUp <- function(x) 10^round(log10(x))

parse_filename <- function(file, indir)
{
  prefix <- str_remove(file, indir)
  
  params <- str_split(prefix, "/")[[1]]
  
  experiment <- params[1]
  bc.details <- str_split(params[2], "_")[[1]]
  bed.details <- str_split(params[3], "_")[[1]]
  
  barcode <- bc.details[2]
  channel <- bc.details[4]
  type <- bed.details[1]
  
  array <- c(barcode, channel, type, experiment)
  
  array 
}


plot.gaps <- FALSE
indir.top = "./assemblies/Inspector/"

bed.dirs <- Sys.glob(paste(indir.top, "*", sep = ""))

# downsample bed.dirs
#bed.dirs <- bed.dirs[11:length(bed.dirs)]

#indir <- "./assemblies/Inspector/CPS_v_WGS_WGS_only_minimap2_sup"

l <- 1
for (l in 1:length(bed.dirs)){
  #iterate over bed files, identify number of barcodes
  indir <- bed.dirs[l]
  bed_files <- Sys.glob(paste(indir, "*/*/*_ref.bed", sep = ""))
  barcode.list <- list()
  
  experiment <- str_split(indir, "/")[[1]]
  experiment <- experiment[length(experiment)]
  
  for (file in bed_files)
  {
    parse.file <- parse_filename(file, indir)
    barcode.list <- append(barcode.list, parse.file[1])
  }
  
  unique.barcodes <- unique(barcode.list)
  num.barcodes <- length(unique.barcodes)
  
  # generate holder list
  bed.list <- vector("list", num.barcodes * 8)
  base.df <- data.frame(chr = c("NA", "NA"), start = c(-1000,-1000), end = c(-1000,-1000), Numeric = c(0,0), Type = c("NA", "NA"), Var = c("NA","NA"))
  
  for (i in 1:length(bed.list))
  {
    bed.list[[i]] <- base.df
  }
  
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
  
  prev.barcode <- "none"
  list.index <- 1
  i <- 1
  for (i in 1:length(bed_files))
  {
    file <- bed_files[i]
    parse.file <- parse_filename(file, indir)
    curr.barcode <- parse.file[1]
    curr.channel <- parse.file[2]
    curr.type <- parse.file[3]
    
    
    # determine where in list to append
    if (curr.barcode != prev.barcode & prev.barcode != "none")
    {
      list.index <- list.index + 8
    }
    
    prev.barcode <- curr.barcode
    
    # generate placement within bed.list, NAS then control
    # 1,2 are contigs
    # 3,4 are small variants
    # 5,6 are large variants
    
    channel.index <- 0
    
    if (curr.channel == "control")
    {
      channel.index <- 1
    }
    
    type.index <- 0
    
    if (curr.type == "coverage")
    {
      type.index <- 2
    }
    else if (curr.type == "small")
    {
      type.index <- 4
    } else if (curr.type == "structural")
    {
      type.index <- 6
    }
    
    assign.index <- list.index + channel.index + type.index
    
    df <- try(read.table(file, sep = ' ', header = FALSE), silent = TRUE)
    if(inherits(df, "try-error"))
    {
      next
    } else {
      space.count <- max(str_count(df$V1, "\t"))
      
      col.num <- LETTERS[seq(space.count + 1)]
      
      df.tibble <- df %>% separate(V1, col.num, sep="\t")
      
      df <- as.data.frame(df.tibble)
      
      
      # # check if empty
      # if (dim(df)[1] == 0)
      # {
      #   next
      # }
      
      if (type.index == 0)
      {
        bed.type <- "Contig"
        var.details <- "NA"
      } else if (type.index == 4)
      {
        bed.type <- "Small_Var"
        var.details <- df$D
      } else if (type.index == 6)
      {
        bed.type <- "Large_var"
        var.details <- df$D
      }
      
      if (channel.index == 0)
      {
        colour.num <- 0.25
      } else 
      {
        colour.num <- 0.75
      }
      
      # if (type.index == 0)
      # {
      #   #bed.df <- data.frame(chr = df$A[1], start = 0, end = 0, Numeric = 0, Type = "Boundary", Var = "NA")
      #   to.append <- data.frame(chr = df$A, start = as.numeric(df$B), end = as.numeric(df$C), Numeric = 0.25, Type = "Contig", Var = "NA")
      # } else if (type.index == 2)
      # {
      #   to.append <- data.frame(chr = df$A, start = as.numeric(df$B), end = as.numeric(df$C), Numeric = 0.5, Type = "Small_Var", Var = df$D)
      # } else if (type.index == 4)
      # {
      #   to.append <- data.frame(chr = df$A, start = as.numeric(df$B), end = as.numeric(df$C), Numeric = 0.75, Type = "Large_var", Var = df$D)
      # }
      
      # append coverage differently
      if (type.index == 2)
      {
        to.append <- data.frame(chr = df$A, start = as.numeric(df$B), end = as.numeric(df$C), Numeric = as.numeric(df$D), Type = "Coverage", Var = "NA")
      } else {
        to.append <- data.frame(chr = df$A, start = as.numeric(df$B), end = as.numeric(df$C), Numeric = colour.num, Type = bed.type, Var = var.details)
      }
    }
    
    
    
    #bed.df <- rbind(bed.df, to.append)
    bed.list[[assign.index]] <- to.append
  }
  
  # sort reference-based measures
  if(type == "WGS")
  {
    # Spn23F.label.bed <- data.frame(
    #   chr = "NA",
    #   start = c(303559, 123682, 291868, 331533, 1973324, 2062591, 1612429, 1278602, 1487660),
    #   end = c(322212, 125858, 294120, 333692, 1975519, 2065056, 1614471, 1280521, 1526939),
    #   label = c("CBL", "pspA", "pbp2X", "pbp1A", "pbp2A", "pbp1B","pbp2B", "tetM", "MM1")
    # )
    Spn23F.label.bed <- data.frame(
      chr = "NA",
      start = c(303559,1487660, 1207639),
      end = c(322212, 1526939, 1288735),
      label = c("CBL", "\u03C6MM1", "ICE*Sp*23FST81")
    )
  }
  
  # iterate over paired lists to same plot
  barcode.index = 0
  j <- 2
  for (j in 1:num.barcodes)
  {
    references <- unique(bed.list[[barcode.index + 1]]$chr)
    
    
    # skip if no alignments
    if (length(references) == 1)
    {
      if (references == "NA")
      {
        barcode.index <- barcode.index + 8
        next
      }
    }
    
    
    ref.bed <- data.frame(chr=character(), start = integer(), end = integer())
    init.bed <- data.frame(chr=character(), start = integer(), end = integer())
    
    # manually add in references
    if (type == "multiCPS")
    {
      if (unique.barcodes[j] == "barcode01")
      {
        references <- c("23F")
      } else if (unique.barcodes[j] == "barcode02")
      {
        references <- c("19A")
      } else if (unique.barcodes[j] == "barcode03")
      {
        references <- c("19F")
      } else if (unique.barcodes[j] == "barcode04")
      {
        references <- c("03")
      } else if (unique.barcodes[j] == "barcode05")
      {
        references <- c("06B")
      } else if (unique.barcodes[j] == "barcode06")
      {
        references <- c("19F")
      } else if (unique.barcodes[j] == "barcode07")
      {
        references <- c("23F", "19A")
      } else if (unique.barcodes[j] == "barcode08")
      {
        references <- c("23F", "19F")
      } else if (unique.barcodes[j] == "barcode09")
      {
        references <- c("23F", "03")
      } else if (unique.barcodes[j] == "barcode10")
      {
        references <- c("23F", "06B")
      } else if (unique.barcodes[j] == "barcode11")
      {
        references <- c("23F", "19F")
      } else if (unique.barcodes[j] == "barcode12")
      {
        references <- c("23F")
      }
    }
    
    for (ref in references)
    {
      start.pos <- 0
      
      if (ref == "23F")
      {
        end.pos <- 18654
      } else if (ref == "19A")
      {
        end.pos <- 14949
      } else if (ref == "19F")
      {
        end.pos <- 14820
      } else if (ref == "03")
      {
        end.pos <- 4586
      } else if (ref == "06B")
      {
        end.pos <- 15935
      } else if (ref == "FM211187.1")
      {
        end.pos <- 2221315
        Spn23F.label.bed$chr <- ref
      }
      
      to.append <- data.frame(chr=ref, start = start.pos, end = end.pos, Numeric = 0, Type = "NA", Var = "NA")
      ref.bed <- rbind(ref.bed, to.append)
      to.append <- data.frame(chr=ref, start = 0, end = 0, Numeric = 0, Type = "NA", Var = "NA")
      init.bed <- rbind(init.bed, to.append)
    }
    
    contig.bed.df <- list(rbind(bed.list[[barcode.index + 1]], init.bed), rbind(bed.list[[barcode.index + 2]], init.bed))
    coverage.bed.df <- list(rbind(bed.list[[barcode.index + 3]]), rbind(bed.list[[barcode.index + 4]]))
    small.var.bed.df <- list(rbind(bed.list[[barcode.index + 5]], init.bed), rbind(bed.list[[barcode.index + 6]], init.bed))
    large.var.bed.df <- list(rbind(bed.list[[barcode.index + 7]], init.bed), rbind(bed.list[[barcode.index + 8]], init.bed))
    
    # interate over each list, ensuring has some information included
    # TODO need to add in missing references if multiple references present
    for (i in 1:length(contig.bed.df))
    {
      df <- contig.bed.df[[i]]
      
      # check if empty
      if (df$chr[1] == "NA")
      {
        for (k in 1:length(references)){
          ref <- references[k]
          contig.bed.df[[i]]$chr[k] = ref
        }
      }
    }
    for (i in 1:length(small.var.bed.df))
    {
      df <- small.var.bed.df[[i]]
      
      # check if empty
      if (df$chr[1] == "NA")
      {
        for (k in 1:length(references)){
          ref <- references[k]
          small.var.bed.df[[i]]$chr[k] = ref
        }
      }
    }
    for (i in 1:length(large.var.bed.df))
    {
      df <- large.var.bed.df[[i]]
      
      # check if empty
      if (df$chr[1] == "NA")
      {
        for (k in 1:length(references)){
          ref <- references[k]
          large.var.bed.df[[i]]$chr[k] = ref
        }
      }
    }
    
    # normalise coverage
    max.cov <- 1
    
    for (i in 1:length(coverage.bed.df))
    {
      coverage.bed.df[[i]]$Numeric <- coverage.bed.df[[i]]$Numeric / max.cov
    }
    
    coverage.bed.df[[1]]$Channel <- "NAS"
    coverage.bed.df[[2]]$Channel <- "Control"
    coverage.bed.df[[3]] <- subset(rbind(coverage.bed.df[[1]], coverage.bed.df[[2]]), Type != "NA")
    coverage.bed.df[[3]]$pos <- round((coverage.bed.df[[3]]$start + coverage.bed.df[[3]]$end) / 2)
    
    coverage.bed.df[[3]]$Channel <- factor(coverage.bed.df[[3]]$Channel, levels = (c("NAS", "Control")))
    
    coverage.plot <- ggplot(data = coverage.bed.df[[3]], aes(x = pos, y = Numeric, colour = Channel)) + theme_light() + xlab("") + ylab("Read coverage") + stat_smooth(geom='line', alpha=0.7, se=FALSE, linewidth=1, span=0.1) + scale_x_continuous(limits = c(0, end.pos)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18)) + scale_y_continuous(limit=c(0,NA),oob=squish) + theme(plot.margin = unit(c(1,0,0,1), 'lines'))
    coverage.plot
    
    coverage.plot.line <- ggplot(data = coverage.bed.df[[3]], aes(x = pos, y = Numeric, colour = Channel)) + theme_light() + xlab("") + ylab("Read coverage") + geom_line(alpha=0.7, linewidth=1) + scale_x_continuous(limits = c(0, end.pos)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18)) + scale_y_continuous(limit=c(0,NA),oob=squish) + theme(plot.margin = unit(c(1,0,0,1), 'lines'))
    coverage.plot.line
    
    contig.bed.df[[1]]$Channel <- "NAS"
    contig.bed.df[[2]]$Channel <- "Control"
    contig.bed.df[[3]] <- subset(rbind(contig.bed.df[[1]], contig.bed.df[[2]]), Type != "NA")
    
    contig.bed.df[[3]]$Channel <- factor(contig.bed.df[[3]]$Channel, levels = (c("NAS", "Control")))
    
    contig.plot <- ggplot(subset(contig.bed.df[[3]], Numeric !=0), aes(x = as.factor(Numeric), xend = as.factor(Numeric), y = start, yend = end, colour = Channel)) + geom_segment(linewidth=30) + xlab("Assembly Coverage") + ylab("") + coord_flip() + scale_y_continuous(limits = c(0, end.pos)) + theme_light() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,0,1), 'lines'))
    contig.plot
    
    small.var.bed.df[[1]]$Channel <- "NAS"
    small.var.bed.df[[2]]$Channel <- "Control"
    small.var.bed.df[[3]] <- subset(rbind(small.var.bed.df[[1]], small.var.bed.df[[2]]), Type != "NA")
    
    small.var.bed.df[[3]]$Channel <- factor(small.var.bed.df[[3]]$Channel, levels = (c("NAS", "Control")))
    
    # small.var.plot <- ggplot(data = small.var.bed.df[[3]], aes(x = start, colour = Channel, y = after_stat(n * count))) + geom_density(adjust=1/2,alpha = 0.7, linewidth=1) + ylab("No. small errors") + xlab("") + theme(legend.position = "none") + scale_x_continuous(limits = c(0, end.pos)) + theme_light() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,0,1), 'lines'))
    # small.var.plot
    
    small.var.plot <- ggplot(data = small.var.bed.df[[3]], aes(x = start, colour = Channel)) + geom_freqpoly(bins=round(100), alpha = 0.7, linewidth=1) + ylab("No. small errors") + xlab("") + theme(legend.position = "none") + scale_x_continuous(limits = c(0, end.pos)) + theme_light() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,0,1), 'lines'))
    small.var.plot
    
    # access geom_freq_poly data for smoothing
    h_plotdata <- ggplot_build(small.var.plot)$data[[1]]
    
    if (nrow(h_plotdata) == 0)
    {
      x <- seq(0, end.pos, by = end.pos / 99)
      count <- rep(0, 100)
      Channel <- rep("NAS", 100)
      h_plotdata1 <- data.frame(x, count)
      h_plotdata <- h_plotdata1
      h_plotdata$Channel <- "NAS"
      h_plotdata1$Channel <- "Control"
      h_plotdata <- rbind(h_plotdata, h_plotdata1)
    } else if (length(unique(small.var.bed.df[[3]]$Channel)) == 1) {
      # check if only one item present
      colour <- unique(small.var.bed.df[[3]]$Channel)[[1]]
      h_plotdata$Channel <- colour
      
      h_plotdata1 <- h_plotdata
      
      h_plotdata1$Channel <- if (colour == "Control") "NAS" else "Control"
      h_plotdata1$count <- 0
      
      h_plotdata <- rbind(h_plotdata, h_plotdata1)
      
    } else {
      h_plotdata$Channel <- "NAS"
      h_plotdata$Channel[h_plotdata$colour == "#00BFC4"] <- "Control"
    }
    
    small.var.plot2 <- ggplot(data = h_plotdata, aes(x = x, y = count, colour = Channel)) + stat_smooth(geom='line', alpha=0.7, se=FALSE, linewidth=1, span=0.25) + ylab("No. small errors") + xlab("") + theme(legend.position = "none") + scale_x_continuous(limits = c(0, end.pos)) + theme_light() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,0,1), 'lines')) + rremove("legend") + scale_y_continuous(limit=c(0,NA),oob=squish) #+ ylim(limits = c(0, NA))
    small.var.plot2 
    
    scaleFactor <- max(max(h_plotdata$count) / max(coverage.bed.df[[3]]$Numeric), 1)
    scaleFactor <- roundUp(scaleFactor)
    
    cov.small.var.plot <- ggplot() + theme_light() + xlab("") + ylab("Read coverage") + stat_smooth(data = coverage.bed.df[[3]], aes(x = pos, y = Numeric, colour = Channel, linetype="solid"), geom='line', alpha=0.7, se=FALSE, linewidth=1.5, span=0.25, inherit.aes = FALSE) + stat_smooth(data = h_plotdata, aes(x = x, y = count / scaleFactor, colour = Channel, linetype="dashed"), geom='line', alpha=0.7, se=FALSE, linewidth=0.75, span=0.25, inherit.aes = FALSE) + scale_x_continuous(limits = c(0, end.pos)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=14), axis.title.y = element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18), legend.spacing = unit(0, "lines"), legend.box.background = element_rect(fill = "transparent", colour = NA), legend.key = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent")) + scale_y_continuous(limit=c(0,NA),oob=squish, sec.axis = sec_axis(~ . * scaleFactor, name = "No. small errors")) + theme(plot.margin = unit(c(1,0,0,1), 'lines')) + scale_linetype_manual("Statistic", values = c(4, 1), labels = c("Small errors", "Coverage"))
    cov.small.var.plot
    
    large.var.bed.df[[1]]$Channel <- "NAS"
    large.var.bed.df[[2]]$Channel <- "Control"
    large.var.bed.df[[3]] <- subset(rbind(large.var.bed.df[[1]], large.var.bed.df[[2]]), Type != "NA")
    
    large.var.bed.df[[3]]$Channel <- factor(large.var.bed.df[[3]]$Channel, levels = (c("NAS", "Control")))
    
    large.var.plot <- ggplot(subset(large.var.bed.df[[3]], Numeric !=0), aes(x = as.factor(Numeric), xend = as.factor(Numeric), y = start, yend = end, colour = Channel)) + geom_segment(linewidth=30) + xlab("Large errors") + ylab("Locus coordinate (bp)")  + coord_flip() + scale_y_continuous(limits = c(0, end.pos)) + theme_light() + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=14), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,1,1), 'lines'))
    large.var.plot
    
    # plot assembly gaps
    assembly.gaps.df <- large.var.bed.df[[3]][0,]
    
    k <- 1
    
    # issue here, when only single contig
    subset.contig.bed.df <- subset(contig.bed.df[[3]], Numeric !=0 & Channel == "NAS")
    if (nrow(subset.contig.bed.df) == 0)
    {
      assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- c("NA", 0, end.pos, 0.25, "Gap", "NA", "NAS")
    } else if (nrow(subset.contig.bed.df) == 1)
    {
      curr.row <- subset.contig.bed.df[1,]
      if (curr.row$start != 0) {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- 0
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- curr.row$start - 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
      
      if (curr.row$end != end.pos)
      {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- curr.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- end.pos
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
    } else {
      for (k in 1:(nrow(subset.contig.bed.df) - 1))
      {
        curr.row <- subset.contig.bed.df[k,]
        if (k == 1 & curr.row$start != 0)
        {
          assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
          assembly.gaps.df[nrow(assembly.gaps.df),]$start <- 0
          assembly.gaps.df[nrow(assembly.gaps.df),]$end <- curr.row$start - 1
          assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
        }
        
        next.row <- subset.contig.bed.df[k + 1,]
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- curr.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- next.row$start - 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
      # add gap to end
      if (next.row$end != end.pos)
      {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- next.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- next.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- end.pos
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
    }
    
    subset.contig.bed.df <- subset(contig.bed.df[[3]], Numeric !=0 & Channel == "Control")
    if (nrow(subset.contig.bed.df) == 0)
    {
      assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- c("NA", 0, end.pos, 0.75, "Gap", "NA", "Control")
    } else if (nrow(subset.contig.bed.df) == 1)
    {
      curr.row <- subset.contig.bed.df[1,]
      if (curr.row$start != 0) {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- 0
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- curr.row$start - 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
      
      if (curr.row$end != end.pos)
      {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- curr.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- end.pos
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
    } else {
      for (k in 1:(nrow(subset.contig.bed.df) - 1))
      {
        curr.row <- subset.contig.bed.df[k,]
        if (k == 1 & curr.row$start != 0)
        {
          assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
          assembly.gaps.df[nrow(assembly.gaps.df),]$start <- 0
          assembly.gaps.df[nrow(assembly.gaps.df),]$end <- curr.row$start - 1
          assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
        }
        
        next.row <- subset.contig.bed.df[k + 1,]
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- curr.row
        
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- curr.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- next.row$start - 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
      # add gap to end
      if (next.row$end != end.pos)
      {
        assembly.gaps.df[nrow(assembly.gaps.df) + 1,] <- next.row
        assembly.gaps.df[nrow(assembly.gaps.df),]$start <- next.row$end + 1
        assembly.gaps.df[nrow(assembly.gaps.df),]$end <- end.pos
        assembly.gaps.df[nrow(assembly.gaps.df),]$Type <- "Gap"
      }
    }
    
    contig.large.var.df <- rbind(large.var.bed.df[[3]], contig.bed.df[[3]], assembly.gaps.df)
    # add border for geom_segment
    contig.large.var.df[nrow(contig.large.var.df) + 1,] <- c("NA", 0, end.pos, 0.25, "Gap", "NA", "NAS")
    contig.large.var.df[nrow(contig.large.var.df),]$start <- 0
    contig.large.var.df[nrow(contig.large.var.df),]$end <- end.pos
    contig.large.var.df[nrow(contig.large.var.df),]$Numeric <- 0.25
    contig.large.var.df[nrow(contig.large.var.df),]$Type <- "Border"
    contig.large.var.df[nrow(contig.large.var.df),]$Channel <- "NAS"
    contig.large.var.df[nrow(contig.large.var.df) + 1,] <- contig.large.var.df[nrow(contig.large.var.df),]
    contig.large.var.df[nrow(contig.large.var.df),]$Channel <- "Control"
    contig.large.var.df[nrow(contig.large.var.df),]$Numeric <- 0.75
    
    contig.large.var.df <- subset(contig.large.var.df, Numeric !=0)
    contig.large.var.df$Numeric <- as.numeric(contig.large.var.df$Numeric)
    contig.large.var.df$start <- as.numeric(contig.large.var.df$start)
    contig.large.var.df$end <- as.numeric(contig.large.var.df$end)
    
    contig.large.var.df$rect_start <- contig.large.var.df$Numeric - 0.1
    contig.large.var.df$rect_end <- contig.large.var.df$Numeric + 0.3
    
    if (plot.gaps == TRUE)
    {
      contig.large.var.plot <- ggplot() + geom_rect(data = subset(contig.large.var.df, Type == "Gap"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end, fill = Channel)) + geom_rect(data = subset(contig.large.var.df, Type == "Large_var"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end), fill = "black") + xlab("Assembly gaps/errors") + ylab("Locus coordinate (bp)")  + coord_flip() + scale_y_continuous(limits = c(0, end.pos)) + scale_x_continuous(limits = c(0, 1.15)) + theme_light() + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(size=14), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,1,1), 'lines')) + geom_rect(data = subset(contig.large.var.df, Type == "Border"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end), colour = "#989898", linewidth=0.8, fill = NA)
    } else {
      contig.large.var.plot <- ggplot() + geom_rect(data = subset(contig.large.var.df, Type == "Contig"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end, fill = Channel)) + geom_rect(data = subset(contig.large.var.df, Type == "Large_var"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end), fill = "black") + xlab("Aligned contigs") + ylab("Locus coordinate (bp)")  + coord_flip() + scale_y_continuous(limits = c(0, end.pos)) + scale_x_continuous(limits = c(0, 1.15)) + theme_light() + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(size=14), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) + theme(plot.margin = unit(c(1,0,1,1), 'lines')) + geom_rect(data = subset(contig.large.var.df, Type == "Border"), aes(xmin = rect_start, xmax = rect_end, ymin = start, ymax = end), colour = "#989898", linewidth=0.8, fill = NA)
    }
        
    contig.large.var.plot 
    
    # combined.plot <- ggarrange(coverage.plot + rremove("x.text") + rremove("x.title"), small.var.plot2 + rremove("x.text") + rremove("x.title"), contig.large.var.plot, nrow = 3, ncol = 1, common.legend = TRUE, legend = "right", align = "v") 
    # combined.plot
    # 
    # ggsave(file = paste(experiment, unique.barcodes[j], "linear_plot.jpg", sep = "_"), plot=combined.plot, height = 9, width = 12)
    # 
    
    if (type == "WGS")
    {
      cov.small.var.plot2 <- cov.small.var.plot + geom_rect(data = Spn23F.label.bed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.6)
      cov.small.var.plot2
      contig.large.var.plot2 <- contig.large.var.plot + geom_rect(data = Spn23F.label.bed, aes(ymin = start, ymax = end, xmin = -Inf, xmax = Inf), fill = "grey", alpha = 0.6) + geom_richtext(data = Spn23F.label.bed, aes(x = 1.125, y = start, label = label), hjust = 1, size = 4, fill = NA, label.color = NA,)
      contig.large.var.plot2
      
      combined.plot <- ggarrange(cov.small.var.plot2 + rremove("x.text") + rremove("x.title"), contig.large.var.plot2, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right", align = "v") 
      combined.plot
    } else
    {
      combined.plot <- ggarrange(cov.small.var.plot + rremove("x.text") + rremove("x.title"), contig.large.var.plot, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right", align = "v") 
      combined.plot
    }
    
    ggsave(file = paste(experiment, unique.barcodes[j], "linear_plot.svg", sep = "_"), plot=combined.plot, height = 6, width = 10)
    
    barcode.index <- barcode.index + 8
  }
}

